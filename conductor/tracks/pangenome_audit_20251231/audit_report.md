# Audit Report: Pangenome Features in `sarek_first_attempt`

## 1. Module Inventory
The following modules were found in the legacy codebase and identified as critical for pangenome and STR analysis.

### Pangenome Alignment (`modules/local/vg`)
*   **`vg/giraffe`**: Main alignment module. Maps reads to the pangenome graph.
*   **`vg/surject`**: Converts graph alignments (GAM) to linear reference coordinates (BAM).
*   **`vg/index`**, `vg/autoindex`, `vg/pack`: Graph indexing and utility modules.
*   **`vg/deconstruct`**: Variant calling from graph (likely out of scope but noted).
*   **`vg/convert`**, `vg/stats`: Utilities.

### Advanced Variant Calling (`modules/local`)
*   **`pangenie`**: Genotype imputation and variant calling using a pangenome panel.
    *   Submodules: `calc_mendelian_violations`, `filter_pangenie_variants`.
*   **`expansionhunter`**: STR caller.
*   **`danbing-tk`**: STR genotyping tool.
    *   Submodules: `align`, `predict`.
*   **`deepvariant`**: Standard DeepVariant (also present in nf-core, need to check for pangenome-specific configs).

### "Glue" & Utility Modules (`modules/local`)
These modules handle data preparation, especially for "Personalized Genome" workflows:
*   `make_personalized_gbz.nf`, `make_personalized_gfa.nf`
*   `add_paths_to_gfa.nf`
*   `split_by_contig.nf`, `split_ref_fasta.nf`
*   `split_sample_from_vcf.nf`, `split_vcf_by_contig.nf`
*   `unify_insertion_fastas.nf`
*   `shapeit4`, `shapeit5` (Phasing)

---

## 2. Command-Line Parameters

### `vg giraffe`
Located in: `modules/local/vg/giraffe/main.nf`
*   **Inputs:** `gbz_index` (-Z), `xg_index` (-x), `dist_index` (-d), `min_index` (-m), `ref_chrom_names` (--ref-paths).
*   **Command Structure:**
    ```bash
    vg giraffe --sample ${meta.id} \
        $args \
        --ref-paths $ref_chrom_names --prune-low-cplx \
        --read-group $readgroup \
        --output-format BAM \
        -Z $gbz_index -x $xg_index \
        --max-multimaps 2 \
        -d $dist_index -m $min_index \
        -t $runtime_cpus \
        --progress \
        -P \
        --report-name giraffe_${prefix}.tsv \
        -f $reads ...
    ```
*   **Key Observations:**
    *   Forces output format to BAM immediately (no intermediate GAM file handling visible in the script, though `surject` module exists).
    *   Explicitly handles single-end vs paired-end via `-f` flags.
    *   Uses `--ref-paths` and `--prune-low-cplx` by default for surjection.
    *   **Crucial:** The input signature accepts a tuple of indices `tuple path(gbz_index), path(dist_index), path(min_index), path(xg_index)`. This confirms the need for a robust index management strategy in the `meta` map.

### `pangenie`
Located in: `modules/local/pangenie/pangenie/main.nf`
*   **Inputs:** `jf_db` (Jellyfish database), `ref_fasta`, `panel_vcf`.
*   **Command Structure:**
    ```bash
    pangenie -i $jf_db \
             -v $panel_vcf \
             -r $ref_fasta \
             -o ${prefix} \
             -s ${meta.id} \
             ${args} \
             -t ${task.cpus} -g -d -u
    ```
*   **Key Observations:**
    *   Requires a pre-built Jellyfish database (`-i`).
    *   Requires a Panel VCF (`-v`).
    *   Runs with `-g` (generate graph), `-d` (genotype), `-u` (update).
    *   Post-processes output with `bcftools view` to compress to `.vcf.gz`.

---

## 3. Workflow Audit (`subworkflows/local`)

### `VG_GIRAFFE_MAP` (`giraffe_mapping.nf`)
*   **Purpose:** Handles alignment and surjection.
*   **Branching Logic:** Contains a major fork based on `params.surject_to_insertions`.
    *   **Standard Path (`else`):**
        *   Joins reads with global indices (`gbz`, `dist`, `min`, `xg`).
        *   Joins with global reference FASTA.
        *   Runs `VG_GIRAFFE` -> `VG_SURJECT`.
    *   **Personalized Path (`if params.surject_to_insertions`):**
        *   **Dynamic Index Generation:** Uses `KMC` to count kmers from reads.
        *   **`MAKE_PERSONALIZED_GBZ`:** Creates a sample-specific graph using the kmers and the global graph.
        *   **`VG_INDEX` & `VG_MINIMIZER`:** Builds `dist` and `min` indices on the fly for the personalized graph.
        *   **Result:** `vg giraffe` runs against this *new, sample-specific* set of indices.
*   **Data Flow Issue:** The subworkflow expects `VG_GIRAFFE.out.gam` to pass to `VG_SURJECT`, but the module audit showed `vg giraffe` configured to output BAM directly. This discrepancy suggests the legacy code is in a "broken/transitional" state and must be harmonized (likely enforce GAM output for surjection, or skip surjection if BAM is produced).

### `PANGENIE_CALL_STRUCTURAL_VARIANTS` (`pangenie_map_and_call.nf`)
*   **Purpose:** Genotype imputation and SV calling.
*   **Workflow:**
    1.  `JELLYFISH_COUNT`: Generates Kmer profile from reads.
    2.  `PANGENIE`: Runs using Kmers + Panel VCF + Reference.
    3.  `CONVERT_TO_VCF_GZ` (BCFTools): Standardization.
    4.  **Cohort Merging:** Collects *all* sample VCFs and merges them (`MERGE_BCF`).
    5.  `FILTER_PANGENIE_VARIANTS`: Applies custom filtering (likely regression-based models).
*   **Key Insight:** This is designed as a "Joint Genotyping" workflow, not just per-sample calling.

---

## 4. Architecture Design Implications

### Reference Strategy ("Meta-Map")
The audit confirms that pangenome indices are not just static assets but can be dynamic (per-sample).
*   **Proposed `meta.pangenome` Structure:**
    Instead of passing loose files, we should attach a map to the sample `meta`:
    ```groovy
    meta.pangenome = [
        graph: "path/to/sample.gbz", // or global.gbz
        index_dist: "path/to/sample.dist",
        index_min: "path/to/sample.min",
        is_personalized: true/false
    ]
    ```
*   **Subworkflow decoupling:** `PANGENOME_ALIGN` should accept this `meta` structure. The logic to *create* the personalized indices (if needed) should happen in a separate `PREPARE_PANGENOME_INDICES` subworkflow upstream.

### Personalized Logic Trace
The "Personalized Genome" logic is currently embedded deeply within the alignment subworkflow.
*   **Refactor Plan:** Extract `KMC` -> `MAKE_PERSONALIZED_GBZ` -> `VG_INDEX` into a distinct `PREPARE_PERSONALIZED_GRAPH` subworkflow.
*   **Benefit:** This allows the alignment logic to be "dumb" (just take graph + reads), regardless of whether the graph is global or personalized.
