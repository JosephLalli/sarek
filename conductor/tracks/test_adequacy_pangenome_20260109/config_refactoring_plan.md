# Configuration Refactoring Plan

## Objective
Standardize `VG_ALIGN` (Giraffe) output organization and file naming to match Sarek pipeline conventions. This ensures seamless integration with the pipeline's downstream processing (preprocessing/mapped) and verification tests.

## Findings
*   **BWA-MEM Behavior:** The standard pipeline aligners (BWA-MEM) output **BAM** files, not CRAM. CRAM conversion happens downstream in `FASTQ_PREPROCESS_GATK`.
*   **GATK Compatibility:** GATK MarkDuplicates requires BAM input.
*   **Giraffe Dedup:** `VG_ALIGN` can perform efficient deduplication (`samtools markdup`) internally. If this is done, we must skip the redundant and expensive `BAM_MARKDUPLICATES` step in `FASTQ_PREPROCESS_GATK`.

## Proposed Changes

### 1. `conf/modules/aligner.config`

**Goal:** Treat `VG_GIRAFFE` and `VG_SURJECT` as standard aligners that output BAMs.

**A. Collective Block Updates:**
Add `VG_GIRAFFE` and `VG_SURJECT` to the standard aligner blocks. This inherits the complex logic for publishing BAMs to `preprocessing/`.

```groovy
withName: 'BWAMEM.*_MEM|DRAGMAP_ALIGN|SENTIEON_BWAMEM|VG_GIRAFFE|VG_SURJECT' {
    publishDir = [
        mode: params.publish_dir_mode,
        path: { "${params.outdir}/preprocessing/" },
        pattern: "*bam",
        saveAs: {
            if (params.save_output_as_bam &&
                (
                    params.save_mapped ||
                    (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) &&
                    !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
                ) && (meta.size * meta.num_lanes == 1)
            ) { "mapped/${meta.id}/${it}" }
            else { null }
        }
    ]
}

withName: 'BWAMEM.*_MEM|DRAGMAP_ALIGN|VG_GIRAFFE|VG_SURJECT' {
    // Default to .sorted prefix.
    // NOTE: Even if internal markdup runs, we output as .sorted so downstream tools treat it as the "mapped" input.
    // The prefix change to .md will happen during CRAM conversion if we skip GATK MD.
    ext.prefix = { params.split_fastq > 1 ? "${meta.id}".concat('.').concat(reads.get(0).name.tokenize('.')[0]) : "${meta.id}.sorted" }
}
```

**B. VG Specific Block Updates:**
Simplify specific blocks.

```groovy
// GIRAFFE (pangenome alignment)
withName: 'VG_GIRAFFE' {
    ext.when   = { params.aligner == 'giraffe' }
    ext.args   = "--output-format GAM"
    ext.args2  = "-O cram" 
    
    // NOTE: BAM/CRAM publishing is handled by the collective BWAMEM block above.
    publishDir = [
        // Standardize GAM path
        [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/preprocessing/gam/${meta.id}/" },
            pattern: "*.gam",
            saveAs: { params.save_mapped ? it : null }
        ],
        // Standardize QC path (tsv, log)
        [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/reports/vg_stats/${meta.id}/" },
            pattern: "*.{tsv,log}"
        ],
        // Publish MarkDup logs if generated
        [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/reports/markduplicates/${meta.id}/" },
            pattern: "*.markdup.log"
        ]
    ]
}
```

### 2. `subworkflows/local/fastq_preprocess_gatk/main.nf`

**Goal:** Logic to skip GATK MarkDuplicates if Giraffe has already done it.

**Logic:**
```groovy
// Check if markduplicates is skipped globally
def skip_global_md = params.skip_tools && params.skip_tools.split(',').contains('markduplicates')

// Check if Giraffe performed internal markduplicates (implied if aligner is giraffe and MD not skipped globally)
def giraffe_did_md = params.aligner == 'giraffe' && !skip_global_md

// Combined condition to skip downstream MD step
// We skip if:
// 1. Globally skipped
// 2. Giraffe did it (so we treat it as "done" and proceed to conversion/QC)
def skip_downstream_md = skip_global_md || giraffe_did_md

// Update the if-condition that decides whether to run MD or simple conversion
if (
    params.save_mapped ||
    (
        skip_downstream_md && 
        !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
    )
) {
    // Run BAM_TO_CRAM_MAPPING (or merge then convert)
    // ...
}

// ... later in the workflow ...
if (skip_downstream_md && ...) {
    // CRAM QC NO MD path
} else {
    // BAM_MARKDUPLICATES path
}
```

### 3. `conf/modules/markduplicates.config`

**Goal:** Ensure `BAM_TO_CRAM_MAPPING` produces files named `.md.cram` if Giraffe did the deduplication, so the final output looks standard.

```groovy
withName: 'BAM_TO_CRAM_MAPPING' {
    // If aligner is giraffe and MD was not skipped (so giraffe did it), name it .md.sorted (or just .md)
    // Else use .sorted (standard raw mapping)
    ext.prefix = { 
        if (params.aligner == 'giraffe' && !(params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) {
            "${meta.id}.md" 
        } else {
            "${meta.id}.sorted" 
        }
    }
    
    // ... logic for publishDir remains same, handled by existing rules ...
}
```