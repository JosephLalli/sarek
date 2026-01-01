# Feature Backlog: Legacy Gap Analysis

Generated: 2024-12-31
Status: Review Complete

---

## Configuration Requirements (HIGH PRIORITY)

### T2T Reference Support
- [ ] Restore T2T-specific configuration settings/defaults from legacy configs
- [ ] snpEff CHM13v2.110 genome configuration
- [ ] T2T-specific interval files and resources

### Pangenome Configuration
- [ ] Support for **pangenome v1.1** (current/legacy)
- [ ] Support for **pangenome v2** (new release)
- [ ] Pangenome-specific default parameters
- [ ] Reference path configurations for both versions

### Parabricks Integration
- [ ] Review legacy Parabricks implementation vs current Sarek
- [ ] Determine if custom implementation preferences should be retained

---

## KEEP: Modules to Integrate

### Pangenome Tools (11 modules) - *Need updates for vg version changes*
| Module | Path | Notes |
|--------|------|-------|
| pangenie | `pangenie/pangenie/main.nf` | Core genotyper |
| vg/giraffe | `vg/giraffe/main.nf` | Pangenome aligner |
| vg/surject | `vg/surject/main.nf` | GAM to BAM conversion |
| vg/haplotypes | `vg/haplotypes/main.nf` | Haplotype extraction |
| vg/stats | `vg/stats/main.nf` | Graph statistics |
| vg/deconstruct | `vg/deconstruct/deconstruct_and_filter.nf` | VCF from graph |
| vg/gbwt | `vg/gbwt/main.nf` | GBWT index |
| vg/minimizer | `vg/minimizer/main.nf` | Minimizer index |
| vg/index | `vg/index/main.nf` | Graph indexing |
| vg/paths | `vg/paths/main.nf` | Path operations |
| vg/convert | `vg/convert/main.nf` | Format conversion |

### Phasing Tools (7 modules + 1 new) - *Use nf-core implementations where available*
| Module | Path | Notes |
|--------|------|-------|
| shapeit5/phase_common | `shapeit5/phase_common/main.nf` | Common variant phasing |
| shapeit5/phase_rare | `shapeit5/phase_rare/main.nf` | Rare variant phasing |
| shapeit5/ligate | `shapeit5/ligate/main.nf` | Chunk ligation |
| shapeit5/switch | `shapeit5/switch/main.nf` | Switch error calculation |
| shapeit4/phase_common | `shapeit4/phase_common/main.nf` | Legacy phasing |
| whatshap/phase | `whatshap/phase/main.nf` | Read-backed phasing |
| whatshap/stats | `whatshap/stats/main.nf` | Phasing statistics |
| **hapcut** | *NEW* | To be added - not in legacy |

### Specialty Tools (5 legacy + 3 new)
| Module | Path | Notes |
|--------|------|-------|
| jellyfish/count | `jellyfish/count/main.nf` | K-mer counting |
| kmc/kmc | `kmc/kmc/main.nf` | K-mer counting |
| kmc/kmc_dump | `kmc/kmc_dump/main.nf` | K-mer dump |
| expansionhunter | `expansionhunter/main.nf` | STR caller |
| illumina/hap.py | `illumina/hap.py/main.nf` | Variant benchmarking |
| deepvariant/convert_haploid_regions | `deepvariant/convert_haploid_regions/main.nf` | Haploid conversion |
| **STRling** | *NEW* | STR caller - see https://github.com/quinlan-lab/STRling-nf |
| **GangSTR** | *NEW* | STR caller |

### Subworkflows (3 workflows)
| Workflow | Path | Notes |
|----------|------|-------|
| phase_and_impute | `phase_and_impute.nf` | Main phasing orchestrator |
| pangenie_map_and_call | `pangenie_map_and_call.nf` | PanGenie workflow |
| giraffe_mapping | `giraffe_mapping.nf` | vg giraffe alignment workflow |

---

## REDESIGN: Functionality Needed, New Implementation Required

### Pangenome Preparation (9 modules)
These modules are too "home-grown" and need nf-core-style reimplementation:

| Functionality | Legacy Module | Notes |
|---------------|---------------|-------|
| Personalized GBZ creation | `make_personalized_gbz.nf` | Graph building |
| Personalized GFA creation | `make_personalized_gfa.nf` | Graph building |
| Path addition to GFA | `add_paths_to_gfa.nf` | Graph manipulation |
| Contig splitting | `split_by_contig.nf` | Parallelization |
| Sample extraction from VCF | `split_sample_from_vcf.nf` | VCF processing |
| PAR region handling | `split_PARs_off_X.nf` | X chromosome handling |
| Reference splitting | `split_ref_fasta.nf` | Reference processing |
| VCF contig splitting | `split_vcf_by_contig.nf` | VCF processing |
| Insertion FASTA unification | `unify_insertion_fastas.nf` | SV handling |

---

## IGNORE: Not Needed

### Pangenome Tools (Dropped)
- `vg/pack`, `vg/call`, `vg/snarls`, `vg/autoindex`
- `pangenie/filter_variants`, `pangenie/calc_mendelian_violations`
- `danbing-tk/*` (predict, align)

### Phasing Tools (Dropped)
- `shapeit5/merge_switch_reports`

### Utility Modules (Use nf-core or rebuild as needed)
- All `samtools/*` extensions (use nf-core)
- All `bcftools/*` extensions (rebuild as needed)
- All interval/bed tools (`bedtools/slop`, `d4tools/slop`, `build_intervals/*`, etc.)
- `freebayes/bamleftalign`, `freebayes/left_align_and_sort`
- `calc_stats/*`
- `abra/*`
- `gatk/realignertargetcreator`
- `prep_for_gatk_metrics`

### Subworkflows (Dropped)
- All CSV generators (`mapping_csv.nf`, `markduplicates_csv.nf`, etc.) - current pipeline has equivalents
- `realign_indels.nf` - GATK4 handles internally
- Legacy phasing variants (`phase_and_impute_old.nf`, `*_panel_merging.nf`, `*_read_informed*.nf`)
- Standard variant calling workflows (use current sarek_pangenome)
- `split_out_haploid_vcfs.nf`

### Modified nf-core Modules
- Use current nf-core versions for all ~40+ modified modules
- Exception: Retain snpEff T2T (CHM13v2) configuration approach

---

## Implementation Priority

### Phase 1: Configuration
1. T2T reference support in configs
2. Pangenome v1.1 + v2 support in configs
3. Review Parabricks implementation

### Phase 2: Core Pangenome
1. vg tools integration (update for current vg version)
2. PanGenie integration
3. giraffe_mapping subworkflow

### Phase 3: Phasing
1. SHAPEIT5 modules
2. WhatsHap modules
3. HapCUT (new)
4. phase_and_impute subworkflow

### Phase 4: Specialty Tools
1. ExpansionHunter (legacy)
2. STRling (new) - https://github.com/quinlan-lab/STRling-nf
3. GangSTR (new)
4. hap.py benchmarking
5. K-mer tools (jellyfish, kmc)

### Phase 5: Redesign
1. Pangenome preparation modules (nf-core style)

---

## Notes

- vg tools require updates due to tool API changes in recent versions
- Use nf-core module implementations where available (phasing tools)
- Pangenome v2 may require different index formats than v1.1
- T2T support should be parameterized, not hardcoded
- STR callers (ExpansionHunter, STRling, GangSTR) form a cohesive analysis group
