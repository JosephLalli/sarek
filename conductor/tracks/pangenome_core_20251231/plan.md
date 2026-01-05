# Track Plan: Core Pangenome Tools

## Phase 1: Assessment & Preparation
- [x] Task: Review current vg version and API changes
    - [x] Legacy uses vg 1.54.0 (jlalli/vg:1.54.0 container)
    - [x] Conda references are outdated (1.41.0, 1.51.0)
    - [x] Need to verify current vg version and update containers
- [x] Task: Review legacy module implementations
    - [x] 15 vg modules found in sarek_first_attempt/modules/local/vg/
    - [x] Common issues: params.* usage, broken stubs, deprecated syntax
    - [x] VG_HAPLOTYPES module incorrectly named as VG_GIRAFFE
- [x] Task: Check for existing nf-core vg modules
    - [x] Only vg/construct, vg/deconstruct, vg/index exist in nf-core
    - [x] Likely outdated - will implement fresh modules

### Legacy Module Audit Summary

| Module | Status | Issues |
|--------|--------|--------|
| vg/giraffe | Review | Broken shell script (orphan `fi`), params.* |
| vg/giraffe_to_gam | Review | Variant of giraffe outputting GAM |
| vg/surject | Review | Complex embedded samtools pipeline |
| vg/stats | OK | Simple, needs cleanup |
| vg/convert | Review | Stub has wrong variable references |
| vg/index | Review | Needs strict syntax |
| vg/gbwt | Review | Multiple use cases (aliased) |
| vg/minimizer | Review | Needs strict syntax |
| vg/haplotypes | Review | Incorrectly named as VG_GIRAFFE |
| vg/paths | Review | Needs strict syntax |
| vg/pack | SKIP | Not needed per backlog |
| vg/call | SKIP | Not needed per backlog |
| vg/snarls | SKIP | Not needed per backlog |
| vg/autoindex | SKIP | Not needed per backlog |
| pangenie | Review | Uses jlalli/pangenie:2.1.1 |

## Phase 2: vg Module Implementation
- [x] Task: Implement vg/giraffe (vg 1.70.0, quay.io/vgteam/vg:v1.70.0)
- [x] Task: Implement vg/surject
- [x] Task: Implement vg/stats
- [x] Task: Implement vg/convert
- [x] Task: Implement vg/index
- [x] Task: Implement vg/gbwt
- [x] Task: Implement vg/minimizer
- [x] Task: Implement vg/haplotypes
- [x] Task: Implement vg/paths
- [x] Task: Implement vg/deconstruct
- [x] Task: Optimize Giraffe Mapping (Direct Surjection)
    - [x] Research `vg giraffe --output-format` vs `VG_SURJECT` module.
    - [x] Implement optional direct alignment to CRAM/BAM (skipping GAM on disk).
    - [x] Ensure reheader/sort/fixmate logic is preserved.
    - [x] Use `--supplemental` flag for direct CRAM output.
- [x] Task: Update `vg/surject` configuration
    - [x] Add `--supplemental` flag to default arguments.
- [x] Task: Implement KMC module
    - [x] Generate `.kff.gz` from reads.
    - [x] Default: min-count 2 (`-ci2`).
    - [x] Stream output through gzip.
- [x] Task: Implement VG_HAPLOTYPES module
    - [x] Standalone module (cleaner than legacy wrap).
    - [x] Support gzipped k-mer input (decompression/streaming).

All modules include: main.nf, environment.yml, meta.yml
Location: modules/local/vg/

## Phase 3: PanGenie Implementation
- [x] Task: Review legacy pangenie module
- [x] Task: Implement pangenie module (4ccbe4f)
    - [x] Strict syntax compliance
    - [x] meta.yml and environment.yml
- [x] Task: Implement JELLYFISH_COUNT module (3e17cd8)
    - [x] Required for PanGenie k-mer genotyping.
- [x] Task: Test PanGenie with sample data (4ccbe4f)
- [x] Task: Implement PanGenie Preprocessing Subworkflow (4ccbe4f)
    - [x] Integrate JELLYFISH_COUNT.
    - [x] Handle VCF conversion/indexing.

Location: modules/local/pangenie/

## Phase 4: Subworkflow Integration
- [~] Task: Implement giraffe_mapping subworkflow
    - [x] Input channel design (reads, gbz, dist, min, ref_paths)
    - [x] VG_GIRAFFE -> VG_SURJECT -> SAMTOOLS_SORT/INDEX pipeline
    - [x] Output: gam, bam, bai, bam_bai, cram, cram_crai, reports, versions
    - [ ] Task: Integrate Modular Personalized Flow:
        - [ ] Optional `KMC` -> `VG_HAPLOTYPES` -> `VG_GIRAFFE`.
    - [ ] Task: Make VG_STATS conditional on `--tools vg_stats` (or similar).
- [x] Task: Implement pangenie_genotyping subworkflow (4ccbe4f)
    - [x] Input from reads + reference + panel VCF
    - [x] Integrate `JELLYFISH_COUNT` preprocessing.
    - [x] Handle VCF output with index and chromosome merging if parallelized.
    - [x] (Note: `FILTER_PANGENIE_VARIANTS` removed per user request).

Location: subworkflows/local/giraffe_mapping/, subworkflows/local/pangenie_genotyping/

## Phase 5: Main Workflow Integration
- [x] Task: Add pangenome aligner option to main.nf
    - [x] Parameter: --aligner giraffe
    - [x] Conditional workflow paths (parabricks/giraffe/GATK branches)
    - [x] Added FASTQ_PREPROCESS_GIRAFFE subworkflow
    - [x] Added pangenome parameters to nextflow.config
    - [x] Added VG_GIRAFFE/SURJECT/STATS configs to aligner.config
- [x] Task: Add pangenie calling option
    - [x] Parameter: --tools pangenie
    - [x] Parameter: --aligner none (skip alignment, pangenie only)
    - [x] Parallel execution: pangenie + alignment when both specified
    - [x] Validation: aligner=none requires pangenie in tools
    - [x] Validation: pangenie requires pangenie_panel_vcf
- [x] Task: Update nextflow_schema.json with pangenome parameters
    - [x] Added pangenome_gbz, pangenome_dist, pangenome_min, pangenome_ref_paths
    - [x] Added pangenie_panel_vcf
    - [x] Added 'giraffe' and 'none' to aligner enum

## Phase 6: Testing & Optimization
- [x] Task: Create test data subset
- [x] Task: Write nf-test case for giraffe aligner
    - [x] Created tests/aligner-giraffe.nf.test
    - [x] Verified execution with local test data
- [ ] Task: Fix Giraffe Integration Test
    - [ ] Update `tests/aligner-giraffe.nf.test` to verify CRAM/CRAI existence.
    - [ ] Ensure `save_mapped` is handled correctly in tests.
- [~] Task: Write nf-test cases for each module
- [ ] Task: Integration test with full workflow
