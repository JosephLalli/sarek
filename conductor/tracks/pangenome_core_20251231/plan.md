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

All modules include: main.nf, environment.yml, meta.yml
Location: modules/local/vg/

## Phase 3: PanGenie Implementation
- [x] Task: Review legacy pangenie module
- [x] Task: Implement pangenie module (pangenie 3.0.2)
    - [x] Strict syntax compliance
    - [x] meta.yml and environment.yml
- [ ] Task: Test with sample data

Location: modules/local/pangenie/

## Phase 4: Subworkflow Integration
- [x] Task: Implement giraffe_mapping subworkflow
    - [x] Input channel design (reads, gbz, dist, min, ref_paths)
    - [x] VG_GIRAFFE -> VG_SURJECT -> SAMTOOLS_SORT/INDEX pipeline
    - [x] Output: gam, bam, bai, bam_bai, reports, versions
- [x] Task: Implement pangenie_genotyping subworkflow
    - [x] Input from reads + reference + panel VCF
    - [x] VCF output with index

Location: subworkflows/local/giraffe_mapping/, subworkflows/local/pangenie_genotyping/

## Phase 5: Main Workflow Integration
- [ ] Task: Add pangenome aligner option to main.nf
    - [ ] Parameter: --aligner giraffe
    - [ ] Conditional workflow paths
- [ ] Task: Add pangenie calling option
    - [ ] Parameter: --tools pangenie
- [ ] Task: Update nextflow_schema.json

## Phase 6: Testing
- [ ] Task: Create test data subset
- [ ] Task: Write nf-test cases for each module
- [ ] Task: Integration test with full workflow
