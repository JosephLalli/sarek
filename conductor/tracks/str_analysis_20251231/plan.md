# Track Plan: STR Analysis Tools

## Phase 1: Assessment
- [x] Task: Review existing STR tool landscape [checkpoint: 3337732]
- [x] Task: Catalog file requirements [checkpoint: 3337732]
- [x] Task: Review legacy ExpansionHunter module [checkpoint: 3337732]

## Phase 2: ExpansionHunter Implementation
- [x] Task: Update ExpansionHunter module [c18f087]
    - [x] Strict syntax compliance
    - [x] meta.yml and environment.yml
    - [x] Support --sex flag (multi-catalog dropped)
- [x] Task: Create/obtain variant catalogs [3950271]
    - [x] GRCh38 catalog (PlatinumTRs v1.0)
    - [x] T2T catalog (PlatinumTRs v1.0)
    - [x] Clinical vs research catalogs (Pending specific clinical sets)
    - [x] Created `assets/str_catalogs/` and test subsets

## Phase 3: STRling Implementation (NEW)
- [ ] Task: Study STRling-nf reference implementation
    - [ ] Understand workflow structure
    - [ ] Identify module requirements
- [ ] Task: Implement strling/extract module
    - [ ] Per-sample STR signal extraction
- [ ] Task: Implement strling/merge module
    - [ ] Cohort-level merging
- [ ] Task: Implement strling/call module
    - [ ] Outlier detection
- [ ] Task: Create STRling subworkflow
    - [ ] Full pipeline integration

## Phase 4: GangSTR Implementation (NEW)
- [ ] Task: Research GangSTR requirements
    - [ ] Input format requirements
    - [ ] Reference panel needs
- [ ] Task: Implement gangstr module
    - [ ] Strict syntax compliance
    - [ ] Support for different references
- [ ] Task: Test with sample data

## Phase 5: Supporting Tools
- [ ] Task: Implement jellyfish/count
    - [ ] K-mer counting module
- [ ] Task: Implement kmc/kmc and kmc/kmc_dump
    - [ ] Alternative k-mer counter
- [ ] Task: Implement illumina/hap.py
    - [ ] Variant benchmarking
    - [ ] Useful for STR validation
- [ ] Task: Review deepvariant/convert_haploid_regions
    - [ ] May not be needed with current DeepVariant

## Phase 6: Integration
- [ ] Task: Create STR analysis subworkflow
    - [ ] Support multiple STR callers
    - [ ] Merged output format
- [ ] Task: Add to main workflow
    - [ ] Parameter: --tools str
    - [ ] Parameter: --str_caller [expansionhunter|strling|gangstr|all]
- [ ] Task: Update documentation

## Phase 7: Testing
- [ ] Task: Obtain test samples with known STR expansions
- [ ] Task: Test each STR caller independently
- [ ] Task: Validate against known calls
- [ ] Task: Integration test with full workflow
