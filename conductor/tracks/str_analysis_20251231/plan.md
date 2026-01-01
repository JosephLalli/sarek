# Track Plan: STR Analysis Tools

## Phase 1: Assessment
- [ ] Task: Review existing STR tool landscape
    - [ ] Check nf-core/modules for implementations
    - [ ] Review STRling-nf pipeline structure
    - [ ] Document input requirements per tool
- [ ] Task: Catalog file requirements
    - [ ] Identify catalogs for GRCh38
    - [ ] Identify/create catalogs for T2T
    - [ ] Document catalog format differences
- [ ] Task: Review legacy ExpansionHunter module
    - [ ] Audit sarek_first_attempt implementation
    - [ ] Note configuration patterns

## Phase 2: ExpansionHunter Implementation
- [ ] Task: Update ExpansionHunter module
    - [ ] Strict syntax compliance
    - [ ] meta.yml and environment.yml
    - [ ] Support multiple catalog files
- [ ] Task: Create/obtain variant catalogs
    - [ ] GRCh38 catalog
    - [ ] T2T catalog
    - [ ] Clinical vs research catalogs

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
