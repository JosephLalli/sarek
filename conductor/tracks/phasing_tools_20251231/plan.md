# Track Plan: Phasing Tools

## Phase 1: Assessment
- [ ] Task: Check nf-core/modules for existing implementations
    - [ ] Search for shapeit5, shapeit4, whatshap, hapcut
    - [ ] Evaluate quality and compatibility
- [ ] Task: Review legacy phasing modules
    - [ ] Audit sarek_first_attempt implementations
    - [ ] Document input/output patterns
- [ ] Task: Define phasing strategy options
    - [ ] Population-based (SHAPEIT5 + panel)
    - [ ] Read-backed (WhatsHap/HapCUT)
    - [ ] Hybrid approach

## Phase 2: SHAPEIT5 Implementation
- [ ] Task: Implement shapeit5/phase_common
    - [ ] Strict syntax compliance
    - [ ] Chunk-based parallelization
    - [ ] Reference panel input handling
- [ ] Task: Implement shapeit5/phase_rare
    - [ ] Scaffold from phase_common
    - [ ] Rare variant handling
- [ ] Task: Implement shapeit5/ligate
    - [ ] Chunk ligation logic
- [ ] Task: Implement shapeit5/switch
    - [ ] Switch error calculation
    - [ ] QC metrics output

## Phase 3: Read-Backed Phasing
- [ ] Task: Implement whatshap/phase
    - [ ] BAM + VCF input
    - [ ] Haplotagging output option
- [ ] Task: Implement whatshap/stats
    - [ ] Phasing N50 calculation
    - [ ] Block statistics
- [ ] Task: Implement hapcut2 (NEW)
    - [ ] Research current best practices
    - [ ] Create module from scratch
    - [ ] Test with sample data

## Phase 4: Legacy Support
- [ ] Task: Implement shapeit4/phase_common
    - [ ] For backwards compatibility
    - [ ] May deprecate in future

## Phase 5: Subworkflow Integration
- [ ] Task: Implement phase_and_impute subworkflow
    - [ ] Chromosome chunking logic
    - [ ] Parallel phasing per chunk
    - [ ] Ligation and merging
    - [ ] Optional imputation step
- [ ] Task: Integrate with main workflow
    - [ ] Parameter: --tools phasing
    - [ ] Parameter: --phasing_method [shapeit5|whatshap|hapcut]

## Phase 6: Testing
- [ ] Task: Create test VCF with known phase
- [ ] Task: Test each phasing module
- [ ] Task: Validate switch error rates
- [ ] Task: Integration test with variant calling
