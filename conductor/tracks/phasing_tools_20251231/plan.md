# Track Plan: Phasing Tools

## Phase 1: Assessment
- [x] Task: Check nf-core/modules for existing implementations [5555555]
    - [x] Search for shapeit5, shapeit4, whatshap, hapcut
    - [x] Evaluate quality and compatibility
- [x] Task: Review legacy phasing modules [5555555]
    - [x] Audit sarek_first_attempt implementations
    - [x] Document input/output patterns
- [x] Task: Define phasing strategy options [5555555]
    - [x] Population-based (SHAPEIT5 + panel)
    - [x] Read-backed (WhatsHap/HapCUT)
    - [x] Read-corrected (SAPPHIRE)
    - [x] Hybrid approach (WhatsHap -> SHAPEIT4 -> SHAPEIT5)

## Phase 2: SHAPEIT5 Implementation
- [x] Task: Implement shapeit5/phase_common [12c1277]
    - [x] Strict syntax compliance
    - [x] Chunk-based parallelization
    - [x] Reference panel input handling
- [x] Task: Implement shapeit5/phase_rare [9e18a95]
    - [x] Scaffold from phase_common
    - [x] Rare variant handling
- [x] Task: Implement shapeit5/ligate [95bb62a]
    - [x] Chunk ligation logic
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

## Phase 4: Hybrid & Correction Support
- [ ] Task: Implement shapeit4/phase_common
    - [ ] Required for Hybrid flow (scaffolding read-backed blocks)
    - [ ] Implement from legacy or scratch
- [ ] Task: Implement sapphire/phase (NEW)
    - [ ] Research input requirements
    - [ ] Create module for phase correction

## Phase 5: Subworkflow Integration
- [ ] Task: Implement phase_and_impute subworkflow
    - [ ] Logic for `--phasing_method` selection
    - [ ] Path 1: Population (SHAPEIT5 suite)
    - [ ] Path 2: Read-backed (WhatsHap/HapCUT2)
    - [ ] Path 3: Read-corrected (SHAPEIT5 -> SAPPHIRE)
    - [ ] Path 4: Hybrid (WhatsHap -> SHAPEIT4 -> SHAPEIT5_rare)
- [ ] Task: Integrate with main workflow
    - [ ] Add new parameters to schema
    - [ ] Default to `population`

## Phase 6: Testing
- [ ] Task: Create test VCF with known phase
- [ ] Task: Test each phasing module
- [ ] Task: Validate switch error rates
- [ ] Task: Integration test with variant calling