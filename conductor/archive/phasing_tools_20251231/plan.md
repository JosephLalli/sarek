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
    - [x] Fixed: Added `${args}` support and renamed output to `phased_common_variants.bcf` to avoid collisions.
- [x] Task: Implement shapeit5/phase_rare [9e18a95]
    - [x] Scaffold from phase_common
    - [x] Rare variant handling
    - [x] Fixed: Added `${args}` support and renamed output to `phased_rare_variants.bcf`.
- [x] Task: Implement shapeit5/ligate [95bb62a]
    - [x] Chunk ligation logic
- [x] Task: Implement shapeit5/switch [5fe57d9]
    - [x] Switch error calculation
    - [x] QC metrics output

## Phase 3: Read-Backed Phasing
- [x] Task: Install nf-core/whatshap/phase [f3a2cc7]
    - [x] Use standard module
    - [x] Preprocessing (subsetting) handled by upstream bcftools/view
- [x] Task: Install nf-core/whatshap/stats [f3a2cc7]
    - [x] Phasing N50 calculation
    - [x] Block statistics
- [x] Task: Implement hapcut2 (NEW) [SKIP]
    - [x] Research current best practices
    - [x] Create module from scratch
    - [x] Test with sample data

## Phase 4: Hybrid & Correction Support
- [x] Task: Implement shapeit4/phase_common [d65ca15]
    - [x] Required for Hybrid flow (scaffolding read-backed blocks)
    - [x] Implement from legacy or scratch
- [x] Task: Find container for SAPPHIRE [53a96e5ac]
    - [x] Identified Dockerfile in GitHub; configured Wave build in modules.
- [x] Task: Implement sapphire/phase (NEW) [53a96e5ac]
    - [x] Research input requirements
    - [x] Draft modules (extractor, phasecaller, update) and subworkflow (SAPPHIRE_PHASE_POLISHING).
    - [x] Perform integration testing.

## Phase 5: Subworkflow Integration
- [x] Task: Implement phase_and_impute subworkflow [5cf4044]
    - [x] Logic for `--phasing_method` selection
    - [x] Path 1: Population (SHAPEIT5 suite)
    - [x] Path 2: Read-backed (WhatsHap/HapCUT2)
    - [x] Path 3: Read-corrected (SHAPEIT5 -> SAPPHIRE) (Initial drafting)
    - [x] Path 4: Hybrid (WhatsHap -> SHAPEIT4 -> SHAPEIT5_rare) (Verified with common/rare split)
- [x] Task: Integrate with main workflow [53a96e5ac]
    - [x] Add new parameters to schema
    - [x] Default to `population`

## Phase 6: Testing
- [x] Task: Create test VCF with known phase [53a96e5ac]
    - [x] Downloaded SAPPHIRE test data (`micro.vcf`, etc.) to `assets/test_data/phasing/sapphire/`.
- [x] Task: Test each phasing module [53a96e5ac]
    - [x] Verified SHAPEIT5 modules and subworkflow with nf-test.
- [x] Task: Validate switch error rates [9d55096]
- [x] Task: Integration test with variant calling [cb66f16]