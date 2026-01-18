# Plan: STR Enhancement

## Phase 1: Testing & Baseline
- [x] Task: Create integration test case running ExpansionHunter, STRling, and GangSTR on a test CRAM.
- [x] Task: Verify current output format and compatibility. [checkpoint: 3969747]

## Phase 2: Merging Strategy
- [ ] Task: Research merging tools (Truvari, custom) for STR VCFs.
- [ ] Task: Implement merging process in `STR_ANALYSIS` subworkflow.

## Phase 3: Joint Calling
- [x] Task: Implement STRling joint calling logic. fab1593
- [x] Task: Verify joint calling output. fab1593

## Phase 4: Validation
- [ ] Task: Verify STRling test data (investigate why allele1_est is 0.00).