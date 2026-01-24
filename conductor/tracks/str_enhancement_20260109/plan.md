# Plan: STR Enhancement

## Phase 1: Testing & Baseline
- [x] Task: Create integration test case running ExpansionHunter, STRling, and GangSTR on a test CRAM.
- [x] Task: Verify current output format and compatibility. [checkpoint: 3969747]

## Phase 2: Merging Strategy
- [x] Task: Research merging tools (Truvari, custom) for STR VCFs. (Superseded by `str_consensus_merge` track)
- [x] Task: Implement merging process in `STR_ANALYSIS` subworkflow. (Superseded by `str_consensus_merge` track)

## Phase 3: Joint Calling
- [x] Task: Implement STRling joint calling logic. fab1593
- [x] Task: Verify joint calling output. fab1593

## Phase 4: Validation
- [x] Task: Verify STRling test data (investigate why allele1_est is 0.00). (Resolved: Confirmed as artifact of small test datasets in `str_consensus_merge` track)
