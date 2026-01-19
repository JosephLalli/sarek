# Implementation Plan: STR Consensus Merging

This plan outlines the development of a production-grade STR consensus merging strategy, including format adapters and the priority-based merging tool.

## Phase 1: STRling-to-VCF Adapter Development [checkpoint: bb0eb73]
- [x] Task: Create Red-Phase test for `STRLING_TO_VCF` using a mock `-genotype.txt` and STR catalog.
- [x] Task: Implement `bin/strling_to_vcf.py` to pass the tests.
- [x] Task: Create `modules/local/strling/to_vcf/main.nf` and corresponding `nf-test`.
- [x] Task: Conductor - User Manual Verification 'STRling Adapter' (Protocol in workflow.md)

## Phase 2: Consensus Merger Development
- [ ] Task: Create Red-Phase test for `STR_CONSENSUS_MERGE` with sample VCFs from EH, GS, and SL.
- [ ] Task: Implement `bin/str_consensus_merge.py` with priority and motif-length logic.
- [ ] Task: Create `modules/local/str/consensus_merge/main.nf` and corresponding `nf-test`.
- [ ] Task: Conductor - User Manual Verification 'Consensus Merger' (Protocol in workflow.md)

## Phase 3: Pipeline Integration & Subworkflow
- [ ] Task: Integrate new modules into the `STR_ANALYSIS` subworkflow.
- [ ] Task: Update the pipeline integration test to verify the full flow (EH + GS + SL -> Consensus VCF).
- [ ] Task: Verify that `FORMAT` fields for individual callers (`EH_GT`, etc.) are correctly populated.
- [ ] Task: Conductor - User Manual Verification 'Pipeline Integration' (Protocol in workflow.md)

## Phase 4: Final Validation & Documentation
- [ ] Task: Update `docs/output.md` to describe the new STR consensus VCF fields.
- [ ] Task: Perform a final regression test using the HPRC test dataset.
- [ ] Task: Conductor - User Manual Verification 'Final Validation' (Protocol in workflow.md)
