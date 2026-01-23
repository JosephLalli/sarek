# Implementation Plan: STR Consensus Merging

This plan outlines the development of a production-grade STR consensus merging strategy, including format adapters and the priority-based merging tool.

## Phase 1: STRling-to-VCF Adapter Development [checkpoint: bb0eb73]
- [x] Task: Create Red-Phase test for `STRLING_TO_VCF` using a mock `-genotype.txt` and STR catalog.
- [x] Task: Implement `bin/strling_to_vcf.py` to pass the tests.
- [x] Task: Create `modules/local/strling/to_vcf/main.nf` and corresponding `nf-test`.
- [x] Task: Conductor - User Manual Verification 'STRling Adapter' (Protocol in workflow.md)

## Phase 2: Consensus Merger Development
- [x] Task: Create Red-Phase test for `STR_CONSENSUS_MERGE` with sample VCFs from EH, GS, and SL.
- [x] Task: Implement `bin/str_consensus_merge.py` with priority and motif-length logic.
- [x] Task: Create `modules/local/str/consensus_merge/main.nf` and corresponding `nf-test`.
- [x] Task: Conductor - User Manual Verification 'Consensus Merger' (Protocol in workflow.md)

## Phase 2b: EnsembleTR Compatibility Refinement
- [x] Task: Update `conductor/tracks/str_consensus_merge_20260118/spec.md` with EnsembleTR compatibility, de novo STRling support, and nf-test content assertions.
- [x] Task: Select a public EnsembleTR container image and pin the version; update `conductor/tech-stack.md` before implementation.
- [x] Task: Create or extend STRling fixtures to guarantee at least one novel call and one ALT variant.
- [x] Task: Create Red-Phase nf-test for `STRLING_TO_VCF` asserting content (HipSTR-style REF/ALT, `##command=hipstr`, START/END/PERIOD, REPCN/Q).
- [x] Task: Implement `bin/strling_to_vcf.py` updates for novel calls and HipSTR-compatible alleles to pass tests.
- [x] Task: Create Red-Phase nf-test for `STR_CONSENSUS_MERGE` asserting coordinate overlap and singleton novel retention.
- [x] Task: Implement coordinate-aware merge and normalization to pass tests.
- [x] Task: Add optional `ENSEMBLETR` Nextflow module with a public container and nf-test (if fixtures allow).
- [x] Task: Emit STRling VCFs as `.vcf.gz` with tabix indexes; update STRLING_TO_VCF tests/snapshots.
- [x] Task: Post-process EnsembleTR output labels to replace HipSTR with STRling for clarity.
- [x] Task: Conductor - User Manual Verification 'EnsembleTR Compatibility' (Protocol in workflow.md)

## Phase 3: Pipeline Integration & Subworkflow
- [x] Task: Integrate new modules into the `STR_ANALYSIS` subworkflow.
- [x] Task: Add integration test for EH + GS -> EnsembleTR using STR enhancement chr4 data (placeholder SL VCF).
- [x] Task: Extend EnsembleTR integration test to include STRling via STRLING_TO_VCF wiring.
- [x] Task: Update the pipeline integration test to verify the full flow (EH + GS + SL -> Consensus VCF).
- [x] Task: Verify that `FORMAT` fields for individual callers (`EH_GT`, etc.) are correctly populated.
- [x] Task: Conductor - User Manual Verification 'Pipeline Integration' (Protocol in workflow.md)

## Phase 4: Final Validation & Documentation [checkpoint: e034e02]
- [x] Task: Update `docs/output.md` to describe the new STR consensus VCF fields.
- [x] Task: Perform a final regression test using the HPRC test dataset. (Resolved: STRling allele1_est=0.00 anomaly confirmed as artifact of small test datasets)
- [x] Task: Conductor - User Manual Verification 'Final Validation' (Protocol in workflow.md)
