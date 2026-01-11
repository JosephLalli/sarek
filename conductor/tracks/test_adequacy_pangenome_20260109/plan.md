# Plan: Test Adequacy Evaluation (Pangenome & DeepVariant)

## Phase 1: Test Data & Baseline Setup
This phase ensures we have the necessary downsampled datasets and a clean baseline for comparison.

- [x] Task: Prepare downsampled test data (e.g., subset FASTQs and reference graphs to chr21) if not already fully staged.
- [x] Task: Create or update `nf-test` snapshots for a baseline pangenome run to capture current behavior.
- [x] Task: Conductor - User Manual Verification 'Test Data & Baseline Setup' (Protocol in workflow.md)

## Phase 2: Architecture Refactoring (VG Integration)
Refactor the alignment architecture to treat VG Giraffe as a first-class aligner within the standard preprocessing workflow.

- [x] Task: Refactor `PREPARE_GENOME`: Bundle Giraffe assets (GBZ, Dist, Min) into a standard `index_alignment` structure.
- [x] Task: **TDD:** Write module tests for `VG_GIRAFFE` verifying the command string for Reheader/Fixmate/Sort/MarkDup pipes. Ensure read count preservation logic is tested (input reads = output alignments).
- [x] Task: Refactor `VG_GIRAFFE` Module: 
    - **Reference:** Use `sarek_first_attempt/modules/local/vg/surject/main.nf` as the template for the pipe logic.
    - Implement direct BAM/CRAM output (no internal `vg surject`).
    - Implement configurable post-processing pipe: `Reheader` -> `Fixmate` (optional) -> `Sort` (optional) -> `MarkDup` (optional).
    - Ensure input is a unified tuple `[meta, reads, gbz, dist, min]`.
- [x] Task: **TDD:** Write module tests for `VG_SURJECT` verifying the command string matches `VG_GIRAFFE`'s pipe logic. Ensure gam record count preservation logic is tested.
- [x] Task: Refactor `VG_SURJECT` Module:
    - Implement the same configurable post-processing pipe: `Reheader` -> `Fixmate` -> `Sort` -> `MarkDup`.
- [x] Task: Implement `VG_ALIGN` Subworkflow (new subworkflow, replaces `GIRAFFE_MAPPING` logic):
    - **Logic:** `Input [meta, reads, global_index]` -> `Personalized Check`.
    - **Personalized:** `KMC` -> `VG_HAPLOTYPES` -> `Join [reads + sample_index]`.
    - **Standard:** `Combine [reads + global_index]`.
    - **Alignment:** `VG_GIRAFFE`.
    - **Surjection Switch:** If `save_gam` (or Classic mode): `VG_GIRAFFE` (GAM) -> `VG_SURJECT` (BAM/CRAM). Run `VG_STATS` on GAM if 'vg_stats' in tools. Merge reports. Else: `VG_GIRAFFE` (BAM/CRAM).
    - **Output:** Standardized `bam`, `bai`, `reports`, `versions` channels matching `FASTQ_ALIGN` interface. (683d739)
- [ ] Task: Refactor `FASTQ_ALIGN`: 
    - Add `VG_ALIGN` as a supported aligner (via `params.aligner == 'giraffe'`).
    - Route `index_alignment` (Giraffe bundle) to `VG_ALIGN`.
    - Mix `VG_ALIGN.out.bam` and `VG_ALIGN.out.reports` into the main output channels.
- [ ] Task: Refactor `FASTQ_PREPROCESS_GATK`: Add logic to detect if `VG_ALIGN` has performed deduplication (similar to `sentieon_dedup` check) and skip `BAM_MARKDUPLICATES` if so.
- [ ] Task: **Testing:** Write Comprehensive Integration Tests for `FASTQ_PREPROCESS_GATK` with `aligner='giraffe'`:
    - **Case 1: Direct Surjection (Multi-Lane & CRAM):** Run with 2 lanes, verify merging, verify CRAM output flows through.
    - **Case 2: Explicit Surjection (Classic):** Verify GAM intermediate, `vg stats` execution, MultiQC report presence, and final CRAM.
    - **Case 3: Haplotype Sampling:** Verify successful execution of KMC->GBWT->Haplotypes->Giraffe chain.
    - **Case 4: BAM Output Regression:** Minimal run configured to output BAMs from `VG_ALIGN` to ensure compatibility.
- [ ] Task: **Cleanup:** Remove `FASTQ_PREPROCESS_GIRAFFE` and `GIRAFFE_MAPPING` once `FASTQ_ALIGN` integration is verified.
- [ ] Task: Conductor - User Manual Verification 'Architecture Refactoring' (Protocol in workflow.md)

## Phase 3: Logic & Configuration Verification (Dry-Run/Stub)
Verify that the correct tools are invoked with the correct parameters without requiring full execution.

- [ ] Task: Write `nf-test` unit tests (using stubs) to verify that `vg giraffe` receives correct inputs and CLI flags when `--aligner giraffe` is set.
- [ ] Task: Write `nf-test` unit tests (using stubs) to verify that the DeepVariant subworkflow is correctly triggered and configured in pangenome mode.
- [ ] Task: Verify that the pipeline logic prevents illegal combinations (e.g., pangenome alignment with incompatible variant callers).
- [ ] Task: Conductor - User Manual Verification 'Logic & Configuration Verification' (Protocol in workflow.md)

## Phase 4: Integration Testing & Output Audit
Execute the pipeline with small datasets to verify end-to-end functionality and file publishing.

- [ ] Task: Run full integration tests with standard `nf-core` test datasets to ensure no regressions in pangenome mode.
- [ ] Task: Run integration tests using the downsampled real-world data (chr21) and verify exit code 0.
- [ ] Task: Audit the `results/` directory to ensure alignment (CRAM/BAM) and variant calling (VCF) files are published to the expected locations with correct naming.
- [ ] Task: Conductor - User Manual Verification 'Integration Testing & Output Audit' (Protocol in workflow.md)

## Phase 5: Final Evaluation & Documentation
Summarize the findings and ensure the testing framework is repeatable.

- [ ] Task: Review all test outputs and logs to confirm tool execution order and parameter passing matches the design intent.
- [ ] Task: Update the project's testing documentation to include instructions for running these specific pangenome evaluation tests.
- [ ] Task: Conductor - User Manual Verification 'Final Evaluation & Documentation' (Protocol in workflow.md)
