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
- [x] Task: Implement `FASTQ_PREPROCESS_PANGENOME`: 
    - Duplicate `FASTQ_PREPROCESS_GATK` structure but use `VG_ALIGN`.
    - Implement logic to use `samtools markdup` inside `VG_ALIGN` pipeline if configured (`skip_markduplicates` + `tools:samtools_markdup`), bypassing `BAM_MARKDUPLICATES`.
    - Remove BQSR steps.
- [x] Task: Integrate `FASTQ_PREPROCESS_PANGENOME` into `workflows/sarek/main.nf`:
    - Replace `FASTQ_PREPROCESS_GIRAFFE` with `FASTQ_PREPROCESS_PANGENOME`.
    - Ensure argument lists match.
- [x] Task: **Cleanup:** Remove `FASTQ_PREPROCESS_GIRAFFE` and `GIRAFFE_MAPPING`.
- [x] Task: **Testing:** Write Comprehensive Integration Tests for `FASTQ_PREPROCESS_PANGENOME`:
    - **Case 1: Direct Surjection (Multi-Lane & CRAM):** Run with 2 lanes, verify merging, verify CRAM output flows through.
    - **Case 2: Explicit Surjection (Classic):** Verify GAM intermediate, `vg stats` execution, MultiQC report presence, and final CRAM.
    - **Case 3: Haplotype Sampling:** Verify successful execution of KMC->GBWT->Haplotypes->Giraffe chain.
    - **Case 4: BAM Output Regression:** Minimal run configured to output BAMs from `VG_ALIGN` to ensure compatibility.
- [x] Task: Conductor - User Manual Verification 'Architecture Refactoring' (Protocol in workflow.md)

## Phase 3: Standardization of Publishing & Configuration
Align VG_ALIGN (Giraffe) output organization and file naming with the standard Sarek pipeline conventions.

- [x] Task: **Refactor Config:** Update `conf/modules/aligner.config` to integrate VG modules into standard publishing blocks.
    - Configured `VG_GIRAFFE` and `VG_SURJECT` to output CRAM if `samtools_markdup` logic is active, else BAM.
    - Set defaults for `samtools markdup` arguments (`-S -d 2500 ...`).
- [x] Task: **Update Tests:** Replaced `tests/aligner-giraffe.nf.test` with comprehensive `tests/fastq_preprocess_pangenome.nf.test` covering output paths.
- [x] Task: **Verify:** Run `tests/fastq_preprocess_pangenome.nf.test` to confirm the pipeline completes and files are found in the standardized locations.
- [x] Task: Conductor - User Manual Verification 'Standardization of Publishing & Configuration' (Protocol in workflow.md)

## Phase 4: Logic & Configuration Verification (Dry-Run/Stub)
Verify that the correct tools are invoked with the correct parameters without requiring full execution.

- [x] Task: Write `nf-test` unit tests (using stubs) to verify that `vg giraffe` receives correct inputs and CLI flags when `--aligner giraffe` is set. (Covered by `tests/fastq_preprocess_pangenome.nf.test`)
- [x] Task: Write `nf-test` unit tests (using stubs) to verify that the DeepVariant subworkflow is correctly triggered and configured in pangenome mode.
- [x] Task: Verify that the pipeline logic prevents illegal combinations (e.g., pangenome alignment with incompatible variant callers). (Covered by Invalid Config test)
- [x] Task: Add PAR bedfiles to iGenomes registry and provide them to DeepVariant.
- [x] Task: Implement Chromosome Y fallback to standard DeepVariant when using v1.1 pangenome.
- [x] Task: Conductor - User Manual Verification 'Logic & Configuration Verification' (Protocol in workflow.md)

## Phase 5: Integration Testing & Output Audit
Execute the pipeline with small datasets to verify end-to-end functionality and file publishing.

- [x] Task: Run full integration tests with standard `nf-core` test datasets to ensure no regressions in pangenome mode. (Covered by `tests/fastq_preprocess_pangenome.nf.test`)
- [x] Task: Run integration tests using the downsampled real-world data (chr21) and verify exit code 0. (Used `micb-kir3dl1` test data)
- [x] Task: Audit the `results/` directory to ensure alignment (CRAM/BAM) and variant calling (VCF) files are published to the expected locations with correct naming.
- [x] Task: Conductor - User Manual Verification 'Integration Testing & Output Audit' (Protocol in workflow.md)

## Phase 6: Final Evaluation & Documentation
Summarize the findings and ensure the testing framework is repeatable.

- [x] Task: Review all test outputs and logs to confirm tool execution order and parameter passing matches the design intent.
- [x] Task: Update the project's testing documentation to include instructions for running these specific pangenome evaluation tests.
- [ ] Task: Conductor - User Manual Verification 'Final Evaluation & Documentation' (Protocol in workflow.md)
