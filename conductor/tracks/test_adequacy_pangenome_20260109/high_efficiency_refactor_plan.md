# Plan: High-Efficiency CRAM & MarkDup Architecture

## Objective
Refactor the alignment and preprocessing logic to support direct CRAM emission and efficient `samtools markdup` piping during alignment. This optimizes IO/storage by avoiding intermediate BAMs and redundant sort/read operations, specifically prioritizing `VG_GIRAFFE` while enabling this architecture for other aligners.

## Requirements Implementation

### 1. New Parameters & Logic
-   **New Param:** `params.use_samtools_markduplicates` (default: false).
-   **Logic in `FASTQ_PREPROCESS_GATK` / `FASTQ_ALIGN`:**
    -   Determine if **Internal MarkDup** should run:
        -   `run_internal_md = params.use_samtools_markduplicates` (Primary trigger).
        -   Also trigger if `params.aligner == 'giraffe'` (optimizing for VG defaults).
    -   Determine **Output Format**:
        -   If `run_internal_md` is TRUE AND `params.skip_tools` contains 'markduplicates' (GATK MD skipped): **Output CRAM**.
        -   If `params.save_output_as_bam`: **Output BAM**.
        -   Default (GATK MD needed downstream): **Output BAM**.

### 2. Refactor `FASTQ_ALIGN` Workflow
-   **Inputs:** Add `run_markdup`, `output_format` (bam/cram) to inputs.
-   **VG_ALIGN Integration:**
    -   Pass `run_markdup` and `output_format` to `VG_ALIGN`.
    -   `VG_ALIGN` passes these to `VG_GIRAFFE`/`VG_SURJECT`.
-   **Standard Aligners (BWAMEM/DRAGMAP):**
    -   **Constraint:** Existing modules may not support piped `markdup`.
    -   **Action:** Modify/Patch local `BWAMEM` and `DRAGMAP` modules to accept `run_markdup` and build the pipe: `| samtools sort | samtools markdup ... | samtools view -O [fmt]`.
-   **Outputs:** Emit `cram` and `crai` channels alongside `bam`/`bai`.

### 3. Refactor `FASTQ_PREPROCESS_GATK`
-   **Handle CRAM Input:** Update logic to accept CRAM output from `FASTQ_ALIGN`.
-   **Skip GATK MarkDuplicates:**
    -   If `FASTQ_ALIGN` performed markdup (based on `run_internal_md` flag logic), skip `BAM_MARKDUPLICATES`.
    -   Route directly to `CRAM_QC` or Recalibration.

### 4. Configuration & Publishing (`conf/modules/`)
-   **Naming Convention:**
    -   If `run_internal_md` ran: Prefix file with `.md`.
    -   Else: Prefix file with `.sorted`.
-   **Publishing Paths:**
    -   `VG_GIRAFFE` / `VG_SURJECT` / `BWAMEM`:
        -   If `.md.cram`: Publish to `preprocessing/markduplicates/`.
        -   If `.sorted.cram` (mapped): Publish to `preprocessing/mapped/`.
-   **Disable Downstream Publishing:** Ensure `BAM_TO_CRAM_MAPPING` does not re-publish or overwrite if the file is already finalized.

## Implementation Steps

1.  **Refactor `FASTQ_ALIGN`:** Update logic to decide `run_markdup` and `output_format`. Update module calls.
2.  **Patch Modules:** Update `VG_GIRAFFE`, `VG_SURJECT`, and `BWAMEM` local modules to implement the piping logic and filename prefixing.
3.  **Refactor `FASTQ_PREPROCESS_GATK`:** Handle CRAM flow and skipping logic.
4.  **Update Configs:** `aligner.config` for publishing rules.
5.  **Tests:** Verify with `tests/aligner-giraffe.nf.test` (expecting `.md.cram` in `markduplicates/` if configured).

