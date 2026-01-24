# Track Specification: Optimize BCFTOOLS Indexing with -W Argument

## 1. Overview
The goal of this track is to optimize the Sarek pipeline by leveraging the internal index generation capability (`-W` / `--write-index`) of `bcftools` commands. This will streamline the pipeline by removing redundant, separate indexing processes (e.g., `TABIX`, `BCFTOOLS_INDEX`) and reducing the total number of jobs submitted.

## 2. Functional Requirements

### 2.1. Module Updates
*   **Target Modules:** All `bcftools` modules used in the pipeline *except* `BCFTOOLS_REHEADER`.
*   **Action:**
    *   Update module configurations (via `modules.config` or `ext.args`) to include the `--write-index` (or `-W`) argument.
    *   **Crucial:** Verify and update the module's `main.nf` output definition to capture the generated index file (e.g., `*.tbi` or `*.csi`) if it is not already captured.

### 2.2. Workflow/Subworkflow Refactoring
*   **Logic Update:** Identify all workflows and subworkflows where a `bcftools` command is immediately followed by an indexing step (like `TABIX` or `BCFTOOLS_INDEX`).
*   **Hard Removal:** Delete the separate indexing process invocations.
*   **Channel Rewiring:** Update the workflow logic to use the index file emitted directly by the upstream `bcftools` process instead of the removed indexing process.

## 3. Non-Functional Requirements
*   **Performance:** Reduction in total job count and potentially reduced wall-clock time due to elimination of queuing/setup time for small indexing jobs.
*   **Compatibility:** Generated indices must be compatible with downstream tools (standard `.tbi` for VCFs, `.csi` where appropriate).

## 4. Acceptance Criteria
*   **Index Generation:** Targeted `bcftools` processes successfully generate an index file alongside the variant file.
*   **Pipeline Cleanliness:** No redundant `TABIX` or `BCFTOOLS_INDEX` jobs are spawned for the outputs of the modified `bcftools` commands.
*   **Regression Success:** Existing `nf-test` snapshots and pipeline integration tests pass, confirming that the changes haven't broken the pipeline's logic or data integrity.

## 5. Out of Scope
*   `BCFTOOLS_REHEADER` module updates.
*   Creation of new test cases (reliance on existing regression suite).
