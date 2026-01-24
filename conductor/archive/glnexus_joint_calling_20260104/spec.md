# Specification: GLnexus Joint Calling for DeepVariant

## 1. Overview
This track implements GLnexus as the designated joint variant caller for DeepVariant gVCFs within the Sarek pipeline. This feature mirrors the existing GATK joint genotyping workflow, allowing users to perform population-scale variant calling from DeepVariant gVCFs triggered by the `--joint_germline` flag.

## 2. Functional Requirements

### 2.1. Integration Logic
*   **Trigger:** The workflow MUST be triggered when `--tools deepvariant` (or just `deepvariant` in the tool list) AND `--joint_germline` are specified.
*   **Input Handling:**
    *   If running from the beginning (FASTQ/BAM), DeepVariant MUST be configured to output gVCFs (`.g.vcf.gz`) when `--joint_germline` is active.
    *   If running from a `joint_germline` step (start from gVCFs), the pipeline MUST recognize DeepVariant gVCFs and route them to GLnexus.
*   **Subworkflow Structure:** A new subworkflow `BAM_JOINT_CALLING_GERMLINE_GLNEXUS` (or similar) will be created to encapsulate the logic.

### 2.2. Component Configuration
*   **DeepVariant:**
    *   Ensure proper gVCF generation parameters are set when `joint_germline` is true.
*   **GLnexus:**
    *   Implement the `nf-core/glnexus` module (or a local adaptation if strictly necessary, but prefer nf-core).
    *   Support configuration of GLnexus presets (e.g., `--config DeepVariantWGS`).
    *   Handle memory requirements dynamically or via config.

### 2.3. Output Processing
*   **Raw Output:** Generate the raw multi-sample VCF/BCF from GLnexus.
*   **Normalization:** Apply `bcftools norm` to split multiallelic sites and left-align variants.
*   **Filtering (Optional):**
    *   Provide hooks or logic for optional VQSR (if applicable/desired, though noted as difficult with DeepVariant scores) or `bcftools filter` based on quality/PASS status.
    *   Ensure the final output follows standard Sarek naming: `joint_germline.vcf.gz`.

### 2.4. Multi-sample Merging
*   The system MUST support merging gVCFs from multiple samples into a single cohort VCF.

## 3. Non-Functional Requirements
*   **Performance:** GLnexus should be optimized for memory usage, potentially using its "scratch" directory features if supported/exposed.
*   **Consistency:** The user experience (flags, file names) must align with the existing GATK and Sentieon joint calling workflows.

## 4. Acceptance Criteria
1.  **End-to-End Run:** A command like `nextflow run . --tools deepvariant --joint_germline --input samplesheet.csv` successfully produces a `joint_germline.vcf.gz`.
2.  **Start-from-gVCF:** A command like `nextflow run . --step joint_germline --tools deepvariant --input gvcf_samplesheet.csv` successfully produces a `joint_germline.vcf.gz`.
3.  **Correct Output:** The output VCF contains calls for all input samples and has been normalized.
4.  **Tests:** nf-test cases cover both single-sample (trivial joint call) and multi-sample scenarios.

## 5. Out of Scope
*   Implementation of VQSR models specifically trained for DeepVariant (unless standard GATK VQSR is sufficient).
*   Integration with other callers (e.g., using GLnexus for GATK gVCFs).
