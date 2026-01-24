# Specification: Long-read PacBio Support

## 1. Overview
This track implements comprehensive support for PacBio long-read data (specifically HiFi/CCS) within the Sarek pipeline. It covers preprocessing, alignment, and variant calling for small variants, structural variants, and tandem repeats.

## 2. Functional Requirements

### 2.1. Input Handling
*   The pipeline MUST accept PacBio data (likely HiFi FASTQ/BAM) via the samplesheet.
*   The system MUST recognize PacBio data types to trigger appropriate default parameters and DeepVariant models.

### 2.2. Preprocessing
*   The pipeline MUST support the following tools, enabled via the `--tools` parameter or similar logic:
    *   **Lima:** For demultiplexing.
    *   **Filtering:** Basic length/quality filtering (e.g., using `chopper`).
*   **Out of Scope:** `ccs` generation (assume input is already HiFi/CCS).

### 2.3. Alignment
*   **Aligner Selection:** Users MUST be able to select the aligner using `--aligner`.
*   **Supported Aligners:**
    *   `minimap2`: Standard long-read aligner (configured for PacBio HiFi).
    *   `giraffe`: Support for `vg giraffe`'s long-read mode.
*   **Exclusion:** `pbmm2` is explicitly excluded for now.

### 2.4. Variant Calling
*   **Small Variants:**
    *   **DeepVariant:** MUST be configured with the specific PacBio model.
*   **Structural Variants:**
    *   **PBSV:** Integrated for SV calling.
*   **Tandem Repeats:**
    *   **TRGT:** Integrated for tandem repeat genotyping.
*   **Caller Selection:** Users MUST be able to select these via `--tools`.

## 3. Non-Functional Requirements
*   **Modularity:** Each tool must be implemented as a separate module/subworkflow to allow mix-and-match usage.
*   **Performance:** Long-read alignment and calling can be resource-intensive; default resource configurations (CPUs/Memory) must be appropriate for these tasks.

## 4. Acceptance Criteria
1.  **Preprocessing Run:** `nextflow run . --input pacbio_sample.bam --tools lima` produces demultiplexed output.
2.  **Alignment Run:** `nextflow run . --input pacbio_sample.fastq --aligner minimap2` produces a valid BAM.
3.  **Variant Calling Run:** `nextflow run . --input pacbio_sample.fastq --aligner minimap2 --tools deepvariant,pbsv,trgt` produces VCFs for small variants, SVs, and TRs.
4.  **Tests:** nf-test cases for each new module and subworkflow.

## 5. Out of Scope
*   `ccs` generation.
*   `isoseq` analysis.
*   `sniffles2` (for now).
*   `pbmm2` aligner.
