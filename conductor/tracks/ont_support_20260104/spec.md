# Specification: Long-read ONT Support

## 1. Overview
This track implements comprehensive support for Oxford Nanopore Technologies (ONT) long-read data within the Sarek pipeline. It covers the full analysis lifecycle from preprocessing and quality control to alignment and variant calling, adhering to Sarek's existing modular design.

## 2. Functional Requirements

### 2.1. Input Handling
*   The pipeline MUST accept ONT FASTQ files via the samplesheet.
*   The system MUST recognize ONT data types (likely via a new `data_type` column or inference) to trigger appropriate default parameters.
*   **Out of Scope:** Conversion from pod5 to FASTQ is explicitly excluded.

### 2.2. Preprocessing & QC
*   The pipeline MUST support the following tools, enabled via the `--tools` parameter:
    *   **NanoPlot:** For QC visualization (analogous to FastQC).
    *   **Porechop:** For adapter trimming.
    *   **Chopper/Filtlong:** For read length and quality filtering.
*   These steps should be optional and configurable via parameters (e.g., `--min_read_length`).

### 2.3. Alignment
*   **Aligner Selection:** Users MUST be able to select the aligner using `--aligner`.
*   **Supported Aligners:**
    *   `minimap2`: Standard long-read aligner.
    *   `giraffe`: Support for `vg giraffe`'s long-read mode.
*   **Output:** Sorted and indexed BAM/CRAM files.

### 2.4. Variant Calling
*   **Caller Selection:** Users MUST be able to select variant callers using `--tools`.
*   **Supported Callers:**
    *   `pepper_deepvariant`: For integrated small variant calling.
    *   `clair3`: Alternative deep-learning based caller.
*   **Configuration:** Callers MUST be configured with models appropriate for the specific ONT chemistry (e.g., R9.4.1 vs R10.4.1) if provided/inferable.

## 3. Non-Functional Requirements
*   **Modularity:** Each tool must be implemented as a separate module/subworkflow to allow mix-and-match usage.
*   **Performance:** Long-read alignment and calling can be resource-intensive; default resource configurations (CPUs/Memory) must be appropriate for these tasks.

## 4. Acceptance Criteria
1.  **QC Run:** `nextflow run . --input ont_sample.csv --tools nanoplot` produces QC reports.
2.  **Alignment Run:** `nextflow run . --input ont_sample.csv --aligner minimap2` produces a valid BAM.
3.  **Variant Calling Run:** `nextflow run . --input ont_sample.csv --aligner minimap2 --tools pepper_deepvariant` produces a VCF.
4.  **Giraffe Integration:** `nextflow run . --input ont_sample.csv --aligner giraffe` produces a valid BAM/GAM.
5.  **Tests:** nf-test cases for each new module and subworkflow.

## 5. Out of Scope
*   Methylation calling.
*   De novo assembly.
*   Hybrid assembly (short + long read).
*   pod5 -> fastq processing.
