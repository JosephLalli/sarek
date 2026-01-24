# Specification: Test Adequacy Evaluation (Pangenome & DeepVariant)

## Overview
This track focuses on evaluating and strengthening the testing infrastructure for the newly integrated pangenome-aware features in the Sarek pipeline. The goal is to ensure that the logic for conditional subworkflow execution, tool configurations (vg giraffe, DeepVariant), and output publishing are robust and correctly implemented.

## Functional Requirements
- **Logic Verification:** Confirm that the pipeline correctly branches into pangenome-aware subworkflows based on input parameters (e.g., `--aligner giraffe`).
- **Configuration Validation:** Verify that correct command-line arguments and resource allocations are passed to `vg giraffe` and `DeepVariant`.
- **Integration Testing:** Execute the pipeline end-to-end using both standard nf-core test datasets and downsampled real-world data (e.g., chr21) to ensure valid output generation.
- **Output Audit:** Validate that the directory structure and published files match expectations for pangenome-aware runs (e.g., BAM/CRAM files, VCFs).

## Non-Functional Requirements
- **Test Efficiency:** Utilize dry-runs and stubs where possible to minimize compute usage during logic verification.
- **Reproducibility:** Ensure test datasets and configurations are documented and easily re-runnable.

## Acceptance Criteria
1.  **Dry-run Success:** Nextflow dry-runs (using `-preview` or examining command logs) confirm that `vg giraffe` and `DeepVariant` receive correct inputs/settings when pangenome mode is enabled.
2.  **Integration Pass:** A full integration test run (nf-test or local execution) with small test data completes successfully with an exit code of 0.
3.  **Output Verification:** Physical inspection of the results directory confirms the presence of expected alignment and variant calling artifacts.
4.  **Test Coverage:** Documentation or nf-test reports show coverage of the pangenome-specific logic branches.

## Out of Scope
- Performance benchmarking or optimization of the tools themselves.
- Adding completely new bioinformatics tools (this is about testing existing ones).
- Production-scale validation (full-sized datasets).
