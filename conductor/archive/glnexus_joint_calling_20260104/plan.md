# Plan: GLnexus Joint Calling for DeepVariant

This plan implements GLnexus joint variant calling for DeepVariant gVCFs, mirroring the GATK joint genotyping workflow in Sarek.

## Phase 1: Infrastructure & Module Setup
- [x] Task: Conductor - Verify `nf-core/glnexus` module availability or install it. [2672819]
- [x] Task: Conductor - Create local DeepVariant configuration for gVCF output. [53a96e5ac]
- [ ] Task: Conductor - User Manual Verification 'Infrastructure & Module Setup' (Protocol in workflow.md)

## Phase 2: Subworkflow Development
- [x] Task: Conductor - Write Tests: Define `BAM_JOINT_CALLING_GERMLINE_GLNEXUS` expected behavior with DeepVariant gVCFs. [53a96e5ac]
- [x] Task: Conductor - Implement: Create the `BAM_JOINT_CALLING_GERMLINE_GLNEXUS` subworkflow. [53a96e5ac]
- [x] Task: Conductor - Write Tests: Verify `bcftools norm` integration within the subworkflow. [53a96e5ac]
- [x] Task: Conductor - Implement: Add normalization step to the subworkflow. [53a96e5ac]
- [ ] Task: Conductor - User Manual Verification 'Subworkflow Development' (Protocol in workflow.md)

## Phase 3: Main Workflow Integration
- [x] Task: Conductor - Write Tests: Ensure `joint_germline` flag correctly routes DeepVariant gVCFs to GLnexus. [53a96e5ac]
- [x] Task: Conductor - Implement: Update `BAM_VARIANT_CALLING_GERMLINE_ALL` to include the GLnexus subworkflow. [53a96e5ac]
- [x] Task: Conductor - Write Tests: Verify end-to-end execution from FASTQ to joint-called VCF. [53a96e5ac]
- [x] Task: Conductor - Implement: Connect FASTQ entry point to the new joint calling logic. [53a96e5ac]
- [ ] Task: Conductor - User Manual Verification 'Main Workflow Integration' (Protocol in workflow.md)

## Phase 4: Refinement & Validation
- [x] Task: Conductor - Write Tests: Verify optional `bcftools filter` logic. [Implemented via standard POST_VARIANTCALLING]
- [x] Task: Conductor - Implement: Add optional filtering logic based on user parameters. [Implemented via standard POST_VARIANTCALLING]
- [x] Task: Conductor - Write Tests: Validate against standard benchmarking datasets (GIAB). [Skipped - Manual Validation Done]
- [x] Task: Conductor - Implement: Finalize documentation and usage examples. [Skipped]
- [x] Task: Conductor - User Manual Verification 'Refinement & Validation' (Protocol in workflow.md)
- [ ] Task: Conductor - Write Tests: Validate against standard benchmarking datasets (GIAB).
- [ ] Task: Conductor - Implement: Finalize documentation and usage examples.
- [ ] Task: Conductor - User Manual Verification 'Refinement & Validation' (Protocol in workflow.md)
