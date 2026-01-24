# Plan: Long-read ONT Support

This plan implements support for ONT data, covering preprocessing, alignment (Minimap2/Giraffe), and variant calling (Pepper-DeepVariant/Clair3).

## Phase 1: Preprocessing & QC
- [ ] Task: Conductor - Verify or install `nf-core/nanoplot` and `nf-core/porechop`.
- [ ] Task: Conductor - Write Tests: Verify ONT QC output with NanoPlot.
- [ ] Task: Conductor - Implement: Integrate NanoPlot into the pipeline.
- [ ] Task: Conductor - Write Tests: Verify adapter trimming with Porechop.
- [ ] Task: Conductor - Implement: Integrate Porechop into the pipeline.
- [ ] Task: Conductor - User Manual Verification 'Preprocessing & QC' (Protocol in workflow.md)

## Phase 2: Alignment Integration
- [ ] Task: Conductor - Verify or install `nf-core/minimap2`.
- [ ] Task: Conductor - Write Tests: Verify ONT alignment with Minimap2.
- [ ] Task: Conductor - Implement: Add Minimap2 as a supported aligner.
- [ ] Task: Conductor - Write Tests: Verify ONT alignment with `vg giraffe` (long-read mode).
- [ ] Task: Conductor - Implement: Add `vg giraffe` long-read support to the alignment subworkflow.
- [ ] Task: Conductor - User Manual Verification 'Alignment Integration' (Protocol in workflow.md)

## Phase 3: Variant Calling Integration
- [ ] Task: Conductor - Verify or install `nf-core/pepper_deepvariant` and `nf-core/clair3`.
- [ ] Task: Conductor - Write Tests: Verify small variant calling with Pepper-DeepVariant.
- [ ] Task: Conductor - Implement: Integrate Pepper-DeepVariant into the variant calling subworkflow.
- [ ] Task: Conductor - Write Tests: Verify small variant calling with Clair3.
- [ ] Task: Conductor - Implement: Integrate Clair3 into the variant calling subworkflow.
- [ ] Task: Conductor - User Manual Verification 'Variant Calling Integration' (Protocol in workflow.md)

## Phase 4: Workflow Orchestration & Validation
- [ ] Task: Conductor - Write Tests: Ensure data type detection correctly triggers ONT-specific parameters.
- [ ] Task: Conductor - Implement: Update `sarek.nf` to orchestrate ONT flows.
- [ ] Task: Conductor - Write Tests: End-to-end ONT test run from FASTQ to VCF.
- [ ] Task: Conductor - Implement: Finalize documentation and usage examples for ONT analysis.
- [ ] Task: Conductor - User Manual Verification 'Workflow Orchestration & Validation' (Protocol in workflow.md)
