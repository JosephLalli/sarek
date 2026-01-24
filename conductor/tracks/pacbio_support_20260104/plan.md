# Plan: Long-read PacBio Support

This plan implements support for PacBio HiFi data, covering preprocessing, alignment (Minimap2/Giraffe), and variant calling (DeepVariant/PBSV/TRGT).

## Phase 1: Preprocessing & QC
- [ ] Task: Conductor - Verify or install `nf-core/lima`.
- [ ] Task: Conductor - Write Tests: Verify PacBio demultiplexing with `lima`.
- [ ] Task: Conductor - Implement: Integrate `lima` into the pipeline.
- [ ] Task: Conductor - Write Tests: Verify basic length/quality filtering for PacBio.
- [ ] Task: Conductor - Implement: Add basic filtering logic (e.g., using `chopper`).
- [ ] Task: Conductor - User Manual Verification 'Preprocessing & QC' (Protocol in workflow.md)

## Phase 2: Alignment Integration
- [ ] Task: Conductor - Write Tests: Verify PacBio HiFi alignment with Minimap2.
- [ ] Task: Conductor - Implement: Configure Minimap2 for PacBio HiFi.
- [ ] Task: Conductor - Write Tests: Verify PacBio alignment with `vg giraffe`.
- [ ] Task: Conductor - Implement: Ensure `vg giraffe` support for PacBio in the alignment subworkflow.
- [ ] Task: Conductor - User Manual Verification 'Alignment Integration' (Protocol in workflow.md)

## Phase 3: Variant Calling Integration
- [ ] Task: Conductor - Verify or install `nf-core/pbsv` and `nf-core/trgt`.
- [ ] Task: Conductor - Write Tests: Verify small variant calling with DeepVariant (PacBio model).
- [ ] Task: Conductor - Implement: Update DeepVariant integration to support the PacBio model.
- [ ] Task: Conductor - Write Tests: Verify structural variant calling with `pbsv`.
- [ ] Task: Conductor - Implement: Integrate `pbsv` into the variant calling subworkflow.
- [ ] Task: Conductor - Write Tests: Verify tandem repeat calling with `trgt`.
- [ ] Task: Conductor - Implement: Integrate `trgt` into the variant calling subworkflow.
- [ ] Task: Conductor - User Manual Verification 'Variant Calling Integration' (Protocol in workflow.md)

## Phase 4: Workflow Orchestration & Validation
- [ ] Task: Conductor - Write Tests: Ensure data type detection correctly triggers PacBio-specific parameters.
- [ ] Task: Conductor - Implement: Update `sarek.nf` to orchestrate PacBio flows.
- [ ] Task: Conductor - Write Tests: End-to-end PacBio test run from FASTQ to VCFs.
- [ ] Task: Conductor - Implement: Finalize documentation and usage examples for PacBio analysis.
- [ ] Task: Conductor - User Manual Verification 'Workflow Orchestration & Validation' (Protocol in workflow.md)
