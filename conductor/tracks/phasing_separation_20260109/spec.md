# Specification: Phasing Pipeline Separation

## Overview
Decouple the `PHASE_AND_IMPUTE` subworkflow from the main Sarek pipeline into a standalone pipeline. Update Sarek to produce the necessary input samplesheet for this new pipeline.

## Functional Requirements
- **Decoupling:** Extract phasing logic into a new independent workflow structure.
- **Input Generation:** Sarek must generate a CSV containing VCFs and BAMs/CRAMs suitable for the phasing pipeline.

## Acceptance Criteria
1.  **Standalone Execution:** The new phasing pipeline runs independently given a samplesheet.
2.  **Sarek Output:** Sarek produces a valid samplesheet for the phasing pipeline.
