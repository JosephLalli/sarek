# Track Plan: Pangenome Feature Audit & Architecture

## Phase 1: Investigation & Audit
- [x] Task: Audit `sarek_first_attempt` Modules
    - [x] Task: List all unique modules in `modules/local` and `modules/nf-core` added/modified for pangenome support.
    - [x] Task: Document the command-line parameters used for `vg giraffe` and `pangenie`.
- [x] Task: Audit `sarek_first_attempt` Subworkflows & Workflows
    - [x] Task: Map the channel logic for `vg giraffe` alignment.
    - [x] Task: Trace the "personalized genome" logic (how per-sample files are associated with samples).
- [x] Task: Conductor - User Manual Verification 'Investigation & Audit' (Protocol in workflow.md)

## Phase 2: Architecture Design
- [x] Task: Design "Meta-Map" Reference Strategy
    - [x] Task: Propose a structure for the `meta` map to handle optional `graph`, `index`, and `reference` paths per sample.
- [x] Task: Draft Pseudocode Skeleton
    - [x] Task: Create pseudocode for a simplified `PANGENOME_ALIGN` subworkflow.
    - [x] Task: Create pseudocode for `PANGENIE` and `STR` integration points.
- [x] Task: Conductor - User Manual Verification 'Architecture Design' (Protocol in workflow.md)

## Phase 3: Final Documentation & Transition
- [x] Task: Compile Audit Report
    - [x] Task: Summarize findings and list "spaghetti code" patterns to avoid.
- [x] Task: Final Review
    - [x] Task: Present the skeleton and report to the user for approval.
- [x] Task: Conductor - User Manual Verification 'Final Documentation & Transition' (Protocol in workflow.md)
