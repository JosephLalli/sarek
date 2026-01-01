# Track Plan: Legacy Feature Gap Analysis

## Phase 1: Automated Inventory & Diff
- [x] Task: Inventory `sarek_first_attempt`
    - [x] Task: List and hash all relevant `.nf` files in `modules` and `subworkflows`.
    - [x] Task: List and hash all relevant `.nf` files in `modules` and `subworkflows`.
- [x] Task: Inventory `sarek_JLL`
    - [x] Task: Attempted to list and hash `.nf` files in `modules` and `subworkflows` of `/mnt/ssd/lalli/nf_stage/sarek_JLL`. No such directories found. `sarek_JLL` appears to contain reference data (`local_references/`).

- [~] Task: Comparative Analysis
    - [x] Task: Compare legacy inventories against the current `sarek_pangenome` baseline.
    - [x] Task: Read content of "Modified" candidates to confirm logical differences.
- [x] Task: Conductor - User Manual Verification 'Automated Inventory & Diff' (Protocol in workflow.md)

## Phase 2: Report & Selection
- [x] Task: Generate Diff Report
    - [x] Task: Compile findings into `gap_analysis_report.md`.
- [x] Task: Interactive Review Session
    - [x] Task: Present findings to the user.
    - [x] Task: Record "Keep/Ignore" decisions for each item in a `feature_backlog.md` file.
- [x] Task: Conductor - User Manual Verification 'Report & Selection' (Protocol in workflow.md)

## Summary of Decisions
- **KEEP**: 11 vg modules, 7 phasing modules (+hapcut new), 6 specialty tools (+STRling, +GangSTR new), 3 subworkflows
- **REDESIGN**: 9 pangenome preparation modules (need nf-core style)
- **IGNORE**: danbing-tk, legacy utilities, CSV generators, modified nf-core modules
- **CONFIG**: T2T support, pangenome v1.1+v2 support, Parabricks review needed
