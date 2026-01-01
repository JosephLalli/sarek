# Track Specification: Legacy Feature Gap Analysis

## Overview
This track aims to perform a comprehensive, unbiased gap analysis between the legacy codebases (`sarek_first_attempt` and `sarek_JLL`) and the current standard Sarek pipeline. The goal is to uncover *all* custom features, modified logic, and unique workflows that may have been overlooked. This will result in a detailed "Diff Report" used for an interactive "Keep/Ignore" review session with the user.

## Functional Requirements

### 1. Comprehensive Inventory
*   **Scope:** Analyze both `/mnt/ssd/lalli/nf_stage/sarek_first_attempt` and `/mnt/ssd/lalli/nf_stage/sarek_JLL`.
*   **Targets:** `modules/local`, `modules/nf-core`, `subworkflows/local`, `subworkflows/nf-core`, `workflows`.

### 2. Deep Logic Comparison
*   **Unique Components:** List all files present in legacy paths that do not exist in the current project.
*   **Modified Components:** For files that exist in both, analyze the `main.nf` content to detect significant logical deviations (e.g., changed inputs, added process steps, modified commands).
    *   *Note:* Ignore trivial formatting or version string differences.

### 3. Reporting
*   **Categorized Output:** The report must group findings by:
    *   **New Modules:** Completely custom tools.
    *   **Modified nf-core Modules:** Standard tools with custom logic.
    *   **New/Modified Subworkflows:** Custom orchestration logic.
    *   **Workflow-Level Changes:** Modifications to the main entry point (`sarek.nf`).

## Non-Functional Requirements
*   **Independence:** The agent must search blindly for *all* differences, not just pangenome-related ones.
*   **Clarity:** The report should be concise enough for rapid review (e.g., "Module X: Added params A, B").

## Acceptance Criteria
*   **Diff Report Generated:** A markdown file listing all identified gaps.
*   **Interactive Review:** The user has reviewed the list and assigned a "Keep" or "Ignore" status to each item.

## Out of Scope
*   Implementation or porting of any features (this is purely analysis).
