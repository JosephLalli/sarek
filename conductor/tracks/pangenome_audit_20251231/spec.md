# Track Specification: Pangenome Feature Audit & Architecture

## Overview
This track is the "Investigation and Design" phase for integrating pangenome features from the legacy `sarek_first_attempt` into the current codebase. Instead of immediate implementation, this track focuses on auditing the existing prototype to understand its logic (especially for personalized genomes) and designing a modern, simplified architecture. The final deliverable is a detailed audit report and a lightweight, documented pseudocode skeleton of the proposed implementation.

## Functional Requirements

### 1. Codebase Audit (`sarek_first_attempt`)
*   **Module Inventory:** Identify all new or modified modules (e.g., `vg/giraffe`, `pangenie`, `str_callers`).
*   **Workflow Mapping:** Map the data flow in `workflows/sarek.nf` (or equivalent) in the prototype.
*   **Personalized Reference Analysis:** Document exactly how per-sample references (graphs, indices, VCFs) are passed through the pipeline.

### 2. Architecture Design
*   **Channel Strategy:** Design a simplified, strict-syntax strategy for joining samples with their correct references (global or per-sample).
*   **Metadata Management:** Define how the `meta` map should be structured to handle optional pangenome indices.
*   **Module Interfaces:** Define the inputs/outputs for the new modules to ensure they follow nf-core/sarek standards.

### 3. Deliverables
*   **Audit Report:** A summary of tools, logic, and "lessons learned" from the prototype.
*   **Pseudocode Skeleton:** A lightweight, documented skeleton (Nextflow-like pseudocode) for the proposed `PANGENOME_ALIGN` and `VARIANT_CALLING_ADVANCED` subworkflows.

## Non-Functional Requirements
*   **Simplification:** The proposed design must explicitly aim to reduce complexity compared to the legacy "spaghetti code."
*   **Best Practices:** The design must adhere to strict DSL2 syntax and nf-core modularity principles.

## Acceptance Criteria
*   **Comprehensive Audit:** All relevant pangenome logic in `sarek_first_attempt` is documented.
*   **Approved Design:** The user approves the pseudocode skeleton and architecture for the subsequent implementation track.

## Out of Scope
*   Actual implementation of functional code or tool execution.
*   Fixing bugs in the legacy `sarek_first_attempt` codebase.
