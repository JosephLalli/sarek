# Specification: STR Enhancement

## Overview
This track focuses on improving the Short Tandem Repeat (STR) analysis capabilities of the pipeline. It aims to implement call merging from multiple tools, explore joint calling, and ensure robust test coverage.

## Functional Requirements
- **Call Merging:** Implement a method to merge STR calls from ExpansionHunter, STRling, and GangSTR into a unified VCF or report.
- **Joint Calling:** Implement joint calling for STRling (and others if supported).
- **Testing:** Create integration tests that verify all 3 tools run successfully on CRAM inputs.

## Acceptance Criteria
1.  **Unified Output:** A merged VCF or report containing calls from all enabled STR tools.
2.  **Test Coverage:** Integration tests pass for ExpansionHunter, STRling, and GangSTR using pipeline-generated CRAMs.
