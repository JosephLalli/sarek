# Track Specification: Strict Syntax Upgrade

## Overview
This track executes the "Zazzy Rolling Stearns" plan to upgrade the Sarek codebase to strict Nextflow DSL2 standards. The primary focus is eliminating deprecated syntax (capitalized `Channel.`, implicit `it` closures) and standardizing module structure (explicit `ext.args`, `meta.yml` documentation).

## Goals
1.  **Eliminate Deprecations:** Replace all instances of `Channel.` with `channel.` and ensure all closures use explicit parameters.
2.  **Standardize Modules:** Ensure local modules use `task.ext` pattern for arguments and have comprehensive `meta.yml` documentation.
3.  **Clean Configuration:** Fix formatting issues and ensure nf-core modules have necessary resource labels and environment files.
4.  **Regression Safety:** Ensure no functionality is lost during these refactors by running existing tests.

## Key Changes
*   **Subworkflows:** Refactor `Channel` factories and closure parameters.
*   **Local Modules:** Inject `task.ext.args/prefix` and add `meta.yml`.
*   **Config:** Minor formatting fixes.
*   **nf-core Modules:** Add missing `environment.yml` and resource labels.

## Verification
*   **Syntax Check:** `nextflow run . -preview` (or similar) should not warn about deprecated syntax.
*   **Functionality:** Existing `test_full` and `test_full_germline` profiles must pass.
