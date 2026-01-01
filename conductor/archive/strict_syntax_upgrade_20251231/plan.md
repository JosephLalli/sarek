# Track Plan: Strict Syntax Upgrade

## Completion Notes
**Status:** Completed
**Date:** 2025-12-31

**Important Note on Syntax Parser:**
Although the codebase has been upgraded to strict syntax, running with `NXF_SYNTAX_PARSER=v2` is currently **broken for test runs**. It is recommended to run tests **without** setting `NXF_SYNTAX_PARSER=v2`. Testing is fully functional with the default parser.

## Phase 1: Critical Syntax Updates
- [x] Task: Refactor `Channel.` to `channel.` in all subworkflows (Scope: ~60+ occurrences)
- [x] Task: Refactor implicit `it` to explicit closure parameters in subworkflows (Scope: ~50+ closures)
- [x] Task: Conductor - User Manual Verification 'Critical Syntax Updates' (Protocol in workflow.md)

## Phase 2: Local Module Standardization
- [x] Task: Add `task.ext.args` and `task.ext.prefix` support to local modules (`add_info_to_vcf`, `reindex_bam`, `create_intervals_bed`)
- [x] Task: Create `meta.yml` documentation for local modules (`add_info_to_vcf`, `reindex_bam`, `create_intervals_bed`)
- [x] Task: Conductor - User Manual Verification 'Local Module Standardization' (Protocol in workflow.md)

## Phase 3: Configuration & Metadata
- [x] Task: Fix minor formatting issues in `base.config` and `reindex_bam` module
- [x] Task: Add `environment.yml` to nf-core modules lacking conda support (`parabricks/fq2bam`, `deepvariant/rundeepvariant`)
- [x] Task: Add resource labels to nf-core index modules (`bwamem2/index`, `bwa/index`)
- [x] Task: Conductor - User Manual Verification 'Configuration & Metadata' (Protocol in workflow.md)