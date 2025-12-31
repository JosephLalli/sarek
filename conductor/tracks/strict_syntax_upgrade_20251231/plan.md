# Track Plan: Strict Syntax Upgrade

## Phase 1: Critical Syntax Updates
- [x] Task: Refactor `Channel.` to `channel.` in all subworkflows (Scope: ~60+ occurrences)
- [~] Task: Refactor implicit `it` to explicit closure parameters in subworkflows (Scope: ~50+ closures)
- [ ] Task: Conductor - User Manual Verification 'Critical Syntax Updates' (Protocol in workflow.md)

## Phase 2: Local Module Standardization
- [ ] Task: Add `task.ext.args` and `task.ext.prefix` support to local modules (`add_info_to_vcf`, `reindex_bam`, `create_intervals_bed`)
- [ ] Task: Create `meta.yml` documentation for local modules (`add_info_to_vcf`, `reindex_bam`, `create_intervals_bed`)
- [ ] Task: Conductor - User Manual Verification 'Local Module Standardization' (Protocol in workflow.md)

## Phase 3: Configuration & Metadata
- [ ] Task: Fix minor formatting issues in `base.config` and `reindex_bam` module
- [ ] Task: Add `environment.yml` to nf-core modules lacking conda support (`parabricks/fq2bam`, `deepvariant/rundeepvariant`)
- [ ] Task: Add resource labels to nf-core index modules (`bwamem2/index`, `bwa/index`)
- [ ] Task: Conductor - User Manual Verification 'Configuration & Metadata' (Protocol in workflow.md)
