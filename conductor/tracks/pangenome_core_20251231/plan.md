# Track Plan: Core Pangenome Tools

## Phase 1: Assessment & Preparation
- [ ] Task: Review current vg version and API changes
    - [ ] Compare legacy vg commands vs current vg CLI
    - [ ] Document breaking changes
- [ ] Task: Review legacy module implementations
    - [ ] Audit each vg module from sarek_first_attempt
    - [ ] Note input/output patterns
- [ ] Task: Check for existing nf-core vg modules
    - [ ] Search nf-core/modules for vg tools
    - [ ] Evaluate adoption vs custom implementation

## Phase 2: vg Module Implementation
- [ ] Task: Implement vg/giraffe
    - [ ] Create module with strict syntax
    - [ ] Add meta.yml documentation
    - [ ] Add environment.yml
- [ ] Task: Implement vg/surject
- [ ] Task: Implement vg/stats
- [ ] Task: Implement vg/convert
- [ ] Task: Implement vg/index
- [ ] Task: Implement vg/gbwt
- [ ] Task: Implement vg/minimizer
- [ ] Task: Implement vg/haplotypes
- [ ] Task: Implement vg/paths
- [ ] Task: Implement vg/deconstruct

## Phase 3: PanGenie Implementation
- [ ] Task: Review legacy pangenie module
- [ ] Task: Implement pangenie module
    - [ ] Strict syntax compliance
    - [ ] meta.yml and environment.yml
- [ ] Task: Test with sample data

## Phase 4: Subworkflow Integration
- [ ] Task: Implement giraffe_mapping subworkflow
    - [ ] Input channel design
    - [ ] Integration with existing alignment paths
    - [ ] Output to variant calling
- [ ] Task: Implement pangenie_map_and_call subworkflow
    - [ ] Input from aligned reads or raw FASTQ
    - [ ] VCF output formatting

## Phase 5: Main Workflow Integration
- [ ] Task: Add pangenome aligner option to main.nf
    - [ ] Parameter: --aligner giraffe
    - [ ] Conditional workflow paths
- [ ] Task: Add pangenie calling option
    - [ ] Parameter: --tools pangenie
- [ ] Task: Update nextflow_schema.json

## Phase 6: Testing
- [ ] Task: Create test data subset
- [ ] Task: Write nf-test cases for each module
- [ ] Task: Integration test with full workflow
