# Track Plan: Pangenome Preparation Redesign

## Phase 1: Analysis & Design
- [x] Task: Audit legacy module implementations [53a96e5ac]
    - [x] Read each module's code
    - [x] Document actual tool commands used
    - [x] Identify external dependencies
    - [x] **Crucial:** Confirm if modules generate inputs for phasing or joint calling pipelines.
- [x] Task: Check nf-core for existing equivalents [53a96e5ac]
    - [x] bcftools/view for sample extraction
    - [x] samtools/faidx for FASTA splitting
    - [x] tabix for VCF operations
- [x] Task: Design module architecture [53a96e5ac]
    - [x] Define input/output channel structures
    - [x] Identify shared patterns
    - [x] Plan module naming convention

## Phase 2: Splitting Utilities
- [x] Task: Implement split_vcf_by_region [ARCHIVED]
    - [x] Generic VCF splitting by BED/interval
    - [x] Use bcftools view or tabix
    - [x] Handle index files
- [ ] Task: Implement split_fasta_by_region [SKIP]
    - [ ] Generic FASTA splitting
    - [ ] Use samtools faidx
- [ ] Task: Implement extract_sample_from_vcf [SKIP]
    - [ ] Single sample extraction
    - [ ] Use bcftools view -s
- [ ] Task: Implement handle_par_regions [SKIP]
    - [ ] Reference-agnostic PAR handling
    - [ ] Configurable PAR coordinates

## Phase 3: Graph Preparation
- [ ] Task: Implement build_personalized_graph [SKIP]
    - [ ] GBZ output format
    - [ ] Wrap vg construct/autoindex
    - [ ] Sample-specific variant integration
- [ ] Task: Implement add_paths_to_graph [SKIP]
    - [ ] Add reference/sample paths
    - [ ] Wrap vg paths commands
- [ ] Task: Implement merge_insertion_sequences [SKIP]
    - [ ] SV insertion handling
    - [ ] FASTA merging logic

## Phase 4: Subworkflow Design
- [ ] Task: Design pangenome_index_build subworkflow [SKIP]
    - [ ] Input: Reference FASTA + population VCF
    - [ ] Output: Complete pangenome index set
    - [ ] Support v1.1 and v2 formats
- [ ] Task: Design personalized_graph_build subworkflow [SKIP]
    - [ ] Per-sample graph customization
    - [ ] Integration with giraffe alignment

## Phase 5: Implementation
- [ ] Task: Implement all modules with strict syntax [SKIP]
- [ ] Task: Create meta.yml for each module [SKIP]
- [ ] Task: Create environment.yml for each module [SKIP]
- [ ] Task: Add process labels (process_low/medium/high) [SKIP]

## Phase 6: Testing
- [ ] Task: Create minimal test reference [SKIP]
    - [ ] Small chromosome subset
    - [ ] Sample VCF with variants
- [ ] Task: Write nf-test for each module [SKIP]
- [ ] Task: Integration test for full index build [SKIP]
- [ ] Task: Validate output formats [SKIP]

## Phase 7: Documentation
- [ ] Task: Document pangenome index requirements [SKIP]
- [ ] Task: Usage examples for each module [SKIP]
- [ ] Task: Workflow diagrams [SKIP]
