# Track Plan: Pangenome Preparation Redesign

## Phase 1: Analysis & Design
- [ ] Task: Audit legacy module implementations
    - [ ] Read each module's code
    - [ ] Document actual tool commands used
    - [ ] Identify external dependencies
- [ ] Task: Check nf-core for existing equivalents
    - [ ] bcftools/view for sample extraction
    - [ ] samtools/faidx for FASTA splitting
    - [ ] tabix for VCF operations
- [ ] Task: Design module architecture
    - [ ] Define input/output channel structures
    - [ ] Identify shared patterns
    - [ ] Plan module naming convention

## Phase 2: Splitting Utilities
- [ ] Task: Implement split_vcf_by_region
    - [ ] Generic VCF splitting by BED/interval
    - [ ] Use bcftools view or tabix
    - [ ] Handle index files
- [ ] Task: Implement split_fasta_by_region
    - [ ] Generic FASTA splitting
    - [ ] Use samtools faidx
- [ ] Task: Implement extract_sample_from_vcf
    - [ ] Single sample extraction
    - [ ] Use bcftools view -s
- [ ] Task: Implement handle_par_regions
    - [ ] Reference-agnostic PAR handling
    - [ ] Configurable PAR coordinates

## Phase 3: Graph Preparation
- [ ] Task: Implement build_personalized_graph
    - [ ] GBZ output format
    - [ ] Wrap vg construct/autoindex
    - [ ] Sample-specific variant integration
- [ ] Task: Implement add_paths_to_graph
    - [ ] Add reference/sample paths
    - [ ] Wrap vg paths commands
- [ ] Task: Implement merge_insertion_sequences
    - [ ] SV insertion handling
    - [ ] FASTA merging logic

## Phase 4: Subworkflow Design
- [ ] Task: Design pangenome_index_build subworkflow
    - [ ] Input: Reference FASTA + population VCF
    - [ ] Output: Complete pangenome index set
    - [ ] Support v1.1 and v2 formats
- [ ] Task: Design personalized_graph_build subworkflow
    - [ ] Per-sample graph customization
    - [ ] Integration with giraffe alignment

## Phase 5: Implementation
- [ ] Task: Implement all modules with strict syntax
- [ ] Task: Create meta.yml for each module
- [ ] Task: Create environment.yml for each module
- [ ] Task: Add process labels (process_low/medium/high)

## Phase 6: Testing
- [ ] Task: Create minimal test reference
    - [ ] Small chromosome subset
    - [ ] Sample VCF with variants
- [ ] Task: Write nf-test for each module
- [ ] Task: Integration test for full index build
- [ ] Task: Validate output formats

## Phase 7: Documentation
- [ ] Task: Document pangenome index requirements
- [ ] Task: Usage examples for each module
- [ ] Task: Workflow diagrams
