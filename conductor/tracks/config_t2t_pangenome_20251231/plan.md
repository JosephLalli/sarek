# Track Plan: T2T and Pangenome Configuration

## Phase 1: Audit Existing Configuration
- [x] Task: Review current `conf/` structure
    - [x] Identify where genome-specific configs live
    - [x] Review igenomes.config pattern
- [x] Task: Locate legacy T2T configuration
    - [x] Find T2T settings in sarek_first_attempt
    - [x] Document what was configured
- [x] Task: Locate legacy pangenome configuration
    - [x] Find pangenome settings in sarek_first_attempt
    - [x] Document v1.1 requirements

## Phase 2: T2T Reference Configuration
- [x] Task: Add T2T genome to igenomes config
    - [x] Reference FASTA path (placeholder)
    - [x] BWA/BWA-MEM2 index paths (placeholder)
    - [x] Known sites (dbSNP, Mills, etc.) (placeholder)
- [x] Task: Configure snpEff for CHM13v2
    - [x] snpeff_db = 'CHM13v2.0'
- [x] Task: Configure VEP for T2T
    - [x] vep_genome = 'T2T-CHM13v2.0'

## Phase 3: Pangenome Configuration
- [x] Task: Create pangenome.config
    - [x] v1.1 index paths (GBZ, minimizer, distance, etc.) - placeholder
    - [x] v2 index paths - placeholder
    - [x] v1.1-grch38 for GRCh38 compatibility
    - [x] Version selection parameter (pangenome_version)
    - [x] PanGenie panel resources
    - [x] Phasing resources (T2T and GRCh38)
- [x] Task: Include pangenome.config in nextflow.config
- [x] Task: Add giraffe to aligner options

## Phase 4: Parabricks Review
- [ ] Task: Compare implementations
    - [ ] Review legacy Parabricks module
    - [ ] Review current Sarek Parabricks
    - [ ] Document differences
- [ ] Task: Decision on approach
    - [ ] Keep current / adopt legacy / hybrid

## Phase 5: Validation
- [ ] Task: Test T2T configuration
- [ ] Task: Test pangenome v1.1 configuration
- [ ] Task: Test pangenome v2 configuration

## Notes
- All paths use TODO placeholders - need to fill in actual S3 URLs
- Pangenome paths follow igenomes pattern: ${params.igenomes_base}/Homo_sapiens/Pangenome/...
- T2T paths follow igenomes pattern: ${params.igenomes_base}/Homo_sapiens/T2T/CHM13v2/...
