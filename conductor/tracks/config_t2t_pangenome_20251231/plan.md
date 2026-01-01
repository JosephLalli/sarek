# Track Plan: T2T and Pangenome Configuration

## Phase 1: Audit Existing Configuration
- [ ] Task: Review current `conf/` structure
    - [ ] Identify where genome-specific configs live
    - [ ] Review igenomes.config pattern
- [ ] Task: Locate legacy T2T configuration
    - [ ] Find T2T settings in sarek_first_attempt
    - [ ] Document what was configured
- [ ] Task: Locate legacy pangenome configuration
    - [ ] Find pangenome settings in sarek_first_attempt
    - [ ] Document v1.1 requirements

## Phase 2: T2T Reference Configuration
- [ ] Task: Add T2T genome to igenomes config
    - [ ] Reference FASTA path
    - [ ] BWA/BWA-MEM2 index paths
    - [ ] Known sites (dbSNP, Mills, etc.)
- [ ] Task: Configure snpEff for CHM13v2
    - [ ] Database setup
    - [ ] Config file modifications
- [ ] Task: Configure VEP for T2T
    - [ ] Cache paths
    - [ ] Species/assembly settings

## Phase 3: Pangenome Configuration
- [ ] Task: Create pangenome.config
    - [ ] v1.1 index paths (GBZ, minimizer, distance, etc.)
    - [ ] v2 index paths
    - [ ] Version selection parameter
- [ ] Task: Add pangenome parameters to schema
    - [ ] Document new parameters
    - [ ] Add validation

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
