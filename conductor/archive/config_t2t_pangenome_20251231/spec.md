# Track Specification: T2T and Pangenome Configuration

## Objective
Add configuration support for T2T (CHM13) reference genome and pangenome (v1.1 + v2) to sarek_pangenome pipeline.

## Scope
1. **T2T Reference Support**
   - Reference genome paths and resources
   - snpEff CHM13v2.110 genome configuration
   - T2T-specific interval files
   - Annotation resources (VEP, dbSNP, etc.)

2. **Pangenome Configuration**
   - Pangenome v1.1 index paths and parameters
   - Pangenome v2 index paths and parameters
   - Version detection/selection mechanism
   - Default parameters for pangenome alignment

3. **Parabricks Review**
   - Compare legacy implementation vs current Sarek
   - Document differences
   - Decide on approach

## Deliverables
- Updated `conf/igenomes.config` or equivalent with T2T paths
- New `conf/pangenome.config` for pangenome-specific settings
- Parameter documentation in `nextflow_schema.json`
- Test configurations for validation

## Dependencies
- Requires access to T2T reference files
- Requires pangenome v1.1 and v2 index locations
- Builds on strict-syntax compliant codebase
