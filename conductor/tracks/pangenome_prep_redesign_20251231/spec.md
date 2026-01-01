# Track Specification: Pangenome Preparation Redesign

## Objective
Redesign the "home-grown" pangenome preparation modules using nf-core conventions and best practices.

## Background
The legacy sarek_first_attempt contained 9 custom modules for preparing pangenome inputs. These were functional but did not follow nf-core module standards. This track will redesign these as proper nf-core-style modules.

## Scope

### Modules to Redesign

| Legacy Module | Function | Redesign Approach |
|---------------|----------|-------------------|
| make_personalized_gbz.nf | Build personalized GBZ graph | Use vg tools, standardize I/O |
| make_personalized_gfa.nf | Build personalized GFA graph | Use vg tools, standardize I/O |
| add_paths_to_gfa.nf | Add sample paths to GFA | Wrap vg paths command |
| split_by_contig.nf | Split files by chromosome | Generic utility module |
| split_sample_from_vcf.nf | Extract sample from multi-sample VCF | Use bcftools view |
| split_PARs_off_X.nf | Handle PAR regions on X chromosome | Specialized interval handling |
| split_ref_fasta.nf | Split reference by contig | Use samtools faidx |
| split_vcf_by_contig.nf | Split VCF by chromosome | Use bcftools or tabix |
| unify_insertion_fastas.nf | Merge insertion sequences | Custom logic needed |

## Design Principles
1. **Single responsibility**: One tool per module
2. **Standard I/O**: Follow nf-core channel patterns
3. **Parameterized**: No hardcoded paths or values
4. **Documented**: meta.yml for each module
5. **Containerized**: environment.yml with conda deps
6. **Testable**: nf-test cases for each module

## Technical Considerations
- Some functionality may already exist in nf-core/modules
- Consider using existing bcftools/samtools modules where possible
- Graph operations should wrap vg commands consistently
- PAR handling is reference-specific (GRCh38 vs T2T)

## Dependencies
- Requires: Core pangenome tools (vg modules)
- Requires: Understanding of pangenome index building workflow

## Deliverables
- 9 redesigned modules following nf-core conventions
- Pangenome index building subworkflow
- Documentation for custom graph building
- Test cases with small reference subset
