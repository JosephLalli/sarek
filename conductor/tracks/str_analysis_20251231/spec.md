# Track Specification: STR Analysis Tools

## Objective
Integrate short tandem repeat (STR) analysis tools into sarek_pangenome pipeline for detection and genotyping of repeat expansions.

## Scope

### ExpansionHunter (legacy - 1 module)
| Tool | Function |
|------|----------|
| expansionhunter | Targeted STR genotyping with variant catalog |

### STRling (NEW - 1 module)
| Tool | Function |
|------|----------|
| strling | Genome-wide STR detection and outlier calling |

Reference: https://github.com/quinlan-lab/STRling-nf

### GangSTR (NEW - 1 module)
| Tool | Function |
|------|----------|
| gangstr | Genome-wide STR profiling |

### Additional Tools (5 legacy modules)
| Tool | Function |
|------|----------|
| jellyfish/count | K-mer counting for STR analysis |
| kmc/kmc | K-mer counting (alternative) |
| kmc/kmc_dump | Export k-mer database |
| illumina/hap.py | Variant benchmarking |
| deepvariant/convert_haploid_regions | Haploid region handling |

## Technical Considerations
- STR tools require repeat catalog files (BED format)
- Different tools have different catalog formats
- Consider T2T-specific STR catalogs
- Sex chromosome handling (chrX/chrY)
- Integration with phasing for haplotype-resolved STRs

## Use Cases
1. **Clinical**: Detect pathogenic repeat expansions (Huntington's, fragile X, etc.)
2. **Population**: Genome-wide STR profiling
3. **QC**: K-mer-based quality metrics

## Dependencies
- Requires: Aligned BAM/CRAM files
- Requires: STR catalog files per reference genome
- Optional: Phased variants for haplotype-resolved calls

## Deliverables
- 3 STR caller modules (ExpansionHunter, STRling, GangSTR)
- 3 k-mer modules (jellyfish, kmc)
- 2 utility modules (hap.py, convert_haploid_regions)
- STR analysis subworkflow
- Catalog files for GRCh38 and T2T
