# Track Specification: Phasing Tools

## Objective
Integrate haplotype phasing tools (SHAPEIT5, SHAPEIT4, WhatsHap, HapCUT) into sarek_pangenome pipeline.

## Scope

### SHAPEIT5 (4 modules)
| Tool | Function |
|------|----------|
| shapeit5/phase_common | Phase common variants (MAF > 0.1%) |
| shapeit5/phase_rare | Phase rare variants using scaffold |
| shapeit5/ligate | Ligate phased chunks |
| shapeit5/switch | Calculate switch error rates |

### SHAPEIT4 (1 module)
| Tool | Function |
|------|----------|
| shapeit4/phase_common | Legacy common variant phasing |

### WhatsHap (2 modules)
| Tool | Function |
|------|----------|
| whatshap/phase | Read-backed phasing |
| whatshap/stats | Phasing statistics and N50 |

### HapCUT (NEW - 1 module)
| Tool | Function |
|------|----------|
| hapcut2 | Read-backed haplotype assembly |

### Subworkflows (1 workflow)
| Workflow | Function |
|----------|----------|
| phase_and_impute | Orchestrate full phasing pipeline |

## Technical Considerations
- Check nf-core/modules for existing implementations
- SHAPEIT5 requires reference panel for population phasing
- WhatsHap/HapCUT use read information (BAM input)
- Support phasing with and without reference panel
- Chunk-based parallelization for large chromosomes

## Dependencies
- Requires: Variant calls (VCF) from upstream callers
- Optional: Reference panel for SHAPEIT5
- Optional: BAM files for read-backed phasing

## Deliverables
- 8 phasing modules in `modules/local/` or `modules/nf-core/`
- 1 phase_and_impute subworkflow
- Integration with variant calling output
- Support for trio-aware phasing (optional)
