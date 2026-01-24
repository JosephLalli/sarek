# Track Specification: Phasing Tools

## Objective
Integrate comprehensive phasing tools (SHAPEIT5, SHAPEIT4, WhatsHap, HapCUT2, SAPPHIRE) into sarek_pangenome pipeline.

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
| shapeit4/phase_common | Ligate/Scaffold read-backed blocks (Hybrid flow) |

### Read-Backed Tools (3 modules)
| Tool | Function |
|------|----------|
| whatshap/phase | Read-backed phasing |
| whatshap/stats | Phasing statistics and N50 |
| hapcut2 | Read-backed haplotype assembly |

### Correction Tools (1 module)
| Tool | Function |
|------|----------|
| sapphire | Read-based phase correction for population phasing |

### Subworkflows (1 workflow)
| Workflow | Function |
|----------|----------|
| phase_and_impute | Orchestrate population, read-backed, corrected, and hybrid phasing |

## Technical Considerations
- **Strategy:** Support 4 distinct workflows via `--phasing_method`:
    1. `population`: SHAPEIT5 (common -> ligate -> rare)
    2. `read_backed`: WhatsHap or HapCUT2
    3. `read_corrected`: SHAPEIT5 -> SAPPHIRE
    4. `hybrid`: WhatsHap/HapCUT2 -> SHAPEIT4 -> SHAPEIT5_rare
- Check nf-core/modules for existing implementations
- SAPPHIRE requires new local module
- Chunk-based parallelization for large chromosomes

## Dependencies
- Requires: Variant calls (VCF) from upstream callers
- Optional: Reference panel (SHAPEIT workflows)
- Optional: BAM files (Read-backed/Corrected workflows)

## Deliverables
- 9 phasing modules
- 1 subworkflow with 4 logic paths
- Integration with variant calling output
- Test cases for all supported tools