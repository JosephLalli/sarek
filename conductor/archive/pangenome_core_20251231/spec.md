# Track Specification: Core Pangenome Tools

## Objective
Integrate core pangenome alignment and genotyping tools (vg, PanGenie) into sarek_pangenome pipeline.

## Scope

### vg Tools (11 modules)
| Tool | Function |
|------|----------|
| vg/giraffe | Pangenome-aware short read alignment |
| vg/surject | Convert GAM alignments to BAM |
| vg/haplotypes | Extract haplotypes from graph |
| vg/stats | Graph and alignment statistics |
| vg/deconstruct | Generate VCF from graph paths |
| vg/gbwt | GBWT haplotype index operations |
| vg/minimizer | Minimizer index for giraffe |
| vg/index | General graph indexing |
| vg/paths | Path manipulation |
| vg/convert | Format conversion (GFA/VG/GBZ) |

### PanGenie (1 module)
| Tool | Function |
|------|----------|
| pangenie | Genotype variants from pangenome graph |

### Subworkflows (2 workflows)
| Workflow | Function |
|----------|----------|
| giraffe_mapping | vg giraffe alignment pipeline |
| pangenie_map_and_call | PanGenie genotyping workflow |

## Technical Considerations
- vg tool APIs have changed significantly - modules need updating
- Support both pangenome v1.1 and v2 index formats
- Integrate with existing Sarek channel structures
- Maintain strict DSL2 syntax compliance

## Dependencies
- Requires: config_t2t_pangenome_20251231 (pangenome paths configured)
- Blocked by: None after config track complete

## Deliverables
- 11 updated vg modules in `modules/local/vg/`
- 1 pangenie module in `modules/local/pangenie/`
- 2 subworkflows in `subworkflows/local/`
- Integration with main workflow
- Test cases for each module
