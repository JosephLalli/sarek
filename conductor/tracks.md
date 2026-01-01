# Sarek Pangenome Development Tracks

This file tracks all major development tracks. Each track has its own detailed plan in its respective folder.

---

## Track Dependency Graph

```
strict_syntax_upgrade (COMPLETED)
         │
         ▼
legacy_gap_analysis (COMPLETED)
         │
         ▼
config_t2t_pangenome ◄─────────────────────────────┐
         │                                          │
         ├──────────────┬──────────────┐           │
         ▼              ▼              ▼           │
 pangenome_core    phasing_tools   str_analysis    │
         │              │              │           │
         ▼              │              │           │
pangenome_prep_redesign │              │           │
         │              │              │           │
         └──────────────┴──────────────┴───────────┘
                        │
                        ▼
                   INTEGRATION
```

---

## Completed Tracks

### [x] Track: Strict Syntax Upgrade
*Link: [./tracks/strict_syntax_upgrade_20251231/](./tracks/strict_syntax_upgrade_20251231/)*

DSL2 strict syntax compliance (NXF_SYNTAX_PARSER=v2). Replaced deprecated Channel. with channel., added def keywords, fixed variable shadowing, converted if/else to ternary.

### [x] Track: Legacy Feature Gap Analysis
*Link: [./tracks/legacy_gap_analysis_20251231/](./tracks/legacy_gap_analysis_20251231/)*

Inventoried ~90 files from sarek_first_attempt. Created feature_backlog.md with Keep/Ignore/Redesign decisions.

---

## Active Tracks

### [ ] Track: T2T and Pangenome Configuration
*Link: [./tracks/config_t2t_pangenome_20251231/](./tracks/config_t2t_pangenome_20251231/)*
*Priority: HIGH*

- T2T (CHM13) reference configuration
- Pangenome v1.1 + v2 index paths
- snpEff CHM13v2 setup
- Parabricks implementation review

---

## Planned Tracks

### [ ] Track: Core Pangenome Tools
*Link: [./tracks/pangenome_core_20251231/](./tracks/pangenome_core_20251231/)*
*Priority: HIGH | Depends on: config_t2t_pangenome*

11 vg modules + 1 PanGenie module + 2 subworkflows:
- vg/giraffe, surject, haplotypes, stats, deconstruct, gbwt, minimizer, index, paths, convert
- pangenie/pangenie
- giraffe_mapping, pangenie_map_and_call subworkflows

### [ ] Track: Phasing Tools
*Link: [./tracks/phasing_tools_20251231/](./tracks/phasing_tools_20251231/)*
*Priority: HIGH | Depends on: config_t2t_pangenome*

8 phasing modules + 1 subworkflow:
- shapeit5/phase_common, phase_rare, ligate, switch
- shapeit4/phase_common
- whatshap/phase, stats
- hapcut2 (NEW)
- phase_and_impute subworkflow

### [ ] Track: STR Analysis
*Link: [./tracks/str_analysis_20251231/](./tracks/str_analysis_20251231/)*
*Priority: MEDIUM | Depends on: config_t2t_pangenome*

9 STR/utility modules:
- expansionhunter, strling (NEW), gangstr (NEW)
- jellyfish/count, kmc/kmc, kmc/kmc_dump
- illumina/hap.py, deepvariant/convert_haploid_regions

### [ ] Track: Pangenome Prep Redesign
*Link: [./tracks/pangenome_prep_redesign_20251231/](./tracks/pangenome_prep_redesign_20251231/)*
*Priority: MEDIUM | Depends on: pangenome_core*

Redesign 9 "home-grown" modules to nf-core standards:
- VCF/FASTA splitting utilities
- PAR region handling
- Personalized graph building
- Insertion sequence merging

---

## Summary

| Track | Status | Modules | Subworkflows |
|-------|--------|---------|--------------|
| strict_syntax_upgrade | DONE | - | - |
| legacy_gap_analysis | DONE | - | - |
| config_t2t_pangenome | ACTIVE | - | - |
| pangenome_core | PLANNED | 12 | 2 |
| phasing_tools | PLANNED | 8 | 1 |
| str_analysis | PLANNED | 9 | 1 |
| pangenome_prep_redesign | PLANNED | 7 | 1 |
| **TOTAL** | | **36** | **5** |
