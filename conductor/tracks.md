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

### [x] Track: Legacy Feature Gap Analysis
*Link: [./tracks/legacy_gap_analysis_20251231/](./tracks/legacy_gap_analysis_20251231/)*

Inventoried ~90 files from sarek_first_attempt. Created feature_backlog.md with Keep/Ignore/Redesign decisions.

### [x] Track: Test Adequacy Evaluation (Pangenome & DeepVariant)
*Link: [./conductor/tracks/test_adequacy_pangenome_20260109/](./conductor/tracks/test_adequacy_pangenome_20260109/)*

Successfully implemented and verified pangenome-aware alignment (VG Giraffe) and variant calling (DeepVariant) using rebased CHM13 test data.

### [x] Track: STR Consensus Merging (EnsembleTR-style)
*Link: [./conductor/tracks/str_consensus_merge_20260118/](./conductor/tracks/str_consensus_merge_20260118/)*

Implemented production-grade STR consensus merging with format adapters for STRling, coordinate-aware normalization, and EnsembleTR integration.

---

## Active & Planned Tracks

## [~] Track: STR Enhancement
*Link: [./conductor/tracks/str_enhancement_20260109/](./conductor/tracks/str_enhancement_20260109/)*

---

## [ ] Track: Phasing Pipeline Separation
*Link: [./conductor/tracks/phasing_separation_20260109/](./conductor/tracks/phasing_separation_20260109/)*

---

### [ ] Track: Long-read ONT Support
*Link: [./conductor/tracks/ont_support_20260104/](./conductor/tracks/ont_support_20260104/)*

Comprehensive support for ONT data: preprocessing, alignment (Minimap2/Giraffe), and variant calling (Pepper-DeepVariant/Clair3).

---

## [ ] Track: Long-read PacBio Support
*Link: [./conductor/tracks/pacbio_support_20260104/](./conductor/tracks/pacbio_support_20260104/)*
*Priority: LOW | Depends on: strict_syntax_upgrade*

Long-read PacBio support including preprocessing (lima), alignment (Minimap2/Giraffe), and variant calling (DeepVariant/PBSV/TRGT).

---

## Summary

| Track | Status | Modules | Subworkflows |
|-------|--------|---------|--------------|
| legacy_gap_analysis | DONE | - | - |
| test_adequacy_pangenome | DONE | - | - |
| str_enhancement | PLANNED | - | - |
| str_consensus_merge | DONE | 3 | 1 |
| phasing_separation | PLANNED | - | - |
| ont_support | PLANNED | - | - |
| pacbio_support | PLANNED | - | - |
| **TOTAL** | | **3** | **1** |
