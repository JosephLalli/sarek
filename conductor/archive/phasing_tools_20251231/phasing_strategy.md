# Phasing Strategy for Sarek Pangenome

This document defines the supported phasing strategies for the pipeline.

## 1. Population-based Phasing (SHAPEIT5)
Best for large cohorts or when a high-quality reference panel (e.g., 1000G, HPRC) is available.
- **Tools:** SHAPEIT5_common (or SHAPEIT4 for read-backed phasing support), and SHAPEIT5_rare.
- **Inputs:** VCF, Reference Panel (VCF/BCF), Genetic Map.
- **Pros:** High accuracy for common variants; can impute missing genotypes.
- **Cons:** Dependent on reference panel quality

## 2. Read-backed Phasing (WhatsHap / HapCUT2)
Best for high-depth short reads or long-read data (ONT/PacBio).
- **Tools:** WhatsHap, HapCUT2.
- **Inputs:** VCF, BAM/CRAM.
- **Pros:** Phases variants based on direct read evidence; excellent for rare variants. Can be used as input to SHAPEIT4 for read-backed statistical phasing.
- **Cons:** Phase blocks are limited by read length; cannot phase across gaps.

## 3. Read-backed phase correction (SAPPHIRE)
Best for high-depth short reads or long-read data (ONT/PacBio).
- **Tools:** SAPPHIRE (https://github.com/rwk-unil/sapphire)
- **Inputs:** Statistically phased VCF, BAM/CRAM.
- **Pros:** Adjusts phasing of rare variants based on read support. 
- **Cons:** Can be finicky, not commonly used

## 3. Hybrid Phasing (Recommended)
Combines read-backed phasing with population-based ligation.
- **Tools:** WhatsHap or HapCUT2 followed by SHAPEIT4_common and SHAPEIT5_rare.
- **Pros:** Maximizes phase block size and accuracy.
- **Cons:** Very computationally expensive
- **Workflow:** 
    1. Phase locally using `WhatsHap` and read evidence.
    2. Use `WhatsHap` output as a scaffold for `SHAPEIT4_common` to ligate blocks and phase remaining variants using population data.
    3. Use SHAPEIT5_rare to phase rare variants.

## 4. Population phasing with read based correction
Combines population phasing with read-backed phase correction.
- **Tools:** SHAPEIT5_common -> SHAPEIT5_rare -> SAPPHIRE.
- **Pros:** Only performs computationally difficult read-backed phasing in regions where population phasing fails. More tractable.
- **Cons:** New software with unclear long term support. Not commonly used.
- **Workflow:** 
    1. Phase using SHAPEIT5 toolset.
    2. Use SAPPHIRE to examine reads and adjust phaasing around variants with little phasing support from population panels.



## Implementation Plan
The pipeline will introduce the following parameters:
- `--phasing`: Enable phasing.
- `--phasing_method`: Select `population`, `read_backed`, `read_corrected`, or `hybrid` (default: `population`).
- `--phasing_panel`: Path to reference panel.
- `--phasing_map`: Path to genetic map.
