# Specification: STR Consensus Merging (EnsembleTR-compatible)

## Overview
Implement a robust, production-grade Short Tandem Repeat (STR) consensus merging strategy within the Sarek pipeline. This track integrates outputs from ExpansionHunter (EH), GangSTR (GS), and STRling (SL) into a single, high-confidence "Consensus VCF" using an EnsembleTR-compatible approach. The goal is to normalize SL output so it can be consumed alongside default EH/GS inputs by EnsembleTR while preserving novel STRling calls.

## Compatibility Target
- **EnsembleTR:** Target the latest EnsembleTR release and pin its version when implemented.
- **Inputs Known to Work:** EnsembleTR accepts default GangSTR, HipSTR, and ExpansionHunter outputs.
- **Objective:** Produce SL and consensus VCFs that follow EnsembleTR expectations (HipSTR-compatible alleles, required INFO/FORMAT fields, and stable locus definitions).
- **Labeling Note:** EnsembleTR hardcodes method labels as `AdVNTR, EH, HipSTR, GangSTR`. STRling output is encoded in HipSTR-compatible format, then the EnsembleTR output is post-processed to relabel `HipSTR` to `STRling` for clarity without changing bit positions.

## Functional Requirements

### 1. STRling-to-VCF Adapter
- **Goal:** Convert STRling's native `-genotype.txt` into a VCF 4.2 compliant file that is EnsembleTR compatible.
- **Language:** Python.
- **Key Logic:**
    - Parse the `repeatunit` column for repeat motif and period.
    - Do not discard non-catalog calls; preserve novel STRling loci.
    - **Catalog calls:** Use the shared STR catalog to resolve absolute coordinates, locus IDs, and reference copy numbers.
    - **Novel calls:** Generate stable IDs (e.g., `STRling_Novel_<chrom>_<pos>`) and derive `END` and reference copy numbers from the STRling locus length and `repeatunit`.
    - **HipSTR-style alleles:** Emit sequence-based alleles and mark the file as HipSTR-compatible:
        - `##command=hipstr` header line for TRTools/EnsembleTR inference.
        - `INFO/START`, `INFO/END`, `INFO/PERIOD` (required for HipSTR parsing).
        - `INFO/RU` for repeat unit (used by downstream consensus merge).
        - `FORMAT/REPCN` (absolute repeat copy numbers, float) retained for consensus merge.
        - `FORMAT/Q` (derived quality score from STRling evidence).
    - **Header completeness:** Add VCF header lines for all INFO/FORMAT fields above.
    - **Contig lengths:** Emit `##contig=<ID=...,length=...>` lines derived from the reference FASTA index (`.fai`) so headers are VCF‑spec compliant.
    - **Compression:** Emit `*.vcf.gz` and `*.vcf.gz.tbi` by default (bgzip + tabix). The uncompressed VCF is not retained.
    - **Quality mapping:** Define a deterministic mapping for `Q` based on STRling read evidence (baseline: `Q = spanning + anchored`, capped at 99). Document the formula in code and in this spec if it changes.

### 2. Multi-Caller Consensus Merger
- **Goal:** Join EH, GS, and SL VCFs using coordinate-aware locus definitions and normalize across slight shifts.
- **Language:** Python (use `cyvcf2` or `pandas`).
- **Merging Logic:**
    - **Locus definition:** Use `POS` and `INFO/END` (fallback: `POS + len(REF) - 1`) to build intervals.
    - **Overlap rule:** Same contig and interval intersection after applying a configurable slop (`max_bp_shift`, default 1 bp).
    - **Normalization checks:** Require `RU` and `PERIOD` match; require `REPCN` difference within a configurable tolerance (`max_copy_delta`, default 4). If checks fail, keep calls separate.
    - **Priority Rule (Baseline):** EH > GS > SL, with motif-length toggle retained and configurable.
    - **Filter awareness:** Only consider `PASS` calls for consensus `GT`.
    - **Singleton handling:** Preserve SL novel calls that do not overlap catalog loci as standalone variants.

### 3. Rich Provenance Output
- **Consensus VCF Structure:**
    - `GT`: The "winning" genotype based on the priority logic.
    - `FORMAT/EH_GT`, `FORMAT/GS_GT`, `FORMAT/SL_GT`: Original genotypes for auditability.
    - `INFO/MERGE_SRC`: Tool providing the consensus `GT`.
    - `INFO/MERGE_RULE`: Rule applied (e.g., `Priority_EH`, `MotifLen_GS`).

### 4. EnsembleTR Validation Module
- **Goal:** Provide a standalone Nextflow module (`main.nf`) to run EnsembleTR from a public container image.
- **Containerization:** Use `community.wave.seqera.io/library/pysam_samtools_pip_ensembletr:343f01fd8b77ac0a` and document it in `conductor/tech-stack.md` before implementation.
- **Inputs:** EH/GS/SL VCFs per EnsembleTR defaults; HipSTR inputs can be passed via `task.ext.args` when needed.
- **Outputs:** EnsembleTR merged VCF and logs suitable for validation.
- **Post-processing:** Replace `hipstr=` with `strling=` in `INPUTS` and update the `METHODS` description string in the VCF header to reference STRling.

### 5. Nextflow Integration
- **Process 1:** `STRLING_TO_VCF` - Standalone modular process.
- **Process 2:** `STR_CONSENSUS_MERGE` - Accepts a tuple of VCFs (EH, GS, SL) and produces the final merged VCF.
- **Process 3 (Optional):** `ENSEMBLETR` - Runs the external tool for compatibility validation.

## Non-Functional Requirements
- **Standardization:** Adhere to nf-core module standards (meta.yml, `task.ext.args`, `task.ext.prefix`).
- **Testability:** Provide nf-test integration tests for the adapter and merger with content assertions (not only file presence). These content assertions count as coverage for this track.
- **Code Quality:** Python scripts must include docstrings and target >20% code coverage.
- **Compatibility:** All VCF outputs must be parseable by `bcftools`.

## Acceptance Criteria
- [ ] `STRLING_TO_VCF` outputs a bgzipped VCF parseable by `bcftools` with HipSTR-style REF/ALT sequences, `##command=hipstr`, and header definitions for `START`, `END`, `PERIOD`, `RU`, `REPCN`, `Q`.
- [ ] `STRLING_TO_VCF` emits a `.vcf.gz.tbi` index alongside the VCF.
- [ ] A novel STRling locus is preserved in output with a generated ID and consistent `END` and `REPCN`.
- [ ] `STR_CONSENSUS_MERGE` merges overlapping calls using coordinate overlap plus normalization, and preserves singleton SL novel calls.
- [ ] nf-test asserts output content (not only file presence) for both modules, including at least one ALT variant.
- [ ] EnsembleTR validation module can run on test fixtures and produces an output VCF (when fixtures are available).
- [ ] EnsembleTR output VCF labels STRling in `METHODS` and `INPUTS` after post-processing.

## Out of Scope
- Re-aligning reads or re-calling variants.
- Adding HipSTR to the pipeline (use if already available).
