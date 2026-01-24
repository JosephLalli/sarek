# Product Guide: Sarek Enhanced (Strict Syntax & Advanced Variant Calling)

## # Initial Concept
A modernized iteration of the `nf-core/sarek` pipeline that first upgrades the codebase to strict Nextflow DSL2 syntax to eliminate deprecations, and then builds upon this robust foundation to introduce advanced variant calling features (Graph alignment, T2T reference, STR calling, Pangenie SV imputation) for BrainVar analysis.

## Project Vision
To engineer a future-proof, high-precision variant calling pipeline. The immediate priority is to refactor the existing Sarek codebase to adhere to strict Nextflow syntax, eliminating deprecated patterns to ensure long-term maintainability. On this modernized architecture, we will integrate cutting-edge genomic tools to benchmark improvements against standard Trio datasets and deploy for production-scale BrainVar data analysis.

## Target Audience
1.  **Bioinformatics Developers:** Contributors needing a clean, strict-syntax codebase for easier maintenance and extension.
2.  **Research & Clinical Analysts:** Users requiring advanced features (Graph alignment, STRs) with clear documentation on usage and interpretation.
3.  **nf-core Community:** Beneficiaries of a modernized, "strict-mode" compliant Sarek implementation.

## Core Goals
1.  **Modernize Architecture (Strict Syntax Upgrade):** Execute the "Zazzy Rolling Stearns" plan to ensure full compliance with modern Nextflow style guides.
    *   *Complete:* Replace deprecated capitalized `Channel.` with `channel.`.
    *   *Complete:* Refactor implicit `it` closures to explicit parameters.
    *   *Complete:* Standardize local modules with `task.ext.args`, `prefix`, and `meta.yml` documentation.
    *   *Complete:* Address missing resource labels and environments in nf-core modules.
2.  **Regression Validation:** Rigorous verification using existing Sarek test datasets (`test_full`, `test_full_germline`) to ensure the syntax upgrade preserves exact functionality.
3.  **Advanced Feature Integration:** Implement Graph alignment, T2T reference support, STR calling, Pangenie SV imputation, and comprehensive Haplotype Phasing.
4.  **Benchmarking:** Evaluate the new features using standard Trio datasets (GIAB) to quantify precision/recall gains.
5.  **Production Deployment:** Validate end-to-end execution on BrainVar pilot data.

## Key Features
*   **Strict Syntax Codebase:** A refactored Nextflow implementation compliant with modern standards (no implicit closures, fully documented local modules).
*   **Graph & T2T Alignment:** Full configuration support for T2T-CHM13 v2.0 and HPRC Pangenome (v1.1 and v2.0) with automated index generation. Integrated `vg giraffe` for pangenome-aware short-read alignment.
*   **Comprehensive Phasing Suite:** Integrated phasing workflow supporting population-based (SHAPEIT5), read-backed (WhatsHap), and hybrid/corrected (SHAPEIT5 + SAPPHIRE) strategies for high-accuracy haplotype reconstruction.
*   **Complex Variant Calling:** Native support for Short Tandem Repeats (STRs) and Genotype Imputation via PanGenie.
*   **Phased Validation Suite:** A structured testing pipeline (Regression -> Benchmarking -> Pilot).
