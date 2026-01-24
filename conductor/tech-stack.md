# Technology Stack

## Core Workflow Engine
*   **Platform:** Nextflow (DSL2)
*   **Language:** Groovy (Nextflow DSL)
*   **Note:** strict-syntax upgrade completed. However, running with `NXF_SYNTAX_PARSER=v2` is currently **broken for test runs**. Recommend running tests without this environment variable.

## Environment & Containerization
*   **Container Runtime:** Docker, Singularity (Apptainer)
*   **Dependency Management:** nf-core modules (Conda/Bioconda supported but not required; Docker/Singularity preferred)

## Languages & Scripts
*   **Scripting:** Bash, Python
*   **Configuration:** Config files (`.config`)

## Infrastructure
*   **Execution:** Local, HPC (Slurm, etc.), Cloud (AWS for CI)
*   **Testing:** nf-test, GitHub Actions

## Bioinformatics Tools (Key Dependencies)
*   **Alignment:** BWA-MEM, BWA-MEM2, Dragmap, Sentieon, vg giraffe
*   **Variant Calling:** GATK (HaplotypeCaller, Mutect2), DeepVariant, Strelka, Manta, FreeBayes, PanGenie, ExpansionHunter, STRling, GangSTR, etc.
*   **STR Consensus:** EnsembleTR (container: `community.wave.seqera.io/library/pysam_samtools_pip_ensembletr:343f01fd8b77ac0a`)
*   **Phasing:** SHAPEIT5, WhatsHap, SHAPEIT4, SAPPHIRE
*   **QC:** FastQC, MultiQC, Samtools
