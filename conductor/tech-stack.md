# Technology Stack

## Core Workflow Engine
*   **Platform:** Nextflow (DSL2)
*   **Language:** Groovy (Nextflow DSL)

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
*   **Alignment:** BWA-MEM, BWA-MEM2, Dragmap, Sentieon
*   **Variant Calling:** GATK (HaplotypeCaller, Mutect2), DeepVariant, Strelka, Manta, FreeBayes, etc.
*   **QC:** FastQC, MultiQC, Samtools
*   **Upcoming:** Pangenie, Graph aligners (minigraph-cactus/vg)