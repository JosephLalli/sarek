process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::expansionhunter=4.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunter:4.0.2--he785bd8_0' :
        'jlalli/expansionhunter:5.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path (fasta)
    path (variant_catalog)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.json"), emit: json_report
    tuple val(meta), path("*.bed"), emit: str_realigned_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = (meta.sex == 'male' || meta.sex == 1 || meta.sex == 'XY') ? "male" : "female"
    """
    /opt/expansionhunter/bin/ExpansionHunter \\
        $args \\
        --reads $bam \\
        --output-prefix $prefix \\
        --reference $fasta \\
        --variant-catalog $variant_catalog \\
        --sex $sex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(/opt/expansionhunter/bin/ExpansionHunter --version 2>&1) | sed 's/^.*ExpansionHunter v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.json
    touch ${prefix}.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(/opt/expansionhunter/bin/ExpansionHunter --version 2>&1) | sed 's/^.*ExpansionHunter v//')
    END_VERSIONS
    """
}
