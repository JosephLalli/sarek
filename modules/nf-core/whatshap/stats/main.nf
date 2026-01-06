process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.3--py310h30d963c_0' :
        'quay.io/biocontainers/whatshap:2.3--py310h30d963c_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.stats.tsv"), emit: stats_tsv
    tuple val(meta), path("*.blocks.gtf"), emit: blocks_gtf, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    whatshap stats \
        --tsv=${prefix}.stats.tsv \
        --block-list=${prefix}.blocks.gtf \
        ${args} \
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats.tsv
    touch ${prefix}.blocks.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 
    END_VERSIONS
    """
}
