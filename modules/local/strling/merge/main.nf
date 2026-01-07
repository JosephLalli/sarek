process STRLING_MERGE {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strling:0.6.0--h7b50bb2_0' :
        'quay.io/biocontainers/strling:0.6.0--h7b50bb2_0' }"

    input:
    path(bins)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    path("*-bounds.txt"), emit: bounds
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "strling_joint"
    """
    strling merge \
        $args \
        -f $fasta \
        -o $prefix \
        $bins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "strling_joint"
    """
    touch ${prefix}-bounds.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: 
    END_VERSIONS
    """
}
