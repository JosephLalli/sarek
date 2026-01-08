process STRLING_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strling:0.6.0--h7b50bb2_0' :
        'quay.io/biocontainers/strling:0.6.0--h7b50bb2_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(str_index)

    output:
    tuple val(meta), path("*.bin"), emit: bin
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_arg = str_index ? "-g ${str_index}" : ""
    """
    strling extract \
        $args \
        $index_arg \
        -f $fasta \
        $bam \
        ${prefix}.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: \$(strling 2>&1 | grep 'version:' | sed 's/strling version: // ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: 0.6.0
    END_VERSIONS
    """
}