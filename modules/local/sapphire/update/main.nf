process SAPPHIRE_UPDATE {
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/jlalli/sapphire:53a96e5ac"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bin)

    output:
    tuple val(meta), path("*.phased_corrected.bcf"), emit: bcf
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /usr/src/sapphire/bin/pp_update \\
        -f ${vcf} \\
        -b ${bin} \\
        -o ${prefix}.phased_corrected.bcf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sapphire: \$(/usr/src/sapphire/bin/pp_extract --version 2>&1 | head -n 1 | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.phased_corrected.bcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sapphire: 1.0.0
    END_VERSIONS
    """
}
