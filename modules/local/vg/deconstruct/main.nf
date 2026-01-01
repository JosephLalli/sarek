process VG_DECONSTRUCT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(graph)
    path ref_paths

    output:
    tuple val(meta), path("*.vcf"),    emit: vcf,    optional: true
    tuple val(meta), path("*.vcf.gz"), emit: vcf_gz, optional: true
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paths_arg = ref_paths ? "-P ${ref_paths}" : ''
    def graph_arg = graph.name.endsWith('.gbz') ? "-g ${graph}" :
                    graph.name.endsWith('.xg') ? "-x ${graph}" :
                    "-g ${graph}"
    """
    vg deconstruct \\
        ${args} \\
        ${graph_arg} \\
        ${paths_arg} \\
        -t ${task.cpus} \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
