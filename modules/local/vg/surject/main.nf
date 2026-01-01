process VG_SURJECT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(gam)
    tuple val(meta2), path(graph)
    path ref_paths

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_format = task.ext.output_format ?: 'bam'
    def output_flag = output_format == 'bam' ? '-b' : ''
    def paths_arg = ref_paths ? "--into-paths ${ref_paths}" : ''
    def graph_arg = graph.name.endsWith('.gbz') ? "-Z ${graph}" :
                    graph.name.endsWith('.xg') ? "-x ${graph}" :
                    "-g ${graph}"
    """
    vg surject \\
        ${args} \\
        ${graph_arg} \\
        ${paths_arg} \\
        -t ${task.cpus} \\
        --sample ${meta.sample ?: meta.id} \\
        ${output_flag} \\
        ${gam} \\
        > ${prefix}.${output_format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_format = task.ext.output_format ?: 'bam'
    """
    touch ${prefix}.${output_format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
