process VG_CONVERT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(input_graph)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa, optional: true
    tuple val(meta), path("*.gbz"), emit: gbz, optional: true
    tuple val(meta), path("*.vg"),  emit: vg,  optional: true
    tuple val(meta), path("*.pg"),  emit: pg,  optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-f'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_ext = args.contains('-f') || args.contains('--gfa-out') ? 'gfa' :
                     args.contains('-G') || args.contains('--gbz-out') ? 'gbz' :
                     args.contains('-v') || args.contains('--vg-out') ? 'vg' :
                     args.contains('-p') || args.contains('--packed-out') ? 'pg' : 'gfa'
    """
    vg convert \\
        ${args} \\
        -t ${task.cpus} \\
        ${input_graph} \\
        > ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '-f'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_ext = args.contains('-f') || args.contains('--gfa-out') ? 'gfa' :
                     args.contains('-G') || args.contains('--gbz-out') ? 'gbz' :
                     args.contains('-v') || args.contains('--vg-out') ? 'vg' :
                     args.contains('-p') || args.contains('--packed-out') ? 'pg' : 'gfa'
    """
    touch ${prefix}.${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
