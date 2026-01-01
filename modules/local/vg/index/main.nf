process VG_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.dist"), emit: dist, optional: true
    tuple val(meta), path("*.xg"),   emit: xg,   optional: true
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dist_cmd = args.contains('-j') || args.contains('--dist-name') ? '' : "-j ${prefix}.dist"
    def xg_cmd = args.contains('-x') || args.contains('--xg-name') ? '' : "-x ${prefix}.xg"
    """
    vg index \\
        ${args} \\
        ${dist_cmd} \\
        ${xg_cmd} \\
        -t ${task.cpus} \\
        ${graph}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dist
    touch ${prefix}.xg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
