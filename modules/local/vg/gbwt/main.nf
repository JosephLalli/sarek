process VG_GBWT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.gbwt"), emit: gbwt, optional: true
    tuple val(meta), path("*.gbz"),  emit: gbz,  optional: true
    tuple val(meta), path("*.gg"),   emit: gg,   optional: true
    tuple val(meta), path("*.txt"),  emit: translation, optional: true
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_arg = args.contains('-o') || args.contains('--output') ? '' : "-o ${prefix}.gbwt"
    """
    vg gbwt \\
        ${args} \\
        ${output_arg} \\
        -t ${task.cpus} \\
        ${input_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gbwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
