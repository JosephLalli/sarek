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
    tuple val(meta), path("*.ri"),   emit: ri,   optional: true
    path "versions.yml",             emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Determine input flag if not provided in args
    def input_arg = ""
    if (!(args.contains('-Z') || args.contains('--gbz-input') || args.contains('-g') || args.contains('--graph-input'))) {
        if (input_file.name.endsWith('.gbz')) {
            input_arg = "-Z ${input_file}"
        } else {
            input_arg = "${input_file}"
        }
    } else {
        input_arg = "${input_file}"
    }
    
    // Determine main output argument if not -r (r-index) or -Z (GBZ) or -g (graph)
    // If just -r is used, we don't need -o
    def output_arg = ""
    if (!(args.contains('-r') || args.contains('--r-index') || args.contains('-Z') || args.contains('--gbz-format'))) {
       output_arg = args.contains('-o') || args.contains('--output') ? '' : "-o ${prefix}.gbwt"
    }

    """
    vg gbwt \\
        ${args} \\
        ${output_arg} \\
        ${input_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gbwt
    touch ${prefix}.ri

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
