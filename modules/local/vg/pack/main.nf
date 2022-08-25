process SV_VG_PACK {
    tag "$meta.id"
    label 'process_med'

    conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.41.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(gam)
    path (ref_pangenome_index) // gbz or xg

    output:
    tuple val(meta), path("*.pack"), emit: pack
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg pack \\
        -x $ref_pangenome_index \\
        -g $gam \\
        -Q 5 \\
        -t $task.cpus \\
        -o ${prefix}.pack

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.pack
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
