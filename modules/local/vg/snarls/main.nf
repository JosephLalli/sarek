process AUTOINDEX_GFA {
    tag "$gbz"
    label 'process_med'

    conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.41.0--h9ee0642_0' }"

    input:
    path (ref_pangenome_index) // gbz or xg

    output:
    path ('*.pb'), emit: snarls
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vg snarls \\
        -a \\
        -t $task.cpus
        $ref_pangenome_index > reference_pangenome.pb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch reference_pangenome.pb
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
