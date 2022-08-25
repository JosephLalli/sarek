process VG_SURJECT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.41.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(gam)
    path (ref_pangenome_index) // xg or gbz
    path (ref_chrom_names)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg surject \\
          $args
          -F $ref_chrom_names \\
          -x $ref_pangenome_index \\
          -t $task.cpus \\
          --bam-output \\
          --sample ${meta.id} \\
          --prune-low-cplx \\
          --interleaved \\
          $gam > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ""
    def max_frag_len = 3000
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
