process FREEBAYES_BAM_LEFT_ALIGN {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3' :
        'quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3' }"

    input:
    tuple val(meta), path(bam)
    path (ref_fasta)
    path (ref_fasta_fai)

    output:
    tuple val(meta), path("*.bam") , emit: bam
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bamleftalign \\
        < $bam \\
        > ${prefix}.bam \\
        --fasta-reference $ref_fasta \\
        --compressed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
