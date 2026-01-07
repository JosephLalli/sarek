process STRLING_CALL {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strling:0.6.0--h7b50bb2_0' :
        'quay.io/biocontainers/strling:0.6.0--h7b50bb2_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bin)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(bounds)

    output:
    tuple val(meta), path("*-bounds.txt")  , emit: bounds
    tuple val(meta), path("*-genotype.txt"), emit: genotype
    tuple val(meta), path("*.vcf.gz")      , emit: vcf, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bounds_arg = bounds ? "-b $bounds" : ""
    """
    strling call \\
        $args \\
        $bounds_arg \\
        -f $fasta \\
        -o $prefix \\
        $bam \\
        $bin

    if [ -f ${prefix}.vcf ]; then
        gzip ${prefix}.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: \$(strling --version | sed "s/^.*strling //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-bounds.txt
    touch ${prefix}-genotype.txt
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling: \$(strling --version | sed "s/^.*strling //")
    END_VERSIONS
    """
}