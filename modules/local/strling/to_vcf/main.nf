process STRLING_TO_VCF {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/biocontainers/python:3.10.4"

    input:
    tuple val(meta), path(genotype)
    path(catalog)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    strling_to_vcf.py \
        $args \
        --strling $genotype \
        --catalog $catalog \
        --output ${prefix}.vcf \
        --sample ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling_to_vcf: v0.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strling_to_vcf: v0.1
    END_VERSIONS
    """
}
