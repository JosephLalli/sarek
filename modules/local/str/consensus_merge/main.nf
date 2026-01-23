process STR_CONSENSUS_MERGE {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/biocontainers/python:3.10.4"

    input:
    tuple val(meta), \
        path(eh_vcf, stageAs: "eh/*"), \
        path(gs_vcf, stageAs: "gs/*"), \
        path(sl_vcf, stageAs: "sl/*")

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    str_consensus_merge.py \
        $args \
        --eh $eh_vcf \
        --gs $gs_vcf \
        --sl $sl_vcf \
        --output ${prefix}.consensus.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        str_consensus_merge: v0.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        str_consensus_merge: v0.1
    END_VERSIONS
    """
}
