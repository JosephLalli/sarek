process ENSEMBLETR {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/pysam_samtools_pip_ensembletr:343f01fd8b77ac0a"

    input:
    tuple val(meta), \
        path(eh_vcf, stageAs: "eh/*"), \
        path(gs_vcf, stageAs: "gs/*"), \
        path(sl_vcf, stageAs: "sl/*"), \
        path(ref_fasta), \
        path(ref_fai)

    output:
    tuple val(meta), path("*.ensembletr.vcf"), emit: vcf
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def base_args = args ?: "--vcfs ${eh_vcf},${gs_vcf},${sl_vcf} --ref ${ref_fasta}"
    """
    EnsembleTR \
        ${base_args} \
        --out ${prefix}.ensembletr.vcf

    perl -pe 's/\\bhipstr=/strling=/g; s/AdVNTR, EH, HipSTR, GangSTR/AdVNTR, EH, STRling, GangSTR/' \
        ${prefix}.ensembletr.vcf > ${prefix}.ensembletr.vcf.tmp
    mv ${prefix}.ensembletr.vcf.tmp ${prefix}.ensembletr.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembletr: latest
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ensembletr.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembletr: latest
    END_VERSIONS
    """
}
