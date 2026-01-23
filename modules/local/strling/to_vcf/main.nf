process STRLING_TO_VCF {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/pysam_samtools_pip_ensembletr:343f01fd8b77ac0a"

    input:
    tuple val(meta), path(genotype)
    path(catalog)
    path(ref_fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
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
        --ref-fai $ref_fai \
        --output ${prefix}.vcf \
        --sample ${meta.id}

    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz
    rm ${prefix}.vcf

    printf '"%s":\n    strling_to_vcf: v0.1\n' "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 -c 'import gzip; f=gzip.open("'"${prefix}.vcf.gz"'","wt"); f.write(""); f.close()'
    touch ${prefix}.vcf.gz.tbi

    printf '"%s":\n    strling_to_vcf: v0.1\n' "${task.process}" > versions.yml
    """
}
