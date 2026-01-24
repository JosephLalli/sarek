process DEEPVARIANT_PANGENOME_AWARE {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepvariant:pangenome_aware_deepvariant-1.9.0"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(par_bed)
    tuple val(meta6), path(pangenome_gbz)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")             , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.{tbi,csi}")   , emit: vcf_index
    tuple val(meta), path("${prefix}.g.vcf.gz")           , emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.{tbi,csi}") , emit: gvcf_index, optional: true
    tuple val(meta), path("${prefix}.visual_report.html") , emit: report, optional: true
    path "versions.yml"                                   , emit: versions


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""
    def par_regions = par_bed ? "--par_regions_bed=${par_bed}" : ""

    """
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \\
        --ref=${fasta} \\
        --reads=${input} \\
        --pangenome=${pangenome_gbz} \\
        --output_vcf=${prefix}.vcf.gz \\
        ${regions} \\
        ${par_regions} \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=tmp \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: 1.9.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: 1.9.0
    END_VERSIONS
    """
}
