nextflow.enable.dsl=2

process WHATSHAP_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.8--py312hf731ba3_0' :
        'quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(bam), path(bai)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: phased_vcf
    tuple val(meta), path("*.tbi")   , emit: phased_tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    whatshap phase \
        --output ${prefix}.phased.vcf.gz \
        --reference ${fasta} \
        --ignore-read-groups \
        ${args} \
        ${vcf} \
        ${bam}

    tabix -p vcf ${prefix}.phased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version | sed "s/WhatsHap //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.phased.vcf.gz
    touch ${prefix}.phased.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version | sed "s/WhatsHap //")
    END_VERSIONS
    """
}
