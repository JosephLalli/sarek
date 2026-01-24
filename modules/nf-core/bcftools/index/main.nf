process BCFTOOLS_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    tuple val(meta), path("*.tbi"), optional:true, emit: tbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools \
        index \
        $args \
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = args.contains("--tbi") || args.contains("-t") ? "tbi" : "csi"
    """
    touch ${vcf}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 
    END_VERSIONS
    """
}
