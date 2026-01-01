process PANGENIE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pangenie:3.0.2--h4ac6f70_0'
        : 'quay.io/biocontainers/pangenie:3.0.2--h4ac6f70_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(panel_vcf)

    output:
    tuple val(meta), path("*_genotyping.vcf.gz"),     emit: vcf
    tuple val(meta), path("*_genotyping.vcf.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(meta), path("*.log"),                   emit: log
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_arg = reads.collect { "-i ${it}" }.join(' ')
    """
    pangenie \\
        ${reads_arg} \\
        -v ${panel_vcf} \\
        -r ${reference} \\
        -o ${prefix} \\
        -s ${meta.sample ?: meta.id} \\
        -t ${task.cpus} \\
        ${args} \\
        2> ${prefix}.log

    bgzip -@ ${task.cpus} ${prefix}_genotyping.vcf
    tabix -p vcf ${prefix}_genotyping.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: \$(pangenie --version 2>&1 | head -n 1 | sed 's/PanGenie version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genotyping.vcf.gz
    touch ${prefix}_genotyping.vcf.gz.tbi
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: \$(pangenie --version 2>&1 | head -n 1 | sed 's/PanGenie version: //')
    END_VERSIONS
    """
}
