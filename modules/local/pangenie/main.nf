process PANGENIE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://mgibio/pangenie:v4.2.1-bookworm'
        : 'docker.io/mgibio/pangenie:v4.2.1-bookworm'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(panel_vcf)
    tuple val(meta4), path(panel_index)

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
    def reads_arg = reads.collect { it.name.endsWith(".gz") ? "-i <(gzip -cd ${it})" : "-i ${it}" }.join(' ')
    
    // PanGenie-index outputs multiple files. We need to find the prefix used for indexing.
    // The prefix is usually the base name of the VCF used in indexing.
    // We assume the index files are passed as a list of paths.
    def index_prefix = panel_index ? panel_index[0].name.replaceFirst(/(_Graph\.cereal|_UniqueKmersMap\.cereal|_path_segments\.fasta|_kmers\.tsv\.gz)$/, '') : ""
    def index_arg = panel_index ? "-P ${index_prefix}" : ""
    """
    PanGenie \\
        ${reads_arg} \\
        -v ${panel_vcf} \\
        ${index_arg} \\
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
        pangenie: \$(PanGenie --version 2>&1 | head -n 1 | sed 's/PanGenie version: //')
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
        pangenie: 4.2.1
    END_VERSIONS
    """
}
