process KMC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmc:3.2.4--hf1761c0_0'
        : 'quay.io/biocontainers/kmc:3.2.4--hf1761c0_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.kff"), emit: kmer_file
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kmer_len = task.ext.kmer_size ?: 29
    def mem_gb = task.memory ? task.memory.toGiga() : 8
    def reads_str = meta.single_end ? "${reads}" : "${reads[0]}\\n${reads[1]}"
    """
    echo -e "${reads_str}" > fastqs.txt

    kmc \\
        -k${kmer_len} \\
        -m${mem_gb} \\
        -okff \\
        -t${task.cpus} \\
        ${args} \\
        @fastqs.txt \\
        ${prefix} \\
        .

    mv ${prefix}.kff ${prefix}.kff || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: \$(kmc 2>&1 | head -n 1 | sed 's/.*ver. \\([^ ]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.kff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: \$(kmc 2>&1 | head -n 1 | sed 's/.*ver. \\([^ ]*\\).*/\\1/')
    END_VERSIONS
    """
}
