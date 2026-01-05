process KMC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmc:3.2.1--hf1761c0_2' :
        'quay.io/biocontainers/kmc:3.2.1--hf1761c0_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.kff.gz"), emit: kff
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kmer_len = task.ext.kmer_size ?: 29
    
    // Handle list of files
    def input_files = reads instanceof List ? reads.join("\n") : "${reads}"
    
    """
    mkdir kmc_tmp
    printf "${input_files}" > fastqs.txt

    kmc \
        -k${kmer_len} \
        -ci2 \
        -okff \
        -t${task.cpus} \
        -m${task.memory.giga} \
        ${args} \
        @fastqs.txt \
        ${prefix} \
        kmc_tmp

    gzip ${prefix}.kff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.kff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: 
    END_VERSIONS
    """
}
