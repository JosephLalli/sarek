process JELLYFISH_COUNT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmer-jellyfish:2.3.1--py311pl5321he264feb_6' :
        'quay.io/biocontainers/kmer-jellyfish:2.3.1--py311pl5321he264feb_6' }"

    input:
    tuple val(meta), path(reads)
    path fasta

    output:
    tuple val(meta), path("*.jf"), emit: kmer_file
    path "versions.yml",           emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_cmd = ""
    if (meta.single_end) {
        input_cmd = "<(gzip -cd ${reads})"
    } else {
        input_cmd = "<(gzip -cd ${reads[0]}) <(gzip -cd ${reads[1]})"
    }

    def if_arg = fasta ? "--if ${fasta}" : ""

    """
    jellyfish count \\
        -m 31 \\
        -s 100M \\
        -t ${task.cpus} \\
        -C \\
        -o ${prefix}.jf \\
        ${if_arg} \\
        ${args} \\
        ${input_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version | sed 's/jellyfish //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.jf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version | sed 's/jellyfish //')
    END_VERSIONS
    """
}
