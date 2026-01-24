process PANGENIE_INDEX {
    tag "${vcf.baseName}"
    label 'process_medium'

    // Add memory requirement for large k-mer hash
    memory { 12.GB * task.attempt }
    
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://mgibio/pangenie:v4.2.1-bookworm'
        : 'docker.io/mgibio/pangenie:v4.2.1-bookworm'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.cereal"), path("*.fasta"), path("*.tsv.gz"), emit: index
    path "versions.yml",                                                 emit: versions

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${vcf.baseName}"
    def input_vcf = vcf.name.endsWith(".gz") ? "input.vcf" : vcf
    def decompress_cmd = vcf.name.endsWith(".gz") ? "gunzip -c ${vcf} > input.vcf" : ""
    """
    ${decompress_cmd}
    
    PanGenie-index \
        -v ${input_vcf} \
        -r ${fasta} \
        -o ${prefix} \
        -t ${task.cpus} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: \$(PanGenie-index --version 2>&1 | head -n 1 | sed 's/PanGenie-index version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${vcf.baseName}"
    """
    touch ${prefix}_Graph.cereal
    touch ${prefix}_UniqueKmersMap.cereal
    touch ${prefix}_path_segments.fasta
    touch ${prefix}_kmers.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: 4.2.1
    END_VERSIONS
    """
}