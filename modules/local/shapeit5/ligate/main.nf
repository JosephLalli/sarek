process SHAPEIT5_LIGATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--h34261f4_2' :
        'quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2' }"

    input:
    tuple val(meta), path(input_list), path(input_list_index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variant
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "%s\\n" ${input_list.join(' ')} > vcf_list.txt

    SHAPEIT5_ligate \\
        --input vcf_list.txt \\
        --output ${prefix}.bcf \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: \$(SHAPEIT5_ligate | head -n 1 | sed 's/^.*v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: 
    END_VERSIONS
    """
}
