process SHAPEIT5_SWITCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--h34261f4_2' :
        'quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2' }"

    input:
    tuple val(meta), path(estimate), path(estimate_index), val(region), path(pedigree), path(truth), path(truth_index), path(freq), path(freq_index)

    output:
    tuple val(meta), path("*.txt.gz"), emit: error_rate
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region_command = region ? "--region ${region}" : ""
    def pedigree_command = pedigree ? "--pedigree ${pedigree}" : ""
    def truth_command = truth ? "--validation ${truth}" : ""
    def freq_command = freq ? "--frequency ${freq}" : ""

    """
    SHAPEIT5_switch \\
        --estimation ${estimate} \\
        ${region_command} \\
        ${pedigree_command} \\
        ${truth_command} \\
        ${freq_command} \\
        --output ${prefix}.error.txt.gz \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        shapeit5: \$(SHAPEIT5_switch 2>&1 | sed -n 's|.*Version[[:space:]]*:[[:space:]]*\\([^/]*[^[:space:]]\\)[[:space:]]*/[[:space:]]*commit[[:space:]]*=[[:space:]]*\\([^/]*\\).*|\\1-\\2|p')
    END_VERSIONS
    """

    

    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.error.txt.gz    

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        shapeit5: \$(SHAPEIT5_switch 2>&1 | sed -n 's|.*Version[[:space:]]*:[[:space:]]*\\([^/]*[^[:space:]]\\)[[:space:]]*/[[:space:]]*commit[[:space:]]*=[[:space:]]*\\([^/]*\\).*|\\1-\\2|p')
    END_VERSIONS
    """
}
