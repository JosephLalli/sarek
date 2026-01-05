process SHAPEIT5_PHASECOMMON {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--h34261f4_2' :
        'quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2' }"

    input:
    tuple val(meta), path(input), path(input_index), path(pedigree), val(region), path(reference), path(reference_index), path(scaffold), path(scaffold_index), path(map)

    output:
    tuple val(meta), path("*.{bcf,graph,bh}"), emit: phased_variant
    tuple val(meta), path("*.log")          , emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pedigree_command = pedigree ? "--pedigree ${pedigree}" : ""
    def reference_command = reference ? "--reference ${reference}" : ""
    def scaffold_command = scaffold ? "--scaffold ${scaffold}" : ""
    def map_command = map ? "--map ${map}" : ""
    def region_command = region ? "--region ${region}" : ""

    """
    SHAPEIT5_phase_common \
        --input ${input} \
        ${pedigree_command} \
        ${reference_command} \
        ${scaffold_command} \
        ${map_command} \
        ${region_command} \
        --output ${prefix}.bcf \
        --thread ${task.cpus} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: 
    END_VERSIONS
    """
}
