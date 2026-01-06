process SHAPEIT5_PHASERARE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--h34261f4_2' :
        'quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2' }"

    input:
    tuple val(meta), path(input), path(input_index), path(pedigree), val(input_region), path(scaffold), path(scaffold_index), val(scaffold_region), path(map)

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
    def map_command = map ? "--map ${map}" : ""
    def input_region_command = input_region ? "--input-region ${input_region}" : ""
    def scaffold_region_command = scaffold_region ? "--scaffold-region ${scaffold_region}" : ""

    """
    SHAPEIT5_phase_rare \\
        --input ${input} \\
        --scaffold ${scaffold} \\
        ${pedigree_command} \\
        ${map_command} \\
        ${input_region_command} \\
        ${scaffold_region_command} \\
        --output ${prefix}.bcf \\
        --thread ${task.cpus} \\
        ${args} \\
        > ${prefix}.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: \$(SHAPEIT5_phase_rare | head -n 1 | sed 's/^.*v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit5: \$(SHAPEIT5_phase_rare | head -n 1 | sed 's/^.*v//; s/ .*//')
    END_VERSIONS
    """
}
