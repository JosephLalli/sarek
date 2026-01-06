nextflow.enable.dsl=2

process SHAPEIT4_PHASECOMMON {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit4:4.2.2--h24bf969_1' :
        'quay.io/biocontainers/shapeit4:4.2.2--h24bf969_1' }"

    input:
    tuple val(meta), path(input), path(input_index), path(reference), path(reference_index), path(map), path(scaffold), path(scaffold_index)
    val region

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: phased_variant
    tuple val(meta), path("*.log")                   , emit: log
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_command = reference ? "--reference ${reference}" : ""
    def scaffold_command = scaffold ? "--scaffold ${scaffold}" : ""
    def map_command = map ? "--map ${map}" : ""
    def region_command = region ? "--region ${region}" : ""

    """
    shapeit4 \
        --input ${input} \
        ${reference_command} \
        ${scaffold_command} \
        ${map_command} \
        ${region_command} \
        --output ${prefix}.bcf \
        --thread ${task.cpus} \
        --log ${prefix}.log \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit4: \$(shapeit4 | head -n 1 | sed "s/^.*v//; s/ .*/")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shapeit4: \$(shapeit4 | head -n 1 | sed "s/^.*v//; s/ .*//")
    END_VERSIONS
    """
}
