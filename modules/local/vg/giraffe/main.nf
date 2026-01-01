process VG_GIRAFFE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(gbz)
    tuple val(meta3), path(dist)
    tuple val(meta4), path(min)

    output:
    tuple val(meta), path("*.gam"), emit: gam
    tuple val(meta), path("*.tsv"), emit: stats, optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-f ${reads}" : "-f ${reads[0]} -f ${reads[1]}"
    def read_group = meta.read_group ? "--read-group '${meta.read_group}'" : ""
    def sample = "--sample ${meta.sample ?: meta.id}"
    """
    vg giraffe \\
        ${args} \\
        -Z ${gbz} \\
        -d ${dist} \\
        -m ${min} \\
        ${sample} \\
        ${read_group} \\
        -t ${task.cpus} \\
        --progress \\
        --report-name ${prefix}.giraffe_stats.tsv \\
        ${reads_command} \\
        > ${prefix}.gam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gam
    touch ${prefix}.giraffe_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
