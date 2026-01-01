process VG_HAPLOTYPES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(gbz), path(dist)
    tuple val(meta2), path(kmer_input)
    path hapl_index
    path r_index

    output:
    tuple val(meta), path("*.gbz"),  emit: personalized_gbz, optional: true
    tuple val(meta), path("*.hapl"), emit: hapl, optional: true
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kmer_arg = kmer_input ? "-k ${kmer_input}" : ''
    def hapl_arg = hapl_index ? "-i ${hapl_index}" : "-H ${prefix}.hapl"
    def r_arg = r_index ? "-r ${r_index}" : ''
    def dist_arg = dist ? "-d ${dist}" : ''
    """
    vg haplotypes \\
        -t ${task.cpus} \\
        -g ${gbz} \\
        ${args} \\
        ${hapl_arg} \\
        ${r_arg} \\
        ${dist_arg} \\
        ${kmer_arg} \\
        -o ${prefix}.gbz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gbz
    touch ${prefix}.hapl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
