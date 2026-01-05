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
    path hapl_input
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

    def kmer_prep = ""
    def kmer_arg = ""
    if (kmer_input) {
        if (kmer_input.name.endsWith('.gz')) {
            kmer_prep = "gzip -cd ${kmer_input} > input.kff"
            kmer_arg = "-k input.kff"
        } else {
            kmer_arg = "-k ${kmer_input}"
        }
    }

    def hapl_arg = (hapl_input && !(hapl_input instanceof List && hapl_input.isEmpty())) ? "-i ${hapl_input}" : "-H ${prefix}.hapl"
    def r_arg = (r_index && !(r_index instanceof List && r_index.isEmpty())) ? "-r ${r_index}" : ''
    def dist_arg = (dist && !(dist instanceof List && dist.isEmpty())) ? "-d ${dist}" : ''

    """
    ${kmer_prep}

    vg haplotypes \\
        -t ${task.cpus} \\
        ${args} \\
        ${hapl_arg} \\
        ${r_arg} \\
        ${dist_arg} \\
        ${kmer_arg} \\
        --include-reference \\
        -g ${prefix}.gbz \\
        ${gbz}

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
