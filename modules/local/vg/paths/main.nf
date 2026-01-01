process VG_PATHS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vg:1.70.0--h9ee0642_0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    input:
    tuple val(meta), path(graph)
    path ref_paths

    output:
    tuple val(meta), path("*.fa"),      emit: fasta,     optional: true
    tuple val(meta), path("*.fa.fai"),  emit: fasta_fai, optional: true
    tuple val(meta), path("*.txt"),     emit: path_list, optional: true
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paths_arg = ref_paths ? "-p ${ref_paths}" : ''
    def extract_fasta = args.contains('-F') || args.contains('--fasta') ? true : false
    """
    vg paths \\
        ${args} \\
        ${paths_arg} \\
        -x ${graph} \\
        -t ${task.cpus} \\
        > ${prefix}.${extract_fasta ? 'fa' : 'txt'}

    if [[ -f "${prefix}.fa" ]]; then
        samtools faidx ${prefix}.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def extract_fasta = args.contains('-F') || args.contains('--fasta') ? true : false
    """
    touch ${prefix}.${extract_fasta ? 'fa' : 'txt'}
    if [[ "${extract_fasta}" == "true" ]]; then
        touch ${prefix}.fa.fai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
