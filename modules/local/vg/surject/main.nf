process VG_SURJECT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-3e7dd4852eb2d35e55d00a00eec396d9bc28e4b3:8771576bd6a6b760e7981c30558d33a68a93ef18-0'
        : 'quay.io/biocontainers/mulled-v2-3e7dd4852eb2d35e55d00a00eec396d9bc28e4b3:8771576bd6a6b760e7981c30558d33a68a93ef18-0'}"

    input:
    tuple val(meta), path(gam)
    tuple val(meta2), path(graph)
    path ref_paths
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)
    tuple val(meta5), path(dict)
    val sort_bam
    val run_fixmate
    val run_markdup

    output:
    tuple val(meta), path("*.bam"),         emit: bam,          optional: true
    tuple val(meta), path("*.bai"),         emit: bai,          optional: true
    tuple val(meta), path("*.cram"),        emit: cram,         optional: true
    tuple val(meta), path("*.crai"),        emit: crai,         optional: true
    tuple val(meta), path("*.markdup.log"), emit: markdup_stats, optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fixmate_cpus = task.ext.fixmate_cpus ?: 4
    def sort_memory = "-m ${Math.round(task.memory.bytes * 0.2 / fixmate_cpus)}"

    // Determine output format
    def extension = args2.contains("--output-fmt cram") ? "cram" :
                    args2.contains("-O cram") ? "cram" :
                    args2.contains("-C") ? "cram" :
                    "bam"
    def index_ext = extension == "cram" ? "crai" : "bai"
    def reference = fasta ? "--reference ${fasta}" : ""

    // Graph argument based on file extension
    def graph_arg = graph.name.endsWith('.gbz') ? "-Z ${graph}" :
                    graph.name.endsWith('.xg') ? "-x ${graph}" :
                    "-g ${graph}"

    // Paths argument
    def paths_arg = ref_paths ? "--into-paths ${ref_paths}" : ''

    // Build the pipeline conditionally
    def reheader_cmd = sort_bam ? "| samtools reheader -P -c 'sed \"s/SN:[^#]*#[^#]*#/SN:/\"' -" : ""
    def fixmate_cmd = run_fixmate ? "| samtools fixmate -O SAM -m --threads ${fixmate_cpus} - -" : ""
    def sort_cmd = sort_bam ? "| samtools sort -@ ${fixmate_cpus} ${sort_memory} ${args2} -u -" : ""
    def markdup_cmd = run_markdup ? "| samtools markdup -u --threads ${fixmate_cpus} -f ${prefix}.markdup.log ${args3} - -" : ""
    def view_cmd = sort_bam ? "| samtools view --threads ${fixmate_cpus} --write-index ${reference} -o ${prefix}.${extension}##idx##${prefix}.${extension}.${index_ext} -" : "> ${prefix}.${extension}"

    // If not sorting, output raw BAM from vg surject
    def surject_output_flag = sort_bam ? "-b" : (extension == "bam" ? "-b" : "")

    """
    export REF_PATH=\$(dirname \$(readlink -f ${fasta}))

    vg surject \\
        ${args} \\
        ${graph_arg} \\
        ${paths_arg} \\
        -t ${task.cpus} \\
        --sample ${meta.sample ?: meta.id} \\
        ${surject_output_flag} \\
        ${gam} \\
    ${reheader_cmd} \\
    ${fixmate_cmd} \\
    ${sort_cmd} \\
    ${markdup_cmd} \\
    ${view_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt cram") ? "cram" :
                    args2.contains("-O cram") ? "cram" :
                    args2.contains("-C") ? "cram" :
                    "bam"
    def index_ext = extension == "cram" ? "crai" : "bai"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.${extension}.${index_ext}
    touch ${prefix}.markdup.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
