process ABRA2 {
    tag "$meta.id"
    // label 'process_high'

    conda (params.enable_conda ? "bioconda::abra2=2.24" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abra2:2.24--h9f5acd7_1' :
        'quay.io/biocontainers/abra2:2.24--h9f5acd7_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path (ref_fasta)
    path (bed)

    output:
    tuple val(meta), path("*.bam"), emit: realigned_bam
    path ("*.abra.log")           , emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    print (prefix)
    def avail_mem = 4
    if (!task.memory) {
        log.info '[ABRA2] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    export JAVA_TOOL_OPTIONS="-Xmx${avail_mem}G"

    abra2 --in $bam \\
          --out ${prefix}.bam \\
          --targets $bed \\
          --ref $ref_fasta \\
          --tmpdir /tmp \\
          --threads ${task.cpus} \\
          --index \\
          > ${prefix}.abra.log 2>&1


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: 2.24
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 4
    if (!task.memory) {
        log.info '[ABRA2] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    touch ${prefix}.bam
    touch ${prefix}.abra.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: 2.24
    END_VERSIONS
    """

}
