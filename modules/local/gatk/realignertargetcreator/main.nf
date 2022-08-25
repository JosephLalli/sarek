process GATK_REALIGNER_TARGET_CREATOR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk=3.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk:3.8--hdfd78af_11':
        'quay.io/biocontainers/gatk:3.8--hdfd78af_11' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.bed"), emit: indel_target_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def indel_window_expansion_length = task.ext.expansion_length ?: 150

    def avail_mem = 4
    if (!task.memory) {
        log.info '[GATK RealignerTargetCreator] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk RealignerTargetCreator \\
        -Xmx${avail_mem}g \\
        -- remove_program_records \\
        -drf DuplicateRead \\
        --disable_bam_indexing
        -nt ${task.cpus} \\
        -R ${fasta} \\
        -L ${meta.id}
        -I ${bam} \\
        $args \\
        -o ${prefix}.intervals \\

    awk -F '[:-]' 'BEGIN { OFS = "\t" } \\
        { if( \$3 == "") { print \$1, \$2-1, \$2 } \\
        else { print \$1, \$2-1, \$3}}' \\
        forIndelRealigner.intervals \\
        > ${prefix}.intervals.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def indel_window_expansion_length = task.ext.expansion_length ?: 150

    def avail_mem = 4
    if (!task.memory) {
        log.info '[GATK RealignerTargetCreator] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    touch ${prefix}.intervals.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """
}