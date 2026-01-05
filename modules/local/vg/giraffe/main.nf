process VG_GIRAFFE {
    tag "${meta.id}"
    label 'process_high'

    // Stage minimizer index as copy (not symlink) so vg can augment it
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://wave.seqera.io/wt/478ac2e4d25c/wave/build:vg-1.70.0_samtools-1.21--ac921c312530e36d'
        : 'wave.seqera.io/wt/478ac2e4d25c/wave/build:vg-1.70.0_samtools-1.21--ac921c312530e36d'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(gbz)
    tuple val(meta3), path(dist)
    tuple val(meta4), path(min)
    path ref_paths
    tuple val(meta5), path(fasta)
    tuple val(meta6), path(fasta_fai)
    tuple val(meta7), path(dict)
    val sort_bam
    val run_fixmate
    val run_markdup

    output:
    tuple val(meta), path("*.gam"),         emit: gam,          optional: true
    tuple val(meta), path("*.tsv"),         emit: stats,        optional: true
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
    def reads_command = meta.single_end ? "-f ${reads}" : "-f ${reads[0]} -f ${reads[1]}"
    def read_group = meta.read_group ? "--read-group '${meta.read_group}'" : ""
    def sample = "--sample ${meta.sample ?: meta.id}"

    // Direct surjection logic
    def fixmate_cpus = task.ext.fixmate_cpus ?: 4
    def sort_memory = "-m ${Math.round(task.memory.bytes * 0.2 / fixmate_cpus)}"

    // Determine output format
    def extension = args2.contains("--output-fmt cram") ? "cram" :
                    args2.contains("-O cram") ? "cram" :
                    args2.contains("-C") ? "cram" :
                    "bam"
    def index_ext = extension == "cram" ? "crai" : "bai"
    def reference = fasta ? "--reference ${fasta}" : ""

    // Surjection parameters
    def paths_arg = ref_paths ? "--into-paths ${ref_paths}" : ''
    
    // Determine if we are doing direct surjection
    // If fasta or ref_paths is provided, we can surject
    def do_surject = fasta || ref_paths

    // Use --supplementary flag if specified or by default for CRAM
    def supplementary_flag = (args.contains("--supplementary") || extension == "cram") ? "--supplementary" : ""
    def surject_format = extension == "cram" ? "--cram-output" : "--bam-output"

    if (do_surject) {
        def reheader_cmd = sort_bam ? "| samtools reheader -P -c 'sed \"s/SN:[^#]*#[^#]*#/SN:/\"' -" : ""
        def fixmate_cmd = run_fixmate ? "| samtools fixmate -O SAM -m --threads ${fixmate_cpus} - -" : ""
        def sort_cmd = sort_bam ? "| samtools sort -@ ${fixmate_cpus} ${sort_memory} ${args2} -u -" : ""
        def markdup_cmd = run_markdup ? "| samtools markdup -u --threads ${fixmate_cpus} -f ${prefix}.markdup.log ${args3} - -" : ""
        def view_cmd = sort_bam ? "| samtools view --threads ${fixmate_cpus} --write-index ${reference} -o ${prefix}.${extension}##idx##${prefix}.${extension}.${index_ext} -" : "> ${prefix}.${extension}"

        """
        vg giraffe \
            ${args} \
            -Z ${gbz} \
            -d ${dist} \
            -m ${min} \
            ${sample} \
            ${read_group} \
            -t ${task.cpus} \
            --progress \
            --report-name ${prefix}.giraffe_stats.tsv \
            ${reads_command} \
        | vg surject \
            -x ${gbz} \
            ${paths_arg} \
            -t ${task.cpus} \
            --sample ${meta.sample ?: meta.id} \
            ${supplementary_flag} \
            ${surject_format} \
            - \
        ${reheader_cmd} \
        ${fixmate_cmd} \
        ${sort_cmd} \
        ${markdup_cmd} \
        ${view_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "*.//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //' | sed 's/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        vg giraffe \
            ${args} \
            -Z ${gbz} \
            -d ${dist} \
            -m ${min} \
            ${sample} \
            ${read_group} \
            -t ${task.cpus} \
            --progress \
            --report-name ${prefix}.giraffe_stats.tsv \
            ${reads_command} \
            > ${prefix}.gam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "*.//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def extension = args2.contains("--output-fmt cram") ? "cram" :
                    args2.contains("-O cram") ? "cram" :
                    args2.contains("-C") ? "cram" :
                    "bam"
    def index_ext = extension == "cram" ? "crai" : "bai"
    def do_surject = fasta || ref_paths
    if (do_surject) {
        """
        touch ${prefix}.${extension}
        touch ${prefix}.${extension}.${index_ext}
        touch ${prefix}.markdup.log
        touch ${prefix}.giraffe_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "*.//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //' | sed 's/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}.gam
        touch ${prefix}.giraffe_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "*.//')
        END_VERSIONS
        """
    }
}
