process VG_GIRAFFE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://wave.seqera.io/wt/478ac2e4d25c/wave/build:vg-1.70.0_samtools-1.21--ac921c312530e36d'
        : 'wave.seqera.io/wt/478ac2e4d25c/wave/build:vg-1.70.0_samtools-1.21--ac921c312530e36d'}"

    input:
    tuple val(meta), path(reads), path(gbz), path(dist), path(min), path(zipcodes)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    path ref_paths
    val sort_bam
    val run_markdup

    output:
    tuple val(meta), path("*.gam"),         emit: gam,          optional: true
    tuple val(meta), path("*.tsv"),         emit: stats,        optional: true
    tuple val(meta), path("*.bam"),         emit: bam,          optional: true
    tuple val(meta), path("*.{bai,csi}"),   emit: bai,          optional: true
    tuple val(meta), path("*.cram"),        emit: cram,         optional: true
    tuple val(meta), path("*.crai"),        emit: crai,         optional: true
    tuple val(meta), path("*.markdup.log"), emit: markdup_stats, optional: true
    path "versions.yml",                    emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-f ${reads}" : "-f ${reads[0]} -f ${reads[1]}"
    def zip_arg = zipcodes ? "--zipcode-name ${zipcodes}" : ''
    def read_group = meta.read_group ? "--read-group '${meta.read_group}'" : ""
    def sample = "--sample ${meta.sample ?: meta.id}"

    // Conditional inputs
    def dist_arg = (dist && dist.name != '[]') ? "-d ${dist}" : ''
    def min_arg = (min && min.name != '[]') ? "-m ${min}" : ''

    // Post-processing logic
    def fixmate_cpus = task.ext.fixmate_cpus ?: 4
    def sort_memory = "-m ${Math.round(task.memory.bytes * 0.2 / fixmate_cpus)}"

    // Determine output format based on args2 (e.g. -O bam) or default to CRAM
    def extension = task.ext.extension
    def index_ext = extension == "cram" ? "crai" : "bai"
    def reference_arg = fasta ? "--reference ${fasta}" : ""
    
    // Surjection/Path args
    def paths_arg = ref_paths ? "--ref-paths ${ref_paths}" : ''
    
    // Command Construction
    // 1. Reheader: Strip pangenome prefixes (e.g. GRCh38#0#chr1 -> chr1)
    def reheader_cmd = "| samtools reheader -P -c 'sed \"s/SN:[^#]*#[^#]*#/SN:/\"' -"

    // 2. Fixmate (Optional) - must happen before sort
    def fixmate_cmd = run_markdup ? "| samtools fixmate -u -m --threads ${fixmate_cpus} - -" : ""

    // 3. Sort (Optional)
    def sort_cmd = sort_bam ? "| samtools sort -@ ${fixmate_cpus} ${sort_memory} -u -" : ""

    // 4. MarkDup (Optional) - must happen after sort
    def markdup_cmd = run_markdup ? "| samtools markdup -u --threads ${fixmate_cpus} -f ${prefix}.markdup.log ${args3} - -" : ""
    
    // 5. Final View (Write output)
    def view_cpus = Math.min(2, task.cpus)
    def output_fmt_flag = task.ext.extension == "cram" ? "-O cram" : "-O bam"
    def view_cmd = "| samtools view --threads ${view_cpus} ${reference_arg} --write-index ${output_fmt_flag} -o ${prefix}.${extension} -"

    if (task.ext.direct_surject) {
        """
        vg giraffe \
            ${args} \
            -o bam \
            -Z ${gbz} \
            -d ${dist} \
            -m ${min} \
            ${zip_arg} \
            ${sample} \
            ${read_group} \
            ${paths_arg} \
            -t ${task.cpus} \
            --progress \
            --report-name ${prefix}.giraffe_stats.tsv \
            ${reads_command} \
        ${reheader_cmd} \
        ${fixmate_cmd} \
        ${sort_cmd} \
        ${markdup_cmd} \
        ${view_cmd} 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        vg giraffe \
            ${args} \
            -Z ${gbz} \
            ${dist_arg} \
            ${min_arg} \
            ${zip_arg} \
            ${sample} \
            ${read_group} \
            -t ${task.cpus} \
            --progress \
            --report-name ${prefix}.giraffe_stats.tsv \
            ${reads_command} \
            > ${prefix}.gam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS

        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def extension = task.ext.extension
    def index_ext = extension == "cram" ? "crai" : "bai"

    if (task.ext.direct_surject) {
        """
        touch ${prefix}.${extension}
        touch ${prefix}.${extension}.${index_ext}
        touch ${prefix}.markdup.log
        touch ${prefix}.giraffe_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}.gam
        touch ${prefix}.giraffe_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
