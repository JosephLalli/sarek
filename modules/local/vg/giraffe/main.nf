process VG_GIRAFFE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.41.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    tuple path(gbz_index), path(dist_index), path(min_index), path(xg_index), path(gbwt), path(ggbwt)
    path (ref_chrom_names)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path ("*.txt"), emit: alignment_report
    // tuple val(meta), path("*.gam"), emit: gam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readgroup = "${meta.read_group}"
    def reads_command = meta.single_end ? "-f $reads" : "-f ${reads[0]} -f ${reads[1]}"
    def index_command = gbz_index ? "-Z $gbz_index -x $xg_index" : "-H $gbwt -g $ggbwt -x $xg_index"
    def surject_options = "--ref-paths $ref_chrom_names --prune-low-cplx"

    """
    vg giraffe --sample ${meta.id} \\
        $args \\
        $surject_options \\
        --read-group $readgroup \\
        --output-format BAM \\
        $index_command \\
        -d $dist_index \\
        -m $min_index \\
        -t $task.cpus \\
        --progress \\
        --ref-paths $ref_chrom_names \\
        -P \\
        --report-name giraffe_${prefix}.tsv \\
        $reads_command > ${prefix}.bam



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readgroup = "${meta.read_group}"
    def reads_command = meta.single_end ? "-f $reads" : "-f ${reads[0]} -f ${reads[1]}"
    def index_command = gbz_index ? "-Z $gbz_index -x $xg_index" : "-H $gbwt -g $ggbwt -x $xg_index"
    def surject_options = "--ref-paths $ref_chrom_names --prune-low-cplx"

    """
    touch ${prefix}.bam
    touch giraffe_${prefix}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
//        # -b fast # SV-workflow-specific argument used - maybe don't here? \\
//        # --hard-hit-cap 500 # SV-workflow-specific argument used - maybe don't here? \\
//--read-group "ID:1 LB:lib1 SM:${meta.id} PL:illumina PU:unit1" \\
