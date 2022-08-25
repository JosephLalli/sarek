process DANBING-TK_PREDICT {
    tag "$meta.id"
    label 'process_high'

    // conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'jlalli/pangenie:latest' }"

    input:
    tuple val(meta), path(reads)
    path (ref_fasta)
    path (panel_vcf)

    output:
    tuple val(meta), path("_genotyping.vcf"), emit: vcf
    tuple val(meta), path("_path_segments.fasta"), emit: path_fastas
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // note: max cpus is number of chromosomes
    """
    pangenie -i $reads      \\
             -v $panel_vcf  \\
             -r $ref_fasta  \\
             -o ${prefix}  \\
             -s ${meta.id}  \\
             -j $task.cpu   \\
             -t $task.cpu -g -d &> {log}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: 1.0.1
    END_VERSIONS
    """
}
//        # -b fast # SV-workflow-specific argument used - maybe don't here? \\
//        # --hard-hit-cap 500 # SV-workflow-specific argument used - maybe don't here? \\
//--read-group "ID:1 LB:lib1 SM:${meta.id} PL:illumina PU:unit1" \\
