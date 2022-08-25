process DANBING-TK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'jlalli/pangenie:latest' }"

    input:
    tuple val(meta), path(reads)
    tuple path()

    output:
    tuple val(meta), path("_genotyping.vcf"), emit: vcf
    tuple val(meta), path("_path_segments.fasta"), emit: path_fastas
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    danbing-tk -gc 80       \\
               -ae          \\
               -kf 4 0      \\
               -cth 45      \\
               -k 21        \\
               -qs $ref_file_prefix \\
               -fai $reads  \\
               -p ${task.cpus}  \\
               -o ${prefix} \\
               $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        danbing-tk 1.3
    END_VERSIONS
    """
}
//        # -b fast # SV-workflow-specific argument used - maybe don't here? \\
//        # --hard-hit-cap 500 # SV-workflow-specific argument used - maybe don't here? \\
//--read-group "ID:1 LB:lib1 SM:${meta.id} PL:illumina PU:unit1" \\
