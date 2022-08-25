process PANGENIE {
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
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.fasta"), emit: path_fastas
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // note: max cpus is number of chromosomes
    // need to unzip and then remove fastqs here;
    // it takes forever, but I will use up 100TB quick if I don't.
    """
    zcat *.fastq.gz > reads.fastq 

    PanGenie -i reads.fastq  \\
             -v $panel_vcf  \\
             -r $ref_fasta  \\
             -o ${prefix}  \\
             -s ${meta.id}  \\
             -j ${task.cpus}   \\
             -t ${task.cpus} -g -d &> ${prefix}.log

    rm *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: 1.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat --version
    touch ${prefix}_genotyping.vcf
    touch ${prefix}_path_segments.fasta
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: 1.0.1
    END_VERSIONS
    """
}
//        # -b fast # SV-workflow-specific argument used - maybe don't here? \\
//        # --hard-hit-cap 500 # SV-workflow-specific argument used - maybe don't here? \\
//--read-group "ID:1 LB:lib1 SM:${meta.id} PL:illumina PU:unit1" \\
