process AUTOINDEX_GFA {
    // tag "reference"
    label 'process_med'

    conda (params.enable_conda ? "bioconda::vg=1.41.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.41.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.41.0--h9ee0642_0' }"

    input:
    path (gfa)
    path (ref_chrom_names)

    output:
    path("*.gbz"), emit: giraffe
    path("*.dist"), emit: giraffe_dist
    path("*.min"), emit: giraffe_min
    path("*.pb"), emit: giraffe_snarls
    path("*.xg"), emit: giraffe_xgindex
    path("*.fa"), emit: fasta
    path("*.txt"), emit: contig_lengths
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    vg autoindex \\
        -p reference_pangenome \\
        -w giraffe \\
        -g $gfa \\
        -t ${task.cpus}
        
    vg snarls \\
        -a \\
        -t ${task.cpus}\\
        reference_pangenome.gbz > reference_pangenome.pb

    vg paths \\
        --extract-fasta \\
        -p $ref_chrom_names \\
        --xg reference_pangenome.gbz > ref.fa

    vg paths \\
        --lengths \\
        -p $ref_chrom_names \\
        --xg reference_pangenome.gbz > contig_lengths.txt

    vg convert -x reference_pangenome.gbz > reference_pangenome.xg


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch reference_pangenome.xg
    touch reference_pangenome.gbz
    touch reference_pangenome.min
    touch reference_pangenome.dist
    touch reference_pangenome.pb
    touch ref.fa
    touch ref_paths.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ ".*//')
    END_VERSIONS
    """
}
