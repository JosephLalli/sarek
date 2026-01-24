process VG_HAPLOTYPES_PERSONALIZED {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/jlalli/vg:1.54.0"

    input:
    tuple val(meta), path(kmers)
    tuple path(gbz), path(dist), path(hapl), path(min)

    output:
    tuple val(meta), path("*.gbz"), emit: gbz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    vg haplotypes \
        -t ${task.cpus} \
        ${args} \
        -i ${hapl} \
        -k ${kmers} \
        -g ${prefix}.gbz \
        ${gbz}

    # Optional: Re-index if needed, or leave as is. 
    # Legacy had: mv ${prefix}.gbz tmp.gbz && vg gbwt --set-reference recombination -Z tmp.gbz -g ${prefix}.gbz --gbz-format && rm tmp.gbz
    # We will assume vg haplotypes output is sufficient unless specified via args or update later.
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: 
$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gbz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: 
$(vg version | head -n 1 | sed 's/vg version //' | sed 's/ "//')
    END_VERSIONS
    """
}
