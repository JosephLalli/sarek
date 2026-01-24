process VCF_SPLIT_BY_REGION {
    tag "${meta.id}-${region_id}"
    label 'process_low'

    container "quay.io/biocontainers/bcftools:1.20--h8b25389_0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(region_id), val(region_string)

    output:
    tuple val(new_meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    new_meta = meta.clone()
    new_meta.region_id = region_id
    
    // Determine if region_string is a file (BED) or a string (CHR:START-END)
    def region_arg = region_string.toString().endsWith('.bed') ? "--regions-file ${region_string}" : "--regions ${region_string}"
    
    """
    bcftools view \\
        ${args} \\
        ${region_arg} \\
        --output ${prefix}.${region_id}.vcf.gz \\
        --output-type z \\
        --threads ${task.cpus} \\
        ${vcf}

    bcftools index --tbi ${prefix}.${region_id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/^.*bcftools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    new_meta = meta.clone()
    new_meta.region_id = region_id
    """
    touch ${prefix}.${region_id}.vcf.gz
    touch ${prefix}.${region_id}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/^.*bcftools //; s/Using.*\$//')
    END_VERSIONS
    """
}
