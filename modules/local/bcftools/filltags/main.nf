process BCFTOOLS_FILLTAGS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.{vcf.gz,bcf}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"),    emit: index, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-t AC,AN,AF'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = vcf.name.endsWith('.bcf') ? 'bcf' : 'vcf.gz'
    def output_type = extension == 'bcf' ? 'b' : 'z'
    def write_index = args.contains("--write-index") || args.contains("-W") ? "" : "--write-index"
    """
    bcftools +fill-tags \\
        ${vcf} \\
        -O${output_type} \\
        ${write_index} \\
        -o ${prefix}.filled.${extension} \\
        -- ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = vcf.name.endsWith('.bcf') ? 'bcf' : 'vcf.gz'
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
                args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
                args.contains("--write-index") || args.contains("-W") ? "csi" :
                "tbi" // default to tbi if we added it in script
    """
    touch ${prefix}.filled.${extension}
    touch ${prefix}.filled.${extension}.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
