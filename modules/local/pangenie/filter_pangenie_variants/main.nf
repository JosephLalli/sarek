process FILTER_PANGENIE_VARIANTS {
    tag "${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'jlalli/process_pangenie_variants:1.0.2' }"

    input:
    tuple val(meta), path (merged_sv_calls)
    path (bi_panel_vcf)
    path (model_file)

    output:
    path ("*.parquet"), emit: parquet_report
    path ("*.csv"), emit: csv
    path ("*.log"), emit: log
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "filtered_pangenie_calls"
    def cutoff = task.ext.cutoff ?: '0.5'
    // produces filtered_pangenie_calls.csv, filtered_pangenie_calls.parquet
    """
    postprocess_vcf.py -g $merged_sv_calls \\
                       -p $bi_panel_vcf \\
                       --cutoff $cutoff \\
                       --model $model_file \\
                       --SVs \\
                       --medium_indels \\
                       -o filtered_pangenie_calls &> filtered_pangenie_calls.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python postprocess_vcf.py --version
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "filtered_pangenie_calls"
    def cutoff = task.ext.cutoff ?: '0.5'

    """
    touch filtered_pangenie_calls.csv
    touch filtered_pangenie_calls.parquet
    touch filtered_pangenie_calls.log
    postprocess_vcf.py --version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_pangenie_variants: \$( echo \$(postprocess_vcf.py --version 2>&1))
    END_VERSIONS
    """

}