process SAPPHIRE_PHASECALLER {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/jlalli/sapphire:53a96e5ac"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bin)
    path bams
    path bais

    output:
    tuple val(meta), path("*.polished.bin"), emit: bin
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Create the samples file for Sapphire
    # Format: SampleID,SampleName,CRAMPath
    # Using a simple loop to create the mapping from the input BAMs
    
    rm -f samples.csv
    
    # Get ordered list of samples from VCF
    bcftools query -l ${vcf} > vcf_samples.txt
    
    for bam in ${bams}; do
        # Extract sample name from BAM header or filename
        SAMPLENAME=\$(samtools view -H \$bam | grep '^@RG' | sed 's/.*SM:\\([^\\t]*\\).*/\\1/' | head -n 1)
        if [ -z "\$SAMPLENAME" ]; then
            SAMPLENAME=\$(basename \$bam .bam)
        fi
        
        # Find index of sample in VCF (0-based)
        SAMPLE_IDX=\$(grep -n "^\${SAMPLENAME}\$" vcf_samples.txt | cut -d: -f1)
        if [ -z "\$SAMPLE_IDX" ]; then
            echo "Sample \${SAMPLENAME} not found in VCF"
            exit 1
        fi
        # Convert 1-based grep line number to 0-based index
        SAMPLE_IDX=\$((SAMPLE_IDX - 1))
        
        echo "\$SAMPLE_IDX,\$SAMPLENAME,\$bam" >> samples.csv
    done

    # Copy input bin to avoid modifying the input file in the work dir directly
    cp ${bin} ${prefix}.polished.bin

    /usr/src/sapphire/bin/phase_caller \\
        -f ${vcf} \\
        -S samples.csv \\
        --cram-path-from-samples-file \\
        -b ${prefix}.polished.bin \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sapphire: \$(/usr/src/sapphire/bin/pp_extract --version 2>&1 | head -n 1 | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.polished.bin
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sapphire: 1.0.0
    END_VERSIONS
    """
}
