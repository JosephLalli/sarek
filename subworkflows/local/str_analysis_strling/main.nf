//
// STRling analysis subworkflow
//

include { STRLING_EXTRACT } from '../../../modules/local/strling/extract/main'
include { STRLING_MERGE   } from '../../../modules/local/strling/merge/main'
include { STRLING_CALL    } from '../../../modules/local/strling/call/main'

workflow STR_ANALYSIS_STRLING {
    take:
    ch_bam          // channel: [ val(meta), bam, bai ]
    ch_fasta        // channel: [ val(meta), fasta ]
    ch_fasta_fai    // channel: [ val(meta), fasta_fai ]
    ch_str_index    // channel: [ val(meta), str_index ]
    val_joint       // boolean: true/false

    main:
    ch_versions = channel.empty()

    // 1. Extract
    STRLING_EXTRACT(ch_bam, ch_fasta, ch_fasta_fai, ch_str_index.map{ it[1] }.collect().ifEmpty([]))
    ch_versions = ch_versions.mix(STRLING_EXTRACT.out.versions)

    // 2. Prepare inputs for calling
    // Join on meta.id
    ch_bam_bin = ch_bam.join(STRLING_EXTRACT.out.bin) 

    if (val_joint) {
        // Collect all bins for merging
        ch_bins_collected = STRLING_EXTRACT.out.bin.map{ it[1] }.collect()

        STRLING_MERGE(ch_bins_collected, ch_fasta, ch_fasta_fai)
        ch_versions = ch_versions.mix(STRLING_MERGE.out.versions)
        
        // Pass merged bounds to call
        STRLING_CALL(
            ch_bam_bin,
            ch_fasta,
            ch_fasta_fai,
            STRLING_MERGE.out.bounds.collect() // Use collect to convert to value channel for reuse
        )
    } else {
        // Single sample calling (no bounds provided)
        STRLING_CALL(
            ch_bam_bin,
            ch_fasta,
            ch_fasta_fai,
            []
        )
    }
    
    ch_versions = ch_versions.mix(STRLING_CALL.out.versions)

    emit:
    vcf      = STRLING_CALL.out.vcf
    genotype = STRLING_CALL.out.genotype
    bounds   = STRLING_CALL.out.bounds
    versions = ch_versions
}