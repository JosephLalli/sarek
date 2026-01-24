//
// STRling analysis subworkflow
//

include { STRLING_EXTRACT } from '../../../modules/local/strling/extract/main'
include { STRLING_MERGE   } from '../../../modules/local/strling/merge/main'
include { STRLING_CALL    } from '../../../modules/local/strling/call/main'

workflow STR_ANALYSIS_STRLING {
    take:
    ch_alignment    // channel: [ val(meta), bam_cram, bai_crai ]
    ch_fasta        // channel: [ val(meta), fasta ]
    ch_fasta_fai    // channel: [ val(meta), fasta_fai ]
    ch_str_index    // channel: [ val(meta), str_index ]
    ch_strling_loci // channel: [ val(meta), strling_loci ]
    val_joint       // boolean: true/false

    main:
    ch_versions = channel.empty()

    // 1. Extract
    // Combine str_index and strling_loci for extraction
    ch_extract_loci = ch_str_index.map{ it[1] }
        .mix(ch_strling_loci.map{ it[1] })
        .collect()
        .ifEmpty([])

    STRLING_EXTRACT(ch_alignment, ch_fasta.collect(), ch_fasta_fai.collect(), ch_extract_loci)
    ch_versions = ch_versions.mix(STRLING_EXTRACT.out.versions)

    // 2. Prepare inputs for calling
    // Join on meta.id
    ch_bam_bin = ch_alignment.join(STRLING_EXTRACT.out.bin) 

    if (val_joint) {
        // Collect all bins for merging
        ch_bins_collected = STRLING_EXTRACT.out.bin.map{ it[1] }.collect()

        STRLING_MERGE(ch_bins_collected, ch_fasta.collect(), ch_fasta_fai.collect())
        ch_versions = ch_versions.mix(STRLING_MERGE.out.versions)
        
        // Pass merged bounds to call
        STRLING_CALL(
            ch_bam_bin,
            ch_fasta.collect(),
            ch_fasta_fai.collect(),
            STRLING_MERGE.out.bounds.collect() // Use collect to convert to value channel for reuse
        )
    } else {
        // Single sample calling (no bounds provided)
        STRLING_CALL(
            ch_bam_bin,
            ch_fasta.collect(),
            ch_fasta_fai.collect(),
            []
        )
    }
    
    ch_versions = ch_versions.mix(STRLING_CALL.out.versions)

    emit:
    genotype = STRLING_CALL.out.genotype
    bounds   = STRLING_CALL.out.bounds
    versions = ch_versions
}