//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// NOTE: This is a copy of GATK4 mapping. Will transition to giraffe if necessary.

include { VG_GIRAFFE     } from '../../modules/local/vg/giraffe/main'
include { VG_SURJECT     } from '../../modules/local/vg/surject/main'
include { SAMTOOLS_SORT  } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_REHEADER  } from '../../modules/local/samtools/reheader/main'

workflow VG_GIRAFFE_MAP {
    take:
        ch_reads          // channel: [mandatory] meta, reads
        ch_map_index      // channel: [mandatory] GBZ index files
        giraffe_chr_names // [mandatory] txt file of reference graph paths to call against
        sort              // boolean: [mandatory] true -> sort, false -> don't sort

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    // First map reads to VG_GIRAFFE

    VG_GIRAFFE(ch_reads, ch_map_index, giraffe_chr_names)
    ch_bam = VG_GIRAFFE.out.bam
    
    // Then surject packed reads to bam
    // ch_gam = VG_GIRAFFE.out.gam
    // VG_SURJECT(ch_gam, ch_map_index.map { it[3] }, giraffe_chr_names)
    // ch_bam = VG_SURJECT.out.bam
    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions.first())
    ch_logs = ch_logs.mix(VG_GIRAFFE.out.alignment_report)
    // Remove decoy contigs, etc
    SAMTOOLS_REHEADER(ch_bam, giraffe_chr_names)
    ch_bam = SAMTOOLS_REHEADER.out.bam

    // Get the bam files from the aligner, sort if specified
    if (sort){
        SAMTOOLS_SORT(ch_bam)
        ch_bam = SAMTOOLS_SORT.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    }

    emit:
        // gam      = ch_gam // channel: [ [meta], gam ]
        bam      = ch_bam // channel: [ [meta], bam ]
        versions = ch_versions   // channel: [ versions.yml ]
        logs     = ch_logs
}
