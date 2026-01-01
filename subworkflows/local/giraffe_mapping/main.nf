//
// GIRAFFE MAPPING
//
// Pangenome-aware alignment using vg giraffe
// Produces BAM files compatible with downstream variant calling

include { VG_GIRAFFE  } from '../../../modules/local/vg/giraffe/main'
include { VG_SURJECT  } from '../../../modules/local/vg/surject/main'
include { VG_STATS    } from '../../../modules/local/vg/stats/main'
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow GIRAFFE_MAPPING {
    take:
    ch_reads           // channel: [mandatory] [ meta, reads ]
    ch_gbz             // channel: [mandatory] [ meta, gbz ]
    ch_dist            // channel: [mandatory] [ meta, dist ]
    ch_min             // channel: [mandatory] [ meta, min ]
    ch_ref_paths       // channel: [optional]  ref_paths.txt for surjection
    ch_fasta           // channel: [optional]  [ meta, fasta ] for CRAM output
    ch_fasta_fai       // channel: [optional]  [ meta, fasta_fai ]

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    // Align reads to pangenome graph
    VG_GIRAFFE(
        ch_reads,
        ch_gbz,
        ch_dist,
        ch_min
    )

    ch_gam = VG_GIRAFFE.out.gam
    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions.first())

    // Collect giraffe alignment stats
    ch_reports = ch_reports.mix(VG_GIRAFFE.out.stats)

    // Compute GAM statistics
    VG_STATS(ch_gam)
    ch_reports = ch_reports.mix(VG_STATS.out.stats)
    ch_versions = ch_versions.mix(VG_STATS.out.versions.first())

    // Surject GAM to BAM (project onto reference paths)
    VG_SURJECT(
        ch_gam,
        ch_gbz,
        ch_ref_paths.collect()
    )

    ch_versions = ch_versions.mix(VG_SURJECT.out.versions.first())

    // Sort BAM
    SAMTOOLS_SORT(
        VG_SURJECT.out.bam,
        ch_fasta
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Index BAM
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Combine BAM and BAI
    ch_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true)

    emit:
    gam      = ch_gam                    // channel: [ meta, gam ]
    bam      = SAMTOOLS_SORT.out.bam     // channel: [ meta, bam ]
    bai      = SAMTOOLS_INDEX.out.bai    // channel: [ meta, bai ]
    bam_bai  = ch_bam_bai                // channel: [ meta, bam, bai ]
    reports  = ch_reports                // channel: [ meta, stats ]
    versions = ch_versions               // channel: [ versions.yml ]
}
