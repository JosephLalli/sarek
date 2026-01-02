//
// GIRAFFE MAPPING
//
// Pangenome-aware alignment using vg giraffe
// Produces BAM/CRAM files compatible with downstream variant calling

include { VG_GIRAFFE  } from '../../../modules/local/vg/giraffe/main'
include { VG_SURJECT  } from '../../../modules/local/vg/surject/main'
include { VG_STATS    } from '../../../modules/local/vg/stats/main'

workflow GIRAFFE_MAPPING {
    take:
    ch_reads           // channel: [mandatory] [ meta, reads ]
    ch_gbz             // channel: [mandatory] [ meta, gbz ]
    ch_dist            // channel: [mandatory] [ meta, dist ]
    ch_min             // channel: [mandatory] [ meta, min ]
    ch_ref_paths       // channel: [optional]  ref_paths.txt for surjection
    ch_fasta           // channel: [mandatory] [ meta, fasta ] for sorting/CRAM
    ch_fasta_fai       // channel: [mandatory] [ meta, fasta_fai ]
    ch_dict            // channel: [mandatory] [ meta, dict ]
    sort_bam           // val: whether to sort BAM
    run_fixmate        // val: whether to run fixmate
    run_markdup        // val: whether to run markdup

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

    // Surject GAM to BAM/CRAM (project onto reference paths)
    // Includes reheader, fixmate, sort, markdup, and indexing
    VG_SURJECT(
        ch_gam,
        ch_gbz,
        ch_ref_paths.collect(),
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        sort_bam,
        run_fixmate,
        run_markdup
    )

    ch_versions = ch_versions.mix(VG_SURJECT.out.versions.first())

    // Collect markdup stats if available
    ch_reports = ch_reports.mix(VG_SURJECT.out.markdup_stats)

    // Combine BAM and BAI (or CRAM and CRAI)
    ch_bam_bai = VG_SURJECT.out.bam
        .join(VG_SURJECT.out.bai, failOnDuplicate: true, failOnMismatch: true)

    ch_cram_crai = VG_SURJECT.out.cram
        .join(VG_SURJECT.out.crai, failOnDuplicate: true, failOnMismatch: true)

    emit:
    gam       = ch_gam                    // channel: [ meta, gam ]
    bam       = VG_SURJECT.out.bam        // channel: [ meta, bam ]
    bai       = VG_SURJECT.out.bai        // channel: [ meta, bai ]
    bam_bai   = ch_bam_bai                // channel: [ meta, bam, bai ]
    cram      = VG_SURJECT.out.cram       // channel: [ meta, cram ]
    crai      = VG_SURJECT.out.crai       // channel: [ meta, crai ]
    cram_crai = ch_cram_crai              // channel: [ meta, cram, crai ]
    reports   = ch_reports                // channel: [ meta, stats ]
    versions  = ch_versions               // channel: [ versions.yml ]
}
