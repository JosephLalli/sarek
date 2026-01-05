//
// GIRAFFE MAPPING
//
// Pangenome-aware alignment using vg giraffe
// Produces BAM/CRAM files compatible with downstream variant calling

include { VG_GIRAFFE    } from '../../../modules/local/vg/giraffe/main'
include { VG_SURJECT    } from '../../../modules/local/vg/surject/main'
include { VG_STATS      } from '../../../modules/local/vg/stats/main'
include { KMC           } from '../../../modules/local/kmc/main'
include { VG_HAPLOTYPES } from '../../../modules/local/vg/haplotypes/main'
include { VG_GBWT       } from '../../../modules/local/vg/gbwt/main'

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

    ch_giraffe_gbz = ch_gbz

    if (params.pangenome_personalized_flow) {
        // Run KMC to count kmers in reads
        KMC(ch_reads)
        ch_versions = ch_versions.mix(KMC.out.versions.first())

        // Run VG_GBWT to generate r-index from GBZ (needed for haplotypes)
        // We use the common GBZ for all samples.
        // We need to ensure VG_GBWT is configured with "-r" in ext.args in the config or passed here.
        // Since ext.args is module specific, we should rely on config or a new param.
        // For now, let's assume we can configure it.
        VG_GBWT(ch_gbz.map { meta, gbz -> [meta, gbz] }.first()) // Run once for the common graph
        ch_versions = ch_versions.mix(VG_GBWT.out.versions)

        // Run VG_HAPLOTYPES to create personalized GBZ
        // VG_HAPLOTYPES takes [meta, gbz, dist], [meta2, kmer_input], hapl_input, r_index
        // We use the common GBZ and DIST for all samples.
        VG_HAPLOTYPES(
            ch_gbz.combine(ch_dist).map { meta, gbz, meta2, dist -> [meta, gbz, dist] }.collect(),
            KMC.out.kff,
            [], // hapl_input
            VG_GBWT.out.ri.map{ it[1] }.collect()  // r_index from VG_GBWT
        )
        ch_versions = ch_versions.mix(VG_HAPLOTYPES.out.versions.first())

        // Use personalized GBZ for giraffe mapping
        ch_giraffe_gbz = VG_HAPLOTYPES.out.personalized_gbz
    }

    // Align reads using vg giraffe
    VG_GIRAFFE(
        ch_reads,
        ch_giraffe_gbz,
        ch_dist,
        ch_min,
        ch_ref_paths.collect(),
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        sort_bam,
        run_fixmate,
        run_markdup
    )

    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions.first())

    // Collect giraffe alignment stats
    ch_reports = ch_reports.mix(VG_GIRAFFE.out.stats)

    // In direct surjection mode, GAM is not produced by default
    // If we want GAM for stats, we'd need to run a separate giraffe or use tee
    // For optimization, we'll skip VG_STATS if GAM is not available.
    
    if (VG_GIRAFFE.out.gam && (params.tools && params.tools.split(',').contains('vg_stats'))) {
        VG_STATS(VG_GIRAFFE.out.gam)
        ch_reports = ch_reports.mix(VG_STATS.out.stats)
        ch_versions = ch_versions.mix(VG_STATS.out.versions.first())
    }

    // Combine BAM and BAI (or CRAM and CRAI)
    ch_bam_bai = VG_GIRAFFE.out.bam
        .join(VG_GIRAFFE.out.bai, failOnDuplicate: true, failOnMismatch: true)

    ch_cram_crai = VG_GIRAFFE.out.cram
        .join(VG_GIRAFFE.out.crai, failOnDuplicate: true, failOnMismatch: true)

    emit:
    gam       = VG_GIRAFFE.out.gam        // channel: [ meta, gam ]
    bam       = VG_GIRAFFE.out.bam        // channel: [ meta, bam ]
    bai       = VG_GIRAFFE.out.bai        // channel: [ meta, bai ]
    bam_bai   = ch_bam_bai                // channel: [ meta, bam, bai ]
    cram      = VG_GIRAFFE.out.cram       // channel: [ meta, cram ]
    crai      = VG_GIRAFFE.out.crai       // channel: [ meta, crai ]
    cram_crai = ch_cram_crai              // channel: [ meta, cram, crai ]
    reports   = ch_reports                // channel: [ meta, stats ]
    versions  = ch_versions               // channel: [ versions.yml ]
}