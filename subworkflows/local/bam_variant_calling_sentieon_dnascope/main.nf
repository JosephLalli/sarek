//
// SENTIEON HAPLOTYPER germline variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS            as MERGE_SENTIEON_DNASCOPE_GVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS            as MERGE_SENTIEON_DNASCOPE_VCFS  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_DNASCOPE                                           } from '../../../modules/nf-core/sentieon/dnascope/main'

workflow BAM_VARIANT_CALLING_SENTIEON_DNASCOPE {
    take:
    cram                              // channel: [mandatory] [ meta, cram, crai, interval.bed ]
    fasta                             // channel: [mandatory]
    fasta_fai                         // channel: [mandatory]
    dict                              // channel: [mandatory]
    dbsnp                             // channel: [optional]
    dbsnp_tbi                         // channel: [optional]
    dbsnp_vqsr                        // channel: [optional]
    intervals                         // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_germline                    // boolean: [mandatory] [default: false] joint calling of germline variants
    sentieon_dnascope_emit_mode       // string
    sentieon_dnascope_pcr_indel_model // string
    sentieon_dnascope_model           // channel

    main:
    versions = channel.empty()

    gvcf               = channel.empty()
    vcf                = channel.empty()
    genotype_intervals = channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals_for_sentieon = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, _cram, _crai, _intervals, num_intervals -> [
            meta + [
                num_intervals:num_intervals,
                intervals_name:_intervals.simpleName,
                variantcaller:'sentieon_dnascope'],
            _cram,
            _crai,
            _intervals
            ]
        }

    emit_mode_items = sentieon_dnascope_emit_mode.split(',').each{ mode -> mode.toLowerCase().trim() }
    lst = emit_mode_items - 'gvcf'
    emit_vcf = lst.size() > 0 ? lst[0] : ''

    SENTIEON_DNASCOPE(
        cram_intervals_for_sentieon,
        fasta,
        fasta_fai,
        dbsnp.map{ _vcf -> [[:], _vcf]},
        dbsnp_tbi.map{ _tbi -> [[:], _tbi]},
        sentieon_dnascope_model.map{ _model -> [[:], _model]},
        sentieon_dnascope_pcr_indel_model,
        emit_vcf,
        emit_mode_items.any{ it.equals('gvcf') })

    if (joint_germline) {
        genotype_intervals = SENTIEON_DNASCOPE.out.gvcf
            .join(SENTIEON_DNASCOPE.out.gvcf_tbi, failOnMismatch: true)
            .join(cram_intervals_for_sentieon, failOnMismatch: true)
            .map{ meta, _gvcf, _tbi, _cram, _crai, _intervals -> [ meta, _gvcf, _tbi, _intervals ] }
    }

    // Figure out if using intervals or no_intervals
    dnascope_vcf_branch = SENTIEON_DNASCOPE.out.vcf.map{
            meta, _vcf -> [ meta - meta.subMap('interval_name'), _vcf]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    dnascope_vcf_tbi_branch = SENTIEON_DNASCOPE.out.vcf_tbi.map{
            meta, _vcf_tbi -> [ meta - meta.subMap('interval_name'), _vcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_branch = SENTIEON_DNASCOPE.out.gvcf.map{
            meta, _gvcf -> [ meta - meta.subMap('interval_name'), _gvcf]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    haplotyper_gvcf_tbi_branch = SENTIEON_DNASCOPE.out.gvcf_tbi.map{
            meta, _gvcf_tbi -> [ meta - meta.subMap('interval_name'), _gvcf_tbi]
        }
        .branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    vcfs_for_merging = dnascope_vcf_branch.intervals.map{
        meta, _vcf -> [ groupKey(meta, meta.num_intervals), _vcf ]}

    vcfs_for_merging = vcfs_for_merging.map{
        meta, _vcf -> [
            meta - meta.subMap('intervals_name'),
            _vcf]}.groupTuple()

    // VCFs
    // Only when using intervals
    MERGE_SENTIEON_DNASCOPE_VCFS(vcfs_for_merging, dict)

    dnascope_vcf = channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_VCFS.out.vcf,
        dnascope_vcf_branch.no_intervals)

    haplotyper_tbi = channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_VCFS.out.tbi,
        dnascope_vcf_tbi_branch.no_intervals)

    // Remove no longer necessary field: num_intervals
    vcf = dnascope_vcf.map{ meta, _vcf -> [ meta - meta.subMap('num_intervals'), _vcf ] }
    vcf_tbi = haplotyper_tbi.map{ meta, _tbi -> [ meta - meta.subMap('num_intervals'), _tbi ] }

    // GVCFs
    // Only when using intervals
    gvcfs_for_merging = haplotyper_gvcf_branch.intervals.map{
        meta, _gvcf -> [groupKey(meta, meta.num_intervals), _gvcf]}

    gvcfs_for_merging = gvcfs_for_merging.map{
        meta, _gvcf -> [ meta - meta.subMap('intervals_name'), _gvcf ]
    }.groupTuple()

    MERGE_SENTIEON_DNASCOPE_GVCFS(gvcfs_for_merging, dict)

    gvcf = channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_GVCFS.out.vcf,
        haplotyper_gvcf_branch.no_intervals)

    gvcf_tbi = channel.empty().mix(
        MERGE_SENTIEON_DNASCOPE_GVCFS.out.tbi,
        haplotyper_gvcf_tbi_branch.no_intervals)

    versions = versions.mix(SENTIEON_DNASCOPE.out.versions)
    versions = versions.mix(MERGE_SENTIEON_DNASCOPE_VCFS.out.versions)
    versions = versions.mix(MERGE_SENTIEON_DNASCOPE_GVCFS.out.versions)

    emit:
    versions
    vcf
    vcf_tbi
    gvcf
    gvcf_tbi
    genotype_intervals // For joint genotyping

}
