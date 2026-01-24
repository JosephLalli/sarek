//
// DEEPVARIANT germline calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPVARIANT_RUNDEEPVARIANT                } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { DEEPVARIANT_PANGENOME_AWARE               } from '../../../modules/local/deepvariant/pangenome_aware/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow BAM_VARIANT_CALLING_DEEPVARIANT {
    take:
    input         // channel: [mandatory] [ meta, input, index ] (BAM or CRAM)
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    pangenome_gbz // channel: [optional]  [ meta, gbz ]

    main:
    versions = channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = input.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, _file, _index, _intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], _file, _index, _intervals ]}

    vcf_from_dv = channel.empty()
    gvcf_from_dv = channel.empty()
    tbi_from_dv = channel.empty()

    if (params.tools && params.tools.split(',').contains('deepvariant_pangenome')) {
        DEEPVARIANT_PANGENOME_AWARE(cram_intervals, fasta, fasta_fai, [ [ id:'null' ], [] ], [ [ id:'null' ], [] ], pangenome_gbz.first())
        
        vcf_from_dv = DEEPVARIANT_PANGENOME_AWARE.out.vcf
        gvcf_from_dv = DEEPVARIANT_PANGENOME_AWARE.out.gvcf
        tbi_from_dv = DEEPVARIANT_PANGENOME_AWARE.out.vcf_index
        versions = versions.mix(DEEPVARIANT_PANGENOME_AWARE.out.versions)
    } else {
        DEEPVARIANT_RUNDEEPVARIANT(cram_intervals, fasta, fasta_fai, [ [ id:'null' ], [] ], [ [ id:'null' ], [] ])
        
        vcf_from_dv = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
        gvcf_from_dv = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        tbi_from_dv = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index
        versions = versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = vcf_from_dv.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = gvcf_from_dv.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict)

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_out = tbi_from_dv.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Mix intervals and no_intervals channels together
    gvcf = channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    tbi = channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.tbi, tbi_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], tbi ] }

    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    emit:
    gvcf
    vcf
    tbi

    versions
}