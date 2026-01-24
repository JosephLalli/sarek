//
// GLNEXUS Joint Calling for DeepVariant
//

include { GLNEXUS        } from '../../../modules/nf-core/glnexus/main'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM   } from '../../../modules/nf-core/bcftools/norm/main'

workflow BAM_JOINT_CALLING_GERMLINE_GLNEXUS {
    take:
    ch_gvcfs // channel: [ val(meta), path(gvcfs), path(custom_config) ]
    ch_bed   // channel: [ val(meta2), path(bed) ]
    fasta    // channel: [ val(meta3), path(fasta) ]

    main:
    ch_versions = channel.empty()

    // 1. Run GLnexus
    GLNEXUS(ch_gvcfs, ch_bed)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // 2. Index the BCF output
    // GLnexus outputs BCF. Index it so BCFTOOLS_NORM can use it
    BCFTOOLS_INDEX(GLNEXUS.out.bcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    // 3. Normalize output VCF
    // Join BCF and its index
    ch_norm_input = GLNEXUS.out.bcf
        .join(BCFTOOLS_INDEX.out.csi)
        .map { meta, bcf, csi -> [ meta, bcf, csi ] }

    BCFTOOLS_NORM(ch_norm_input, fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    emit:
    vcf      = BCFTOOLS_NORM.out.vcf
    tbi      = BCFTOOLS_NORM.out.tbi
    versions = ch_versions
}