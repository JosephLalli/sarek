//
// SAPPHIRE_REPHASING: Polishing phased haplotypes using reads
//

include { SAPPHIRE_EXTRACTOR   } from '../../../modules/local/sapphire/extractor/main'
include { SAPPHIRE_PHASECALLER } from '../../../modules/local/sapphire/phasecaller/main'
include { SAPPHIRE_UPDATE      } from '../../../modules/local/sapphire/update/main'

workflow SAPPHIRE_PHASE_POLISHING {
    take:
    ch_vcf // channel: [ meta, vcf, tbi, region ]
    ch_bam // channel: [ meta, bam, bai ]

    main:
    ch_versions = channel.empty()

    // 1. Extract variants for polishing
    SAPPHIRE_EXTRACTOR(ch_vcf)
    ch_versions = ch_versions.mix(SAPPHIRE_EXTRACTOR.out.versions)

    // 2. Prepare Phasecaller input
    // Join VCF and Extracted Bin by meta
    ch_vcf_bin = ch_vcf
        .join(SAPPHIRE_EXTRACTOR.out.bin)
        .map { meta, vcf, tbi, region, bin -> [ meta, vcf, tbi, bin ] }

    SAPPHIRE_PHASECALLER(
        ch_vcf_bin,
        ch_bam.map{ it[1] }.collect(),
        ch_bam.map{ it[2] }.collect()
    )
    ch_versions = ch_versions.mix(SAPPHIRE_PHASECALLER.out.versions)

    // 3. Update VCF with polished genotypes
    ch_update_input = ch_vcf
        .join(SAPPHIRE_PHASECALLER.out.bin)
        .map { meta, vcf, tbi, region, polished_bin ->
            [ meta, vcf, tbi, polished_bin ]
        }

    SAPPHIRE_UPDATE(ch_update_input)
    ch_versions = ch_versions.mix(SAPPHIRE_UPDATE.out.versions)

    emit:
    vcf      = SAPPHIRE_UPDATE.out.bcf
    versions = ch_versions
}
