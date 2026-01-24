//
// PHASE AND IMPUTE
//
// Phasing and imputation using SHAPEIT5, WhatsHap, and Hybrid methods

include { SHAPEIT5_PHASECOMMON } from '../../../modules/local/shapeit5/phase_common/main'
include { SHAPEIT5_PHASERARE   } from '../../../modules/local/shapeit5/phase_rare/main'
include { SHAPEIT5_LIGATE      } from '../../../modules/local/shapeit5/ligate/main'
include { SHAPEIT4_PHASECOMMON } from '../../../modules/local/shapeit4/phase_common/main'
include { WHATSHAP_PHASE       } from '../../../modules/nf-core/whatshap/phase/main'
include { BCFTOOLS_VIEW                             } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER_MAF } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_COMMON   } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_WHATSHAP } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SHAPEIT4 } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_FILLTAGS                         } from '../../../modules/local/bcftools/filltags/main'

workflow PHASE_AND_IMPUTE {
    take:
    ch_input // channel: [ meta, vcf, tbi, region, reference, reference_index, map, scaffold, scaffold_index, bam, bai ]
    ch_fasta // channel: [ meta, fasta ]
    ch_fai   // channel: [ meta, fai ]

    main:
    ch_versions = channel.empty()
    ch_phased_vcf = channel.empty()

    // Normalize input for easier mapping
    ch_normalized = ch_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai ->
        [ meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai ]
    }

    // PATH 1: Population Phasing (SHAPEIT5)
    if (params.phasing_method == 'population') {
        
        // Phase Common
        SHAPEIT5_PHASECOMMON(
            ch_normalized.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai ->
                [ meta, vcf, tbi, [], region, ref, ref_idx, scaff, scaff_idx, map ]
            }
        )
        ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions)

        // Index common phased output (BCF needs CSI)
        BCFTOOLS_INDEX_COMMON(SHAPEIT5_PHASECOMMON.out.phased_variant)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_COMMON.out.versions)

        ch_index = BCFTOOLS_INDEX_COMMON.out.csi.mix(BCFTOOLS_INDEX_COMMON.out.tbi)

        // Phase Rare (Optional - only if requested or if rare variants expected)
        if (params.run_phasing_rare) {
            ch_rare_input = ch_normalized
                .join(SHAPEIT5_PHASECOMMON.out.phased_variant)
                .join(ch_index)
                .map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, common_vcf, common_idx ->
                    [ meta, vcf, tbi, [], region, common_vcf, common_idx, region, map ]
                }

            SHAPEIT5_PHASERARE(ch_rare_input)
            ch_versions = ch_versions.mix(SHAPEIT5_PHASERARE.out.versions)
            ch_phased_vcf = SHAPEIT5_PHASERARE.out.phased_variant
        } else {
            ch_phased_vcf = SHAPEIT5_PHASECOMMON.out.phased_variant
        }
    }
    // PATH 2: Read-backed (WhatsHap)
    else if (params.phasing_method == 'read_backed') {
        
        ch_whatshap_input = ch_normalized
            .join(ch_fasta)
            .join(ch_fai)
        
        WHATSHAP_PHASE(
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, vcf, tbi] },
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, bam, bai] },
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, fasta, fai] }
        )
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)
        ch_phased_vcf = WHATSHAP_PHASE.out.phased_vcf
    }
    // PATH 4: Hybrid (WhatsHap -> SHAPEIT4 -> SHAPEIT5_rare)
    else if (params.phasing_method == 'hybrid') {
        
        // 1. WhatsHap
        ch_whatshap_input = ch_normalized
            .join(ch_fasta)
            .join(ch_fai)
            
        WHATSHAP_PHASE(
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, vcf, tbi] },
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, bam, bai] },
            ch_whatshap_input.map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, fasta, fai -> [meta, fasta, fai] }
        )
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)
        
        // Index WhatsHap output
        BCFTOOLS_INDEX_WHATSHAP(WHATSHAP_PHASE.out.phased_vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_WHATSHAP.out.versions)
        ch_whatshap_index = BCFTOOLS_INDEX_WHATSHAP.out.csi.mix(BCFTOOLS_INDEX_WHATSHAP.out.tbi)

        // Filter input VCF to keep only common variants (MAF >= threshold) for SHAPEIT4
        // This leaves rare variants unphased so SHAPEIT5_PHASERARE can phase them
        ch_maf_filter_input = ch_normalized
            .map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai ->
                [ meta, vcf, tbi ]
            }
        BCFTOOLS_VIEW_FILTER_MAF(ch_maf_filter_input, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_FILTER_MAF.out.versions)
        ch_maf_filtered_index = BCFTOOLS_VIEW_FILTER_MAF.out.csi.mix(BCFTOOLS_VIEW_FILTER_MAF.out.tbi)

        // 2. SHAPEIT4 (Input = MAF-filtered VCF, Scaffold = WhatsHap)
        ch_shapeit4_input = ch_normalized
            .join(BCFTOOLS_VIEW_FILTER_MAF.out.vcf)
            .join(ch_maf_filtered_index)
            .join(WHATSHAP_PHASE.out.phased_vcf)
            .join(ch_whatshap_index)
            .map { meta, orig_vcf, orig_tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, filt_vcf, filt_idx, wh_vcf, wh_idx ->
                [ meta, filt_vcf, filt_idx, ref, ref_idx, map, wh_vcf, wh_idx, region ]
            }

        SHAPEIT4_PHASECOMMON(
            ch_shapeit4_input.map { it[0..7] }, // tuple inputs
            ch_shapeit4_input.map { it[8] }     // region
        )
        ch_versions = ch_versions.mix(SHAPEIT4_PHASECOMMON.out.versions)
        
        // Index SHAPEIT4 output
        BCFTOOLS_INDEX_SHAPEIT4(SHAPEIT4_PHASECOMMON.out.phased_variant)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_SHAPEIT4.out.versions)
        ch_shapeit4_index = BCFTOOLS_INDEX_SHAPEIT4.out.csi.mix(BCFTOOLS_INDEX_SHAPEIT4.out.tbi)

        // Add AC/AN/AF tags to SHAPEIT4 output (required by SHAPEIT5_PHASERARE)
        BCFTOOLS_FILLTAGS(
            SHAPEIT4_PHASECOMMON.out.phased_variant
                .join(ch_shapeit4_index)
        )
        ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

        ch_filled_index = BCFTOOLS_FILLTAGS.out.index

        // 3. SHAPEIT5_RARE (Scaffold = SHAPEIT4 with tags)
        ch_rare_input = ch_normalized
            .join(BCFTOOLS_FILLTAGS.out.vcf)
            .join(ch_filled_index)
            .map { meta, vcf, tbi, region, ref, ref_idx, map, scaff, scaff_idx, bam, bai, s4_vcf, s4_idx ->
                [ meta, vcf, tbi, [], region, s4_vcf, s4_idx, region, map ]
            }

        SHAPEIT5_PHASERARE(ch_rare_input)
        ch_versions = ch_versions.mix(SHAPEIT5_PHASERARE.out.versions)

        ch_phased_vcf = SHAPEIT5_PHASERARE.out.phased_variant
    }

    emit:
    phased_vcf = ch_phased_vcf
    versions   = ch_versions
}
