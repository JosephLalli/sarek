//
// Unified STR Analysis Subworkflow
//

include { EXPANSIONHUNTER      } from '../../../modules/nf-core/expansionhunter/main'
include { STRLING_TO_VCF      } from '../../../modules/local/strling/to_vcf/main'
include { STR_ANALYSIS_STRLING } from '../str_analysis_strling/main'
include { GANGSTR              } from '../../../modules/nf-core/gangstr/main'
include { STR_CONSENSUS_MERGE  } from '../../../modules/local/str/consensus_merge/main'

workflow STR_ANALYSIS {
    take:
    ch_alignment                // channel: [ val(meta), bam_cram, bai_crai ]
    ch_fasta                    // channel: [ val(meta), fasta ]
    ch_fasta_fai                // channel: [ val(meta), fasta_fai ]
    ch_expansionhunter_catalog  // channel: [ val(meta), catalog ]
    ch_gangstr_catalog          // channel: [ val(meta), catalog ]
    ch_strling_catalog          // channel: [ val(meta), catalog ]
    ch_str_index                // channel: [ val(meta), str_index ]
    ch_strling_loci             // channel: [ val(meta), strling_loci ]
    val_tools                   // list: ['expansionhunter', 'strling', 'gangstr']
    val_joint_strling           // boolean: true/false for STRling joint calling

    main:
    ch_versions = channel.empty()
    
    ch_vcf_eh           = channel.empty()
    ch_vcf_gangstr      = channel.empty()
    ch_vcf_strling      = channel.empty()
    ch_vcf_consensus    = channel.empty()
    ch_strling_genotype = channel.empty()
    ch_strling_bounds   = channel.empty()

    // ExpansionHunter
    if (val_tools.contains('expansionhunter')) {
        ch_bam_eh = ch_alignment.map { meta, bam, bai ->
            def new_meta = meta.clone()
            if (meta.sex == 'XX') new_meta.sex = 'female'
            else if (meta.sex == 'XY') new_meta.sex = 'male'
            else new_meta.sex = null // Remove invalid sex so no --sex flag is passed
            [ new_meta, bam, bai ]
        }
        EXPANSIONHUNTER(
            ch_bam_eh,
            ch_fasta.collect(),
            ch_fasta_fai.collect(),
            ch_expansionhunter_catalog.collect()
        )
        ch_vcf_eh = EXPANSIONHUNTER.out.vcf
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)
    }

    // STRling
    if (val_tools.contains('strling')) {
        STR_ANALYSIS_STRLING(
            ch_alignment,
            ch_fasta,
            ch_fasta_fai,
            ch_str_index,
            ch_strling_loci,
            val_joint_strling
        )
        ch_strling_genotype = STR_ANALYSIS_STRLING.out.genotype
        ch_strling_bounds   = STR_ANALYSIS_STRLING.out.bounds
        ch_versions = ch_versions.mix(STR_ANALYSIS_STRLING.out.versions)

        STRLING_TO_VCF(
            ch_strling_genotype,
            ch_strling_catalog.map { it[1] }.collect(),
            ch_fasta_fai.map { it[1] }.collect()
        )
        ch_vcf_strling = STRLING_TO_VCF.out.vcf
        ch_versions = ch_versions.mix(STRLING_TO_VCF.out.versions)
    }

    // GangSTR
    if (val_tools.contains('gangstr')) {
        // GANGSTR input: tuple val(meta), path(alignment_files), path(alignment_indices), path(ref_regions)
        // Combine ch_alignment with ch_gangstr_catalog
        ch_gangstr_input = ch_alignment.combine(ch_gangstr_catalog.collect())
            .map { meta, bam, bai, cat_meta, catalog ->
                [ meta, [bam], [bai], catalog ]
            }
        
        GANGSTR(
            ch_gangstr_input,
            ch_fasta.map{ it[1] }.collect(),
            ch_fasta_fai.map{ it[1] }.collect()
        )
        ch_vcf_gangstr = GANGSTR.out.vcf
        ch_versions = ch_versions.mix(GANGSTR.out.versions)
    }

    if (val_tools.contains('expansionhunter') && val_tools.contains('gangstr') && val_tools.contains('strling')) {
        ch_eh_keyed = ch_vcf_eh.map { meta, vcf -> [ meta.id, meta, file(vcf, stageAs: "${meta.id}.eh.vcf.gz") ] }
        ch_gs_keyed = ch_vcf_gangstr.map { meta, vcf -> [ meta.id, meta, file(vcf, stageAs: "${meta.id}.gs.vcf.gz") ] }
        ch_sl_keyed = ch_vcf_strling.map { meta, vcf -> [ meta.id, meta, file(vcf, stageAs: "${meta.id}.sl.vcf.gz") ] }

        ch_consensus_input = ch_eh_keyed
            .combine(ch_gs_keyed, by: 0)
            .combine(ch_sl_keyed, by: 0)
            .map { _id, meta_eh, eh_vcf, meta_gs, gs_vcf, meta_sl, sl_vcf ->
                [ meta_eh, eh_vcf, gs_vcf, sl_vcf ]
            }

        STR_CONSENSUS_MERGE(ch_consensus_input)
        ch_vcf_consensus = STR_CONSENSUS_MERGE.out.vcf
        ch_versions = ch_versions.mix(STR_CONSENSUS_MERGE.out.versions)
    }

    emit:
    vcf_eh           = ch_vcf_eh
    vcf_gangstr      = ch_vcf_gangstr
    vcf_strling      = ch_vcf_strling
    vcf_consensus    = ch_vcf_consensus
    strling_genotype = ch_strling_genotype
    strling_bounds   = ch_strling_bounds
    versions         = ch_versions
}
