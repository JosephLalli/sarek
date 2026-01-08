//
// Unified STR Analysis Subworkflow
//

include { EXPANSIONHUNTER      } from '../../../modules/nf-core/expansionhunter/main'
include { STR_ANALYSIS_STRLING } from '../str_analysis_strling/main'
include { GANGSTR              } from '../../../modules/nf-core/gangstr/main'

workflow STR_ANALYSIS {
    take:
    ch_bam                      // channel: [ val(meta), bam, bai ]
    ch_fasta                    // channel: [ val(meta), fasta ]
    ch_fasta_fai                // channel: [ val(meta), fasta_fai ]
    ch_expansionhunter_catalog  // channel: [ val(meta), catalog ]
    ch_gangstr_catalog          // channel: [ val(meta), catalog ]
    ch_str_index                // channel: [ val(meta), str_index ]
    val_tools                   // list: ['expansionhunter', 'strling', 'gangstr']
    val_joint_strling           // boolean: true/false for STRling joint calling

    main:
    ch_versions = channel.empty()
    
    ch_vcf_eh       = channel.empty()
    ch_vcf_strling  = channel.empty()
    ch_vcf_gangstr  = channel.empty()

    // ExpansionHunter
    if (val_tools.contains('expansionhunter')) {
        ch_bam_eh = ch_bam.map { meta, bam, bai ->
            def new_meta = meta.clone()
            if (meta.sex == 'XX') new_meta.sex = 'female'
            else if (meta.sex == 'XY') new_meta.sex = 'male'
            else new_meta.sex = null // Remove invalid sex so no --sex flag is passed
            [ new_meta, bam, bai ]
        }
        EXPANSIONHUNTER(
            ch_bam_eh,
            ch_fasta.first(),
            ch_fasta_fai.first(),
            ch_expansionhunter_catalog.first()
        )
        ch_vcf_eh = EXPANSIONHUNTER.out.vcf
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)
    }

    // STRling
    if (val_tools.contains('strling')) {
        STR_ANALYSIS_STRLING(
            ch_bam,
            ch_fasta,
            ch_fasta_fai,
            ch_str_index,
            val_joint_strling
        )
        ch_vcf_strling = STR_ANALYSIS_STRLING.out.vcf
        ch_versions = ch_versions.mix(STR_ANALYSIS_STRLING.out.versions)
    }

    // GangSTR
    if (val_tools.contains('gangstr')) {
        // GANGSTR input: tuple val(meta), path(alignment_files), path(alignment_indices), path(ref_regions)
        // Combine ch_bam with ch_gangstr_catalog
        ch_gangstr_input = ch_bam.combine(ch_gangstr_catalog.first())
            .map { meta, bam, bai, cat_meta, catalog ->
                [ meta, [bam], [bai], catalog ]
            }
        
        GANGSTR(
            ch_gangstr_input,
            ch_fasta.map{ it[1] }.first(),
            ch_fasta_fai.map{ it[1] }.first()
        )
        ch_vcf_gangstr = GANGSTR.out.vcf
        ch_versions = ch_versions.mix(GANGSTR.out.versions)
    }

    emit:
    vcf_eh      = ch_vcf_eh
    vcf_strling = ch_vcf_strling
    vcf_gangstr = ch_vcf_gangstr
    versions    = ch_versions
}