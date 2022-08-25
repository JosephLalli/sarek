//
// REALIGN INDELS
//

include { FREEBAYES_BAM_LEFT_ALIGN                   } from '../../modules/local/freebayes/bamleftalign/main'
include { GATK_REALIGNER_TARGET_CREATOR              } from '../../modules/local/gatk/realignertargetcreator/main'
include { BEDTOOLS_SLOP                              } from '../../modules/local/bedtools/slop/main'
include { ABRA2                                      } from '../../modules/local/abra/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_LEFTALIGN } from '../../modules/nf-core/modules/samtools/index/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ABRA2     } from '../../modules/nf-core/modules/samtools/index/main'

workflow REALIGN_INDELS {
    take:
    bam                   // channel: [ val(meta), bam ]
    ref_fasta
    ref_fasta_fai
    ref_fasta_dict        
    ref_contig_lengths
    
    main:
    ch_logs   = Channel.empty()
    ch_versions  = Channel.empty()


    FREEBAYES_BAM_LEFT_ALIGN(bam, ref_fasta, ref_fasta_fai)
    ch_bam = FREEBAYES_BAM_LEFT_ALIGN.out.bam

    bai = SAMTOOLS_INDEX_LEFTALIGN(bam)
    bam_bai = bam.join(SAMTOOLS_INDEX_LEFTALIGN.out.bai)

    GATK_REALIGNER_TARGET_CREATOR(bam_bai, ref_fasta, ref_fasta_fai, ref_fasta_dict)
    indel_target_bed = GATK_REALIGNER_TARGET_CREATOR.out.indel_target_bed
    
    BEDTOOLS_SLOP(indel_target_bed, ref_contig_lengths)
    widened_indel_target_bed = BEDTOOLS_SLOP.out.widened_indel_target_bed

    // bam_bai.view()
    // ref_fasta.view()
    // widened_indel_target_bed.view()
    ABRA2(bam_bai, ref_fasta, widened_indel_target_bed.map { it[1] })
    ch_bam_out = ABRA2.out.realigned_bam
    // ch_bam_out = ch_bam_out.join(SAMTOOLS_INDEX_ABRA2(ch_bam_out).out.bai)
    

    ch_versions = ch_versions.mix(FREEBAYES_BAM_LEFT_ALIGN.out.versions)
    ch_versions = ch_versions.mix(GATK_REALIGNER_TARGET_CREATOR.out.versions)
    ch_versions = ch_versions.mix(BEDTOOLS_SLOP.out.versions)
    ch_versions = ch_versions.mix(ABRA2.out.versions)
    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX_ABRA2.out.versions)

    ch_logs  = ch_logs.mix(ABRA2.out.log)

    emit:
    bam       = ch_bam_out   //    channel: [ val(meta), bam, bai ]
    reports   = ch_logs      //    path: *.log
    versions  = ch_versions  //    path: versions.yml
}
