//
// PANGENIE GENOTYPING
//
// Genotype structural variants from pangenome graph using PanGenie

include { PANGENIE } from '../../../modules/local/pangenie/main'

workflow PANGENIE_GENOTYPING {
    take:
    ch_reads           // channel: [mandatory] [ meta, reads ] - FASTQ or k-mer counts
    ch_reference       // channel: [mandatory] [ meta, reference ] - reference FASTA
    ch_panel_vcf       // channel: [mandatory] [ meta, panel_vcf ] - panel VCF with variants

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    // Genotype variants using PanGenie
    PANGENIE(
        ch_reads,
        ch_reference,
        ch_panel_vcf
    )

    ch_versions = ch_versions.mix(PANGENIE.out.versions.first())
    ch_reports = ch_reports.mix(PANGENIE.out.log)

    emit:
    vcf      = PANGENIE.out.vcf          // channel: [ meta, vcf.gz ]
    vcf_tbi  = PANGENIE.out.vcf_tbi      // channel: [ meta, vcf.gz.tbi ]
    reports  = ch_reports                 // channel: [ meta, log ]
    versions = ch_versions                // channel: [ versions.yml ]
}
