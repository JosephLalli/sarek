//
// PANGENIE GENOTYPING
//
// Genotype structural variants from pangenome graph using PanGenie

include { PANGENIE       } from '../../../modules/local/pangenie/main'
include { PANGENIE_INDEX } from '../../../modules/local/pangenie/index/main'
include { JELLYFISH_COUNT } from '../../../modules/local/jellyfish/count/main'

workflow PANGENIE_GENOTYPING {
    take:
    ch_reads           // channel: [mandatory] [ meta, reads ] - FASTQ or k-mer counts
    ch_reference       // channel: [mandatory] [ meta, reference ] - reference FASTA
    ch_panel_vcf       // channel: [mandatory] [ meta, panel_vcf ] - panel VCF with variants

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    // Create PanGenie index (runs once per panel VCF)
    PANGENIE_INDEX(ch_panel_vcf, ch_reference)
    ch_versions = ch_versions.mix(PANGENIE_INDEX.out.versions)

    // Optional JELLYFISH_COUNT if we have FASTQ reads and want to pre-count kmers
    // PanGenie can also count them internally, but pre-counting allows using --if with segments FASTA.
    // We branch based on whether input reads are already .jf files or not.
    ch_reads.branch {
        jf: it[1].any { it.name.endsWith('.jf') }
        fastq: true
    }.set { ch_reads_branched }

    // Run Jellyfish on FASTQ reads
    // We use the segments FASTA from PANGENIE_INDEX as the --if input for k-mer counting
    JELLYFISH_COUNT(
        ch_reads_branched.fastq,
        PANGENIE_INDEX.out.index.map { it[2] }.first() // The *.fasta file from index
    )
    ch_versions = ch_versions.mix(JELLYFISH_COUNT.out.versions.first())

    // Combine pre-counted .jf files and existing .jf files
    ch_ready_reads = ch_reads_branched.jf.mix(JELLYFISH_COUNT.out.kmer_file)

    // Genotype variants using PanGenie with index
    PANGENIE(
        ch_ready_reads,
        ch_reference,
        ch_panel_vcf,
        PANGENIE_INDEX.out.index
    )

    ch_versions = ch_versions.mix(PANGENIE.out.versions.first())
    ch_reports = ch_reports.mix(PANGENIE.out.log)

    emit:
    vcf      = PANGENIE.out.vcf          // channel: [ meta, vcf.gz ]
    vcf_tbi  = PANGENIE.out.vcf_tbi      // channel: [ meta, vcf.gz.tbi ]
    reports  = ch_reports                 // channel: [ meta, log ]
    versions = ch_versions                // channel: [ versions.yml ]
}
