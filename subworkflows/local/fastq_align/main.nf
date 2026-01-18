//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'
include { SENTIEON_BWAMEM        } from '../../../modules/nf-core/sentieon/bwamem/main'
include { VG_ALIGN               } from '../vg_align/main'

workflow FASTQ_ALIGN {
    take:
    reads // channel: [mandatory] meta, reads
    index // channel: [mandatory] index
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
    fasta
    fasta_fai
    dict
    run_markdup
    personalized
    save_gam
    run_vg_stats

    main:

    versions = channel.empty()
    reports = channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(reads, index, [[id:'no_fasta'], []], sort) // If aligner is dragmap
    // The sentieon-bwamem-module does sorting as part of the conversion from sam to bam.
    SENTIEON_BWAMEM(reads, index, fasta, fasta_fai) // If aligner is sentieon-bwamem

    // If aligner is giraffe (pangenome)
    ch_vg_bam = channel.empty()
    ch_vg_bai = channel.empty()
    ch_vg_cram = channel.empty()
    ch_vg_crai = channel.empty()
    ch_vg_markdup_stats = channel.empty()
    ch_vg_reports = channel.empty()
    ch_vg_versions = channel.empty()

    if (params.aligner == 'giraffe') {
        VG_ALIGN(
            reads,
            index,
            fasta,
            fasta_fai,
            dict,
            sort,
            run_markdup,
            personalized,
            save_gam,
            run_vg_stats
        )
        ch_vg_bam = VG_ALIGN.out.bam
        ch_vg_bai = VG_ALIGN.out.bai
        ch_vg_cram = VG_ALIGN.out.cram
        ch_vg_crai = VG_ALIGN.out.crai
        ch_vg_markdup_stats = VG_ALIGN.out.markdup_stats
        ch_vg_reports = VG_ALIGN.out.reports
        ch_vg_versions = VG_ALIGN.out.versions
    }

    // Get the bam files from the aligner
    // Only one aligner is run
    bam = channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)
    bam = bam.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, _bam, _bai -> [ meta, _bam ] })
    bam = bam.mix(ch_vg_bam)

    bai = channel.empty()
    bai = bai.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, _bam, _bai -> [ meta, _bai ] })
    bai = bai.mix(ch_vg_bai)

    cram = ch_vg_cram
    crai = ch_vg_crai
    markdup_stats = ch_vg_markdup_stats

    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)
    reports = reports.mix(ch_vg_reports)

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    versions = versions.mix(SENTIEON_BWAMEM.out.versions)
    versions = versions.mix(ch_vg_versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    cram     // channel: [ [meta], cram ]
    crai     // channel: [ [meta], crai ]
    markdup_stats // channel: [ [meta], markdup_stats ]
    reports
    versions // channel: [ versions.yml ]
}
