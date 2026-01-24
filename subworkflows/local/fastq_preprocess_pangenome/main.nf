// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                          } from '../../../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV                 } from '../../../subworkflows/local/channel_markduplicates_create_csv/main'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI         } from '../../../subworkflows/local/bam_convert_samtools/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                             } from '../../../modules/nf-core/fastp/main'

// remove genomic contaminants with bbsplit
include { BBMAP_BBSPLIT                                     } from '../../../modules/nf-core/bbmap/bbsplit'

// Create umi consensus bams from fastq
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO                  } from '../../../subworkflows/local/fastq_create_umi_consensus_fgbio/main'

// Map input reads to reference genome - PANGENOME SPECIFIC
include { VG_ALIGN                                          } from '../../../subworkflows/local/vg_align/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                          } from '../../../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING           } from '../../../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                   } from '../../../modules/nf-core/samtools/convert/main'

// Copy UMIs from read name to RX tag
include { FGBIO_COPYUMIFROMREADNAME                         } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                                } from '../../../subworkflows/local/bam_markduplicates/main'
include { BAM_MARKDUPLICATES_SPARK                          } from '../../../subworkflows/local/bam_markduplicates_spark/main'
include { BAM_SENTIEON_DEDUP                                } from '../../../subworkflows/local/bam_sentieon_dedup/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD        } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

workflow FASTQ_PREPROCESS_PANGENOME {
    take:
        input_fastq
        input_sample
        dict
        fasta
        fasta_fai
        index_alignment
        intervals_and_num_intervals
        intervals_for_preprocessing
        known_sites_indels
        known_sites_indels_tbi
        bbsplit_index

    main:

    // To gather all QC reports for MultiQC
    reports          = channel.empty()
    versions         = channel.empty()
    ch_cram_variant_calling = channel.empty()

    // Configuration Logic for MarkDuplicates
    def tools_list = params.tools ? params.tools.split(',') : []
    def skip_list = params.skip_tools ? params.skip_tools.split(',') : []
    
    def use_samtools_markdup = tools_list.contains('samtools_markdup')
    def skip_markduplicates = skip_list.contains('markduplicates')
    
    if (use_samtools_markdup && !skip_markduplicates) {
        error("Invalid configuration: 'samtools_markdup' cannot be in 'tools' if 'markduplicates' is not in 'skip_tools'.")
    }

    def run_markdup_in_vg = use_samtools_markdup && skip_markduplicates

    // PREPROCESSING

    if (params.step == 'mapping') {

        // STEP 0: QC & TRIM
        if (params.umi_read_structure) {
            FASTQ_CREATE_UMI_CONSENSUS_FGBIO(input_fastq, fasta, fasta_fai, dict, index_alignment, params.group_by_umi_strategy)
            bam_converted_from_fastq = FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.consensusbam.map{ meta, bam -> [ meta, bam, [] ] }
            CONVERT_FASTQ_UMI(bam_converted_from_fastq, [ [ id:"fasta" ], [] ], [ [ id:'null' ], [] ], false)
            reads_for_fastp = CONVERT_FASTQ_UMI.out.reads
            versions = versions.mix(CONVERT_FASTQ_UMI.out.versions, FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.versions)
        } else {
            reads_for_fastp = input_fastq
        }

        if (params.trim_fastq || params.split_fastq > 0 || params.umi_location) {
            FASTP(reads_for_fastp, [], false, false, false)
            reports = reports.mix(FASTP.out.json.collect{ it[1] }, FASTP.out.html.collect{ it[1] })
            reads_for_bbsplit = params.split_fastq ? FASTP.out.reads.map{ m, r -> [m + [n_fastq: r.size()/2], r] }.transpose() : FASTP.out.reads
            versions = versions.mix(FASTP.out.versions)
        } else {
            reads_for_bbsplit = reads_for_fastp
        }

        if (tools_list.contains('bbsplit')) {
            BBMAP_BBSPLIT (reads_for_bbsplit, bbsplit_index, [], [ [], [] ], false)
            reads_for_alignment = BBMAP_BBSPLIT.out.primary_fastq
            reports = reports.mix(BBMAP_BBSPLIT.out.stats.collect{ it[1] })
            versions = versions.mix(BBMAP_BBSPLIT.out.versions.first())
        } else {
            reads_for_alignment = reads_for_bbsplit
        }

        // STEP 1: MAPPING
        reads_for_alignment.map { meta, reads -> [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ] }
            .groupTuple().map { meta, reads -> meta + [ n_fastq: reads.size() ] }.set { reads_grouping_key }

        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        VG_ALIGN(reads_for_alignment, index_alignment, fasta, fasta_fai, dict, true, run_markdup_in_vg, params.pangenome_personalized_flow, params.save_mapped, true)
        versions = versions.mix(VG_ALIGN.out.versions)
        reports = reports.mix(VG_ALIGN.out.reports)

        // MIX outputs
        aligned_bam = VG_ALIGN.out.bam.mix(VG_ALIGN.out.cram)
        aligned_bai = VG_ALIGN.out.bai.mix(VG_ALIGN.out.crai)

        if (params.umi_in_read_header || params.umi_location) {
            FGBIO_COPYUMIFROMREADNAME(aligned_bam.map{m, f -> [m, f, []]})
            aligned_bam = FGBIO_COPYUMIFROMREADNAME.out.bam
            aligned_bai = FGBIO_COPYUMIFROMREADNAME.out.bai
            versions = versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)
        }

        bam_mapped = aligned_bam.join(aligned_bai)
            .combine(reads_grouping_key)
            .filter { m1, f, i, m2 -> m1.sample == m2.sample }
            .map { m1, f, i, m2 -> [ m1 + m2 - m1.subMap('id', 'read_group', 'data_type', 'num_lanes', 'size', 'sample_lane_id') + [ id: m1.sample ], f, i ] }
            .map { m, f, i -> [ groupKey( m, m.n_fastq), f, i ] }
            .groupTuple()

        if (params.save_mapped || (skip_markduplicates && !use_samtools_markdup)) {
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped.map{ m, f, i -> [m, f] })
            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions, BAM_TO_CRAM_MAPPING.out.versions)
        }
    }

    // STEP 2: MARKDUPLICATES
    if (params.step in ['mapping', 'markduplicates']) {
        if (run_markdup_in_vg) {
            ch_cram_variant_calling = bam_mapped.map { m, f, i -> [ m + [data_type: f[0].name.endsWith('.cram') ? 'cram' : 'bam'], f[0], i[0] ] }
        } else {
            cram_for_md = params.step == 'mapping' ? bam_mapped.map{ m, f, i -> [m, f] } : input_sample.map{ m, f, i -> [m, f] }
            BAM_MARKDUPLICATES(cram_for_md, fasta, fasta_fai, intervals_for_preprocessing)
            ch_cram_variant_calling = BAM_MARKDUPLICATES.out.cram
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ it[1] })
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }
    }
    
    emit:
    cram_variant_calling = ch_cram_variant_calling
    reports
    versions
}
