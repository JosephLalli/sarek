//
// PREPROCESSING WITH GIRAFFE
//
// Pangenome-aware alignment using vg giraffe, then standard GATK preprocessing

// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV        } from '../../../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV } from '../../../subworkflows/local/channel_markduplicates_create_csv/main'

// Pangenome alignment
include { GIRAFFE_MAPPING                 } from '../../../subworkflows/local/giraffe_mapping/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS        } from '../../../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING } from '../../../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES              } from '../../../subworkflows/local/bam_markduplicates/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

workflow FASTQ_PREPROCESS_GIRAFFE {
    take:
    ch_reads                 // channel: [mandatory] meta, reads
    ch_input_sample          // channel: input samples for restart
    ch_fasta                 // channel: [mandatory] meta, fasta
    ch_fasta_fai             // channel: [mandatory] meta, fasta_fai
    ch_dict                  // channel: [mandatory] meta, dict
    ch_gbz                   // channel: [mandatory] meta, gbz
    ch_dist                  // channel: [mandatory] meta, dist
    ch_min                   // channel: [mandatory] meta, min
    ch_ref_paths             // channel: [optional]  ref_paths.txt
    ch_intervals_for_preprocessing // channel: intervals

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    if (params.step == 'mapping') {

        // Calculate number of lanes for each sample for grouping
        ch_reads.map { meta, reads ->
                [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]
            }
            .groupTuple()
            .map { meta, reads ->
                meta + [ n_fastq: reads.size() ]
            }
            .set { reads_grouping_key }

        ch_reads_for_alignment = ch_reads.map { meta, reads ->
            // Update meta.id to meta.sample when single lane
            if (meta.size * meta.num_lanes == 1) [ meta + [ id: meta.sample ], reads ]
            else [ meta, reads ]
        }

        // STEP 1: GIRAFFE ALIGNMENT
        // sort_bam=true, run_fixmate=true, run_markdup=false (markdup handled separately)
        GIRAFFE_MAPPING(
            ch_reads_for_alignment,
            ch_gbz,
            ch_dist,
            ch_min,
            ch_ref_paths,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            true,   // sort_bam
            true,   // run_fixmate
            false   // run_markdup - handled by BAM_MARKDUPLICATES later
        )

        ch_versions = ch_versions.mix(GIRAFFE_MAPPING.out.versions)
        ch_reports = ch_reports.mix(GIRAFFE_MAPPING.out.reports)

        // Group bams from the same sample
        bam_mapped = GIRAFFE_MAPPING.out.bam
            .combine(reads_grouping_key)
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }
            .map { meta1, bam, meta2 -> [ meta1 + meta2, bam ] }
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'size', 'sample_lane_id') + [ data_type: 'bam', id: meta.sample ], bam ]
            }
            .map { meta, bam -> [ groupKey(meta, meta.n_fastq), bam ] }
            .groupTuple()

        // Merge and save mapped BAMs if requested
        if (params.save_mapped || (params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) {
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)
            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, ch_fasta, ch_fasta_fai)

            if (params.save_output_as_bam) {
                CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, params.outdir, params.save_output_as_bam)
            } else {
                CHANNEL_ALIGN_CREATE_CSV(
                    BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true),
                    params.outdir,
                    params.save_output_as_bam
                )
            }

            ch_versions = ch_versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            ch_versions = ch_versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }
    }

    // STEP 2: MARKDUPLICATES
    if (params.step in ['mapping', 'markduplicates']) {

        cram_markduplicates = channel.empty()
        cram_skip_markduplicates = channel.empty()

        ch_cram_for_markduplicates = params.step == 'mapping'
            ? bam_mapped
            : ch_input_sample.map { meta, input, _index -> [ meta, input ] }

        if (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) {
            if (params.step == 'mapping') {
                cram_skip_markduplicates = BAM_TO_CRAM_MAPPING.out.cram
                    .join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true)
            } else {
                cram_skip_markduplicates = ch_input_sample
            }

            CRAM_QC_NO_MD(cram_skip_markduplicates, ch_fasta, ch_intervals_for_preprocessing)
            ch_reports = ch_reports.mix(CRAM_QC_NO_MD.out.reports.collect { _meta, report -> [ report ] })
            ch_versions = ch_versions.mix(CRAM_QC_NO_MD.out.versions)
        } else {
            BAM_MARKDUPLICATES(
                ch_cram_for_markduplicates,
                ch_fasta,
                ch_fasta_fai,
                ch_intervals_for_preprocessing
            )

            cram_markduplicates = BAM_MARKDUPLICATES.out.cram
            ch_reports = ch_reports.mix(BAM_MARKDUPLICATES.out.reports.collect { _meta, report -> [ report ] })
            ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)
        }

        cram_variant_calling = channel.empty()
            .mix(cram_markduplicates, cram_skip_markduplicates)
            .map { meta, cram, crai -> [ meta + [ data_type: 'cram' ], cram, crai ] }

        // Create CSV for restart
        CHANNEL_MARKDUPLICATES_CREATE_CSV(
            cram_variant_calling,
            'markduplicates',
            params.outdir,
            params.save_output_as_bam
        )
    }

    emit:
    cram_variant_calling = cram_variant_calling  // channel: [ meta, cram, crai ]
    reports              = ch_reports            // channel: [ reports ]
    versions             = ch_versions           // channel: [ versions.yml ]
}
