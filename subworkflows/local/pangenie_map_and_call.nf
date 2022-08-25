//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { PANGENIE                             } from '../../modules/local/pangenie/pangenie/main'
include { FILTER_PANGENIE_VARIANTS             } from '../../modules/local/pangenie/filter_pangenie_variants/main'
include { BCFTOOLS_VIEW as CONVERT_TO_BCF_GZ   } from '../../modules/nf-core/modules/bcftools/view/main'
include { BCFTOOLS_VIEW as SELECT_VARIANTS     } from '../../modules/nf-core/modules/bcftools/view/main'
include { BCFTOOLS_MERGE as MERGE_BCF          } from '../../modules/nf-core/modules/bcftools/merge/main'
include { BCFTOOLS_SORT as SORT_BCF            } from '../../modules/nf-core/modules/bcftools/sort/main'
include { TABIX_TABIX as INDEX_PANGENIE_OUTBCF } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as INDEX_MERGED_BCF      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as INDEX_SORTED_BCF      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { CAT_CAT as UNZIP                     } from '../../modules/nf-core/modules/cat/cat/main'

workflow PANGENIE_CALL_STRUCTURAL_VARIANTS {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ref_fasta    // channel: [mandatory] fasta
        panel_vcf    // channel: [mandatory] panel_vcf
        biallelic_panel_vcf // channel [mandatory] biallelic panel vcf
        model_file  // channel: [mandatory] pickled files of models used to regress

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // First map reads to graph using pangenie, converting output to bcf and indexing
    // UNZIP(ch_reads)
    // ***TODO***: How to best unzip these reads for PANGENIE?
    // Worse, I think that I need to combine them into one file while doing it
    PANGENIE(ch_reads, ref_fasta, panel_vcf)
    CONVERT_TO_BCF_GZ(PANGENIE.out.vcf.map {[it[0], it[1], []]}, [], [], [], [])
    ch_pangenie_raw_bcf = CONVERT_TO_BCF_GZ.out.vcf

    INDEX_PANGENIE_OUTBCF(ch_pangenie_raw_bcf)
    ch_pangenie_out_index = INDEX_PANGENIE_OUTBCF.out.tbi

    // Collect all mapped bcfs
    ch_to_merge = ch_pangenie_raw_bcf.join(ch_pangenie_out_index)

    MERGE_BCF( ch_to_merge, [], [], [] )  // to make this work, MERGE_BCF will have to have input path('*.bcf.gz'), and command will have '*.bcf.gz' as input command
    merged_sv_calls = MERGE_BCF.out.merged_variants.map{[[id:"all_pangenie_SVs"], it[1]]}
    INDEX_MERGED_BCF(merged_sv_calls)

    merged_indexed_sv_calls = merged_sv_calls.join(INDEX_MERGED_BCF.out.tbi)

    // Apply post-processing regression-based filtering of variant calls
    FILTER_PANGENIE_VARIANTS ( merged_sv_calls, biallelic_panel_vcf, model_file )
    SELECT_VARIANTS ( merged_indexed_sv_calls, FILTER_PANGENIE_VARIANTS.out.csv, [], [], [] )

    // Sort and index selected variants for output, publishing dataframe of results
    SORT_BCF( SELECT_VARIANTS.out.vcf )
    ch_bcf_out = SORT_BCF.out.vcf

    // INDEX_SORTED_BCF(sorted_bcf)
    // sorted_bcf_index = INDEX_SORTED_BCF.out.tbi

    // ch_bcf_out = sorted_bcf.join(sorted_bcf_index)


    ch_versions = ch_versions.mix(PANGENIE.out.versions)
    ch_versions = ch_versions.mix(FILTER_PANGENIE_VARIANTS.out.versions)
    ch_versions = ch_versions.mix(SORT_BCF.out.versions)
    
    ch_logs = ch_logs.mix(PANGENIE.out.log)
    ch_logs = ch_logs.mix(FILTER_PANGENIE_VARIANTS.out.log)

    emit:
        pangenie_sv_bcf   = ch_bcf_out    // channel: [ [meta], bam, bai ]
        versions          = ch_versions   // channel: [ versions.yml ]
        logs              = ch_logs
}


// # make version without snps (why?)
// bcftools view --exclude-types snps -o {sample}_multi_nosnps.bcf.gz -O b pangenie-{sample}_genotyping.vcf
// tabix output_multi_nosnps.bcf.gz
// bcftools view --exclude-types snps -o {sample}_bi_nosnps.bcf.gz -O b pangenie-{sample}_genotyping.vcf
// tabix output_bi_nosnps.bcf.gz
// bcftools view --include ID==filtered_pangenie_calls.csv $combined_var_calls -O b > filtered_pangenie_variants.bcf.gz

// ## Merge all the result vcfs
// out_multi_vcfs.collect()
// bcftools merge -l {outvcfs} -O b > all-samples_multi_all.bcf.gz
// out_bi_vcfs.collect()
// bcftools merge -l {outvcfs} -O b > all-samples_bi_all.bcf.gz



// python postprocess_vcf.py -p panel.vcf -g all-samples_bi_all.bcf.gz -o filtered_variants.parquet

// bcftools view --include ID==filtered_variants.csv all-samples_bi_all.bcf.gz -O b > filtered_pangenie_variants.bcf.gz
// output:
// filtered_variants.parquet, filtered_pangenie_variants.bcf.gz, all-samples_multi_all.bcf.gz, panel.vcf, logs #(can reuse if pangenie settings are the same)