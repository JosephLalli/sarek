# Legacy Feature Gap Analysis Report

## Unique Modules/Subworkflows in `sarek_first_attempt` (not in `sarek_pangenome`)

### `modules/local` (19 unique files)
*   `make_personalized_gbz.nf`
*   `make_personalized_gfa.nf`
*   `add_paths_to_gfa.nf`
*   `split_by_contig.nf`
*   `split_sample_from_vcf.nf`
*   `split_PARs_off_X.nf`
*   `split_ref_fasta.nf`
*   `split_vcf_by_contig.nf`
*   `unify_insertion_fastas.nf`
*   `abra/abra2_and_markdup.nf`
*   `abra/main.nf`
*   `jellyfish/count/main.nf`
*   `whatshap/stats/main.nf`
*   `whatshap/phase/main.nf`
*   `bedtools/slop/main.nf`
*   `shapeit4/phase_common/main.nf`
*   `build_intervals_from_bed/main.nf`
*   `d4tools/slop/main.nf`
*   `calc_stats/per_base/main.nf`
*   `calc_stats/per_read/main.nf`
*   `build_intervals/main.nf`
*   `gatk/realignertargetcreator/main.nf`
*   `generate_phasing_intervals.nf`
*   `samtools/reheader/main.nf`
*   `samtools/view/main.nf`
*   `samtools/cat/main.nf`
*   `samtools/markdup/main.nf`
*   `samtools/dict/main.nf`
*   `samtools/sort_and_leftalign.nf`
*   `samtools/merge_dedup/main.nf`
*   `samtools/postprocess_1.nf`
*   `illumina/hap.py/main.nf`
*   `danbing-tk/predict/main.nf`
*   `danbing-tk/align/main.nf`
*   `bcftools/split_and_merge/main.nf`
*   `bcftools/merge_and_filter/main.nf`
*   `bcftools/view/main.nf`
*   `bcftools/sample_name_to_patient_name/main.nf`
*   `bcftools/trio_switch_error/main.nf`
*   `bcftools/mendelian/main.nf`
*   `bcftools/sort_and_filter/main.nf`
*   `bcftools/add_variants_to_phased_common/main.nf`
*   `bcftools/calc_ref_bias/main.nf`
*   `bcftools/merge_and_annotate/main.nf`
*   `freebayes/bamleftalign/main.nf`
*   `freebayes/left_align_and_sort/main.nf`
*   `deepvariant/convert_haploid_regions/main.nf`
*   `prep_for_gatk_metrics.nf`
*   `kmc/kmc/main.nf`
*   `kmc/kmc_dump/main.nf`
*   `vg/convert/main.nf`
*   `vg/surject/main.nf`
*   `vg/giraffe_to_gam/main.nf`
*   `vg/stats/main.nf`
*   `vg/vg_giraffe_and_sort/main.nf`
*   `vg/pack/main.nf`
*   `vg/deconstruct/deconstruct_and_filter.nf`
*   `vg/call/main.nf`
*   `vg/snarls/main.nf`
*   `vg/paths/main.nf`
*   `vg/autoindex/main.nf`
*   `vg/giraffe/main.nf`
*   `vg/index/main.nf`
*   `vg/haplotypes/main.nf`
*   `vg/gbwt/main.nf`
*   `vg/minimizer/main.nf`
*   `pangenie/pangenie/main.nf`
*   `pangenie/calc_mendelian_violations/main.nf`
*   `pangenie/filter_pangenie_variants/main.nf`
*   `expansionhunter/main.nf`
*   `shapeit5/merge_switch_reports/main.nf`
*   `shapeit5/phase_rare/main.nf`
*   `shapeit5/phase_common/main.nf`
*   `shapeit5/switch/main.nf`
*   `shapeit5/ligate/main.nf`


### `subworkflows/local` (18 unique files)
*   `phase_and_impute_old.nf`
*   `phase_and_impute_panel_merging.nf`
*   `phase_and_impute_read_informed.nf`
*   `phase_and_impute_read_informed_panel_merging.nf`
*   `phase_and_impute.nf`
*   `concat_vcfs_by_contig.nf`
*   `mapping_csv.nf`
*   `prepare_recalibration_csv.nf`
*   `giraffe_mapping.nf`
*   `pair_variant_calling.nf`
*   `germline_variant_calling.nf`
*   `recalibrate_csv.nf`
*   `pangenie_map_and_call.nf`
*   `markduplicates_csv.nf`
*   `annotate.nf`
*   `split_out_haploid_vcfs.nf`
*   `tumor_variant_calling.nf`
*   `realign_indels.nf`
*   `variantcalling_csv.nf`
*   `merge_vcfs.nf`

### `modules/nf-core` (Modified - Same relative path, different hash)
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/snpeff/main.nf` (current `nf-core/snpeff/snpeff/main.nf`)
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/ensemblvep/main.nf` (current `nf-core/ensemblvep/vep/main.nf`)
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/deepvariant/main.nf` (current `nf-core/deepvariant/rundeepvariant/main.nf`)
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/mosdepth/main.nf`
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/ascat/main.nf`
*   `/mnt/ssd/lalli/nf_stage/sarek_pangenome/sarek_first_attempt/modules/nf-core/modules/parabricks/fqtobam/main.nf`
*   Many `gatk4` modules
*   Many `samtools` modules
*   Many `bcftools` modules
*   `msisensorpro/msi_somatic/main.nf`
*   `msisensorpro/scan/main.nf`
*   `cnvkit/*` modules
*   `strelka/somatic/main.nf`
*   `strelka/germline/main.nf`
*   `tabix/bgziptabix/main.nf`
*   `tabix/tabix/main.nf`
*   `bwa/index/main.nf`
*   `bwa/mem/main.nf`
*   `bwamem2/index/main.nf`
*   `bwamem2/mem/main.nf`
*   `freebayes/main.nf`
*   `manta/somatic/main.nf`
*   `manta/germline/main.nf`
*   `manta/tumoronly/main.nf`
*   `tiddit/sv/main.nf`
*   `vcftools/main.nf`
*   `controlfreec/*` modules
*   `svdb/merge/main.nf`
*   `glnexus/main.nf`
*   `fgbio/*` modules
*   `fastp/main.nf`
*   `picard/*` modules
*   `samblaster/main.nf`
*   `cat/*` modules
*   `dragmap/*` modules
*   `biobambam/main.nf`
*   `unzip/main.nf`
*   `untar/main.nf`
*   `fastqc/main.nf`
*   `custom/dumpsoftwareversions/main.nf`
*   `bedtools/maskfasta/main.nf`


### `subworkflows/nf-core` (Modified - Same relative path, different hash)
*   `alignment_to_fastq.nf`
*   `gatk4/tumor_normal_somatic_variant_calling/main.nf`
*   `gatk4/prepare_recalibration_spark/main.nf`
*   `gatk4/markduplicates/main.nf`
*   `gatk4/tumor_only_somatic_variant_calling/main.nf`
*   `gatk4/joint_germline_variant_calling/main.nf`
*   `gatk4/markduplicates_spark/main.nf`
*   `gatk4/mapping/main.nf`
*   `gatk4/recalibrate_spark/main.nf`
*   `gatk4/recalibrate/main.nf`
*   `gatk4/prepare_recalibration/main.nf`
*   `gatk4/single_sample_germline_variant_calling/main.nf`
*   `run_fastqc.nf`
*   `merge_index_cram.nf`
*   `fgbio_create_umi_consensus/main.nf`
*   `cram_qc.nf`
*   `bam_to_cram.nf`
*   `variantcalling/ascat/main.nf`
*   `variantcalling/mpileup/main.nf`
*   `variantcalling/deepvariant/main.nf`
*   `variantcalling/manta/somatic/main.nf`
*   `variantcalling/manta/germline/main.nf`
*   `variantcalling/manta/tumoronly/main.nf`
*   `variantcalling/tiddit/somatic/main.nf`
*   `variantcalling/tiddit/single/main.nf`
*   `variantcalling/haplotypecaller/main.nf`
*   `variantcalling/controlfreec/somatic/main.nf`
*   `variantcalling/controlfreec/tumoronly/main.nf`
*   `variantcalling/freebayes/main.nf`
*   `variantcalling/strelka/somatic/main.nf`
*   `variantcalling/strelka/single/main.nf`
*   `variantcalling/cnvkit/main.nf`
*   `merge_index_bam.nf`
*   `annotation/snpeff/main.nf`
*   `annotation/ensemblvep/main.nf`
*   `vcf_qc.nf`

## Files from `sarek_JLL`
*   The directory `/mnt/ssd/lalli/nf_stage/sarek_JLL` contains primarily reference data (`local_references/1K_fastqs/`, `annotations/`, `deepvariant_models/`, `tandem_repeats/`) and is not structured as a Nextflow pipeline for direct module/subworkflow comparison.

---

**Next Steps for User:** I will present this report in `gap_analysis_report.md` and then initiate the interactive review process. I will ask you to confirm if there are any specific files or categories you would like to investigate further to discern the logical changes, or if you are ready to proceed with reviewing this high-level summary.
