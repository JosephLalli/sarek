include { VG_GIRAFFE }    from '../../../modules/local/vg/giraffe/main'
include { VG_SURJECT }    from '../../../modules/local/vg/surject/main'
include { VG_STATS }      from '../../../modules/local/vg/stats/main'
include { KMC }           from '../../../modules/local/kmc/main'
include { VG_HAPLOTYPES } from '../../../modules/local/vg/haplotypes/main'

workflow VG_ALIGN {

    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_index_alignment  // channel: [ val(meta), [ gbz, dist, min, zipcodes, ref_paths, hapl ] ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fasta_fai        // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    sort_bam            // val: boolean
    run_markdup         // val: boolean
    personalized        // val: boolean
    save_gam            // val: boolean
    run_vg_stats        // val: boolean

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    
    // Convert index to a value channel to allow multiple combinations with reads
    ch_index_val = ch_index_alignment.first()

    ch_giraffe_input = Channel.empty()
    ch_personalized_gbz = Channel.empty()

    if (personalized) {
        // 1. KMC on Reads
        KMC(ch_reads)
        ch_versions = ch_versions.mix(KMC.out.versions)
        
        // 2. VG_HAPLOTYPES (Sampling Mode)
        // Input: [meta, gbz, dist, hapl, r_index]
        // Sampling only requires gbz and hapl.
        ch_reads
            .join(KMC.out.kff)
            .combine(ch_index_val)
            .multiMap { meta, reads, kff, meta_idx, index_list ->
                // index_list: [gbz, dist, min, zipcodes, ref_paths, hapl]
                index: [meta, index_list[0], [], index_list[5], []] // [meta, gbz, dist=[], hapl, r_index=[]]
                kmer:  [meta, kff]
            }
            .set { ch_haplo_inputs }
        
        VG_HAPLOTYPES(
            ch_haplo_inputs.index,
            ch_haplo_inputs.kmer
        )
        ch_versions = ch_versions.mix(VG_HAPLOTYPES.out.versions)
        
        ch_personalized_gbz = VG_HAPLOTYPES.out.personalized_gbz

        // 3. Construct Giraffe Input with Personalized GBZ
        ch_giraffe_input = ch_reads
            .join(ch_personalized_gbz) 
            .combine(ch_index_val)
            .map { meta, reads, gbz, meta_idx, index_list ->
                // index_list[3] is zipcodes
                [meta, reads, gbz, [], [], index_list[3]] // [meta, reads, gbz, dist=[], min=[], zipcodes]
            }

    } else {
        // Standard Flow
        ch_giraffe_input = ch_reads
            .combine(ch_index_val)
            .map { meta, reads, meta_idx, index_list ->
                // index_list: [gbz, dist, min, zipcodes, ref_paths, hapl]
                [meta, reads, index_list[0], index_list[1], index_list[2], index_list[3]]
            }
    }

    // 4. Align (Giraffe)
    VG_GIRAFFE (
        ch_giraffe_input,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_index_val.map { meta_idx, index_list -> index_list[4] }, // ref_paths
        sort_bam,
        run_markdup
    )
    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions)

    // 5. Conditional Surject
    ch_surject_input = Channel.empty()

    if (personalized) {
        ch_surject_input = VG_GIRAFFE.out.gam
            .join(ch_personalized_gbz)
            .map { meta, gam, gbz -> [meta, gam, gbz] }
    } else {
        ch_surject_input = VG_GIRAFFE.out.gam
            .combine(ch_index_val.map { meta_idx, index_list -> index_list[0] }) // global gbz
            .map { meta, gam, gbz -> [meta, gam, gbz] }
    }

    VG_SURJECT (
        ch_surject_input,
        ch_index_val.map { meta_idx, index_list -> index_list[4] }, // ref_paths
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        sort_bam,
        run_markdup
    )
    ch_versions = ch_versions.mix(VG_SURJECT.out.versions)

    // 6. Optional VG_STATS
    ch_reports = VG_GIRAFFE.out.stats
    if (run_vg_stats) {
        VG_STATS(VG_GIRAFFE.out.gam)
        ch_reports = ch_reports.mix(VG_STATS.out.stats)
        ch_versions = ch_versions.mix(VG_STATS.out.versions)
    }

    // 7. Mix Outputs
    bam           = VG_GIRAFFE.out.bam.mix(VG_SURJECT.out.bam)
    bai           = VG_GIRAFFE.out.bai.mix(VG_SURJECT.out.bai)
    cram          = VG_GIRAFFE.out.cram.mix(VG_SURJECT.out.cram)
    crai          = VG_GIRAFFE.out.crai.mix(VG_SURJECT.out.crai)
    markdup_stats = VG_GIRAFFE.out.markdup_stats.mix(VG_SURJECT.out.markdup_stats)
    reports       = ch_reports

    emit:
    bam
    bai
    cram
    crai
    markdup_stats
    reports
    versions = ch_versions
}
