//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr         // channel: [mandatory] meta, cram, crai, table
        tools                   //
        skip_tools              //
        outdir                  //
        save_output_as_bam      //

    main:
        // Creating csv files to restart from this step
        // Determine the output directory based on tools/skip_tools
        cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram_file, crai_file, table_file ->

            def patient = meta.patient
            def sample  = meta.sample
            def sex     = meta.sex
            def status  = meta.status
            def suffix_aligned = save_output_as_bam ? "bam" : "cram"
            def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"

            // Determine output paths based on tools configuration
            def output_subdir = (tools && tools.split(',').contains('sentieon_dedup')) ? "sentieon_dedup" :
                               (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) ? "markduplicates" :
                               "mapped"

            def csv_name = (tools && tools.split(',').contains('sentieon_dedup')) ? "markduplicates.csv" :
                          (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) ? "markduplicates.csv" :
                          "sorted.csv"

            def cram_path = "${outdir}/preprocessing/${output_subdir}/${sample}/${cram_file.baseName}.${suffix_aligned}"
            def crai_path = "${outdir}/preprocessing/${output_subdir}/${sample}/${crai_file.baseName.minus(".cram")}.${suffix_index}"
            def table_path = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

            def type = save_output_as_bam ? "bam" : "cram"
            def type_index = save_output_as_bam ? "bai" : "crai"

            [csv_name, "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${cram_path},${crai_path},${table_path}\n"]
        }
}
