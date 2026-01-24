//
// STR Integration: ExpansionHunter + GangSTR + STRling -> EnsembleTR
//

include { EXPANSIONHUNTER } from '../../../modules/nf-core/expansionhunter/main'
include { GANGSTR } from '../../../modules/nf-core/gangstr/main'
include { STRLING_TO_VCF } from '../../../modules/local/strling/to_vcf/main'
include { ENSEMBLETR } from '../../../modules/local/str/ensembletr/main.nf'

process BCFTOOLS_REHEADER_EH {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.reheader"
    """
    echo "${meta.id}" > sample.txt
    bcftools reheader -s sample.txt -o ${prefix}.vcf.gz ${vcf}
    """
}

process BCFTOOLS_REHEADER_GS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.reheader"
    """
    echo "${meta.id}" > sample.txt
    bcftools reheader -s sample.txt -o ${prefix}.vcf.gz ${vcf}
    """
}

workflow STR_ENSEMBLETR_INTEGRATION {
    take:
    ch_alignment               // channel: [ val(meta), bam, bai ]
    ch_fasta                   // channel: [ val(meta), fasta ]
    ch_fasta_fai               // channel: [ val(meta), fasta_fai ]
    ch_expansionhunter_catalog // channel: [ val(meta), catalog ]
    ch_gangstr_catalog         // channel: [ val(meta), catalog ]
    ch_strling_genotype        // channel: [ val(meta), genotype ]
    ch_strling_catalog         // channel: [ val(meta), catalog ]

    main:
    ch_versions = channel.empty()

    // ExpansionHunter
    ch_bam_eh = ch_alignment.map { meta, bam, bai ->
        def new_meta = meta.clone()
        if (meta.sex == 'XX') new_meta.sex = 'female'
        else if (meta.sex == 'XY') new_meta.sex = 'male'
        else new_meta.sex = null
        [ new_meta, bam, bai ]
    }
    EXPANSIONHUNTER(
        ch_bam_eh,
        ch_fasta.collect(),
        ch_fasta_fai.collect(),
        ch_expansionhunter_catalog.collect()
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)
    ch_vcf_eh_raw = EXPANSIONHUNTER.out.vcf.map { meta, vcf ->
        [ meta, file(vcf, stageAs: "${meta.id}.eh.vcf.gz") ]
    }
    BCFTOOLS_REHEADER_EH(ch_vcf_eh_raw)
    ch_vcf_eh = BCFTOOLS_REHEADER_EH.out.vcf

    // GangSTR
    ch_gangstr_input = ch_alignment.combine(ch_gangstr_catalog.collect())
        .map { meta, bam, bai, cat_meta, catalog ->
            [ meta, [bam], [bai], catalog ]
        }

    GANGSTR(
        ch_gangstr_input,
        ch_fasta.map { it[1] }.collect(),
        ch_fasta_fai.map { it[1] }.collect()
    )
    ch_versions = ch_versions.mix(GANGSTR.out.versions)
    ch_vcf_gangstr_raw = GANGSTR.out.vcf.map { meta, vcf ->
        [ meta, file(vcf, stageAs: "${meta.id}.gs.vcf.gz") ]
    }
    BCFTOOLS_REHEADER_GS(ch_vcf_gangstr_raw)
    ch_vcf_gangstr = BCFTOOLS_REHEADER_GS.out.vcf

    // STRling to VCF
    ch_strling_catalog_path = ch_strling_catalog.map { it[1] }.collect()
    ch_strling_ref_fai = ch_fasta_fai.map { it[1] }.collect()
    STRLING_TO_VCF(
        ch_strling_genotype,
        ch_strling_catalog_path,
        ch_strling_ref_fai
    )
    ch_versions = ch_versions.mix(STRLING_TO_VCF.out.versions)
    ch_vcf_strling = STRLING_TO_VCF.out.vcf.map { meta, vcf ->
        [ meta, file(vcf, stageAs: "${meta.id}.sl.vcf.gz") ]
    }

    // EnsembleTR
    ch_eh_keyed = ch_vcf_eh.map { meta, vcf -> [ meta.id, meta, vcf ] }
    ch_gs_keyed = ch_vcf_gangstr.map { meta, vcf -> [ meta.id, meta, vcf ] }
    ch_sl_keyed = ch_vcf_strling.map { meta, vcf -> [ meta.id, meta, vcf ] }

    ch_ensembletr_base = ch_eh_keyed
        .combine(ch_gs_keyed, by: 0)
        .combine(ch_sl_keyed, by: 0)
        .map { _id, meta_eh, eh_vcf, meta_gs, gs_vcf, meta_sl, sl_vcf ->
            [ meta_eh, eh_vcf, gs_vcf, sl_vcf ]
        }

    ch_ref_fasta = ch_fasta.map { it[1] }.collect()
    ch_ref_fai = ch_fasta_fai.map { it[1] }.collect()

    ch_ensembletr_input = ch_ensembletr_base
        .combine(ch_ref_fasta)
        .combine(ch_ref_fai)
        .map { meta, eh_vcf, gs_vcf, sl_vcf, ref_fasta, ref_fai ->
            def fasta = ref_fasta instanceof List ? ref_fasta[0] : ref_fasta
            def fai = ref_fai instanceof List ? ref_fai[0] : ref_fai
            [ meta, eh_vcf, gs_vcf, sl_vcf, fasta, fai ]
        }

    ENSEMBLETR(ch_ensembletr_input)
    ch_versions = ch_versions.mix(ENSEMBLETR.out.versions)

    emit:
    vcf_eh = ch_vcf_eh
    vcf_gangstr = ch_vcf_gangstr
    vcf_ensembletr = ENSEMBLETR.out.vcf
    versions = ch_versions
}
