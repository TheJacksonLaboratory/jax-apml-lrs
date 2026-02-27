#!/usr/bin/env nextflow

nextflow.enable.dsl=2
nextflow.preview.output = true

/*
========================================================================================
    LRS_ASM_TRIO WORKFLOW
========================================================================================
    Trio-based de novo genome assembly and variant calling pipeline for PacBio
    HiFi long-read sequencing data.

    Steps:
        00  Input validation and samplesheet parsing
        01  Trio-aware de novo assembly (Yak + HiFiasm + samtools)
        02  Structural and small variant calling from assemblies (PAV)
        03  Split PAV output by variant size:
                - Small variants (SNVs + indels < 50 bp)
                - Structural variants (>= 50 bp)
        04  Phenotype-driven SV prioritization (SvAnna + filter_svanna_pav.R)

    Usage:
        nextflow run lrs_asm_trio.nf -profile <profile> --csv_path <samplesheet.csv>

    Required parameters:
        --csv_path                        path to input samplesheet CSV
        --outputDir                       path to output directory
        --reference_file_path             path to reference genome FASTA
        --reference_index_file_path       path to reference genome FASTA index

    Expected samplesheet format:
        sID,proband_hifi_fastq,mat_hifi_fastq,pat_hifi_fastq,HPO

    See nextflow.config for full parameter documentation.
========================================================================================
*/

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (params.csv_path) { ch_input = file(params.csv_path) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                                  } from './subworkflow/00_INPUT_CHECK.nf'
include { YAK;
          ASSEMBLY_HIFIASM_TRIO;
          SAMTOOLS_HIFIASM_TRIO                        } from './subworkflow/01_ASSEMBLE_TRIO_HIFIASM.nf'
include { CALL_SVS_PAV                                 } from './subworkflow/02_CALL_SV_PAV.nf'
include { SPLIT_SMALL_VARIANTS_PAV;
          INDEX_SMALL_VARIANTS_PAV;
          CONCAT_SMALL_VARIANTS_PAV;
          SPLIT_SV_PAV;
          INDEX_SV_PAV;
          CONCAT_SV_PAV                                } from './subworkflow/03_SPLIT_PAV_VARIANT_SIZE.nf'
include { PRIORITIZE_BY_SVANNA;
          RSCRIPT_PROCESS_SVANNA_RESULT                } from './subworkflow/04_PRIORITIZE_SV_SVANNA.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    main:

    // Reference and resource parameters
    reference_file       = params.reference_file_path
    reference_index_file = params.reference_index_file_path
    refID                = params.referenceID
    n_proc               = params.n_proc
    hifiasm_threads      = params.hifiasm.threads

    // SvAnna
    pav_svanna_r         = params.svAnna.pav_svanna_r_script_path

    // ---------------------------------------------------------------------------
    // 00 Input check
    // ---------------------------------------------------------------------------
    input_check_out = INPUT_CHECK(
        ch_input
    )

    // ---------------------------------------------------------------------------
    // 01 Trio assembly
    // ---------------------------------------------------------------------------
    yak_out = YAK(
        input_check_out.bams,
        hifiasm_threads
    )

    assembly_hifiasm_trio_out = ASSEMBLY_HIFIASM_TRIO(
        input_check_out.bams,
        yak_out.out_yak,
        hifiasm_threads
    )

    samtools_hifiasm_trio_out = SAMTOOLS_HIFIASM_TRIO(
        assembly_hifiasm_trio_out.out_assembled_ctg
    )

    // ---------------------------------------------------------------------------
    // 02 Call variants with PAV
    // ---------------------------------------------------------------------------
    call_svs_pav_out = CALL_SVS_PAV(
        assembly_hifiasm_trio_out.out_assembled_ctg,
        samtools_hifiasm_trio_out.out_assembled_utg,
        file(reference_file),
        file(reference_index_file),
        n_proc,
        refID
    )

    // ---------------------------------------------------------------------------
    // 03 Split PAV output by variant size
    // ---------------------------------------------------------------------------

    // Small variant track
    split_small_variants_pav_out = SPLIT_SMALL_VARIANTS_PAV(
        call_svs_pav_out.out_vcf
    )

    index_small_variants_pav_out = INDEX_SMALL_VARIANTS_PAV(
        split_small_variants_pav_out.out_splitted_vcf
    )

    concat_small_variants_pav_out = CONCAT_SMALL_VARIANTS_PAV(
        index_small_variants_pav_out.out_small_variant_vcf_indexed,
        call_svs_pav_out.out_vcf,
        refID
    )

    // SV track
    split_sv_pav_out = SPLIT_SV_PAV(
        call_svs_pav_out.out_vcf
    )

    index_sv_pav_out = INDEX_SV_PAV(
        split_sv_pav_out.out_splitted_vcf
    )

    concat_sv_pav_out = CONCAT_SV_PAV(
        index_sv_pav_out.out_ins_del_vcf_gz,
        refID
    )

    // ---------------------------------------------------------------------------
    // 04 Phenotype-driven SV prioritization
    // ---------------------------------------------------------------------------
    prioritize_by_svanna_out = PRIORITIZE_BY_SVANNA(
        input_check_out.bams,
        concat_small_variants_pav_out.out_small_variant_vcf,
        refID
    )

    rscript_process_svanna_result_out = RSCRIPT_PROCESS_SVANNA_RESULT(
        prioritize_by_svanna_out.out_svanna_csv,
        file(pav_svanna_r)
    )

    /*
    ========================================================================================
        PUBLISH
    ========================================================================================
    */

    publish:
    yak_out.out_yak                                              >> 'yak'
    assembly_hifiasm_trio_out.out_assembled_ctg                  >> 'assembly_hifiasm_trio'
    samtools_hifiasm_trio_out.out_assembled_utg                  >> 'samtools_hifiasm_trio'

    call_svs_pav_out.out_vcf                                     >> 'call_svs_pav'

    split_small_variants_pav_out.out_splitted_vcf                >> 'split_small_variants_pav'
    index_small_variants_pav_out.out_small_variant_vcf_indexed   >> 'index_small_variants_pav'
    concat_small_variants_pav_out.out_small_variant_vcf          >> 'concat_small_variants_pav'

    split_sv_pav_out.out_splitted_vcf                            >> 'split_sv_pav'
    index_sv_pav_out.out_ins_del_vcf_gz                          >> 'index_sv_pav'
    concat_sv_pav_out.out_concat_sv_vcf                          >> 'concat_sv_pav'

    prioritize_by_svanna_out.out_svanna_csv                      >> 'prioritize_by_svanna'
    rscript_process_svanna_result_out.out_csv                    >> 'rscript_process_svanna_result'
}

/*
========================================================================================
    OUTPUT
========================================================================================
*/

output {
    'yak' {
        path { meta, pat_yak, mat_yak -> "${params.outputDir}${meta.id}/yak" }
        mode 'copy'
        overwrite false
    }
    'assembly_hifiasm_trio' {
        path { meta, out_asm_gfa, pat_yak, mat_yak, out_asm_hap1, out_asm_hap2 -> "${params.outputDir}${meta.id}/hifiasm_trio" }
        mode 'copy'
        overwrite false
    }
    'samtools_hifiasm_trio' {
        path { meta, out_utg_fa, out_utg_fa_fai -> "${params.outputDir}${meta.id}/hifiasm_trio" }
        mode 'copy'
        overwrite false
    }
    'call_svs_pav' {
        path { meta, pav_vcf_gz_rid, pav_vcf_gz_tbi_rid -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'split_small_variants_pav' {
        path { meta, pav_snvs_vcf, pav_indels_vcf, pav_indels50bp_vcf -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'index_small_variants_pav' {
        path { meta, pav_snvs_vcf_gz, pav_snvs_vcf_gz_tbi, pav_indels50bp_vcf_gz, pav_indels50bp_vcf_gz_tbi -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'concat_small_variants_pav' {
        path { meta, out_smallvariants_pav_gz, out_pav_vcf -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'split_sv_pav' {
        path { meta, pav_SVins_vcf, pav_SVdel_vcf -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'index_sv_pav' {
        path { meta, pav_SVins_vcf_gz, pav_SVdel_vcf_gz, pav_svins_vcf_gz_tbi, pav_svdel_vcf_gz_tbi -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'concat_sv_pav' {
        path { meta, out_concat_sv_vcf -> "${params.outputDir}${meta.id}/pav" }
        mode 'copy'
        overwrite false
    }
    'prioritize_by_svanna' {
        path { meta, out_pav_svanna_csv -> "${params.outputDir}${meta.id}/svanna" }
        mode 'copy'
        overwrite false
    }
    'rscript_process_svanna_result' {
        path { meta, out_csv -> "${params.outputDir}${meta.id}/svanna" }
        mode 'copy'
        overwrite false
    }
}