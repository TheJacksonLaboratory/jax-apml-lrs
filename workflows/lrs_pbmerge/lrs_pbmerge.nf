#!/usr/bin/env nextflow

nextflow.enable.dsl=2
nextflow.preview.output = true

/*
========================================================================================
    LRS_PBMERGE WORKFLOW
========================================================================================
    Merges multiple HiFi unaligned 5mC BAM files per sample into a single BAM
    using pbmerge, in preparation for downstream alignment and variant calling.

    Steps:
        00  Input validation and samplesheet parsing
        01  Merge HiFi BAM files (pbmerge)

    Usage:
        nextflow run lrs_pbmerge.nf -profile <profile> --csv_path <samplesheet.csv>

    Required parameters:
        --csv_path   path to input samplesheet CSV
        --outputDir  path to output directory

    Expected samplesheet format:
        sID,bam1,bam2[,bam3]
        SAMPLE01,/path/to/run1.bam,/path/to/run2.bam
        SAMPLE02,/path/to/run1.bam,/path/to/run2.bam,/path/to/run3.bam

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

include { INPUT_CHECK } from './subworkflow/00_INPUT_CHECK.nf'
include { PBMERGE     } from './subworkflow/01_PBMERGE_BAMS.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    main:

    // ---------------------------------------------------------------------------
    // 00 Input check
    // ---------------------------------------------------------------------------
    input_check_out = INPUT_CHECK(
        ch_input
    )

    // ---------------------------------------------------------------------------
    // 01 Merge BAMs
    // ---------------------------------------------------------------------------
    pbmerge_out = PBMERGE(
        input_check_out.bams
    )

    /*
    ========================================================================================
        PUBLISH
    ========================================================================================
    */

    publish:
    pbmerge_out.merged_bam >> 'pbmerge'
}

/*
========================================================================================
    OUTPUT
========================================================================================
*/

output {
    'pbmerge' {
        path { meta, out_file_name -> "${params.outputDir}${meta.id}/pbmerge/${params.runDate}" }
        mode 'copy'
        overwrite false
    }
}