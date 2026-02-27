/*
========================================================================================
    01_PBMERGE_BAMS
========================================================================================
    Merges multiple HiFi unaligned 5mC BAM files for a single sample into a
    single merged BAM using pbmerge.

    Processes:
        MERGE_BAMS  -- merge input BAM files with pbmerge

    Input:
        bams  -- channel: [ val(meta), list(bam_files) ]

    Emits:
        merged_bam  -- channel: [ val(meta), path(<sampleID>_pbmerge_merged.5mC.bam) ]

    Notes:
        Input BAM files must be unaligned PacBio HiFi BAMs with 5mC methylation
        tags. The number of input BAMs is flexible — all BAM columns from the
        samplesheet are collected and passed to pbmerge.
========================================================================================
*/

process MERGE_BAMS {
    tag "$meta.id"
    container 'quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0'
    memory '32 GB'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path(out_file_name), emit: merged_bam

    script:
    out_file_name = "${meta.id}_pbmerge_merged.5mC.bam"

    """
    pbmerge -o $out_file_name ${bams.join(' ')}
    """
}

workflow PBMERGE {
    take:
    bams  // channel: [ val(meta), list(bam_files) ]

    main:
    MERGE_BAMS(bams)

    emit:
    merged_bam = MERGE_BAMS.out.merged_bam
}