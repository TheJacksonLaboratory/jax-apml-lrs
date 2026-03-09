/*
========================================================================================
    00_INPUT_CHECK
========================================================================================
    Validates the input samplesheet and creates the input channel for the
    lrs_pbmerge workflow.

    Expected samplesheet format (CSV):
        sID,bam1,bam2[,bam3]
        SAMPLE01,/path/to/run1.bam,/path/to/run2.bam
        SAMPLE02,/path/to/run1.bam,/path/to/run2.bam,/path/to/run3.bam

    Columns:
        sID   -- sample identifier
        bam1  -- path to first HiFi unaligned 5mC BAM file
        bam2  -- path to second HiFi unaligned 5mC BAM file
        bam3  -- path to third HiFi unaligned 5mC BAM file (optional)

    Emits:
        bams      -- channel: [ val(meta), list(bam_files) ]
        versions  -- channel: [ versions.yml ]
========================================================================================
*/

include { SAMPLESHEET_CHECK_PBMERGE } from '../../../modules/samplesheet_check_pbmerge'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK_PBMERGE(samplesheet)
        .csv
        .splitCsv(header: true, sep: ',')
        .map { create_bam_channel(it) }
        .set { bams }

    emit:
    bams                                             // channel: [ val(meta), list(bam_files) ]
    versions = SAMPLESHEET_CHECK_PBMERGE.out.versions
}

def create_bam_channel(LinkedHashMap row) {
    def meta = [ id: row.sID ]

    // Collect all non-empty BAM file paths (all columns except sID)
    // Supports 2 or 3 BAM files per sample
    def bam_paths = []
    row.each { key, value ->
        if (key != "sID" && value?.trim()) {
            bam_paths.add(value)
        }
    }

    if (bam_paths.isEmpty()) {
        exit 1, "ERROR: No BAM files found for sample ${meta.id}"
    }

    bam_paths.each { path ->
        if (!file(path).exists()) {
            exit 1, "ERROR: BAM file does not exist for sample ${meta.id}:\n${path}"
        }
    }

    def input_files = bam_paths.collect { path -> file(path) }

    return [ meta, input_files ]
}