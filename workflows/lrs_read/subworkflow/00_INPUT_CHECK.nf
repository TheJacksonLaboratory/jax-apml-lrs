/*
========================================================================================
    00_INPUT_CHECK
========================================================================================
    Validates the input samplesheet and creates BAM input channels.

    Expected samplesheet columns (CSV, comma-delimited):
        sID                  -- sample identifier
        reads_unaligned_bam  -- path to unaligned PacBio BAM file
        HPO                  -- path to file containing HPO phenotype term IDs (one per line)

    Emits:
        bams  -- channel of [ meta, bam, hpo_file ]
========================================================================================
*/

include { SAMPLESHEET_CHECK_READ } from '../../../modules/samplesheet_check_read'

workflow INPUT_CHECK {

    take:
    samplesheet // file: path to samplesheet CSV

    main:
    SAMPLESHEET_CHECK_READ ( samplesheet )
        .csv
        .splitCsv ( header: true, sep: ',' )
        .map { row -> create_bam_channel(row) }
        .set { bams }

    emit:
    bams                                      // channel: [ val(meta), path(bam), path(hpo) ]
    versions = SAMPLESHEET_CHECK_READ.out.versions // channel: [ versions.yml ]
}

// ---------------------------------------------------------------------------
// Helper: build [ meta, bam, hpo ] tuple from a samplesheet row
// ---------------------------------------------------------------------------
def create_bam_channel(LinkedHashMap row) {

    def meta = [:]
    meta.id = row.sID

    if (!file(row.reads_unaligned_bam).exists()) {
        exit 1, "ERROR: Unaligned BAM file not found for sample '${meta.id}':\n  ${row.reads_unaligned_bam}"
    }
    if (!file(row.HPO).exists()) {
        exit 1, "ERROR: HPO terms file not found for sample '${meta.id}':\n  ${row.HPO}"
    }

    return [ meta, file(row.reads_unaligned_bam), file(row.HPO) ]
}