/*
========================================================================================
    00_INPUT_CHECK
========================================================================================
    Validates the input samplesheet and creates the input channel for the
    lrs_asm_trio workflow.

    Expected samplesheet format (CSV):
        sID,proband_hifi_fastq,mat_hifi_fastq,pat_hifi_fastq,HPO
        SAMPLE01,/path/to/proband.fastq.gz,/path/to/mat.fastq.gz,/path/to/pat.fastq.gz,/path/to/HPO.txt

    Columns:
        sID               -- sample identifier
        proband_hifi_fastq -- path to proband HiFi FASTQ file
        mat_hifi_fastq    -- path to maternal HiFi FASTQ file
        pat_hifi_fastq    -- path to paternal HiFi FASTQ file
        HPO               -- path to HPO phenotype terms file (one term per line)

    Emits:
        bams      -- channel: [ val(meta), path(proband_fastq), path(mat_fastq),
                                path(pat_fastq), path(HPO) ]
        versions  -- channel: [ versions.yml ]
========================================================================================
*/

include { SAMPLESHEET_CHECK_ASM_TRIO } from '../../../modules/samplesheet_check_asm-trio'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK_ASM_TRIO(samplesheet)
        .csv
        .splitCsv(header: true, sep: ',')
        .map { create_fastq_channel(it) }
        .set { bams }

    emit:
    bams                                          // channel: [ val(meta), path(proband_fastq), path(mat_fastq), path(pat_fastq), path(HPO) ]
    versions = SAMPLESHEET_CHECK_ASM_TRIO.out.versions
}

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sID

    if (!file(row.proband_hifi_fastq).exists()) {
        exit 1, "ERROR: Proband HiFi FASTQ file does not exist for sample ${meta.id}:\n${row.proband_hifi_fastq}"
    }
    if (!file(row.mat_hifi_fastq).exists()) {
        exit 1, "ERROR: Maternal HiFi FASTQ file does not exist for sample ${meta.id}:\n${row.mat_hifi_fastq}"
    }
    if (!file(row.pat_hifi_fastq).exists()) {
        exit 1, "ERROR: Paternal HiFi FASTQ file does not exist for sample ${meta.id}:\n${row.pat_hifi_fastq}"
    }
    if (!file(row.HPO).exists()) {
        exit 1, "ERROR: HPO file does not exist for sample ${meta.id}:\n${row.HPO}"
    }

    return [ meta, file(row.proband_hifi_fastq), file(row.mat_hifi_fastq), file(row.pat_hifi_fastq), file(row.HPO) ]
}