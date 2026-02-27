/*
========================================================================================
    00_INPUT_CHECK
========================================================================================
    Validates the input samplesheet and creates the input channel for the
    lrs_asm_single workflow.

    Expected samplesheet format (CSV):
        sID,hifi_fastq,HPO
        SAMPLE01,/path/to/sample.fastq.gz,/path/to/HPO.txt

    Columns:
        sID        -- sample identifier
        hifi_fastq -- path to HiFi FASTQ file
        HPO        -- path to HPO phenotype terms file (one term per line)

    Emits:
        bams      -- channel: [ val(meta), path(hifi_fastq), path(HPO) ]
        versions  -- channel: [ versions.yml ]
========================================================================================
*/

include { SAMPLESHEET_CHECK_ASM_SINGLE } from '../../modules/samplesheet_check_asm-single'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK_ASM_SINGLE(samplesheet)
        .csv
        .splitCsv(header: true, sep: ',')
        .map { create_fastq_channel(it) }
        .set { bams }

    emit:
    bams                                          // channel: [ val(meta), path(hifi_fastq), path(HPO) ]
    versions = SAMPLESHEET_CHECK_ASM_SINGLE.out.versions
}

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sID

    if (!file(row.hifi_fastq).exists()) {
        exit 1, "ERROR: HiFi FASTQ file does not exist for sample ${meta.id}:\n${row.hifi_fastq}"
    }
    if (!file(row.HPO).exists()) {
        exit 1, "ERROR: HPO file does not exist for sample ${meta.id}:\n${row.HPO}"
    }

    return [ meta, file(row.hifi_fastq), file(row.HPO) ]
}