/*
========================================================================================
    SAMPLESHEET_CHECK_ASM_SINGLE
========================================================================================
    Validates the input samplesheet CSV for the lrs_asm_single workflow using
    check_samplesheet_asm-single.py (in bin/).

    Expected samplesheet format:
        sID,hifi_fastq,HPO

    Input:
        samplesheet  -- path to input samplesheet CSV

    Emits:
        csv       -- validated samplesheet CSV
        versions  -- versions.yml
========================================================================================
*/

process SAMPLESHEET_CHECK_ASM_SINGLE {
    tag "$samplesheet"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.11'

    input:
    path samplesheet

    output:
    path '*.csv',        emit: csv
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    """
    check_samplesheet_asm-single.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}