/*
========================================================================================
    SAMPLESHEET_CHECK_ASM_TRIO
========================================================================================
    Validates the input samplesheet CSV for the lrs_asm_trio workflow using
    check_samplesheet_asm-trio.py (in bin/).

    Expected samplesheet format:
        sID,proband_hifi_fastq,mat_hifi_fastq,pat_hifi_fastq,HPO

    Input:
        samplesheet  -- path to input samplesheet CSV

    Emits:
        csv       -- validated samplesheet CSV
        versions  -- versions.yml
========================================================================================
*/

process SAMPLESHEET_CHECK_ASM_TRIO {
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
    check_samplesheet_asm-trio.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}