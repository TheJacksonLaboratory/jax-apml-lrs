/*
========================================================================================
    SAMPLESHEET_CHECK_READ
========================================================================================
    Validates the input samplesheet CSV using check_samplesheet_read.py (in bin/).
    Outputs a validated CSV and a versions.yml for provenance tracking.

    Input:
        samplesheet  -- path to input samplesheet CSV

    Emits:
        csv       -- validated samplesheet CSV
        versions  -- versions.yml
========================================================================================
*/

process SAMPLESHEET_CHECK_READ {
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
    check_samplesheet_read.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}