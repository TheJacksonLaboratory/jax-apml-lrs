/*
========================================================================================
    SAMPLESHEET_CHECK_PBMERGE
========================================================================================
    Validates the input samplesheet CSV for the lrs_pbmerge workflow using
    check_samplesheet_pbmerge.py (in bin/).

    Expected samplesheet format:
        sID,bam1,bam2,bam3

    Input:
        samplesheet  -- path to input samplesheet CSV

    Emits:
        csv       -- validated samplesheet CSV
        versions  -- versions.yml
========================================================================================
*/

process SAMPLESHEET_CHECK_PBMERGE {
    tag "$samplesheet"
    label 'process_single'
    container 'public.ecr.aws/docker/library/python:slim-buster'

    input:
    path samplesheet

    output:
    path '*.csv',        emit: csv
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    """
    check_samplesheet_pbmerge.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}