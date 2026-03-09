/*
========================================================================================
    02_QUALITY_METRICS
========================================================================================
    Computes alignment quality metrics from an indexed BAM file using NanoPlot,
    samtools coverage, and a custom coverage summary script.

    Subworkflows:
        02_QUALITY_METRICS_NANOPLOT      -- generates alignment QC plots and stats
        02_QUALITY_METRICS_SAMTOOLS      -- computes per-chromosome coverage with samtools
        02_QUALITY_METRICS_GETCOVERAGE   -- summarizes mean autosomal depth from samtools output

    Input:
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
        nThread            -- number of threads for NanoPlot
        refID              -- reference genome identifier (e.g. hg38)
        getcoverage_script -- path to coverage.py script

    Emits:
        nanoplot_results  -- channel: [ val(meta), path(nanoplot/*) ]
        samtools_coverage -- channel: [ val(meta), path(coverage.txt) ]
        meandepth         -- channel: [ val(meta), path(meandepth.tsv) ]
========================================================================================
*/

process QUALITY_METRICS_NANOPLOT {
    tag "$meta.id"
    container 'staphb/nanoplot:1.42.0'
    cpus params.quality_metrics.cpus

    input:
    tuple val(meta), path(input_bam), path(bam_index_file_name)
    val   nThread
    val   refID

    output:
    tuple val(meta), path("nanoplot/*"), emit: nanoplot_results

    script:
    """
    NanoPlot -t $nThread --bam $input_bam -o nanoplot
    """
}

process QUALITY_METRICS_SAMTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(input_bam), path(bam_index_file_name)
    val   refID

    output:
    tuple val(meta), path(samtools_out_file_name), emit: samtools_coverage

    script:
    samtools_out_file_name = meta.id + "_" + refID + "_pbmm2_samtools_coverage.txt"

    """
    samtools coverage $input_bam > $samtools_out_file_name
    """
}

process QUALITY_METRICS_GETCOVERAGE {
    tag "$meta.id"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    tuple val(meta), path(samtools_out_file_name)
    path  script_getcoverage
    val   refID

    output:
    tuple val(meta), path(meandepth_out_file_name), emit: meandepth

    script:
    meandepth_out_file_name = meta.id + "_" + refID + "_pbmm2_meandepth.tsv"

    """
    python3 $script_getcoverage --covfile $samtools_out_file_name
    mv meandepth.tsv $meandepth_out_file_name
    """
}