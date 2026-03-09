/*
========================================================================================
    04_CALL_PARALOGS_PARAPHASE
========================================================================================
    Resolves paralogous genes from aligned PacBio reads using Paraphase, then
    extracts copy number and variant information for target genes of interest.

    Processes:
        CALL_PARAPHASE_CALL     -- resolve paralogous gene haplotypes with Paraphase
        CALL_PARAPHASE_EXTRACT  -- extract CN and variant tables from Paraphase JSON
                                   output using convert_paraphase.R
                                   NOTE: this process is currently in development
                                   and its output is not published by default.

    Input:
        aligned_bam_index              -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file                 -- path to reference genome FASTA
        extract_paraphase_json_script  -- path to convert_paraphase.R script
        target_genes_txt               -- path to target genes file (GeneSymbol, Chrom)
        paraphase_threads              -- number of threads for Paraphase

    Emits:
        paraphase_out      -- channel: [ val(meta), path(json), path(bam), path(bai) ]
        paraphase_vcfs_dir -- channel: [ val(meta), path(vcfs/*) ]
========================================================================================
*/

process CALL_PARAPHASE_CALL {
    tag "$meta.id"
    container 'quay.io/pacbio/paraphase:3.1.1_build1'

    input:
    path  reference_file
    tuple val(meta), path(input_bam), path(input_bam_index)
    val   paraphase_threads

    output:
    tuple val(meta), path(out_paraphase_json), path(out_paraphase_bam), path(out_paraphase_bam_bai), emit: paraphase_out
    tuple val(meta1), path(out_dir), emit: paraphase_vcfs_dir

    script:
    out_paraphase_json  = meta.id + ".paraphase.json"
    out_paraphase_bam   = meta.id + ".paraphase.bam"
    out_paraphase_bam_bai = meta.id + ".paraphase.bam.bai"
    out_dir             = meta.id + "_paraphase_vcfs"
    meta1               = meta

    """
    paraphase \\
        -b $input_bam \\
        -o ./ \\
        -r $reference_file \\
        -t $paraphase_threads
    """
}

process CALL_PARAPHASE_EXTRACT {
    tag "$meta.id"
    container 'rocker/tidyverse:4.3'

    input:
    tuple val(meta), path(input_json_file), path(paraphase_bam_file), path(paraphase_bam_bai_file)
    path  extract_paraphase_json_r
    path  target_genes_txt

    output:
    tuple val(meta), path("*.tsv"), emit: paraphase_extract_out

    script:
    """
    Rscript $extract_paraphase_json_r $input_json_file $target_genes_txt
    """
}