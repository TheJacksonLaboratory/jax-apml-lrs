/*
========================================================================================
    04_CALL_CNV_HIFICNV
========================================================================================
    Calls copy number variants (CNVs) from aligned PacBio reads using HiFiCNV.

    Processes:
        CALL_CNVS_HIFICNV  -- call CNVs from aligned BAM, output VCF and depth tracks

    Input:
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file     -- path to reference genome FASTA
        excl_regions_file  -- path to excluded regions BED file (common CNV regions)
        refID              -- reference genome identifier (e.g. hg38)
        nThread            -- number of threads

    Emits:
        hificnv_vcf  -- channel: [ val(meta), path(vcf), path(depth.bw),
                                   path(copynum.bedgraph), path(hificnv.log) ]
========================================================================================
*/

process CALL_CNVS_HIFICNV {
    tag "$meta.id"
    container 'quay.io/pacbio/hificnv:1.0.1_build1'

    input:
    path  excl_regions_file
    path  reference_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)
    val   nThread

    output:
    tuple val(meta), path(hificnv_vcf_out_file_name), path(hificnv_depth_bw_out_file_name), path(hificnv_depth_copynum_bedgraph_out_file_name), path("hificnv.log"), emit: hificnv_vcf

    script:
    output_prefix                              = "hificnv"
    hificnv_vcf_gz_out_file_name               = "hificnv." + meta.id + ".vcf.gz"
    hificnv_vcf_out_file_name                  = "hificnv." + meta.id + ".vcf"
    hificnv_depth_bw_out_file_name             = "hificnv." + meta.id + ".depth.bw"
    hificnv_depth_copynum_bedgraph_out_file_name = "hificnv." + meta.id + ".copynum.bedgraph"

    """
    hificnv \\
        --bam $input_bam \\
        --ref $reference_file \\
        --exclude $excl_regions_file \\
        --threads $nThread \\
        --output-prefix $output_prefix
    gzip -d $hificnv_vcf_gz_out_file_name
    """
}