/*
========================================================================================
    06_GENERATE_BAF_PLOT
========================================================================================
    Generates B-allele frequency (BAF) and copy number (CN) plots per chromosome
    by combining HiFiCNV copy number output with phased DeepVariant variant calls.

    Processes:
        GENERATE_B_ALLELE  -- generate per-chromosome BAF and CN PDF plots
                              using baf_plot.R

    Input:
        hificnv_vcf            -- channel: [ val(meta), path(vcf), path(depth.bw),
                                             path(copynum.bedgraph), path(log) ]
        deepvariant2_filt_vcfs -- channel: [ val(meta), path(FILT.vcf), path(Emedgene.vcf) ]
        b_allele_plot_r        -- path to baf_plot.R script
        refID                  -- reference genome identifier (e.g. hg38)

    Emits:
        report_pdf  -- channel: [ val(meta), path(*_B-Allele_Plot.pdf) ]

    Notes:
        Requires the Bioconductor packages CopyNumberPlots and VariantAnnotation.
        Input VCF must be uncompressed — HiFiCNV output is decompressed in
        04_CALL_CNV_HIFICNV.nf before being passed to this process.
========================================================================================
*/

process GENERATE_B_ALLELE {
    tag "$meta.id"
    container 'quay.io/bioconductor/bioconductor_docker:RELEASE_3_18'

    input:
    path  b_allele_plot_r
    tuple val(meta), path(hificnv_vcf_out_file_name), path(hificnv_depth_bw_out_file_name), path(hificnv_depth_copynum_bedgraph_out_file_name), path(hificnv_log), path(deepvariant2_filt_vcf), path(deepvariant2_Emedgene_vcf)
    val   refID

    output:
    tuple val(meta), path("*.pdf"), emit: report_pdf

    script:
    """
    Rscript $b_allele_plot_r $meta.id $refID
    """
}