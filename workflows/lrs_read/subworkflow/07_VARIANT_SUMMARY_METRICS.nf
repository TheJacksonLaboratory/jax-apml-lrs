/*
========================================================================================
    07_VARIANT_SUMMARY_METRICS
========================================================================================
    Computes summary statistics for small variant and structural variant call sets
    using bcftools stats and SURVIVOR stats.

    Processes:
        VARIANT_SUMMARY_METRICS_BCFTOOLS  -- compute stats on filtered DeepVariant VCF
        VARIANT_SUMMARY_METRICS_SURVIVOR  -- compute stats on PBSV, Sniffles2, and
                                             Delly SV call sets

    Input:
        deepvariant2_filt_vcfs   -- channel: [ val(meta), path(FILT.vcf), path(Emedgene.vcf) ]
        svs_filtered_vcfs        -- pbsv channel: [ val(meta), path(Emedgene.vcf), path(FILT.vcf) ]
        sniffles2_FILT_vcf       -- Sniffles2 channel: [ val(meta), path(FILT.vcf) ]
        delly_out_filt_vcf       -- Delly channel: [ val(meta), path(telomereEX.vcf), path(FILT.vcf) ]
        refID                    -- reference genome identifier (e.g. hg38)

    Emits:
        deepvariant_varmetrics  -- channel: [ val(meta), path(deepvariant_varmetrics.txt) ]
        pbsv_varmetrics         -- channel: [ val(meta), path(pbsv_varmetrics.txt),
                                              path(pbsv_varmetrics_CHR), path(pbsv_varmetrics_support) ]
        sniffles2_varmetrics    -- channel: [ val(meta), path(sniffles2_varmetrics.txt),
                                              path(sniffles2_varmetrics_CHR), path(sniffles2_varmetrics_support) ]
        delly_varmetrics        -- channel: [ val(meta), path(delly_varmetrics.txt),
                                              path(delly_varmetrics_CHR), path(delly_varmetrics_support) ]
========================================================================================
*/

process VARIANT_SUMMARY_METRICS_BCFTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(deepvariant2_FILT_vcf), path(deepvariant2_Emedgene_vcf)

    output:
    tuple val(meta), path(out_deepvariant_varmetrics_file), emit: deepvariant_varmetrics

    script:
    out_deepvariant_varmetrics_file = meta.id + "_" + refID + "_pbmm2_deepvariant_filt_varmetrics.txt"

    """
    bcftools stats $deepvariant2_FILT_vcf > $out_deepvariant_varmetrics_file
    """
}

process VARIANT_SUMMARY_METRICS_SURVIVOR {
    tag "$meta.id"
    container 'mgibio/survivor-cwl:1.0.6.2'

    input:
    tuple val(meta), path(pbsv_filt_output), path(pbsv_Emedgene_output), path(sniffles2_filt_output), path(delly_vcf_out_file_name), path(delly_filt_vcf_out_file_name)
    val   refID

    output:
    tuple val(meta),  path(out_pbsv_varmetrics_file),     path(out_pbsv_varmetrics_CHR_file),     path(out_pbsv_varmetrics_support_file),     emit: pbsv_varmetrics
    tuple val(meta1), path(out_sniffles2_varmetrics_file), path(out_sniffles2_varmetrics_CHR_file), path(out_sniffles2_varmetrics_support_file), emit: sniffles2_varmetrics
    tuple val(meta2), path(out_delly_varmetrics_file),     path(out_delly_varmetrics_CHR_file),     path(out_delly_varmetrics_support_file),     emit: delly_varmetrics

    script:
    out_pbsv_varmetrics_file          = meta.id + "_" + refID + "_pbmm2_pbsv_varmetrics.txt"
    out_pbsv_varmetrics_CHR_file      = meta.id + "_" + refID + "_pbmm2_pbsv_varmetrics.txt_CHR"
    out_pbsv_varmetrics_support_file  = meta.id + "_" + refID + "_pbmm2_pbsv_varmetrics.txtsupport"
    out_sniffles2_varmetrics_file         = meta.id + "_" + refID + "_pbmm2_sniffles2_varmetrics.txt"
    out_sniffles2_varmetrics_CHR_file     = meta.id + "_" + refID + "_pbmm2_sniffles2_varmetrics.txt_CHR"
    out_sniffles2_varmetrics_support_file = meta.id + "_" + refID + "_pbmm2_sniffles2_varmetrics.txtsupport"
    out_delly_varmetrics_file         = meta.id + "_" + refID + "_pbmm2_delly_varmetrics.txt"
    out_delly_varmetrics_CHR_file     = meta.id + "_" + refID + "_pbmm2_delly_varmetrics.txt_CHR"
    out_delly_varmetrics_support_file = meta.id + "_" + refID + "_pbmm2_delly_varmetrics.txtsupport"
    meta1 = meta
    meta2 = meta

    """
    SURVIVOR stats $pbsv_filt_output -1 -1 -1 $out_pbsv_varmetrics_file
    SURVIVOR stats $sniffles2_filt_output -1 -1 -1 $out_sniffles2_varmetrics_file
    SURVIVOR stats $delly_vcf_out_file_name -1 -1 -1 $out_delly_varmetrics_file
    """
}