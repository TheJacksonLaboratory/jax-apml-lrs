/*
========================================================================================
    04_CALL_METHYLATION
========================================================================================
    Computes CpG methylation scores from haplotagged PacBio reads using
    pb-CpG-tools, then generates a methylation region profile using MethBat.

    NOTE: This subworkflow is currently in development and results should be
    considered experimental.

    Processes:
        METHBAT_ALIGNED_BAM  -- compute per-haplotype CpG methylation scores
                                from haplotagged BAM using pb-CpG-tools
        METHBAT_PROFILE      -- generate a methylation region profile using MethBat

    Input:
        whatshap2_bam_bai        -- channel: [ val(meta), path(haplotagged.bam), path(bai) ]
        methbat_input_regions    -- path to methylation region profile model TSV
        refID                    -- reference genome identifier (e.g. hg38)
        nThread                  -- number of threads for pb-CpG-tools

    Emits:
        methbat_bams          -- channel: [ val(meta), path(combined.bed), path(combined.bw),
                                            path(hap1.bed), path(hap1.bw), path(hap2.bed),
                                            path(hap2.bw), path(log), val(prefix) ]
        methbat_region_profile -- channel: [ val(meta), path(methbat.tsv) ]

    Notes:
        pb-CpG-tools (aligned_bam_to_cpg_scores) is bundled within the container
        at /usr/local/bin/. No external binary path is required.
        The pileup calling model is bundled within the pb-CpG-tools container and
        does not need to be supplied separately.
========================================================================================
*/

process METHBAT_ALIGNED_BAM {
    tag "$meta.id"
    container 'quay.io/pacbio/pb-cpg-tools:v2.3.2'
    cpus   params.methbat.cpus
    memory params.methbat.mem

    input:
    val   nThread
    val   refID
    tuple val(meta), path(whatshap2_bam_out_file_name), path(whatshap2_bam_bai_out_file_name)

    output:
    tuple val(meta), path(out_combined_bed), path(out_combined_bw), path(out_hap1_bed), path(out_hap1_bw), path(out_hap2_bed), path(out_hap2_bw), path(out_log), val(out_prefix_out), emit: methbat_bams

    script:
    out_prefix_out   = meta.id + "_" + refID + "_pbmm2"
    out_combined_bed = meta.id + "_" + refID + "_pbmm2.combined.bed"
    out_combined_bw  = meta.id + "_" + refID + "_pbmm2.combined.bw"
    out_hap1_bed     = meta.id + "_" + refID + "_pbmm2.hap1.bed"
    out_hap1_bw      = meta.id + "_" + refID + "_pbmm2.hap1.bw"
    out_hap2_bed     = meta.id + "_" + refID + "_pbmm2.hap2.bed"
    out_hap2_bw      = meta.id + "_" + refID + "_pbmm2.hap2.bw"
    out_log          = meta.id + "_" + refID + "_pbmm2.log"

    """
    /usr/local/bin/aligned_bam_to_cpg_scores \\
        --bam $whatshap2_bam_out_file_name \\
        --output-prefix $out_prefix_out \\
        --model /opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \\
        --threads $nThread
    gunzip ${out_combined_bed}.gz &
    gunzip ${out_hap1_bed}.gz &
    gunzip ${out_hap2_bed}.gz &
    wait
    """
}

process METHBAT_PROFILE {
    tag "$meta.id"
    container 'quay.io/biocontainers/methbat:0.13.2--h9ee0642_0'

    input:
    path  methbat_input_regions
    val   refID
    tuple val(meta), path(out_combined_bed), path(out_combined_bw), path(out_hap1_bed), path(out_hap1_bw), path(out_hap2_bed), path(out_hap2_bw), path(out_log), val(out_prefix_methbat_script)

    output:
    tuple val(meta), path(out_region_profile), emit: methbat_region_profile

    script:
    out_region_profile = meta.id + "_" + refID + "_pbmm2_methbat.tsv"

    """
    mkdir -p methbat
    methbat profile \\
        --input-prefix ./${meta.id}_${refID}_pbmm2 \\
        --input-regions $methbat_input_regions \\
        --output-region-profile $out_region_profile \\
        --profile-label ALL
    """
}