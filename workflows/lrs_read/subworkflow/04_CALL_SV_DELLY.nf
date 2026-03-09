/*
========================================================================================
    04_CALL_SV_DELLY
========================================================================================
    Calls and filters structural variants from aligned PacBio long reads using Delly.

    Processes:
        CALL_SVS_DELLY_LR           -- call SVs in long-read mode, output BCF
        CALL_SVS_DELLY_OUT_AND_FILT -- convert BCF to VCF and filter to PASS variants

    Input:
        aligned_bam_index    -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file       -- path to reference genome FASTA
        hg38_excl_file       -- path to telomere/centromere exclusion BED file
        refID                -- reference genome identifier (e.g. hg38)
        min_mapping_quality  -- minimum mapping quality for read inclusion
        alignment_scoring    -- minimum alignment score
        min_clique_size      -- minimum clique size for SV calling
        nThread              -- number of threads

    Emits:
        delly_out_bcf      -- channel: [ val(meta), path(bcf), path(csi) ]
        delly_out_filt_vcf -- channel: [ val(meta), path(telomereEX.vcf), path(FILT.vcf) ]
========================================================================================
*/

process CALL_SVS_DELLY_LR {
    tag "$meta.id"
    container 'quay.io/biocontainers/delly:1.2.8--hf9970c3_0'
    cpus   params.NTHREAD
    memory '24 GB'

    input:
    path  hg38_excl_file
    path  reference_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)
    val   min_mapping_quality
    val   alignment_scoring
    val   min_clique_size
    val   nThread

    output:
    tuple val(meta), path(delly_bcf_out_file_name), path(delly_bcf_csi_out_file_name), emit: delly_out_bcf

    script:
    delly_bcf_out_file_name     = meta.id + "_" + refID + "_pbmm2_delly_telomereEX.bcf"
    delly_bcf_csi_out_file_name = meta.id + "_" + refID + "_pbmm2_delly_telomereEX.bcf.csi"

    """
    export OMP_NUM_THREADS=$nThread
    delly lr \\
        -y pb \\
        -o $delly_bcf_out_file_name \\
        -g $reference_file \\
        -q $min_mapping_quality \\
        -s $alignment_scoring \\
        -z $min_clique_size \\
        -x $hg38_excl_file \\
        $input_bam
    """
}

process CALL_SVS_DELLY_OUT_AND_FILT {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'
    cpus   params.NTHREAD
    memory '24 GB'

    input:
    val   refID
    tuple val(meta), path(delly_bcf_out_file_name), path(delly_bcf_csi_out_file_name)
    val   nThread

    output:
    tuple val(meta), path(delly_vcf_out_file_name), path(delly_filt_vcf_out_file_name), emit: delly_out_filt_vcf

    script:
    delly_vcf_out_file_name      = meta.id + "_" + refID + "_pbmm2_delly_telomereEX.vcf"
    delly_filt_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_delly_FILT.vcf"

    """
    bcftools view --threads $nThread $delly_bcf_out_file_name > $delly_vcf_out_file_name
    bcftools view --threads $nThread -f PASS $delly_vcf_out_file_name > $delly_filt_vcf_out_file_name
    """
}