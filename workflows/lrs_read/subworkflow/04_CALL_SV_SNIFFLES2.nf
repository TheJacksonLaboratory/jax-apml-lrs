/*
========================================================================================
    04_CALL_SV_SNIFFLES2
========================================================================================
    Calls and filters structural variants from aligned PacBio reads using Sniffles2.

    Processes:
        CALL_SVS_SNIFFLES2_CALL       -- call SVs from aligned BAM
        CALL_SVS_SNIFFLES2_FILTERING  -- filter to PASS variants

    Input:
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
        refID              -- reference genome identifier (e.g. hg38)

    Emits:
        sniffles2_vcf      -- channel: [ val(meta), path(sniffles2.vcf) ]
        sniffles2_FILT_vcf -- channel: [ val(meta), path(sniffles2_FILT.vcf) ]
========================================================================================
*/

process CALL_SVS_SNIFFLES2_CALL {
    tag "$meta.id"
    container 'quay.io/biocontainers/sniffles:2.3.3--pyhdfd78af_0'

    input:
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)

    output:
    tuple val(meta), path(sniffles2_vcf_out_file_name), emit: sniffles2_vcf

    script:
    sniffles2_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_sniffles2.vcf"

    """
    sniffles \\
        -i $input_bam \\
        -v $sniffles2_vcf_out_file_name
    """
}

process CALL_SVS_SNIFFLES2_FILTERING {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(sniffles2_vcf_out_file_name)

    output:
    tuple val(meta), path(sniffles2_FILT_vcf_out_file_name), emit: sniffles2_FILT_vcf

    script:
    sniffles2_FILT_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_sniffles2_FILT.vcf"

    """
    bcftools view -f PASS $sniffles2_vcf_out_file_name > $sniffles2_FILT_vcf_out_file_name
    """
}