/*
========================================================================================
    03_ANNOTATE_SV_SVAFOTATE
========================================================================================
    Annotates and filters structural variants from PAV using SVAFotate population
    allele frequency annotation.

    SVAFotate annotates each SV with population allele frequencies from a reference
    BED file, then variants are filtered into two frequency tiers:
        - RARE_UNIQUE : Max_AF < 0.01
        - LOW_FREQ    : 0.01 <= Max_AF < 0.05

    Processes:
        CALL_SVAFOTATE_ANNOTATE  -- annotate SVs with population allele frequencies
        CALL_SVAFOTATE_FILTER    -- split annotated VCF into RARE_UNIQUE and LOW_FREQ

    Input:
        pav_sv_vcf               -- PAV: channel: [ val(meta), path(pav_SV.vcf) ]
        refID                    -- reference genome identifier (e.g. hg38)
        svafotate_overlap        -- minimum reciprocal overlap for SVAFotate annotation
        svafotate_bed            -- path to SVAFotate population frequency BED file
        nThread                  -- number of threads for SVAFotate

    Emits:
        svafotate_annotate_out   -- channel: [ val(meta), path(SVAFotate-ALL.vcf) ]
        svafotate_filter_out     -- channel: [ val(meta), path(RARE-UNIQUE.vcf), path(LOWFREQ.vcf) ]
========================================================================================
*/

process CALL_SVAFOTATE_ANNOTATE {
    tag "$meta.id"
    container 'jxprismdocker/prism_svafotate:latest'
    cpus   params.svafotate.cpus
    memory params.svafotate.mem

    input:
    val   refID
    tuple val(meta), path(input_vcf)
    val   overlap
    path  svafotate_bed
    val   nThread

    output:
    tuple val(meta), path(vcf_out), emit: svafotate_annotate_out

    script:
    vcf_out = meta.id + "_hifiasm_" + refID + "_pav_SV_Emedgene_SVAFotate-ALL.vcf"

    """
    svafotate annotate \\
        --cpu $nThread \\
        -v $input_vcf \\
        -o $vcf_out \\
        -f $overlap \\
        -b $svafotate_bed
    """
}

process CALL_SVAFOTATE_FILTER {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path(vcf_out_RARE_UNIQUE), path(vcf_out_LOW_FREQ), emit: svafotate_filter_out

    script:
    vcf_out_RARE_UNIQUE = meta.id + "_hifiasm_" + refID + "_pav_SV_Emedgene_SVAFotate-RARE-UNIQUE.vcf"
    vcf_out_LOW_FREQ    = meta.id + "_hifiasm_" + refID + "_pav_SV_Emedgene_SVAFotate-LOWFREQ.vcf"

    """
    sleep 5
    bcftools view -i 'INFO/Max_AF < 0.01' $input_vcf --output $vcf_out_RARE_UNIQUE
    bcftools view -i 'INFO/Max_AF >= 0.01 & INFO/Max_AF < 0.05' $input_vcf --output $vcf_out_LOW_FREQ
    """
}

workflow ANNOTATE_SVS_SVAFOTATE {
    take:
    pav_sv_vcf
    refID
    svafotate_overlap
    svafotate_bed
    nThread

    main:
    CALL_SVAFOTATE_ANNOTATE(
        refID,
        pav_sv_vcf,
        svafotate_overlap,
        svafotate_bed,
        nThread
    )

    CALL_SVAFOTATE_FILTER(
        refID,
        CALL_SVAFOTATE_ANNOTATE.out.svafotate_annotate_out
    )

    emit:
    svafotate_annotate_out = CALL_SVAFOTATE_ANNOTATE.out.svafotate_annotate_out
    svafotate_filter_out   = CALL_SVAFOTATE_FILTER.out.svafotate_filter_out
}
