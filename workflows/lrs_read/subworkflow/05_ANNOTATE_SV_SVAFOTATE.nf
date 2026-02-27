/*
========================================================================================
    05_ANNOTATE_SV_SVAFOTATE
========================================================================================
    Annotates and filters structural variants from PBSV, Sniffles2, and Delly
    using SVAFotate population allele frequency annotation.

    SVAFotate annotates each SV with population allele frequencies from a reference
    BED file, then variants are filtered into two frequency tiers:
        - RARE_UNIQUE : Max_AF < 0.01
        - LOW_FREQ    : 0.01 <= Max_AF < 0.05

    Processes:
        CALL_SVAFOTATE_ANNOTATE  -- annotate SVs with population allele frequencies
                                    (aliased separately for PBSV, Sniffles2, Delly)
        CALL_SVAFOTATE_FILTER    -- split annotated VCF into RARE_UNIQUE and LOW_FREQ
                                    (aliased separately for PBSV, Sniffles2, Delly)

    Input:
        svs_filtered_vcfs        -- pbsv: channel: [ val(meta), path(Emedgene.vcf), path(FILT.vcf) ]
        sniffles2_FILT_vcf       -- Sniffles2: channel: [ val(meta), path(FILT.vcf) ]
        delly_out_filt_vcf       -- Delly: channel: [ val(meta), path(telomereEX.vcf), path(FILT.vcf) ]
        refID                    -- reference genome identifier (e.g. hg38)
        svafotate_overlap        -- minimum reciprocal overlap for SVAFotate annotation
        svafotate_bed            -- path to SVAFotate population frequency BED file
        nThread                  -- number of threads for SVAFotate

    Emits:
        svafotate_annotate_out   -- channel: [ val(meta), path(SVAFotate-ALL.vcf) ]
                                    (one per caller: Emedgene, sniffles2_FILT, delly)
        svafotate_filter_out     -- channel: [ val(meta), path(RARE-UNIQUE.vcf), path(LOWFREQ.vcf) ]
                                    (one per caller: Emedgene, sniffles2_FILT, delly)

    Note on removed processes:
        CALL_PRISM_SVAFOTATE                          -- removed (unused reference process)
        CALL_SVAFOTATE_EMEDGENE_CONVERT               -- removed (superseded)
        CALL_SVAFOTATE_DELLY_EMEDGENE_CONVERT_*       -- removed (Delly is not ingested
                                                          into Emedgene)
========================================================================================
*/

process CALL_SVAFOTATE_ANNOTATE {
    tag "$meta.id"
    container 'quay.io/biocontainers/svafotate:0.2.0--pyhdfd78af_0'
    cpus   params.svafotate.cpus
    memory params.svafotate.mem

    input:
    val   refID
    tuple val(meta), path(input_vcf)
    val   suffix
    val   overlap
    path  svafotate_bed
    val   nThread

    output:
    tuple val(meta), path(vcf_out), emit: svafotate_annotate_out

    script:
    vcf_out = meta.id + "_" + refID + "_" + suffix + "_pbmm2_" + suffix + "_SVAFotate-ALL.vcf"

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
    val   suffix

    output:
    tuple val(meta), path(vcf_out_RARE_UNIQUE), path(vcf_out_LOW_FREQ), emit: svafotate_annotate_out

    script:
    vcf_out_RARE_UNIQUE = meta.id + "_" + refID + "_" + suffix + "_pbmm2_" + suffix + "_SVAFotate-RARE-UNIQUE.vcf"
    vcf_out_LOW_FREQ    = meta.id + "_" + refID + "_" + suffix + "_pbmm2_" + suffix + "_SVAFotate-LOWFREQ.vcf"

    """
    sleep 5
    bcftools view -i 'INFO/Max_AF < 0.01' $input_vcf --output $vcf_out_RARE_UNIQUE
    bcftools view -i 'INFO/Max_AF >= 0.01 & INFO/Max_AF < 0.05' $input_vcf --output $vcf_out_LOW_FREQ
    """
}