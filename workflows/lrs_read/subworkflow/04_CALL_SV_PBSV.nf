/*
========================================================================================
    04_CALL_SV_PBSV
========================================================================================
    Discovers and calls structural variants from aligned PacBio reads using pbsv.

    pbsv runs in two stages:
        1. Discover -- scans aligned BAM for SV signatures, outputs .svsig.gz
        2. Call     -- calls SVs from indexed signatures, outputs VCF

    Processes:
        CALL_SVS_PBSV_DISCOVER  -- discover SV signatures from aligned BAM
        CALL_SVS_PBSV_TABIX     -- index .svsig.gz file
        CALL_SVS_PBSV_CALL      -- call SVs from indexed signatures
        CALL_SVS_PBSV_BCFTOOLS  -- filter to PASS variants and generate Emedgene VCF

    Input:
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file     -- path to reference genome FASTA
        ref_trf            -- path to tandem repeat BED file (for PBSV discover)
        refID              -- reference genome identifier (e.g. hg38)
        nThread            -- number of threads for PBSV call

    Emits:
        svsig           -- channel: [ val(meta), path(svsig.gz) ]
        svsig_indexed   -- channel: [ val(meta), path(svsig.gz), path(tbi) ]
        svs_vcf         -- channel: [ val(meta), path(pbsv.vcf) ]
        svs_filtered_vcfs -- channel: [ val(meta), path(Emedgene.vcf), path(FILT.vcf) ]
========================================================================================
*/

process CALL_SVS_PBSV_DISCOVER {
    tag "$meta.id"
    container 'quay.io/pacbio/pbsv:2.9.0_1.14_build1'

    input:
    path  ref_trf
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)

    output:
    tuple val(meta), path(svsig_out_file_name), emit: svsig

    script:
    svsig_out_file_name = meta.id + "_" + refID + "_pbmm2.svsig.gz"

    """
    pbsv discover \\
        --tandem-repeats $ref_trf \\
        $input_bam \\
        $svsig_out_file_name
    """
}

process CALL_SVS_PBSV_TABIX {
    tag "$meta.id"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    val   refID
    tuple val(meta), path(svsig_out_file_name)

    output:
    tuple val(meta), path(svsig_out_file_name), path(svsig_out_file_name_index), emit: svsig_indexed

    script:
    svsig_out_file_name_index = svsig_out_file_name + ".tbi"

    """
    tabix -c '#' -s 3 -b 4 -e 4 $svsig_out_file_name
    """
}

process CALL_SVS_PBSV_CALL {
    tag "$meta.id"
    container 'quay.io/pacbio/pbsv:2.9.0_1.14_build1'
    cpus   params.NTHREAD
    memory '64 GB'

    input:
    path  reference_file
    val   refID
    val   nThread
    tuple val(meta), path(svsig_out_file_name), path(svsig_out_file_name_index)

    output:
    tuple val(meta), path(pbsv_vcf_out_file_name), emit: svs_vcf

    script:
    pbsv_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_pbsv.vcf"

    """
    pbsv call \\
        --ccs \\
        -j $nThread \\
        $reference_file \\
        $svsig_out_file_name \\
        $pbsv_vcf_out_file_name
    """
}

process CALL_SVS_PBSV_BCFTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(pbsv_vcf_out_file_name)

    output:
    tuple val(meta), path(pbsv_Emedgene_vcf_out_file_name), path(pbsv_FILT_vcf_out_file_name), emit: svs_filtered_vcfs

    script:
    pbsv_FILT_vcf_out_file_name     = meta.id + "_" + refID + "_pbmm2_pbsv_FILT.vcf"
    pbsv_Emedgene_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_pbsv_Emedgene.vcf"

    // Note: the Emedgene VCF is a reformatted copy of the FILT VCF for ingestion
    // into the Emedgene clinical interpretation platform. The PBSV version string
    // and a header line are modified for compatibility. If not using Emedgene,
    // the FILT VCF is the primary output.
    """
    bcftools view -f PASS $pbsv_vcf_out_file_name > $pbsv_FILT_vcf_out_file_name
    cp $pbsv_FILT_vcf_out_file_name $pbsv_Emedgene_vcf_out_file_name
    sed -i '2d' $pbsv_Emedgene_vcf_out_file_name
    sed -i 's/pbsv 2.9.0 (commit v2.9.0-2-gce1559a)/pbsv 2.6.2/g' $pbsv_Emedgene_vcf_out_file_name
    """
}