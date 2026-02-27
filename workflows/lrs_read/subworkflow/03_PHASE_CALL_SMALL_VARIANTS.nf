/*
========================================================================================
    03_PHASE_CALL_SMALL_VARIANTS
========================================================================================
    Calls and phases small variants from aligned PacBio reads using a two-pass
    DeepVariant + WhatsHap strategy:

        Pass 1: DeepVariant (unphased) → WhatsHap phase → haplotag BAM
        Pass 2: DeepVariant (on haplotagged BAM) → bcftools filter

    Processes:
        CALL_SMALL_VARIANTS_DEEPVARIANT1  -- initial DeepVariant SNV/indel calling
        PHASE_SMALL_VARIANTS_WHATSHAP1    -- phase variants with WhatsHap
        PHASE_SMALL_VARIANTS_TABIX        -- index phased VCF
        PHASE_SMALL_VARIANTS_WHATSHAP2    -- haplotag BAM reads using phased VCF
        CALL_SMALL_VARIANTS_SAMTOOLS      -- index haplotagged BAM
        CALL_SMALL_VARIANTS_DEEPVARIANT2  -- second DeepVariant pass on haplotagged BAM
        CALL_SMALL_VARIANTS_BCFTOOLS      -- filter final VCF to PASS variants

    Input:
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file     -- path to reference genome FASTA
        reference_index    -- path to reference genome FASTA index (.fai)
        refID              -- reference genome identifier (e.g. hg38)
        nThread            -- number of threads for samtools

    Emits:
        deepvariant1_vcf      -- channel: [ val(meta), path(vcf.gz), path(tbi), path(html) ]
        whatshap_vcf          -- channel: [ val(meta), path(phased.vcf.gz) ]
        whatshap_vcf_indexed  -- channel: [ val(meta), path(phased.vcf.gz), path(tbi) ]
        whatshap2_bam         -- channel: [ val(meta), path(haplotagged.bam) ]
        whatshap2_bam_bai     -- channel: [ val(meta), path(haplotagged.bam), path(bai) ]
        deepvariant2_vcfs     -- channel: [ val(meta), path(vcf.gz), path(tbi), path(g.vcf.gz), path(g.tbi), path(html) ]
        deepvariant2_filt_vcfs -- channel: [ val(meta), path(FILT.vcf), path(Emedgene.vcf) ]
========================================================================================
*/

process CALL_SMALL_VARIANTS_DEEPVARIANT1 {
    tag "$meta.id"
    container 'google/deepvariant:1.6.1'
    cpus   params.MORECPU
    memory '196 GB'
    time   '24h'

    input:
    path  reference_file
    path  reference_index_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)

    output:
    tuple val(meta), path(deepvariant1_vcf_out_file_name), path(deepvariant1_vcf_tbi_out_file_name), path(deepvariant1_visual_report_out_file_name), emit: deepvariant1_vcf

    script:
    deepvariant1_vcf_out_file_name             = meta.id + "_" + refID + "_pbmm2.vcf.gz"
    deepvariant1_vcf_tbi_out_file_name         = meta.id + "_" + refID + "_pbmm2.vcf.gz.tbi"
    deepvariant1_visual_report_out_file_name   = meta.id + "_" + refID + "_pbmm2.visual_report.html"

    """
    export LD_LIBRARY_PATH=/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/usr/local/lib/python3.8/dist-packages/tensorrt_libs:/usr/local/lib/python3.8/dist-packages/nvidia/cublas/lib/
    /opt/deepvariant/bin/run_deepvariant \\
        --vcf_stats_report \\
        --model_type PACBIO \\
        --ref $reference_file \\
        --reads $input_bam \\
        --output_vcf $deepvariant1_vcf_out_file_name \\
        --num_shards ${params.MORECPU}
    """
}

process PHASE_SMALL_VARIANTS_WHATSHAP1 {
    tag "$meta.id"
    container 'quay.io/biocontainers/whatshap:2.3--py38h2494328_0'
    memory '72 GB'

    input:
    path  reference_file
    path  reference_index_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index), path(input_deepvariant1_vcf), path(input_deepvariant1_vcf_tbi), path(input_deepvariant1_visual_report)

    output:
    tuple val(meta), path(whatshap_vcf_out_file_name), emit: whatshap_vcf

    script:
    whatshap_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2.phased.vcf.gz"

    """
    whatshap phase \\
        --output $whatshap_vcf_out_file_name \\
        --reference $reference_file \\
        $input_deepvariant1_vcf \\
        $input_bam
    """
}

process PHASE_SMALL_VARIANTS_TABIX {
    tag "$meta.id"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index), path(input_deepvariant1_vcf), path(input_deepvariant1_vcf_tbi), path(input_deepvariant1_visual_report), path(whatshap_vcf_out_file_name)

    output:
    tuple val(meta), path(whatshap_vcf_out_file_name), path(whatshap_vcf_index_out_file_name), emit: whatshap_vcf_indexed

    script:
    whatshap_vcf_index_out_file_name = whatshap_vcf_out_file_name + ".tbi"

    """
    tabix -p vcf $whatshap_vcf_out_file_name
    """
}

process PHASE_SMALL_VARIANTS_WHATSHAP2 {
    tag "$meta.id"
    container 'quay.io/biocontainers/whatshap:2.3--py38h2494328_0'
    memory '72 GB'

    input:
    path  reference_file
    path  reference_index_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index), path(input_deepvariant1_vcf), path(input_deepvariant1_vcf_tbi), path(input_deepvariant1_visual_report), path(whatshap_vcf_out_file_name), path(whatshap_vcf_index_out_file_name)

    output:
    tuple val(meta), path(whatshap2_bam_out_file_name), emit: whatshap2_bam

    script:
    whatshap2_bam_out_file_name = meta.id + "_" + refID + "_pbmm2.haplotagged.bam"

    """
    whatshap haplotag \\
        --output $whatshap2_bam_out_file_name \\
        --reference $reference_file \\
        $whatshap_vcf_out_file_name \\
        $input_bam
    """
}

process CALL_SMALL_VARIANTS_SAMTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_1'
    cpus   params.NTHREAD
    memory '12 GB'

    input:
    val   nThread
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index), path(input_deepvariant1_vcf), path(input_deepvariant1_vcf_tbi), path(input_deepvariant1_visual_report), path(whatshap_vcf_out_file_name), path(whatshap_vcf_index_out_file_name), path(whatshap2_bam_out_file_name)

    output:
    tuple val(meta), path(whatshap2_bam_out_file_name), path(whatshap2_bam_bai_out_file_name), emit: whatshap2_bam_bai

    script:
    whatshap2_bam_bai_out_file_name = meta.id + "_" + refID + "_pbmm2.haplotagged.bam.bai"

    """
    samtools index -@ $nThread $whatshap2_bam_out_file_name
    """
}

process CALL_SMALL_VARIANTS_DEEPVARIANT2 {
    tag "$meta.id"
    container 'google/deepvariant:1.6.1'
    cpus   params.MORECPU
    memory '196 GB'
    time   '24h'

    input:
    path  reference_file
    path  reference_index_file
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index), path(input_deepvariant1_vcf), path(input_deepvariant1_vcf_tbi), path(input_deepvariant1_visual_report), path(whatshap_vcf_out_file_name), path(whatshap_vcf_index_out_file_name), path(whatshap2_bam_out_file_name), path(whatshap2_bam_bai_out_file_name)

    output:
    tuple val(meta), path(deepvariant2_vcf_out_file_name), path(deepvariant2_vcf_tbi_out_file_name), path(deepvariant2_g_vcf_out_file_name), path(deepvariant2_g_vcf_tbi_out_file_name), path(deepvariant2_visual_report_file_name), emit: deepvariant2_vcfs

    script:
    deepvariant2_vcf_out_file_name          = meta.id + "_" + refID + "_pbmm2_deepvariant.vcf.gz"
    deepvariant2_vcf_tbi_out_file_name      = meta.id + "_" + refID + "_pbmm2_deepvariant.vcf.gz.tbi"
    deepvariant2_g_vcf_out_file_name        = meta.id + "_" + refID + "_pbmm2_deepvariant.g.vcf.gz"
    deepvariant2_g_vcf_tbi_out_file_name    = meta.id + "_" + refID + "_pbmm2_deepvariant.g.vcf.gz.tbi"
    deepvariant2_visual_report_file_name    = meta.id + "_" + refID + "_pbmm2_deepvariant.visual_report.html"

    """
    export LD_LIBRARY_PATH=/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/usr/local/lib/python3.8/dist-packages/tensorrt_libs:/usr/local/lib/python3.8/dist-packages/nvidia/cublas/lib/
    /opt/deepvariant/bin/run_deepvariant \\
        --vcf_stats_report \\
        --model_type PACBIO \\
        --ref $reference_file \\
        --reads $whatshap2_bam_out_file_name \\
        --output_vcf $deepvariant2_vcf_out_file_name \\
        --output_gvcf $deepvariant2_g_vcf_out_file_name \\
        --num_shards ${params.MORECPU}
    """
}

process CALL_SMALL_VARIANTS_BCFTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(deepvariant2_vcf_out_file_name), path(deepvariant2_vcf_tbi_out_file_name), path(deepvariant2_g_vcf), path(deepvariant2_g_vcf_tbi), path(deepvariant2_visual_report)

    output:
    tuple val(meta), path(deepvariant2_FILT_vcf_out_file_name), path(deepvariant2_Emedgene_vcf_out_file_name), emit: deepvariant2_filt_vcfs

    script:
    deepvariant2_FILT_vcf_out_file_name     = meta.id + "_" + refID + "_pbmm2_deepvariant_FILT.vcf"
    deepvariant2_Emedgene_vcf_out_file_name = meta.id + "_" + refID + "_pbmm2_deepvariant_Emedgene.vcf"

    // Note: the Emedgene VCF is a reformatted copy of the FILT VCF for ingestion
    // into the Emedgene clinical interpretation platform. Specific header lines
    // (lines 5 and 14) and the DeepVariant version string are modified for
    // compatibility. If not using Emedgene, the FILT VCF is the primary output.
    """
    bcftools view -f PASS $deepvariant2_vcf_out_file_name > $deepvariant2_FILT_vcf_out_file_name
    cp $deepvariant2_FILT_vcf_out_file_name $deepvariant2_Emedgene_vcf_out_file_name
    sed -i '5d;14d' $deepvariant2_Emedgene_vcf_out_file_name
    sed -i 's/DeepVariant_version=1.6.0/DeepVariant_version=1.0.0/g' $deepvariant2_Emedgene_vcf_out_file_name
    """
}