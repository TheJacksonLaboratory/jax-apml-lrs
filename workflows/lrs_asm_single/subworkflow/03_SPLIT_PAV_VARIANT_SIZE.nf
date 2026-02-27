/*
========================================================================================
    03_SPLIT_PAV_VARIANT_SIZE
========================================================================================
    Splits, indexes, and merges PAV variant output by size into two tracks:
        - Small variants (SNVs and indels < 50 bp)
        - Structural variants (insertions and deletions >= 50 bp)

    Both tracks take the same PAV VCF as input and run independently.

    Small variant track:
        SPLIT_SMALL_VARIANTS_PAV  -- split into SNV and indel subsets
        INDEX_SMALL_VARIANTS_PAV  -- bgzip and tabix index
        CONCAT_SMALL_VARIANTS_PAV -- merge SNVs + indels < 50 bp into one VCF

    Structural variant track:
        SPLIT_SV_PAV   -- split into insertions (>= 50 bp) and deletions (<= -50 bp)
        INDEX_SV_PAV   -- bgzip and tabix index
        CONCAT_SV_PAV  -- merge SV insertions and deletions into one VCF

    Input:
        out_vcf  -- channel: [ val(meta), path(<sampleID>_hifiasm_<refID>_pav.vcf.gz),
                               path(<sampleID>_hifiasm_<refID>_pav.vcf.gz.tbi) ]
        refID    -- reference genome identifier (e.g. hg38)

    Emits:
        out_small_variant_vcf  -- channel: [ val(meta),
                                             path(<sampleID>_hifiasm_<refID>_pav_smallvariants.vcf.gz),
                                             path(<sampleID>_hifiasm_<refID>_pav.vcf) ]
        out_concat_sv_vcf      -- channel: [ val(meta),
                                             path(<sampleID>_hifiasm_<refID>_pav_SV.vcf) ]

    Note on Emedgene outputs:
        The PAV version string in SV VCFs is modified (PAV 2.2.4 → PAVSV 2.2.4) for
        compatibility with the Emedgene clinical interpretation platform. If not using
        Emedgene, the sed commands in SPLIT_SV_PAV can be removed.
========================================================================================
*/

// ---------------------------------------------------------------------------
// Small variant track (SNVs + indels < 50 bp)
// ---------------------------------------------------------------------------

process SPLIT_SMALL_VARIANTS_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple val(meta), path(pav_vcf_gz_rid), path(pav_vcf_gz_tbi_rid)

    output:
    tuple val(meta), path("pav_snvs.vcf"), path("pav_indels.vcf"), path("pav_indels50bp.vcf"), emit: out_splitted_vcf

    """
    bcftools view -v snps $pav_vcf_gz_rid > pav_snvs.vcf
    bcftools view -v indels $pav_vcf_gz_rid > pav_indels.vcf
    bcftools view -i 'INFO/SVLEN > -50 && INFO/SVLEN < 50' pav_indels.vcf > pav_indels50bp.vcf
    """
}

process INDEX_SMALL_VARIANTS_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    tuple val(meta), path(pav_snvs_vcf), path(pav_indels_vcf), path(pav_indels50bp_vcf)

    output:
    tuple val(meta), path(pav_snvs_vcf_gz), path(pav_snvs_vcf_gz_tbi), path(pav_indels50bp_vcf_gz), path(pav_indels50bp_vcf_gz_tbi), emit: out_small_variant_vcf_indexed

    script:
    pav_snvs_vcf_gz           = "pav_snvs.vcf.gz"
    pav_snvs_vcf_gz_tbi       = "pav_snvs.vcf.gz.tbi"
    pav_indels50bp_vcf_gz     = "pav_indels50bp.vcf.gz"
    pav_indels50bp_vcf_gz_tbi = "pav_indels50bp.vcf.gz.tbi"

    """
    bgzip $pav_snvs_vcf
    bgzip $pav_indels50bp_vcf
    tabix pav_snvs.vcf.gz
    tabix pav_indels50bp.vcf.gz
    """
}

process CONCAT_SMALL_VARIANTS_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple val(meta), path(pav_snvs_vcf_gz), path(pav_snvs_vcf_gz_tbi), path(pav_indels50bp_vcf_gz), path(pav_indels50bp_vcf_gz_tbi)
    tuple val(meta1), path(pav_vcf_gz_rid), path(pav_vcf_gz_tbi_rid)
    val   rID

    output:
    tuple val(meta), path(out_smallvariants_pav_gz), path(out_pav_vcf), emit: out_small_variant_vcf

    script:
    smallvariants_pav        = meta.id + "_hifiasm_" + rID + "_pav_smallvariants.vcf"
    out_smallvariants_pav_gz = meta.id + "_hifiasm_" + rID + "_pav_smallvariants.vcf.gz"
    out_pav_vcf              = meta.id + "_hifiasm_" + rID + "_pav.vcf"

    """
    bcftools concat -a $pav_snvs_vcf_gz $pav_indels50bp_vcf_gz > $smallvariants_pav
    gzip $smallvariants_pav
    gzip -d $pav_vcf_gz_rid
    """
}

// ---------------------------------------------------------------------------
// Structural variant track (insertions + deletions >= 50 bp)
// ---------------------------------------------------------------------------

process SPLIT_SV_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple val(meta), path(pav_vcf_gz_rid), path(pav_vcf_gz_tbi_rid)

    output:
    tuple val(meta), path("pav_SVins.vcf"), path("pav_SVdel.vcf"), emit: out_splitted_vcf

    // Note: PAV version string is modified (PAV 2.2.4 → PAVSV 2.2.4) for
    // Emedgene ingestion compatibility. If not using Emedgene, these sed
    // commands can be removed.
    """
    bcftools view -i 'INFO/SVLEN >= 50' $pav_vcf_gz_rid > pav_SVins.vcf
    bcftools view -i 'INFO/SVLEN <= -50' $pav_vcf_gz_rid > pav_SVdel.vcf
    sed -i 's/PAV 2.2.4/PAVSV 2.2.4/g' pav_SVins.vcf
    sed -i 's/PAV 2.2.4/PAVSV 2.2.4/g' pav_SVdel.vcf
    """
}

process INDEX_SV_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    tuple val(meta), path(pav_svins_vcf), path(pav_svdel_vcf)

    output:
    tuple val(meta), path("pav_SVins.vcf.gz"), path("pav_SVdel.vcf.gz"), path("pav_SVins.vcf.gz.tbi"), path("pav_SVdel.vcf.gz.tbi"), emit: out_ins_del_vcf_gz

    """
    bgzip $pav_svins_vcf
    bgzip $pav_svdel_vcf
    tabix pav_SVins.vcf.gz
    tabix pav_SVdel.vcf.gz
    """
}

process CONCAT_SV_PAV {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple val(meta), path(pav_svins_vcf_gz), path(pav_svdel_vcf_gz), path(pav_svins_vcf_gz_tbi), path(pav_svdel_vcf_gz_tbi)
    val   rID

    output:
    tuple val(meta), path(out_concat_sv_vcf), emit: out_concat_sv_vcf

    script:
    out_concat_sv_vcf = meta.id + "_hifiasm_" + rID + "_pav_SV.vcf"

    """
    bcftools concat -a $pav_svins_vcf_gz $pav_svdel_vcf_gz > $out_concat_sv_vcf
    """
}