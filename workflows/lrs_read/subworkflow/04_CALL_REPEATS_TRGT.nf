/*
========================================================================================
    04_CALL_REPEATS_TRGT
========================================================================================
    Genotypes tandem repeats from aligned PacBio reads using TRGT, then sorts,
    indexes, and reformats outputs for downstream analysis.

    Processes:
        CALL_TRGT_GENOTYPE       -- genotype tandem repeats with TRGT
        CALL_TRGT_BCFTOOLS       -- sort output VCF
        CALL_TRGT_TABIX          -- index sorted VCF
        CALL_TRGT_SAMTOOLS       -- sort and index spanning BAM
        CALL_TRGT_GATK           -- convert VCF to TSV with VariantsToTable
        CALL_TRGT_BCFTOOLS_FINAL -- reformat TSV for readability

    Input:
        aligned_bam_index   -- channel: [ val(meta), path(bam), path(bai) ]
        reference_file      -- path to reference genome FASTA
        reference_index     -- path to reference genome FASTA index (.fai)
        repeats_path        -- path to medically relevant repeats BED file
        refID               -- reference genome identifier (e.g. hg38)
        nThread             -- number of threads for TRGT genotype

    Emits:
        trgt_bam_vcf        -- channel: [ val(meta), path(vcf.gz), path(bam) ]
        trgt_sorted_vcf     -- channel: [ val(meta), path(sorted.vcf.gz) ]
        trgt_vcf_index      -- channel: [ val(meta), path(sorted.vcf.gz), path(tbi) ]
        trgt_bam_sorted     -- channel: [ val(meta), path(sorted.bam), path(bai) ]
        trgt_sorted_tsv     -- channel: [ val(meta), path(vcf.gz), path(tbi), path(tsv) ]
        new_sorted_tsv_trgt_file -- channel: [ val(meta), path(final.tsv) ]
========================================================================================
*/

process CALL_TRGT_GENOTYPE {
    tag "$meta.id"
    container 'quay.io/pacbio/trgt:1.0.0_build1'

    input:
    path  reference_file
    path  reference_index_file
    path  repeats
    val   nThread
    val   refID
    tuple val(meta), path(input_bam), path(input_bam_index)

    output:
    tuple val(meta), path(out_trgt_vcf), path(out_trgt_bam), emit: trgt_bam_vcf

    script:
    out_prefix   = meta.id + "_" + refID + "_pbmm2_TRGT"
    out_trgt_vcf = meta.id + "_" + refID + "_pbmm2_TRGT.vcf.gz"
    out_trgt_bam = meta.id + "_" + refID + "_pbmm2_TRGT.spanning.bam"

    """
    trgt genotype \\
        --genome $reference_file \\
        --repeats $repeats \\
        --reads $input_bam \\
        --output-prefix $out_prefix \\
        --threads $nThread
    """
}

process CALL_TRGT_BCFTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(trgt_vcf_out_file_name), path(trgt_bam_out_file_name)

    output:
    tuple val(meta), path(out_trgt_vcf_sorted_fileName), emit: trgt_sorted_vcf

    script:
    out_trgt_vcf_sorted_fileName = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.vcf.gz"

    """
    bcftools sort -Ob -o $out_trgt_vcf_sorted_fileName $trgt_vcf_out_file_name
    """
}

process CALL_TRGT_TABIX {
    tag "$meta.id"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    val   refID
    tuple val(meta), path(out_trgt_vcf_sorted_fileName)

    output:
    tuple val(meta), path(out_trgt_vcf_sorted_fileName), path(out_trgt_vcf_index_sorted_fileName), emit: trgt_vcf_index

    script:
    out_trgt_vcf_index_sorted_fileName = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.vcf.gz.tbi"

    """
    tabix $out_trgt_vcf_sorted_fileName
    """
}

process CALL_TRGT_SAMTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_1'

    input:
    val   refID
    tuple val(meta), path(trgt_vcf_out_file_name), path(trgt_bam_out_file_name)

    output:
    tuple val(meta), path(out_trgt_bam_sorted_fileName), path(out_trgt_bam_bai_sorted_fileName), emit: trgt_bam_sorted

    script:
    out_trgt_bam_sorted_fileName     = meta.id + "_" + refID + "_pbmm2_TRGT.spanning.sorted.bam"
    out_trgt_bam_bai_sorted_fileName = meta.id + "_" + refID + "_pbmm2_TRGT.spanning.sorted.bam.bai"

    """
    samtools sort -o $out_trgt_bam_sorted_fileName $trgt_bam_out_file_name
    samtools index $out_trgt_bam_sorted_fileName
    """
}

process CALL_TRGT_GATK {
    tag "$meta.id"
    container 'quay.io/biocontainers/gatk4:4.2.6.1--py36hdfd78af_1'

    input:
    val   refID
    tuple val(meta), path(out_trgt_vcf_sorted_fileName), path(out_trgt_vcf_index_sorted_fileName)

    output:
    tuple val(meta), path(out_trgt_vcf_sorted_fileName), path(out_trgt_vcf_index_sorted_fileName), path(out_trgt_tsv_sorted_fileName), emit: trgt_sorted_tsv

    script:
    out_trgt_tsv_sorted_fileName = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.tsv"

    // Note: '-F INF0' in the original was a typo (zero instead of letter O).
    // Corrected to '-F INFO' here.
    """
    gatk VariantsToTable \\
        -V $out_trgt_vcf_sorted_fileName \\
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO \\
        -GF GT -GF AL -GF ALLR -GF SD -GF MC -GF MS -GF AP -GF AM \\
        -O $out_trgt_tsv_sorted_fileName
    """
}

process CALL_TRGT_BCFTOOLS_FINAL {
    tag "$meta.id"
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    val   refID
    tuple val(meta), path(out_trgt_vcf_sorted_fileName), path(out_trgt_vcf_index_sorted_fileName), path(out_trgt_tsv_sorted_fileName)

    output:
    tuple val(meta), path(sorted_tsv_trgt_file_name), emit: new_sorted_tsv_trgt_file

    script:
    out_trgt_vcf_fileName           = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.vcf"
    tempGT_trgt_file_name           = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.tempGT"
    sorted_tsv_temp_trgt_file_name  = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.tsv.temp"
    sorted_tsv_trgt_file_name       = meta.id + "_" + refID + "_pbmm2_TRGT.sorted.tsv"

    """
    bcftools view $out_trgt_vcf_sorted_fileName > $out_trgt_vcf_fileName
    grep -v ^## $out_trgt_vcf_fileName | cut -f10- | sed 's/:\\S*//g' | sed 's/UnnamedSample/GT/g' > tempGT
    awk -v FS='\\t' -v OFS='\\t' 'FNR==NR{a[NR]=\$1;next}{\$9=a[FNR]}1' tempGT $out_trgt_tsv_sorted_fileName > $tempGT_trgt_file_name
    awk '{print \$10}' $tempGT_trgt_file_name | awk -F "," 'BEGIN{OFS="\\t"}NR==1{print "AL_1","AL_2";next}{print \$1,\$2}' > tempAL
    awk -v FS='\\t' -v OFS='\\t' 'FNR==NR{a[NR]=\$1"\\t"\$2;next}{\$10=a[FNR]}1' tempAL $tempGT_trgt_file_name > $sorted_tsv_temp_trgt_file_name
    sed 's/UnnamedSample.//g' $sorted_tsv_temp_trgt_file_name > $sorted_tsv_trgt_file_name
    """
}