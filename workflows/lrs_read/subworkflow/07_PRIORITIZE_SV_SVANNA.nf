/*
========================================================================================
    08_PRIORITIZE_SV_SVANNA
========================================================================================
    Prioritizes structural variants by phenotypic relevance using SvAnna, then
    merges and ranks candidates from pbsv and Sniffles2 callers.

    Two SvAnna runs are performed:
        1. Unfiltered: runs on raw pbsv and Sniffles2 FILT VCFs
        2. Filtered:   runs on SVAFotate RARE+UNIQUE pbsv and Sniffles2 VCFs

    The summary CSV (SvAnna_SV_Candidates.csv) is generated from the filtered run.

    Processes:
        PRIORITIZE_SVS_BY_SVANNA          -- run SvAnna on unfiltered pbsv and Sniffles2 VCFs
        PRIORITIZE_SVS_BY_SVANNA_FILTERED -- run SvAnna on SVAFotate RARE+UNIQUE VCFs
        COMBINE_SVS_BY_SVANNA             -- merge and rank filtered SvAnna candidates
                                             using filter_svanna_pbsv_sniffles2.R

    Input:
        svs_filtered_vcfs          -- pbsv channel: [ val(meta), path(Emedgene.vcf), path(FILT.vcf) ]
        sniffles2_FILT_vcf         -- Sniffles2 channel: [ val(meta), path(FILT.vcf) ]
        svafotate_pbsv_rare_unique -- channel: [ val(meta), path(pbsv_SVAFotate-RARE-UNIQUE.vcf) ]
        sniffles2_rare_unique      -- channel: [ val(meta), path(sniffles2_SVAFotate-RARE-UNIQUE.vcf) ]
        bams                       -- channel: [ val(meta), path(bam), path(hpo) ]
        pbsv_sniffles2_svanna_r    -- path to filter_svanna_pbsv_sniffles2.R script
        refID                      -- reference genome identifier (e.g. hg38)

    Emits:
        csv_files           -- channel: unfiltered svanna outputs (pbsv + sniffles2)
        csv_files_filtered  -- channel: filtered svanna outputs (pbsv + sniffles2 RARE+UNIQUE)
        sv_candidates_file  -- channel: [ val(meta), path(SvAnna_SV_Candidates.csv) ]

    Notes:
        SvAnna requires an HPO phenotype terms file per sample (one term per line).
        This file is passed through from the input samplesheet via the bams channel.
        The bundled SvAnna database (svanna_db_2304_hg38) is expected to be present
        within the container at /app/svanna_db_2304_hg38/.
========================================================================================
*/

process PRIORITIZE_SVS_BY_SVANNA {
    tag "$meta.id"
    container params.svAnna.svanna_sif_path

    input:
    tuple val(meta), path(pbsv_filt_output), path(pbsv_Emedgene_output), path(sniffles2_filt_output), path(bam), path(HPO_file)
    val   refID

    output:
    tuple val(meta), path(out_pbsv_svanna_file), path(out_pbsv_svanna_file_html), path(out_pbsv_svanna_file_vcf_gz), val(meta), path(out_sniffles2_svanna_file), path(out_sniffles2_svanna_file_html), path(out_sniffles2_svanna_file_vcf_gz), emit: csv_files

    script:
    out_pbsv_svanna_gz_file           = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.csv.gz"
    out_pbsv_svanna_file              = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.csv"
    out_pbsv_svanna_file_html         = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.html"
    out_pbsv_svanna_file_vcf_gz       = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.vcf.gz"
    out_sniffles2_svanna_gz_file      = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.csv.gz"
    out_sniffles2_svanna_file         = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.csv"
    out_sniffles2_svanna_file_html    = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.html"
    out_sniffles2_svanna_file_vcf_gz  = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.vcf.gz"
    out_pbsv_svanna_prefix            = meta.id + "_" + refID + "_pbmm2_pbsv_svanna"
    out_sniffles2_svanna_prefix       = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna"

    """
    java -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $pbsv_filt_output \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO_file | tr '\\n' ' ') \\
        --out-dir . \\
        --prefix $out_pbsv_svanna_prefix \\
        --output-format html,vcf,csv
    gzip -d $out_pbsv_svanna_gz_file

    java -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $sniffles2_filt_output \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO_file | tr '\\n' ' ') \\
        --out-dir . \\
        --prefix $out_sniffles2_svanna_prefix \\
        --output-format html,vcf,csv
    gzip -d $out_sniffles2_svanna_gz_file
    """
}

process PRIORITIZE_SVS_BY_SVANNA_FILTERED {
    tag "$meta.id"
    container params.svAnna.svanna_sif_path

    input:
    tuple val(meta), path(pbsv_rare_unique_vcf), path(pbsv_lowfreq_vcf), path(sniffles2_rare_unique_vcf), path(sniffles2_lowfreq_vcf), path(bam), path(HPO_file)
    val   refID

    output:
    tuple val(meta), path(out_pbsv_svanna_file), path(out_pbsv_svanna_file_html), path(out_pbsv_svanna_file_vcf_gz), val(meta), path(out_sniffles2_svanna_file), path(out_sniffles2_svanna_file_html), path(out_sniffles2_svanna_file_vcf_gz), emit: csv_files_filtered

    script:
    out_pbsv_svanna_gz_file           = meta.id + "_" + refID + "_pbmm2_pbsv_SVAFotate-RARE-UNIQUE_svanna.csv.gz"
    out_pbsv_svanna_file              = meta.id + "_" + refID + "_pbmm2_pbsv_SVAFotate-RARE-UNIQUE_svanna.csv"
    out_pbsv_svanna_file_html         = meta.id + "_" + refID + "_pbmm2_pbsv_SVAFotate-RARE-UNIQUE_svanna.html"
    out_pbsv_svanna_file_vcf_gz       = meta.id + "_" + refID + "_pbmm2_pbsv_SVAFotate-RARE-UNIQUE_svanna.vcf.gz"
    out_sniffles2_svanna_gz_file      = meta.id + "_" + refID + "_pbmm2_sniffles2_SVAFotate-RARE-UNIQUE_svanna.csv.gz"
    out_sniffles2_svanna_file         = meta.id + "_" + refID + "_pbmm2_sniffles2_SVAFotate-RARE-UNIQUE_svanna.csv"
    out_sniffles2_svanna_file_html    = meta.id + "_" + refID + "_pbmm2_sniffles2_SVAFotate-RARE-UNIQUE_svanna.html"
    out_sniffles2_svanna_file_vcf_gz  = meta.id + "_" + refID + "_pbmm2_sniffles2_SVAFotate-RARE-UNIQUE_svanna.vcf.gz"
    out_pbsv_svanna_prefix            = meta.id + "_" + refID + "_pbmm2_pbsv_SVAFotate-RARE-UNIQUE_svanna"
    out_sniffles2_svanna_prefix       = meta.id + "_" + refID + "_pbmm2_sniffles2_SVAFotate-RARE-UNIQUE_svanna"

    """
    java -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $pbsv_rare_unique_vcf \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO_file | tr '\\n' ' ') \\
        --out-dir . \\
        --prefix $out_pbsv_svanna_prefix \\
        --output-format html,vcf,csv
    gzip -d $out_pbsv_svanna_gz_file

    java -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $sniffles2_rare_unique_vcf \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO_file | tr '\\n' ' ') \\
        --out-dir . \\
        --prefix $out_sniffles2_svanna_prefix \\
        --output-format html,vcf,csv
    gzip -d $out_sniffles2_svanna_gz_file
    """
}

process COMBINE_SVS_BY_SVANNA {
    tag "$meta.id"
    container 'rocker/tidyverse:4.3'

    input:
    path  pbsv_sniffles2_svanna_r
    val   refID
    tuple val(meta), path(out_pbsv_svanna_file), path(out_pbsv_svanna_file_html), path(out_pbsv_svanna_file_vcf_gz), val(meta1), path(out_sniffles2_svanna_file), path(out_sniffles2_svanna_file_html), path(out_sniffles2_svanna_file_vcf_gz)

    output:
    tuple val(meta), path(out_sv_candidates_file), emit: sv_candidates_file

    script:
    out_sv_candidates_file = "SvAnna_SV_Candidates.csv"

    """
    Rscript $pbsv_sniffles2_svanna_r $out_pbsv_svanna_file $out_sniffles2_svanna_file
    """
}
