/*
========================================================================================
    08_PRIORITIZE_SV_SVANNA
========================================================================================
    Prioritizes structural variants by phenotypic relevance using SvAnna, then
    merges and ranks candidates from pbsv and Sniffles2 callers.

    SvAnna scores SVs against HPO phenotype terms provided per sample, enabling
    phenotype-driven prioritization of clinically relevant variants.

    Processes:
        PRIORITIZE_SVS_BY_SVANNA  -- run SvAnna on PBSV and Sniffles2 FILT VCFs
                                     using per-sample HPO terms
        COMBINE_SVS_BY_SVANNA     -- merge and rank SvAnna candidates from both
                                     callers using filter_svanna_pbsv_sniffles2.R

    Input:
        svs_filtered_vcfs    -- pbsv channel: [ val(meta), path(Emedgene.vcf), path(FILT.vcf) ]
        sniffles2_FILT_vcf   -- Sniffles2 channel: [ val(meta), path(FILT.vcf) ]
        bams                 -- channel: [ val(meta), path(bam), path(hpo) ]
        pbsv_sniffles2_svanna_r -- path to filter_svanna_pbsv_sniffles2.R script
        refID                -- reference genome identifier (e.g. hg38)

    Emits:
        csv_files           -- channel: [ val(meta), path(pbsv_svanna.csv),
                                          path(pbsv_svanna.html), path(pbsv_svanna.vcf.gz),
                                          val(meta), path(sniffles2_svanna.csv),
                                          path(sniffles2_svanna.html), path(sniffles2_svanna.vcf.gz) ]
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
    container 'quay.io/biocontainers/svanna:1.0.4--hdfd78af_0'

    input:
    tuple val(meta), path(pbsv_filt_output), path(pbsv_Emedgene_output), path(sniffles2_filt_output), path(bam), path(HPO_file)
    val   refID

    output:
    tuple val(meta), path(out_pbsv_svanna_file), path(out_pbsv_svanna_file_html), path(out_pbsv_svanna_file_vcf_gz), val(meta), path(out_sniffles2_svanna_file), path(out_sniffles2_svanna_file_html), path(out_sniffles2_svanna_file_vcf_gz), emit: csv_files

    script:
    out_pbsv_svanna_gz_file       = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.csv.gz"
    out_pbsv_svanna_file          = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.csv"
    out_pbsv_svanna_file_html     = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.html"
    out_pbsv_svanna_file_vcf_gz   = meta.id + "_" + refID + "_pbmm2_pbsv_svanna.vcf.gz"
    out_sniffles2_svanna_gz_file      = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.csv.gz"
    out_sniffles2_svanna_file         = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.csv"
    out_sniffles2_svanna_file_html    = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.html"
    out_sniffles2_svanna_file_vcf_gz  = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna.vcf.gz"
    out_pbsv_svanna_prefix        = meta.id + "_" + refID + "_pbmm2_pbsv_svanna"
    out_sniffles2_svanna_prefix   = meta.id + "_" + refID + "_pbmm2_sniffles2_svanna"

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

process COMBINE_SVS_BY_SVANNA {
    tag "$meta.id"
    container 'quay.io/biocontainers/r-dplyr:1.1.4--r43hf9e5401_0'

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