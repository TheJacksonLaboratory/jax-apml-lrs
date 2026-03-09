/*
========================================================================================
    04_PRIORITIZE_SV_SVANNA
========================================================================================
    Prioritizes PAV structural variants by phenotypic relevance using SvAnna,
    then filters and ranks candidates using filter_svanna_pav.R.

    Two SvAnna runs are performed:
        1. Unfiltered: runs on the full PAV SV VCF (pav_SV_Emedgene.vcf)
        2. Filtered:   runs on SVAFotate RARE+UNIQUE variants only

    The summary CSV (SvAnna_PAV_Candidates.csv) is generated from the filtered run.

    Processes:
        PRIORITIZE_BY_SVANNA_UNFILTERED -- score all PAV SVs against HPO phenotype terms
        PRIORITIZE_BY_SVANNA            -- score RARE+UNIQUE PAV SVs against HPO phenotype terms
        RSCRIPT_PROCESS_SVANNA_RESULT   -- filter and rank SvAnna output using
                                           filter_svanna_pav.R

    Input:
        bams                    -- channel: [ val(meta), path(hifi_fastq), path(HPO) ]
        pav_sv_emedgene_vcf     -- channel: [ val(meta), path(pav_SV_Emedgene.vcf) ]
        svafotate_filter_out    -- channel: [ val(meta), path(RARE-UNIQUE.vcf), path(LOWFREQ.vcf) ]
        pav_svanna_r            -- path to filter_svanna_pav.R script
        refID                   -- reference genome identifier (e.g. hg38)

    Emits:
        out_svanna_unfiltered_csv -- channel: [ val(meta), path(<sampleID>_<refID>_pav_svanna.csv),
                                                path(<sampleID>_<refID>_pav_svanna.html),
                                                path(<sampleID>_<refID>_pav_svanna.vcf.gz) ]
        out_svanna_csv            -- channel: [ val(meta), path(<sampleID>_<refID>_pav_SVAFotate-RARE-UNIQUE_svanna.csv),
                                                path(<sampleID>_<refID>_pav_SVAFotate-RARE-UNIQUE_svanna.html),
                                                path(<sampleID>_<refID>_pav_SVAFotate-RARE-UNIQUE_svanna.vcf.gz) ]
        out_csv                   -- channel: [ val(meta), path(SvAnna_PAV_Candidates.csv) ]

    Notes:
        SvAnna requires an HPO phenotype terms file per sample (one term per line).
        The bundled SvAnna database (svanna_db_2304_hg38) is expected to be present
        within the container at /app/svanna_db_2304_hg38/.
========================================================================================
*/

process PRIORITIZE_BY_SVANNA_UNFILTERED {
    tag "$meta.id"
    container params.svAnna.svanna_sif_path
    cpus   params.svAnna.threads
    memory '400 GB'

    input:
    tuple val(meta), path(hifi_fastq), path(HPO)
    tuple val(meta1), path(pav_sv_emedgene_vcf)
    val   rID

    output:
    tuple val(meta), path(out_pav_svanna_csv), path(out_pav_svanna_html), path(out_pav_svanna_vcf_gz), emit: out_svanna_unfiltered_csv

    script:
    pav_svanna_csv_gz  = meta.id + "_" + rID + "_pav_svanna.csv.gz"
    out_pav_svanna_csv = meta.id + "_" + rID + "_pav_svanna.csv"
    out_pav_svanna_html = meta.id + "_" + rID + "_pav_svanna.html"
    out_pav_svanna_vcf_gz = meta.id + "_" + rID + "_pav_svanna.vcf.gz"
    prefix_svanna      = meta.id + "_" + rID + "_pav_svanna"

    """
    java -Xms200G -Xmx300G -XX:+UseZGC -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $pav_sv_emedgene_vcf \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO | tr '\\n' ' ') \\
        --prefix $prefix_svanna \\
        --output-format html,vcf,csv
    gzip -d $pav_svanna_csv_gz
    """
}

process PRIORITIZE_BY_SVANNA {
    tag "$meta.id"
    container params.svAnna.svanna_sif_path
    cpus   params.svAnna.threads
    memory '320 GB'

    input:
    tuple val(meta), path(hifi_fastq), path(HPO)
    tuple val(meta1), path(rare_unique_vcf), path(lowfreq_vcf)
    val   rID

    output:
    tuple val(meta), path(out_pav_svanna_csv), path(out_pav_svanna_html), path(out_pav_svanna_vcf_gz), emit: out_svanna_csv

    script:
    pav_svanna_csv_gz  = meta.id + "_" + rID + "_pav_SVAFotate-RARE-UNIQUE_svanna.csv.gz"
    out_pav_svanna_csv = meta.id + "_" + rID + "_pav_SVAFotate-RARE-UNIQUE_svanna.csv"
    out_pav_svanna_html = meta.id + "_" + rID + "_pav_SVAFotate-RARE-UNIQUE_svanna.html"
    out_pav_svanna_vcf_gz = meta.id + "_" + rID + "_pav_SVAFotate-RARE-UNIQUE_svanna.vcf.gz"
    prefix_svanna      = meta.id + "_" + rID + "_pav_SVAFotate-RARE-UNIQUE_svanna"

    """
    java -Xms200G -Xmx300G -XX:+UseZGC -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $rare_unique_vcf \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO | tr '\\n' ' ') \\
        --prefix $prefix_svanna \\
        --output-format html,vcf,csv
    gzip -d $pav_svanna_csv_gz
    """
}

process RSCRIPT_PROCESS_SVANNA_RESULT {
    tag "$meta.id"
    container 'rocker/tidyverse:4.3'

    input:
    tuple val(meta), path(out_pav_svanna_csv)
    path  pav_svanna_r

    output:
    tuple val(meta), path("SvAnna_PAV_Candidates.csv"), emit: out_csv

    script:
    """
    Rscript $pav_svanna_r $out_pav_svanna_csv
    """
}
