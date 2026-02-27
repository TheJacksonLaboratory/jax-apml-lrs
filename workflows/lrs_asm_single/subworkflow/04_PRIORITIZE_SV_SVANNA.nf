/*
========================================================================================
    04_PRIORITIZE_SV_SVANNA
========================================================================================
    Prioritizes PAV structural variants by phenotypic relevance using SvAnna,
    then filters and ranks candidates using filter_svanna_pav.R.

    Processes:
        PRIORITIZE_BY_SVANNA          -- score PAV SVs against HPO phenotype terms
        RSCRIPT_PROCESS_SVANNA_RESULT -- filter and rank SvAnna output using
                                         filter_svanna_pav.R

    Input:
        bams                  -- channel: [ val(meta), path(hifi_fastq), path(HPO) ]
        out_small_variant_vcf -- channel: [ val(meta), path(smallvariants.vcf.gz),
                                            path(<sampleID>_hifiasm_<refID>_pav.vcf) ]
        pav_svanna_r          -- path to filter_svanna_pav.R script
        refID                 -- reference genome identifier (e.g. hg38)

    Emits:
        out_svanna_csv  -- channel: [ val(meta), path(<sampleID>_<refID>_pav_svanna.csv) ]
        out_csv         -- channel: [ val(meta), path(SvAnna_PAV_Candidates.csv) ]

    Notes:
        SvAnna requires an HPO phenotype terms file per sample (one term per line).
        The bundled SvAnna database (svanna_db_2304_hg38) is expected to be present
        within the container at /app/svanna_db_2304_hg38/.
========================================================================================
*/

process PRIORITIZE_BY_SVANNA {
    tag "$meta.id"
    container 'lrs-svanna:1.0.5'
    cpus   params.svAnna.threads
    memory '320 GB'

    input:
    tuple val(meta), path(hifi_fastq), path(HPO)
    tuple val(meta1), path(out_smallvariants_pav_gz), path(out_pav_vcf)
    val   rID

    output:
    tuple val(meta), path(out_pav_svanna_csv), emit: out_svanna_csv

    script:
    pav_svanna_csv_gz  = meta.id + "_" + rID + "_pav_svanna.csv.gz"
    out_pav_svanna_csv = meta.id + "_" + rID + "_pav_svanna.csv"
    prefix_svanna      = meta.id + "_" + rID + "_pav_svanna"

    """
    java -Xms200G -Xmx300G -XX:+UseZGC -jar /app/svanna-cli.jar prioritize \\
        -d /app/svanna_db_2304_hg38/ \\
        --vcf $out_pav_vcf \\
        \$(while read f; do echo "--phenotype-term \$f"; done < $HPO | tr '\\n' ' ') \\
        --prefix $prefix_svanna \\
        --output-format html,vcf,csv
    gzip -d $pav_svanna_csv_gz
    """
}

process RSCRIPT_PROCESS_SVANNA_RESULT {
    tag "$meta.id"
    container 'quay.io/biocontainers/r-dplyr:1.1.4--r43hf9e5401_0'

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