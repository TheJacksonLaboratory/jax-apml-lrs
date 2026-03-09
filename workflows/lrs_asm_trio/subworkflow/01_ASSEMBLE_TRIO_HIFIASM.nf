/*
========================================================================================
    01_ASSEMBLE_TRIO_HIFIASM
========================================================================================
    Performs trio-aware de novo genome assembly using HiFiasm with Yak k-mer
    counting from parental reads.

    Steps:
        YAK                   -- count k-mers from maternal and paternal HiFi reads
        ASSEMBLY_HIFIASM_TRIO -- assemble proband reads using parental Yak databases
        SAMTOOLS_HIFIASM_TRIO -- convert assembled GFA to FASTA and index with samtools

    Input:
        bams     -- channel: [ val(meta), path(proband_fastq), path(mat_fastq),
                               path(pat_fastq), path(HPO) ]
        threads  -- number of threads for Yak and HiFiasm

    Emits:
        out_yak          -- channel: [ val(meta), path(pat.yak), path(mat.yak) ]
        out_assembled_ctg -- channel: [ val(meta), path(asm.dip.p_utg.gfa),
                                        path(pat.yak), path(mat.yak),
                                        path(hap1.p_ctg.gfa), path(hap2.p_ctg.gfa) ]
        out_assembled_utg -- channel: [ val(meta), path(p_utg.fa), path(p_utg.fa.fai) ]
========================================================================================
*/

process YAK {
    tag "$meta.id"
    container 'quay.io/biocontainers/yak:0.1--h577a1d6_6'
    cpus   params.n_proc
    memory '128 GB'

    input:
    tuple val(meta), path(proband_hifi_fastq), path(mat_hifi_fastq), path(pat_hifi_fastq), path(HPO)
    val   threads

    output:
    tuple val(meta), path("pat.yak"), path("mat.yak"), emit: out_yak

    """
    yak count -k31 -b37 -t${threads} -o pat.yak ${pat_hifi_fastq}
    yak count -k31 -b37 -t${threads} -o mat.yak ${mat_hifi_fastq}
    """
}

process ASSEMBLY_HIFIASM_TRIO {
    tag "$meta.id"
    container 'dnalinux/hifiasm:0.20.0'
    cpus   params.hifiasm.threads
    memory '196 GB'

    input:
    tuple val(meta), path(proband_hifi_fastq), path(mat_hifi_fastq), path(pat_hifi_fastq), path(HPO)
    tuple val(meta1), path(pat_yak), path(mat_yak)
    val   threads

    output:
    tuple val(meta), path(out_asm_gfa), path("pat.yak"), path("mat.yak"), path(out_asm_hap1), path(out_asm_hap2), emit: out_assembled_ctg

    script:
    out_hifiasm_prefix = meta.id + ".asm"
    out_asm_gfa        = meta.id + ".asm.dip.p_utg.gfa"
    out_asm_hap1       = meta.id + ".asm.dip.hap1.p_ctg.gfa"
    out_asm_hap2       = meta.id + ".asm.dip.hap2.p_ctg.gfa"

    """
    hifiasm \\
        -o $out_hifiasm_prefix \\
        -t $threads \\
        -1 pat.yak \\
        -2 mat.yak \\
        $proband_hifi_fastq
    """
}

process SAMTOOLS_HIFIASM_TRIO {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(out_asm_gfa), path(pat_yak), path(mat_yak), path(out_asm_hap1), path(out_asm_hap2)

    output:
    tuple val(meta), path(out_utg_fa), path(out_utg_fa_fai), emit: out_assembled_utg

    script:
    out_utg_fa     = meta.id + ".asm.dip.p_utg.fa"
    out_utg_fa_fai = meta.id + ".asm.dip.p_utg.fa.fai"

    """
    awk '/^S/{print ">"\$2; print \$3}' $out_asm_gfa > $out_utg_fa
    samtools faidx $out_utg_fa
    """
}