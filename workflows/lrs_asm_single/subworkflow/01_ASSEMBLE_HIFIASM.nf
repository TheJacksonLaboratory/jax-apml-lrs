/*
========================================================================================
    01_ASSEMBLE_HIFIASM
========================================================================================
    Single-sample de novo genome assembly using HiFiasm, followed by FASTA
    conversion and indexing with samtools.

    Processes:
        ASSEMBLY_HIFIASM   -- assemble HiFi reads with HiFiasm (single-sample mode)
        SAMTOOLS_HIFIASM   -- index the primary contig assembly with samtools faidx

    Input:
        bams    -- channel: [ val(meta), path(hifi_fastq), path(HPO) ]
        threads -- number of threads for HiFiasm

    Emits:
        assembled_ctg     -- channel: [ val(meta), path(.fa), path(.gfa),
                                        path(hap1.gfa), path(hap2.gfa) ]
        assembled_ctg_fai -- channel: [ val(meta), path(.fa.fai) ]

    Notes:
        HiFiasm is run in single-sample mode (no parental Yak databases).
        Primary contig outputs use the bp.p_ctg prefix.
        Both haplotype GFA files are passed to PAV for phased variant calling.
========================================================================================
*/

process ASSEMBLY_HIFIASM {
    tag "$meta.id"
    container 'dnalinux/hifiasm:0.20.0'

    input:
    tuple val(meta), path(hifi_fastq), path(HPO)
    val   threads

    output:
    tuple val(meta), path(out_file_name), path(out_file_name_gfa), path(out_file_name_hap1_gfa), path(out_file_name_hap2_gfa), emit: assembled_ctg

    script:
    out_hifiasm_file_name  = meta.id + ".asm"
    out_file_name          = meta.id + ".asm.bp.p_ctg.fa"
    out_file_name_gfa      = meta.id + ".asm.bp.p_ctg.gfa"
    out_file_name_hap1_gfa = meta.id + ".asm.bp.hap1.p_ctg.gfa"
    out_file_name_hap2_gfa = meta.id + ".asm.bp.hap2.p_ctg.gfa"

    """
    hifiasm -o $out_hifiasm_file_name -t ${task.cpus} $hifi_fastq
    awk '/^S/{print ">"\$2;print \$3}' $out_file_name_gfa > $out_file_name
    """
}

process SAMTOOLS_HIFIASM {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(assembled_contigs), path(assembled_contigs_gfa), path(assembled_contigs_hap1_gfa), path(assembled_contigs_hap2_gfa)

    output:
    tuple val(meta), path(out_file_name_index), emit: assembled_ctg_fai

    script:
    out_file_name_index = meta.id + ".asm.bp.p_ctg.fa.fai"

    """
    samtools faidx $assembled_contigs
    """
}