/*
========================================================================================
    01_ALIGN_READS_PBMM2
========================================================================================
    Aligns unaligned PacBio BAM file (merged or single) to a reference genome using
    pbmm2, then indexes the resulting aligned BAM with samtools.

    Subworkflows:
        PBMM2_ALIGN    -- aligns unaligned BAM to reference using pbmm2 (CCS preset)
        PBMM2_SAMTOOLS -- indexes the aligned BAM with samtools

    Input:
        bams           -- channel: [ val(meta), path(bam), path(hpo) ]
        reference_file -- path to reference genome FASTA
        refID          -- reference genome identifier (e.g. hg38)
        nThread        -- number of threads for alignment

    Emits:
        aligned_bam        -- channel: [ val(meta), path(bam) ]
        aligned_bam_index  -- channel: [ val(meta), path(bam), path(bai) ]
========================================================================================
*/

process PBMM2_ALIGN {
    tag "$meta.id"
    container 'quay.io/pacbio/pbmm2:1.13.1_build3'
    cpus   params.pbmm2.cpus
    memory params.pbmm2.mem
    time   '48h'

    input:
    tuple val(meta), path(bam), path(hpo)
    path  reference_file
    val   refID
    val   nThread

    output:
    tuple val(meta), path(out_file_name), emit: aligned_bam

    script:
    out_file_name = meta.id + "_" + refID + "_pbmm2.bam"

    """
    pbmm2 align \\
        --log-level INFO \\
        --strip \\
        --num-threads $nThread \\
        --sort \\
        --sort-threads ${params.pbmm2.sort_threads} \\
        --sort-memory ${params.pbmm2.sort_memory} \\
        --sample $meta.id \\
        --preset CCS \\
        $reference_file \\
        $bam \\
        $out_file_name
    """
}

process PBMM2_SAMTOOLS {
    tag "$meta.id"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam_file_name)
    val   refID

    output:
    tuple val(meta), path(bam_file_name), path(bam_index_file_name), emit: aligned_bam_index

    script:
    bam_index_file_name = meta.id + "_" + refID + "_pbmm2.bam.bai"

    """
    samtools index $bam_file_name
    """
}