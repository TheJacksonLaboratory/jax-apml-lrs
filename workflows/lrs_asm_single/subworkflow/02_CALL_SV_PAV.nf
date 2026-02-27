/*
========================================================================================
    02_CALL_SV_PAV
========================================================================================
    Calls structural and small variants from single-sample assembled haplotypes
    using PAV (Phased Assembly Variant caller).

    PAV takes the assembled haplotype GFA files from HiFiasm and aligns them to
    the reference to call variants across the full size spectrum.

    Processes:
        CALL_SVS_PAV  -- run PAV variant calling on assembled haplotypes

    Input:
        assembled_ctg      -- channel: [ val(meta), path(bp.p_ctg.fa),
                                         path(bp.p_ctg.gfa), path(hap1.p_ctg.gfa),
                                         path(hap2.p_ctg.gfa) ]
        assembled_ctg_fai  -- channel: [ val(meta), path(bp.p_ctg.fa.fai) ]
        reference_file     -- path to reference genome FASTA
        reference_index    -- path to reference genome FASTA index (.fai)
        n_proc             -- number of threads for PAV
        refID              -- reference genome identifier (e.g. hg38)

    Emits:
        out_vcf  -- channel: [ val(meta), path(<sampleID>_hifiasm_<refID>_pav.vcf.gz),
                               path(<sampleID>_hifiasm_<refID>_pav.vcf.gz.tbi) ]
========================================================================================
*/

process CALL_SVS_PAV {
    tag "$meta.id"
    container 'becklab/pav:2.4.2'
    cpus   params.n_proc
    memory '128 GB'

    input:
    tuple val(meta),  path(assembled_contigs), path(assembled_contigs_gfa), path(assembled_contigs_hap1_gfa), path(assembled_contigs_hap2_gfa)
    tuple val(meta1), path(assembled_contigs_index)
    path  reference_file
    path  reference_index_file
    val   n_proc
    val   rID

    output:
    tuple val(meta), path(output_pav_vcf), path(output_pav_vcf_tbi), emit: out_vcf

    script:
    input_pav_vcf    = meta.id + ".vcf.gz"
    input_pav_vcf_tbi = meta.id + ".vcf.gz.tbi"
    output_pav_vcf   = meta.id + "_hifiasm_" + rID + "_pav.vcf.gz"
    output_pav_vcf_tbi = meta.id + "_hifiasm_" + rID + "_pav.vcf.gz.tbi"

    """
    echo '{"reference": "$reference_file"}' > config.json
    echo -e "NAME\\tHAP1\\tHAP2\\n${meta.id}\\t${assembled_contigs_hap1_gfa}\\t${assembled_contigs_hap2_gfa}" > assemblies.tsv

    /opt/pav/files/docker/run -c $n_proc --nt

    mv $input_pav_vcf $output_pav_vcf
    mv $input_pav_vcf_tbi $output_pav_vcf_tbi
    """
}