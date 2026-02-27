/*
========================================================================================
    02_CALL_SV_PAV
========================================================================================
    Calls structural and small variants from trio-assembled haplotypes using 
    Phased Assembly Variant (PAV) caller.

    PAV takes the assembled haplotype GFA files from HiFiasm and the unitig FASTA
    as input, aligns them to the reference, and calls variants across the full
    size spectrum.

    Processes:
        CALL_SVS_PAV  -- run PAV variant calling on assembled haplotypes

    Input:
        out_assembled_ctg  -- channel: [ val(meta), path(p_utg.gfa), path(pat.yak),
                                         path(mat.yak), path(hap1.p_ctg.gfa),
                                         path(hap2.p_ctg.gfa) ]
        out_assembled_utg  -- channel: [ val(meta), path(p_utg.fa), path(p_utg.fa.fai) ]
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
    container 'becklab/pav:2.4.0'
    cpus   params.n_proc
    memory '128 GB'

    input:
    tuple val(meta), path(out_asm_gfa), path(pat_yak), path(mat_yak), path(out_asm_hap1), path(out_asm_hap2)
    tuple val(meta1), path(out_utg_fa), path(out_utg_fa_fai)
    path  reference_file
    path  reference_index_file
    val   n_proc
    val   rId

    output:
    tuple val(meta), path(pav_vcf_gz_rid), path(pav_vcf_gz_tbi_rid), emit: out_vcf

    script:
    out_asm_prefix   = meta.id + ".asm"
    pav_vcf_gz       = meta.id + ".vcf.gz"
    pav_vcf_gz_tbi   = meta.id + ".vcf.gz.tbi"
    pav_vcf_gz_rid   = meta.id + "_hifiasm_" + rId + "_pav.vcf.gz"
    pav_vcf_gz_tbi_rid = meta.id + "_hifiasm_" + rId + "_pav.vcf.gz.tbi"

    """
    echo '{"reference": "$reference_file"}' > config.json
    echo -e "NAME\\tHAP1\\tHAP2\\n${meta.id}\\t${out_asm_hap1}\\t${out_asm_hap2}" > assemblies.tsv

    /opt/pav/files/docker/run -c $n_proc --nt

    mv $pav_vcf_gz $pav_vcf_gz_rid
    mv $pav_vcf_gz_tbi $pav_vcf_gz_tbi_rid
    """
}