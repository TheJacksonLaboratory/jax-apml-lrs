#!/usr/bin/env nextflow

nextflow.enable.dsl=2
nextflow.preview.output = true

/*
========================================================================================
    LRS_READ WORKFLOW
========================================================================================
    End-to-end long-read sequencing analysis pipeline for PacBio HiFi data.

    Steps:
        00  Input validation and samplesheet parsing
        01  Read alignment with pbmm2
        02  Alignment quality metrics
        03  Small variant calling and phasing (DeepVariant + WhatsHap)
        04  Parallel variant calling:
                - Structural variants: pbsv, Sniffles2, Delly
                - Copy number variants: HiFiCNV
                - Tandem repeats: TRGT
                - Paralogous genes: Paraphase
                - Methylation: pb-CpG-tools + MethBat (experimental)
        05  SV population frequency annotation: SVAFotate
        06  B-allele frequency and CN plots
        07  Variant summary metrics
        08  Phenotype-driven SV prioritization: SvAnna

    Usage:
        nextflow run lrs_read.nf -profile <profile> --csv_path <samplesheet.csv>

    Required parameters:
        --csv_path            path to input samplesheet CSV
        --outputDir           path to output directory
        --reference_file_path path to reference genome FASTA
        --reference_index_file_path path to reference genome FASTA index

    See nextflow.config for full parameter documentation.
========================================================================================
*/

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (params.csv_path) { ch_input = file(params.csv_path) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                          } from './subworkflow/00_INPUT_CHECK.nf'
include { PBMM2_ALIGN; PBMM2_SAMTOOLS         } from './subworkflow/01_ALIGN_READS_PBMM2.nf'
include { QUALITY_METRICS_NANOPLOT;
          QUALITY_METRICS_SAMTOOLS;
          QUALITY_METRICS_GETCOVERAGE          } from './subworkflow/02_QUALITY_METRICS.nf'
include { CALL_SMALL_VARIANTS_DEEPVARIANT1;
          PHASE_SMALL_VARIANTS_WHATSHAP1;
          PHASE_SMALL_VARIANTS_TABIX;
          PHASE_SMALL_VARIANTS_WHATSHAP2;
          CALL_SMALL_VARIANTS_SAMTOOLS;
          CALL_SMALL_VARIANTS_DEEPVARIANT2;
          CALL_SMALL_VARIANTS_BCFTOOLS         } from './subworkflow/03_PHASE_CALL_SMALL_VARIANTS.nf'
include { CALL_SVS_PBSV_DISCOVER;
          CALL_SVS_PBSV_TABIX;
          CALL_SVS_PBSV_CALL;
          CALL_SVS_PBSV_BCFTOOLS               } from './subworkflow/04_CALL_SV_PBSV.nf'
include { CALL_SVS_SNIFFLES2_CALL;
          CALL_SVS_SNIFFLES2_FILTERING         } from './subworkflow/04_CALL_SV_SNIFFLES2.nf'
include { CALL_SVS_DELLY_LR;
          CALL_SVS_DELLY_OUT_AND_FILT          } from './subworkflow/04_CALL_SV_DELLY.nf'
include { CALL_CNVS_HIFICNV                   } from './subworkflow/04_CALL_CNV_HIFICNV.nf'
include { CALL_TRGT_GENOTYPE;
          CALL_TRGT_BCFTOOLS;
          CALL_TRGT_TABIX;
          CALL_TRGT_SAMTOOLS;
          CALL_TRGT_GATK;
          CALL_TRGT_BCFTOOLS_FINAL             } from './subworkflow/04_CALL_REPEATS_TRGT.nf'
include { CALL_PARAPHASE_CALL;
          CALL_PARAPHASE_EXTRACT               } from './subworkflow/04_CALL_PARALOGS_PARAPHASE.nf'
include { METHBAT_ALIGNED_BAM;
          METHBAT_PROFILE                      } from './subworkflow/04_CALL_METHYLATION.nf'
include { CALL_SVAFOTATE_ANNOTATE;
          CALL_SVAFOTATE_ANNOTATE as CALL_SVAFOTATE_ANNOTATE_1;
          CALL_SVAFOTATE_ANNOTATE as CALL_SVAFOTATE_ANNOTATE_2;
          CALL_SVAFOTATE_FILTER;
          CALL_SVAFOTATE_FILTER as CALL_SVAFOTATE_FILTER_1;
          CALL_SVAFOTATE_FILTER as CALL_SVAFOTATE_FILTER_2   } from './subworkflow/05_ANNOTATE_SV_SVAFOTATE.nf'
include { GENERATE_B_ALLELE                   } from './subworkflow/06_GENERATE_BAF_PLOT.nf'
include { VARIANT_SUMMARY_METRICS_BCFTOOLS;
          VARIANT_SUMMARY_METRICS_SURVIVOR     } from './subworkflow/07_VARIANT_SUMMARY_METRICS.nf'
include { PRIORITIZE_SVS_BY_SVANNA;
          COMBINE_SVS_BY_SVANNA               } from './subworkflow/08_PRIORITIZE_SV_SVANNA.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    main:

    // Reference and resource parameters
    path_to_reference_file       = params.reference_file_path
    path_to_reference_index_file = params.reference_index_file_path
    refID                        = params.referenceID
    nThread                      = params.NTHREAD
    moreThread                   = params.MORECPU

    // Quality metrics
    getCoverageScriptPath        = params.quality_metrics.getcoverage_script_path

    // PBSV
    path_to_ref_trf              = params.ref_trf_path

    // B-allele plot
    b_allele_plot_r              = params.b_allele.b_allele_plot_r

    // Delly
    path_to_hg38_excl            = params.delly.hg38_excl_path
    min_mapping_quality          = params.delly.min_mapping_quality
    alignment_scoring            = params.delly.alignment_scoring
    min_clique_size              = params.delly.min_clique_size

    // HiFiCNV
    excl_region_path             = params.CNVs_HiFiCNV.excl_region_path

    // TRGT
    repeats_path                 = params.repeats_path

    // Paraphase
    paraphase_threads                  = params.Paraphase.paraphase_threads
    extract_paraphase_json_script_path = params.Paraphase.extract_paraphase_json_r_path
    targetGenesTxt                     = params.Paraphase.target_genes_txt_path

    // MethBat
    path_to_methbat_input_regions = params.methbat.meth_input_regions

    // SvAnna
    pbsvSniffles2SvannaRScriptPath = params.svAnna.pbsv_sniffles_2_svanna_r_script_path

    // SVAFotate
    svafotate_overlap            = params.svafotate.overlap
    svafotate_bed                = params.svafotate.svafotate_bed

    // ---------------------------------------------------------------------------
    // 00 Input check
    // ---------------------------------------------------------------------------
    input_check_out = INPUT_CHECK(
        ch_input
    )

    // ---------------------------------------------------------------------------
    // 01 Align reads
    // ---------------------------------------------------------------------------
    pbmm2_align_out = PBMM2_ALIGN(
        input_check_out.bams,
        file(path_to_reference_file),
        refID,
        moreThread
    )

    pbmm2_samtools_out = PBMM2_SAMTOOLS(
        pbmm2_align_out.aligned_bam,
        refID
    )

    // ---------------------------------------------------------------------------
    // 02 Quality metrics
    // ---------------------------------------------------------------------------
    quality_metrics_nanoplot_out = QUALITY_METRICS_NANOPLOT(
        pbmm2_samtools_out.aligned_bam_index,
        nThread,
        refID
    )

    quality_metrics_samtools_out = QUALITY_METRICS_SAMTOOLS(
        pbmm2_samtools_out.aligned_bam_index,
        refID
    )

    quality_metrics_getcoverage_out = QUALITY_METRICS_GETCOVERAGE(
        quality_metrics_samtools_out.samtools_coverage,
        file(getCoverageScriptPath),
        refID
    )

    // ---------------------------------------------------------------------------
    // 03 Phase and call small variants
    // ---------------------------------------------------------------------------
    call_small_variants_deepvariant1_out = CALL_SMALL_VARIANTS_DEEPVARIANT1(
        file(path_to_reference_file),
        file(path_to_reference_index_file),
        refID,
        pbmm2_samtools_out.aligned_bam_index
    )

    ch_phase_small_variants_whatshap1_in = pbmm2_samtools_out.aligned_bam_index
        .join(call_small_variants_deepvariant1_out.deepvariant1_vcf)
    phase_small_variants_whatshap1_out = PHASE_SMALL_VARIANTS_WHATSHAP1(
        file(path_to_reference_file),
        file(path_to_reference_index_file),
        refID,
        ch_phase_small_variants_whatshap1_in
    )

    ch_phase_small_variants_tabix_in = ch_phase_small_variants_whatshap1_in
        .join(phase_small_variants_whatshap1_out.whatshap_vcf)
    phase_small_variants_tabix_out = PHASE_SMALL_VARIANTS_TABIX(
        refID,
        ch_phase_small_variants_tabix_in
    )

    ch_phase_small_variants_whatshap2_in = ch_phase_small_variants_whatshap1_in
        .join(phase_small_variants_tabix_out.whatshap_vcf_indexed)
    phase_small_variants_whatshap2_out = PHASE_SMALL_VARIANTS_WHATSHAP2(
        file(path_to_reference_file),
        file(path_to_reference_index_file),
        refID,
        ch_phase_small_variants_whatshap2_in
    )

    ch_call_small_variants_samtools_in = ch_phase_small_variants_whatshap1_in
        .join(phase_small_variants_tabix_out.whatshap_vcf_indexed)
        .join(phase_small_variants_whatshap2_out.whatshap2_bam)
    call_small_variants_samtools_out = CALL_SMALL_VARIANTS_SAMTOOLS(
        nThread,
        refID,
        ch_call_small_variants_samtools_in
    )

    ch_call_small_variants_deepvariant2_in = ch_phase_small_variants_whatshap1_in
        .join(phase_small_variants_tabix_out.whatshap_vcf_indexed)
        .join(call_small_variants_samtools_out.whatshap2_bam_bai)
    call_small_variants_deepvariant2_out = CALL_SMALL_VARIANTS_DEEPVARIANT2(
        file(path_to_reference_file),
        file(path_to_reference_index_file),
        refID,
        ch_call_small_variants_deepvariant2_in
    )

    call_small_variants_bcftools_out = CALL_SMALL_VARIANTS_BCFTOOLS(
        refID,
        call_small_variants_deepvariant2_out.deepvariant2_vcfs
    )

    // ---------------------------------------------------------------------------
    // 04 Parallel variant calling
    // ---------------------------------------------------------------------------

    // PBSV
    call_svs_pbsv_discover_out = CALL_SVS_PBSV_DISCOVER(
        file(path_to_ref_trf),
        refID,
        pbmm2_samtools_out.aligned_bam_index
    )

    call_svs_pbsv_tabix_out = CALL_SVS_PBSV_TABIX(
        refID,
        call_svs_pbsv_discover_out.svsig
    )

    call_svs_pbsv_call_out = CALL_SVS_PBSV_CALL(
        file(path_to_reference_file),
        refID,
        nThread,
        call_svs_pbsv_tabix_out.svsig_indexed
    )

    call_svs_pbsv_bcftools_out = CALL_SVS_PBSV_BCFTOOLS(
        refID,
        call_svs_pbsv_call_out.svs_vcf
    )

    // Sniffles2
    call_svs_sniffles2_call_out = CALL_SVS_SNIFFLES2_CALL(
        refID,
        pbmm2_samtools_out.aligned_bam_index
    )

    call_svs_sniffles2_filtering_out = CALL_SVS_SNIFFLES2_FILTERING(
        refID,
        call_svs_sniffles2_call_out.sniffles2_vcf
    )

    // Delly
    call_svs_delly_lr_out = CALL_SVS_DELLY_LR(
        file(path_to_hg38_excl),
        file(path_to_reference_file),
        refID,
        pbmm2_samtools_out.aligned_bam_index,
        min_mapping_quality,
        alignment_scoring,
        min_clique_size,
        nThread
    )

    call_svs_delly_out_and_filt_out = CALL_SVS_DELLY_OUT_AND_FILT(
        refID,
        call_svs_delly_lr_out.delly_out_bcf,
        nThread
    )

    // HiFiCNV
    call_cnvs_hificnv_out = CALL_CNVS_HIFICNV(
        file(excl_region_path),
        file(path_to_reference_file),
        refID,
        pbmm2_samtools_out.aligned_bam_index,
        nThread
    )

    // TRGT
    call_trgt_genotype_out = CALL_TRGT_GENOTYPE(
        file(path_to_reference_file),
        file(path_to_reference_index_file),
        repeats_path,
        nThread,
        refID,
        pbmm2_samtools_out.aligned_bam_index
    )

    call_trgt_bcftools_out = CALL_TRGT_BCFTOOLS(
        refID,
        call_trgt_genotype_out.trgt_bam_vcf
    )

    call_trgt_tabix_out = CALL_TRGT_TABIX(
        refID,
        call_trgt_bcftools_out.trgt_sorted_vcf
    )

    call_trgt_samtools_out = CALL_TRGT_SAMTOOLS(
        refID,
        call_trgt_genotype_out.trgt_bam_vcf
    )

    call_trgt_gatk_out = CALL_TRGT_GATK(
        refID,
        call_trgt_tabix_out.trgt_vcf_index
    )

    call_trgt_bcftools_final_out = CALL_TRGT_BCFTOOLS_FINAL(
        refID,
        call_trgt_gatk_out.trgt_sorted_tsv
    )

    // Paraphase
    call_paraphase_call_out = CALL_PARAPHASE_CALL(
        file(path_to_reference_file),
        pbmm2_samtools_out.aligned_bam_index,
        paraphase_threads
    )

    call_paraphase_extract_out = CALL_PARAPHASE_EXTRACT(
        call_paraphase_call_out.paraphase_out,
        file(extract_paraphase_json_script_path),
        file(targetGenesTxt)
    )

    // Methylation (experimental)
    methbat_aligned_bam_out = METHBAT_ALIGNED_BAM(
        nThread,
        refID,
        call_small_variants_samtools_out.whatshap2_bam_bai
    )

    methbat_profile_out = METHBAT_PROFILE(
        file(path_to_methbat_input_regions),
        refID,
        methbat_aligned_bam_out.methbat_bams
    )

    // ---------------------------------------------------------------------------
    // 05 SV annotation with SVAFotate
    // ---------------------------------------------------------------------------
    call_svs_pbsv_bcftools_out.svs_filtered_vcfs
        .map { meta, vcf1, vcf2 -> tuple(meta, vcf1) }
        .set { svafotate_pbsv_input }

    call_svafotate_annotate_Emedgene_vcfs_out = CALL_SVAFOTATE_ANNOTATE(
        refID,
        svafotate_pbsv_input,
        "Emedgene",
        svafotate_overlap,
        svafotate_bed,
        nThread
    )

    call_svafotate_annotate_Sniffles2_vcfs_out = CALL_SVAFOTATE_ANNOTATE_1(
        refID,
        call_svs_sniffles2_filtering_out.sniffles2_FILT_vcf,
        "sniffles2_FILT",
        svafotate_overlap,
        svafotate_bed,
        nThread
    )

    call_svs_delly_out_and_filt_out.delly_out_filt_vcf
        .map { meta, vcf1, vcf2 -> tuple(meta, vcf2) }
        .set { svafotate_delly_input }

    call_svafotate_annotate_delly_vcfs_out = CALL_SVAFOTATE_ANNOTATE_2(
        refID,
        svafotate_delly_input,
        "delly",
        svafotate_overlap,
        svafotate_bed,
        nThread
    )

    call_svafotate_filter_Emedgene_vcf_out = CALL_SVAFOTATE_FILTER(
        refID,
        call_svafotate_annotate_Emedgene_vcfs_out.svafotate_annotate_out,
        "Emedgene"
    )

    call_svafotate_filter_Sniffles2_vcf_out = CALL_SVAFOTATE_FILTER_1(
        refID,
        call_svafotate_annotate_Sniffles2_vcfs_out.svafotate_annotate_out,
        "sniffles2_FILT"
    )

    call_svafotate_filter_delly_vcf_out = CALL_SVAFOTATE_FILTER_2(
        refID,
        call_svafotate_annotate_delly_vcfs_out.svafotate_annotate_out,
        "delly"
    )

    // ---------------------------------------------------------------------------
    // 06 B-allele frequency plot
    // ---------------------------------------------------------------------------
    ch_generate_b_allele_in = call_cnvs_hificnv_out.hificnv_vcf
        .join(call_small_variants_bcftools_out.deepvariant2_filt_vcfs)
    generate_b_allele_out = GENERATE_B_ALLELE(
        file(b_allele_plot_r),
        ch_generate_b_allele_in,
        refID
    )

    // ---------------------------------------------------------------------------
    // 07 Variant summary metrics
    // ---------------------------------------------------------------------------
    variant_summary_metrics_bcftools_out = VARIANT_SUMMARY_METRICS_BCFTOOLS(
        refID,
        call_small_variants_bcftools_out.deepvariant2_filt_vcfs
    )

    ch_variant_summary_metrics_survivor_in = call_svs_pbsv_bcftools_out.svs_filtered_vcfs
        .join(call_svs_sniffles2_filtering_out.sniffles2_FILT_vcf)
        .join(call_svs_delly_out_and_filt_out.delly_out_filt_vcf)
    variant_summary_metrics_survivor_out = VARIANT_SUMMARY_METRICS_SURVIVOR(
        ch_variant_summary_metrics_survivor_in,
        refID
    )

    // ---------------------------------------------------------------------------
    // 08 Phenotype-driven SV prioritization
    // ---------------------------------------------------------------------------
    ch_prioritize_svs_by_svanna_in = call_svs_pbsv_bcftools_out.svs_filtered_vcfs
        .join(call_svs_sniffles2_filtering_out.sniffles2_FILT_vcf)
        .join(input_check_out.bams)
    prioritize_svs_by_svanna_out = PRIORITIZE_SVS_BY_SVANNA(
        ch_prioritize_svs_by_svanna_in,
        refID
    )

    combine_svs_by_svanna_out = COMBINE_SVS_BY_SVANNA(
        file(pbsvSniffles2SvannaRScriptPath),
        refID,
        prioritize_svs_by_svanna_out.csv_files
    )

    /*
    ========================================================================================
        PUBLISH
    ========================================================================================
    */

    publish:
    pbmm2_align_out.aligned_bam                                    >> 'pbmm2_align'
    pbmm2_samtools_out.aligned_bam_index                           >> 'pbmm2_samtools'

    quality_metrics_nanoplot_out.nanoplot_results                  >> 'quality_metrics_nanoplot'
    quality_metrics_samtools_out.samtools_coverage                 >> 'quality_metrics_samtools'
    quality_metrics_getcoverage_out.meandepth                      >> 'quality_metrics_getcoverage'

    call_small_variants_deepvariant1_out.deepvariant1_vcf          >> 'call_small_variants_deepvariant1'
    phase_small_variants_whatshap1_out.whatshap_vcf                >> 'phase_small_variants_whatshap1'
    phase_small_variants_tabix_out.whatshap_vcf_indexed            >> 'phase_small_variants_tabix'
    phase_small_variants_whatshap2_out.whatshap2_bam               >> 'phase_small_variants_whatshap2'
    call_small_variants_samtools_out.whatshap2_bam_bai             >> 'call_small_variants_samtools'
    call_small_variants_deepvariant2_out.deepvariant2_vcfs         >> 'call_small_variants_deepvariant2'
    call_small_variants_bcftools_out.deepvariant2_filt_vcfs        >> 'call_small_variants_bcftools'

    call_svs_pbsv_discover_out.svsig                               >> 'call_svs_pbsv_discover'
    call_svs_pbsv_tabix_out.svsig_indexed                          >> 'call_svs_pbsv_tabix'
    call_svs_pbsv_call_out.svs_vcf                                 >> 'call_svs_pbsv_call'
    call_svs_pbsv_bcftools_out.svs_filtered_vcfs                   >> 'call_svs_pbsv_bcftools'

    call_svs_sniffles2_call_out.sniffles2_vcf                      >> 'call_svs_sniffles2_call'
    call_svs_sniffles2_filtering_out.sniffles2_FILT_vcf            >> 'call_svs_sniffles2_filtering'

    call_svs_delly_lr_out.delly_out_bcf                            >> 'call_svs_delly_lr'
    call_svs_delly_out_and_filt_out.delly_out_filt_vcf             >> 'call_svs_delly_out_and_filt'

    call_cnvs_hificnv_out.hificnv_vcf                             >> 'call_cnvs_hificnv'

    call_trgt_genotype_out.trgt_bam_vcf                            >> 'call_trgt_genotype'
    call_trgt_bcftools_out.trgt_sorted_vcf                         >> 'call_trgt_bcftools'
    call_trgt_tabix_out.trgt_vcf_index                             >> 'call_trgt_tabix'
    call_trgt_samtools_out.trgt_bam_sorted                         >> 'call_trgt_samtools'
    call_trgt_gatk_out.trgt_sorted_tsv                             >> 'call_trgt_gatk'
    call_trgt_bcftools_final_out.new_sorted_tsv_trgt_file          >> 'call_trgt_bcftools_final'

    call_paraphase_call_out.paraphase_out                          >> 'call_paraphase_call'
    call_paraphase_call_out.paraphase_vcfs_dir                     >> 'call_paraphase_vcfs_dir'

    generate_b_allele_out.report_pdf                               >> 'generate_b_allele'

    variant_summary_metrics_bcftools_out.deepvariant_varmetrics    >> 'variant_summary_metrics_bcftools'
    variant_summary_metrics_survivor_out.pbsv_varmetrics           >> 'variant_summary_metrics_survivor_pbsv'
    variant_summary_metrics_survivor_out.sniffles2_varmetrics      >> 'variant_summary_metrics_survivor_sniffles2'
    variant_summary_metrics_survivor_out.delly_varmetrics          >> 'variant_summary_metrics_survivor_delly'

    prioritize_svs_by_svanna_out.csv_files                         >> 'prioritize_svs_by_svanna'
    combine_svs_by_svanna_out.sv_candidates_file                   >> 'combine_svs_by_svanna'

    methbat_aligned_bam_out.methbat_bams                           >> 'methbat_aligned_bam'
    methbat_profile_out.methbat_region_profile                     >> 'methbat_profile'

    call_svafotate_annotate_Emedgene_vcfs_out.svafotate_annotate_out  >> 'svafotate_Emedgene_vcfs_out'
    call_svafotate_annotate_Sniffles2_vcfs_out.svafotate_annotate_out >> 'svafotate_Sniffles2_vcfs_out'
    call_svafotate_annotate_delly_vcfs_out.svafotate_annotate_out     >> 'svafotate_delly_vcfs_out'

    call_svafotate_filter_Emedgene_vcf_out.svafotate_annotate_out  >> 'svafotate_filter_Emedgene_vcfs_out'
    call_svafotate_filter_Sniffles2_vcf_out.svafotate_annotate_out >> 'svafotate_filter_Sniffles2_vcfs_out'
    call_svafotate_filter_delly_vcf_out.svafotate_annotate_out     >> 'svafotate_filter_delly_vcfs_out'
}

/*
========================================================================================
    OUTPUT
========================================================================================
*/

output {
    'pbmm2_align' {
        path { meta, out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/pbmm2" }
        mode 'copy'
        overwrite false
    }
    'pbmm2_samtools' {
        path { meta, bam_file_name, bam_index_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/pbmm2" }
        mode 'copy'
        overwrite false
    }
    'quality_metrics_nanoplot' {
        path { meta, nanoplot_files -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/bam" }
        mode 'copy'
        overwrite false
    }
    'quality_metrics_samtools' {
        path { meta, samtools_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/bam" }
        mode 'copy'
        overwrite false
    }
    'quality_metrics_getcoverage' {
        path { meta, meandepth_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/bam" }
        mode 'copy'
        overwrite false
    }
    'call_small_variants_deepvariant1' {
        path { meta, deepvariant1_vcf_out_file_name, deepvariant1_vcf_tbi_out_file_name, deepvariant1_visual_report_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/deepvariant1" }
        mode 'copy'
        overwrite false
    }
    'phase_small_variants_whatshap1' {
        path { meta, whatshap_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/whatshap" }
        mode 'copy'
        overwrite false
    }
    'phase_small_variants_tabix' {
        path { meta, whatshap_vcf_out_file_name, whatshap_vcf_index_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/whatshap" }
        mode 'copy'
        overwrite false
    }
    'phase_small_variants_whatshap2' {
        path { meta, whatshap2_bam_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/whatshap" }
        mode 'copy'
        overwrite false
    }
    'call_small_variants_samtools' {
        path { meta, whatshap2_bam_out_file_name, whatshap2_bam_bai_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/whatshap" }
        mode 'copy'
        overwrite false
    }
    'call_small_variants_deepvariant2' {
        path { meta, deepvariant2_vcf_out_file_name, deepvariant2_vcf_tbi_out_file_name, deepvariant2_g_vcf_out_file_name, deepvariant2_g_vcf_tbi_out_file_name, deepvariant2_visual_report_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/deepvariant2" }
        mode 'copy'
        overwrite false
    }
    'call_small_variants_bcftools' {
        path { meta, deepvariant2_FILT_vcf_out_file_name, deepvariant2_Emedgene_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/deepvariant2" }
        mode 'copy'
        overwrite false
    }
    'call_svs_pbsv_discover' {
        path { meta, svsig_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/pbsv" }
        mode 'copy'
        overwrite false
    }
    'call_svs_pbsv_tabix' {
        path { meta, svsig_out_file_name, svsig_out_file_name_index -> "${params.outputDir}${meta.id}/${params.referenceID}/pbsv" }
        mode 'copy'
        overwrite false
    }
    'call_svs_pbsv_call' {
        path { meta, pbsv_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/pbsv" }
        mode 'copy'
        overwrite false
    }
    'call_svs_pbsv_bcftools' {
        path { meta, pbsv_FILT_vcf_out_file_name, pbsv_Emedgene_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/pbsv" }
        mode 'copy'
        overwrite false
    }
    'call_svs_sniffles2_call' {
        path { meta, sniffles2_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/sniffles2" }
        mode 'copy'
        overwrite false
    }
    'call_svs_sniffles2_filtering' {
        path { meta, sniffles2_FILT_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/sniffles2" }
        mode 'copy'
        overwrite false
    }
    'call_svs_delly_lr' {
        path { meta, delly_bcf_out_file_name, delly_bcf_csi_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/delly" }
        mode 'copy'
        overwrite false
    }
    'call_svs_delly_out_and_filt' {
        path { meta, delly_vcf_out_file_name, delly_filt_vcf_out_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/delly" }
        mode 'copy'
        overwrite false
    }
    'call_cnvs_hificnv' {
        path { meta, hificnv_vcf_out_file_name, hificnv_depth_bw_out_file_name, hificnv_depth_copynum_bedgraph_out_file_name, hificnv_log -> "${params.outputDir}${meta.id}/${params.referenceID}/hificnv" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_genotype' {
        path { meta, out_trgt_vcf, out_trgt_bam -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_bcftools' {
        path { meta, out_trgt_vcf_sorted_fileName -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_tabix' {
        path { meta, out_trgt_vcf_sorted_fileName, out_trgt_vcf_index_sorted_fileName -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_samtools' {
        path { meta, out_trgt_bam_sorted_fileName, out_trgt_bam_bai_sorted_fileName -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_gatk' {
        path { meta, out_trgt_vcf_sorted_fileName, out_trgt_vcf_index_sorted_fileName, out_trgt_tsv_sorted_fileName -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_trgt_bcftools_final' {
        path { meta, sorted_tsv_trgt_file_name -> "${params.outputDir}${meta.id}/${params.referenceID}/trgt" }
        mode 'copy'
        overwrite false
    }
    'call_paraphase_call' {
        path { meta, out_paraphase_json, out_paraphase_bam, out_paraphase_bam_bai -> "${params.outputDir}${meta.id}/${params.referenceID}/paraphase" }
        mode 'copy'
        overwrite false
    }
    'call_paraphase_vcfs_dir' {
        path { meta, files -> "${params.outputDir}${meta.id}/${params.referenceID}/paraphase/vcfs" }
        mode 'copy'
        overwrite false
    }
    'generate_b_allele' {
        path { meta, pdfs -> "${params.outputDir}${meta.id}/${params.referenceID}/b_allele" }
        mode 'copy'
        overwrite false
    }
    'variant_summary_metrics_bcftools' {
        path { meta, out_deepvariant_varmetrics_file -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/vcf" }
        mode 'copy'
        overwrite false
    }
    'variant_summary_metrics_survivor_pbsv' {
        path { meta, out_pbsv_varmetrics_file, out_pbsv_varmetrics_CHR_file, out_pbsv_varmetrics_support_file -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/vcf" }
        mode 'copy'
        overwrite false
    }
    'variant_summary_metrics_survivor_sniffles2' {
        path { meta, out_sniffles2_varmetrics_file, out_sniffles2_varmetrics_CHR_file, out_sniffles2_varmetrics_support_file -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/vcf" }
        mode 'copy'
        overwrite false
    }
    'variant_summary_metrics_survivor_delly' {
        path { meta, out_delly_varmetrics_file, out_delly_varmetrics_CHR_file, out_delly_varmetrics_support_file -> "${params.outputDir}${meta.id}/${params.referenceID}/metrics/vcf" }
        mode 'copy'
        overwrite false
    }
    'prioritize_svs_by_svanna' {
        path { meta, out_pbsv_svanna_file, out_pbsv_svanna_file_html, out_pbsv_svanna_file_vcf_gz, meta1, out_sniffles2_svanna_file, out_sniffles2_svanna_file_html, out_sniffles2_svanna_file_vcf_gz -> "${params.outputDir}${meta.id}/${params.referenceID}/svanna" }
        mode 'copy'
        overwrite false
    }
    'combine_svs_by_svanna' {
        path { meta, out_sv_candidates_file -> "${params.outputDir}${meta.id}/${params.referenceID}/svanna" }
        mode 'copy'
        overwrite false
    }
    'methbat_aligned_bam' {
        path { meta, out_combined_bed, out_combined_bw, out_hap1_bed, out_hap1_bw, out_hap2_bed, out_hap2_bw, out_log, out_prefix_out -> "${params.outputDir}${meta.id}/${params.referenceID}/methbat" }
        mode 'copy'
        overwrite false
    }
    'methbat_profile' {
        path { meta, out_region_profile -> "${params.outputDir}${meta.id}/${params.referenceID}/methbat" }
        mode 'copy'
        overwrite false
    }
    'svafotate_Emedgene_vcfs_out' {
        path { meta, vcf_out -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
    'svafotate_Sniffles2_vcfs_out' {
        path { meta, vcf_out -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
    'svafotate_delly_vcfs_out' {
        path { meta, vcf_out -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
    'svafotate_filter_Emedgene_vcfs_out' {
        path { meta, vcf_out_rare_unique, vcf_out_low_freq -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
    'svafotate_filter_Sniffles2_vcfs_out' {
        path { meta, vcf_out_rare_unique, vcf_out_low_freq -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
    'svafotate_filter_delly_vcfs_out' {
        path { meta, vcf_out_rare_unique, vcf_out_low_freq -> "${params.outputDir}${meta.id}/${params.referenceID}/svafotate" }
        mode 'copy'
        overwrite false
    }
}