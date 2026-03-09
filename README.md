![JAX-APML-LRS Pipeline](resources/apml-lrs-pipeline.png)

# **JAX-APML-LRS** `jax-internal` is configured for Sumner2

# Overview

**JAX-APML-LRS** was developed for research use by the Advanced Precision Medicine Laboratory (APML) at The Jackson Laboratory (JAX) for Genomic Medicine. 

This is a Nextflow DSL2 pipeline for end-to-end analysis of PacBio HiFi long-read sequencing data against the human [hg38 no-alt](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz) reference genome, supporting:
- HiFi unaligned BAM merging
- *de novo* genome assembly (single-sample or trio-based)
- read-based alignment
- variant calling
- phenotype-driven structural variant prioritization

---

# Requirements

- [Nextflow](https://www.nextflow.io/) version 24.10.5 (bundled bootstrap script included, see [Quick Start](#quick-start))
- [Docker](https://www.docker.com/) or [Apptainer/Singularity](https://apptainer.org/)
- Java version 11 or higher

    *NOTE: All additional tools run inside containers and no local installation is required.*

---

# Workflows

| Workflow | Description |
|---|---|
| `lrs_pbmerge` | Merge unaligned PacBio HiFi BAM files (2 or 3) per sample |
| `lrs_read` | Read alignment, small and structural variant calling, methylation, and phenotype-driven SV prioritization from unaligned PacBio HiFi BAM |
| `lrs_asm_single` | Single-sample *de novo* assembly, variant calling, and phenotype-driven SV prioritization from PacBio HiFi FASTQ |
| `lrs_asm_trio` | Trio-informed *de novo* assembly, variant calling, and phenotype-driven SV prioritization from PacBio HiFi FASTQ |

---

# Execution Profiles

Each workflow's `nextflow.config` defines three profiles.

| Profile | Description | Notes |
|---|---|---|
| `standard` | Local execution | Default |
| `hpc` | SLURM cluster with Singularity | Edit queue & account settings in `nextflow.config`, and uncomment `process.module = 'apptainer'` |
| `gcb` | Google Cloud Batch | Edit project and bucket settings in `nextflow.config` |

<br>

---

<br>

<details>
<summary><b>Repository Structure (Click to expand)</b></summary>

```
jax-apml-lrs/
  nextflow                              Nextflow bootstrap script
  run.sh                                Unified workflow launcher
  download_refs.sh                      Large reference file downloader
  refs/                                 Bundled small reference files
    human.hg38.excl.tsv
    cnv.excluded_regions.common_50.hg38.bed.gz
    medically_relevant_repeats_hg38.bed
    meth_profile_model.tsv
    targetGenes.txt
  scripts/                              Bundled pipeline R and Python scripts
    coverage.py
    baf_plot.R
    convert_paraphase.R
    filter_svanna_pbsv_sniffles2.R
    filter_svanna_pav.R
  workflows/
    lrs_read/
      lrs_read.nf
      nextflow.config
      subworkflow/
        00_INPUT_CHECK.nf
        01_ALIGN_READS_PBMM2.nf
        02_QUALITY_METRICS.nf
        03_PHASE_CALL_SMALL_VARIANTS.nf
        04_CALL_SV_PBSV.nf
        04_CALL_SV_SNIFFLES2.nf
        04_CALL_SV_DELLY.nf
        04_CALL_CNV_HIFICNV.nf
        04_CALL_REPEATS_TRGT.nf
        04_CALL_PARALOGS_PARAPHASE.nf
        04_CALL_METHYLATION.nf
        05_ANNOTATE_SV_SVAFOTATE.nf
        06_VARIANT_SUMMARY_METRICS.nf
        07_PRIORITIZE_SV_SVANNA.nf
      modules/
        samplesheet_check_read.nf
      bin/
        check_samplesheet_read.py
    lrs_asm_single/
      lrs_asm_single.nf
      nextflow.config
      subworkflow/
        00_INPUT_CHECK.nf
        01_ASSEMBLE_HIFIASM.nf
        02_CALL_SV_PAV.nf
        03_SPLIT_PAV_VARIANT_SIZE.nf
	04_ANNOTATE_SV_SVAFOTATE.nf        
	05_PRIORITIZE_SV_SVANNA.nf
      modules/
        samplesheet_check_asm-single.nf
      bin/
        check_samplesheet_asm-single.py
    lrs_asm_trio/
      lrs_asm_trio.nf
      nextflow.config
      subworkflow/
        00_INPUT_CHECK.nf
        01_ASSEMBLE_TRIO_HIFIASM.nf
        02_CALL_SV_PAV.nf
        03_SPLIT_PAV_VARIANT_SIZE.nf
	04_ANNOTATE_SV_SVAFOTATE.nf        
	05_PRIORITIZE_SV_SVANNA.nf
      modules/
        samplesheet_check_asm-trio.nf
      bin/
        check_samplesheet_asm-trio.py
    lrs_pbmerge/
      lrs_pbmerge.nf
      nextflow.config
      subworkflow/
        00_INPUT_CHECK.nf
        01_PBMERGE_BAMS.nf
      modules/
        samplesheet_check_pbmerge.nf
      bin/
        check_samplesheet_pbmerge.py
```

</details>

<br>

---

# Quick Start

## 1. Clone the repository

    git clone https://github.com/TheJacksonLaboratory/jax-apml-lrs.git
    cd jax-apml-lrs

## 2. Install Nextflow

>The repository includes the Nextflow bootstrap script. Install it in the repo root: 

```bash    
curl -s https://get.nextflow.io | bash
```

## 3. Container images

>All containers are pulled automatically from public registries at runtime. No local image builds are required. On first run, images are downloaded and cached locally; subsequent runs will use cached images and skip the download step. On Sumner2, the cache is set to `flashscratch/${USER}/containers/` by default in the `hpc` profile.

## 4. Download large reference files

>The following files are too large to bundle and need to be downloaded. 

| File | Size | Source | Used by | Notes
|---|---|---|---|---|
| `hg38.no_alt.fa` | ~3GB | [1000 Genomes EBI FTP](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/) | pbmm2 and all variant callers | HGSVC reference genome for long-reads |
| `hg38.no_alt.fa.fai` | ~8KB | generated via `samtools faidx` | pbmm2 and all variant callers | Index for `hg38.no_alt.fa` |
| `human_GRCh38_no_alt_analysis_set.trf.bed` | ~7MB | [PacBio pbsv](https://github.com/PacificBiosciences/pbsv/tree/master/annotations) | pbsv | Used to improve SV detection within repeat regions |
| `SVAFotate_core_SV_popAFs.GRCh38.v4.1_LR_ v2.0_INS1bp.bed.gz` | ~451MB | [GitHub Release asset](https://github.com/TheJacksonLaboratory/jax-apml-lrs/releases) | SVAFotate | Compiled SV AFs for filtering, modified from the [SVAFotate core BED](https://github.com/fakedrtom/SVAFotate) to include additional LR population databases (HGSVC2/3, GA4K, CoLoRS); hosted as a release asset on this repository |
| `svanna_v1.0.4.sif` | ~900MB | [GitHub Release asset](https://github.com/TheJacksonLaboratory/jax-apml-lrs/releases) | SvAnna | Pre-built Apptainer/Singularity image of SvAnna v1.0.4 with hg38 database bundled; hosted as a release asset on this repository |

>Use `download_refs.sh` to download all large reference files and copy the bundled small reference files into a single directory. Pass the output directory to `run.sh` with `--refs_path`.

```bash
# download to a specified directory (recommended)
./download_refs.sh --outdir /path/to/refs/

# download to refs/ in the repo root (default)
./download_refs.sh

# to use a different container runtime
./download_refs.sh --outdir /path/to/refs/ --container-runtime singularity
```

>The following small reference files are bundled in the repository under `refs/` and are automatically copied by `download_refs.sh`:

| Bundled File | Used by |
|---|---|
| `human.hg38.excl.tsv` | Delly |
| `cnv.excluded_regions.common_50.hg38.bed.gz` | HiFiCNV |
| `medically_relevant_repeats_hg38.bed` | TRGT |
| `meth_profile_model.tsv` | MethBat |
| `targetGenes.txt` | Paraphase |



## 5. Run a workflow

>Submit the pipeline via sbatch as `/flashscratch` is only mounted on compute nodes. Output will be written to `/flashscratch/${USER}/jax-apml-lrs/` by default

```bash
# Kick off specified workflow
sbatch --job-name=nf-<workflow> \
       --time=48:00:00 \
       --mem=16G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w <workflow> -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/ \
                   --refs_path /path/to/refs" # not required for lrs_pbmerge

# Resume a previous run
sbatch --job-name=nf-<workflow> \
       --time=48:00:00 \
       --mem=16G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w <workflow> -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/ \
                   --refs_path /path/to/refs \ # not required for lrs_pbmerge
                   --resume"

# For full usage info, run
./run.sh --help
```

<br>

---

# Workflow details 

 ## (I) Workflow: `lrs_pbmerge`

Merges two (required) or three (optional) unaligned PacBio HiFi BAM files per sample into a single merged BAM using pbmerge. Intended as a preprocessing step before running `lrs_read`.

#### Sample sheet .csv set-up (requires headers shown below)

```csv 
sID,bam1,bam2,bam3
SampleID,/path/to/unaligned_bam1,/path/to/unaligned_bam2,/path/to/unaligned_bam3[optional]
```

#### Run `lrs_pbmerge`

```bash
sbatch --job-name=nf-lrs_pbmerge \
       --time=48:00:00 \
       --mem=16G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w lrs_pbmerge -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/"
```

#### Output Structure

```
<outputDir>/
└── <sampleID>/
    └── pbmerge/
        └── <YYYYMMDD>/
            └── <sampleID>_pbmerge_merged.5mC.bam
```

#### Tool Versions

| Tool | Version | Container | Reference |
|---|---|---|---|
| pbtk (pbmerge) | 3.1.1 | `quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0` | [github.com/PacificBiosciences/pbtk](https://github.com/PacificBiosciences/pbtk) |

<br>

## (II) Workflow: `lrs_read`

Full variant calling pipeline from unaligned PacBio HiFi BAM files. Produces phased BAMs, small variant calls, structural variants calls, CNV calls, tandem repeat genotypes, paralog-resolved calls, methylation profiles, and phenotype-driven SV prioritization.

#### Sample sheet .csv set-up (requires headers shown below)

```csv
sID,reads_unaligned_bam,HPO
SampleID,/path/to/unaligned_bam,/path/to/HPO.txt
```
>NOTE: HPO.txt is a line-separated list of HPO phenotype terms (e.g. `HP:0001250`)

#### Run `lrs_read`
```bash
sbatch --job-name=nf-<workflow> \
       --time=48:00:00 \
       --mem=16G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w lrs_read -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/ \
		   --refs_path refs/"
```

| Step | Subworkflow | Tools | Notes |
|---|---|---|---|
| 00 | Input validation | — | |
| 01 | Read alignment | pbmm2, samtools | |
| 02 | Alignment QC | NanoPlot, samtools | |
| 03 | Small variant calling and phasing | DeepVariant, WhatsHap, bcftools | DeepVariant processes set `LD_LIBRARY_PATH` for GPU support. Without a GPU, falls back to CPU execution. |
| 04 | Structural variant calling | pbsv, Sniffles2, Delly | |
| 04 | Copy number variant calling | HiFiCNV | |
| 04 | Tandem repeat genotyping | TRGT, GATK4 | |
| 04 | Paralog resolution | Paraphase | |
| 04 | Methylation profiling | pb-CpG-tools, MethBat | *in development* |
| 05 | SV population frequency annotation | SVAFotate | |
| 06 | Variant summary metrics | bcftools stats, SURVIVOR | |
| 07 | Phenotype-driven SV prioritization | SvAnna, R | |

#### Output Structure
```
<outputDir>/
└── <sampleID>/
    └── <read>/
    	└── <YYYYMMDD>/
            ├── pbmm2/              aligned BAM and index
            ├── metrics/
            │   ├── bam/            NanoPlot results, samtools coverage, mean depth TSV
            │   └── vcf/            bcftools stats, SURVIVOR stats
            ├── deepvariant1/       first-pass DeepVariant VCF
            ├── deepvariant2/       second-pass DeepVariant VCF (on haplotagged BAM)
            ├── whatshap/           phased VCF, haplotagged BAM
            ├── pbsv/               pbsv SV calls
            ├── sniffles2/          Sniffles2 SV calls
            ├── delly/              Delly SV calls
            ├── hificnv/            HiFiCNV CNV calls
            ├── trgt/               TRGT tandem repeat genotypes
            ├── paraphase/          Paraphase paralog calls and per-gene VCFs
            ├── svafotate/          SVAFotate-annotated SV VCFs
            ├── svanna/             SvAnna prioritized SVs
            └── methbat/            methylation profiles (in development)
```
>NOTE: **Emedgene outputs:** Several processes produce `*_Emedgene*.vcf` files (DeepVariant, pbsv, SVAFotate). These are reformatted copies for ingestion into the [Emedgene](https://www.emedgene.com/) clinical interpretation platform. If not using Emedgene, use primary VCFs.

>NOTE: **SvAnna outputs:** Full SvAnna results (csv/html/vcf) are available per SV caller (pbsv, Sniffles2) from both unfiltered and filtered (AF < 0.01, `*SVAFotate-RARE-UNIQUE.vcf`) VCFs. The `SvAnna_SV_Candidates.csv` file, however, is specifically generated using `*SVAFotate-RARE-UNIQUE.vcf` with a SvAnna psv score >= 2.0. To adjust the psv threshold or change which SvAnna results are used for the candidates file, modify `resources/scripts/filter_svanna_pbsv_sniffles2.R`. Note that the unfiltered SvAnna run requires more memory than the filtered run; if it fails with OOM, see TROUBLESHOOTING.md section 8.

#### Tool versions

| Tool | Version | Container | Reference |
|---|---|---|---|
| pbmm2 | 1.13.1 | `quay.io/pacbio/pbmm2:1.13.1_build3` | [github.com/PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2) |
| samtools | 1.21 | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` | [github.com/samtools](https://github.com/samtools) |
| NanoPlot | 1.42.0 | `staphb/nanoplot:1.42.0` | [github.com/wdecoster/NanoPlot](https://github.com/wdecoster/NanoPlot) |
| DeepVariant | 1.6.1 | `google/deepvariant:1.6.1` | [github.com/google/deepvariant](https://github.com/google/deepvariant) |
| WhatsHap | 2.3 | `quay.io/biocontainers/whatshap:2.3--py38h2494328_0` | [whatshap.readthedocs.io](https://whatshap.readthedocs.io) |
| bcftools | 1.19 | `quay.io/biocontainers/bcftools:1.19--h8b25389_1` | [github.com/samtools/bcftools](https://github.com/samtools/bcftools) |
| tabix | 1.11 | `quay.io/biocontainers/tabix:1.11--hdfd78af_0` | [github.com/tabixio/tabix](https://github.com/tabixio/tabix)|
| pbsv | 2.9.0 | `quay.io/pacbio/pbsv:2.9.0_1.14_build1` | [github.com/PacificBiosciences/pbsv](https://github.com/PacificBiosciences/pbsv)
| Sniffles2 | 2.3.3 | `quay.io/biocontainers/sniffles:2.3.3--pyhdfd78af_0` | [github.com/fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles) |
| Delly | 1.2.9 | `quay.io/biocontainers/delly:1.2.9--hf9970c3_0` | [github.com/dellytools/delly](https://github.com/dellytools/delly) |
| HiFiCNV | 1.0.1 | `quay.io/pacbio/hificnv:1.0.1_build1` | [github.com/PacificBiosciences/HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) |
| TRGT | 1.0.0 | `quay.io/pacbio/trgt:1.0.0_build1` | [github.com/PacificBiosciences/trgt](https://github.com/PacificBiosciences/trgt) |
| GATK4 | 4.6.2 | `quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1` | [github.com/broadinstitute/gatk](https://github.com/broadinstitute/gatk) |
| Paraphase | 3.1.1 | `quay.io/pacbio/paraphase:3.1.1_build1` | [github.com/PacificBiosciences/paraphase](https://github.com/PacificBiosciences/paraphase) |
| pb-CpG-tools | 2.3.2 | `quay.io/pacbio/pb-cpg-tools:v2.3.2` | [github.com/PacificBiosciences/pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) |
| MethBat | 0.13.2 | `quay.io/biocontainers/methbat:0.13.2--h9ee0642_0` | [github.com/PacificBiosciences/MethBat](https://github.com/PacificBiosciences/MethBat) |
| SVAFotate | 0.2.0 | `jxprismdocker/prism_svafotate:latest` | [github.com/fakedrtom/SVAFotate](https://github.com/fakedrtom/SVAFotate) |
| SvAnna | 1.0.4 | `svanna_v1.0.4.sif`(downloaded via `download_refs.sh` | [github.com/TheJacksonLaboratory/SvAnna](https://github.com/TheJacksonLaboratory/SvAnna) |
| SURVIVOR | 1.0.6.2 | `mgibio/survivor-cwl:1.0.6.2` | [github.com/fritzsedlazeck/SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) |
| R (tidyverse/dplyr) | 4.3 | `rocker/tidyverse:4.3` | [rocker-project.org/](https://rocker-project.org) |
| pandas | 2.2.1 | `quay.io/biocontainers/pandas:2.2.1` | [github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)|

<br>

## (III) Workflow: `lrs_asm_single`

Single-sample *de novo* genome assembly and structural variant calling from PacBio HiFi FASTQ. Assembles with HiFiasm, calls variants with PAV, and prioritizes SVs by phenotype using SvAnna.

#### Sample sheet .csv set-up (requires headers shown below)

```csv
sID,hifi_fastq,HPO
SampleID,/path/to/hifi_fastq,/path/to/HPO.txt
```
>NOTE: HPO.txt is a line-separated list of HPO phenotype terms (e.g. `HP:0001250`)

#### Run `lrs_asm_single`
```bash
sbatch --job-name=nf-<workflow> \
       --time=48:00:00 \
       --mem=16G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w lrs_asm_single -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/ \
                   --refs_path refs/"
```

#### Steps
| Step | Subworkflow | Tools |
|---|---|---|
| 00 | Input validation | — |
| 01 | De novo assembly | HiFiasm, samtools |
| 02 | Variant calling from assembly | PAV |
| 03 | Split PAV output by variant size (`<`50 bp for small variants and >50 bp for SVs) | bcftools, tabix |
| 04 | SV population frequency annotation | SVAFotate | 
| 05 | Phenotype-driven SV prioritization | SvAnna, R (dplyr) |


#### Output Structure
```
<outputDir>/
└── <sampleID>/
    └── asm/
        └── <YYYYMMDD>/
            ├── hifiasm/        assembled contigs (.fa, .gfa) and haplotype GFAs
            ├── pav/            PAV VCFs — full, small variants, and SVs (reformatted for Emedgene ingestion)
            ├── svafotate/          SVAFotate-annotated SV VCFs            
            └── svanna/         SvAnna prioritized SVs and filtered candidates
```
>NOTE: **Emedgene outputs:** VCFs containing `_Emedgene` files are reformatted copies for ingestion into the [Emedgene](https://www.emedgene.com/) clinical interpretation platform. If not using Emedgene, use primary VCFs.

>NOTE: **SvAnna outputs:** Full SvAnna results (csv/html/vcf) are available on PAV from both unfiltered and filtered (AF < 0.01, `*SVAFotate-RARE-UNIQUE.vcf`) SV VCFs. The `SvAnna_PAV_Candidates.csv` file, however, is specifically generated using `*SVAFotate-RARE-UNIQUE.vcf` with a SvAnna psv score >= 4.0. To adjust the psv threshold or change which SvAnna results are used for the candidates file, modify `resources/scripts/filter_svanna_pav.R`. Note that the unfiltered SvAnna run requires more memory than the filtered run; if it fails with OOM, see TROUBLESHOOTING.md section 8.

#### Tool Versions

| Tool | Version | Container | Reference |
|---|---|---|---|
| Yak | 0.1 | `quay.io/biocontainers/yak:0.1--h577a1d6_6` | [github.com/lh3/yak](https://github.com/lh3/yak) |
| HiFiasm | 0.20.0 | `dnalinux/hifiasm:0.20.0` | [github.com/chhylp123/hifiasm](https://github.com/chhylp123/hifiasm) |
| samtools | 1.21 | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` | [github.com/samtools](https://github.com/samtools) |
| PAV | 2.3.4 | `becklab/pav:2.3.4` | [github.com/EichlerLab/pav](https://github.com/EichlerLab/pav) |
| bcftools | 1.19 | `quay.io/biocontainers/bcftools:1.19--h8b25389_1` | [github.com/samtools/bcftools](https://github.com/samtools/bcftools) |
| tabix | 1.11 | `quay.io/biocontainers/tabix:1.11--hdfd78af_0` | [github.com/tabixio/tabix](https://github.com/tabixio/tabix)|
| SVAFotate | 0.2.0 | `jxprismdocker/prism_svafotate:latest` | [github.com/fakedrtom/SVAFotate](https://github.com/fakedrtom/SVAFotate) |
| SvAnna | 1.0.4 | `svanna_v1.0.4.sif`(downloaded via `download_refs.sh` | [github.com/TheJacksonLaboratory/SvAnna](https://github.com/TheJacksonLaboratory/SvAnna) |
| R (tidyverse/dplyr) | 4.3 | `rocker/tidyverse:4.3` | [rocker-project.org](https://rocker-project.org/) |

<br>

## (IV) Workflow: `lrs_asm_trio`

Trio-informed *de novo* genome assembly and structural variant calling from PacBio HiFi FASTQ files for a proband using fastqs from both parents. Uses parental k-mer databases to phase the assembly with HiFiasm.

#### Sample sheet .csv set-up (requires headers shown below)

```csv
sID,proband_hifi_fastq,mat_hifi_fastq,pat_hifi_fastq,HPO
SampleID,/path/to/hifi_fastq,/path/to/maternal_hifi_fastq,/path/to/paternal_hifi_fastq,/path/to/HPO.txt
```
>NOTE: HPO.txt is a line-separated list of HPO phenotype terms (e.g. `HP:0001250`)

#### Run `lrs_asm_trio`
```bash
sbatch --job-name=nf-<workflow> \
       --time=48:00:00 \
       --mem=32G \
       --partition=compute \
       -q batch \
       --wrap="export PATH=/cm/local/apps/apptainer/current/bin:\$PATH && \
               cd /path/to/jax-apml-lrs && \
               ./run.sh -w lrs_asm_trio -p hpc \
                   --csv_path /path/to/samplesheet.csv \
                   --outputDir /flashscratch/\${USER}/jax-apml-lrs/ \
                   --refs_path refs/"
```

#### Steps

| Step | Subworkflow | Tools |
|---|---|---|
| 00 | Input validation | — |
| 01 | Trio-aware de novo assembly | Yak, HiFiasm, samtools |
| 02 | Variant calling from assembly | PAV |
| 03 | Split PAV output by variant size | bcftools, tabix |
| 04 | SV population frequency annotation | SVAFotate |
| 05 | Phenotype-driven SV prioritization | SvAnna, R |

#### Output Structure

```
<outputDir>/
└── <sampleID>/
    └── asm-trio/
        └── <YYYYMMDD>/
            ├── yak/              Yak parental k-mer databases
            ├── hifiasm_trio/     assembled contigs and haplotype FASTAs using trio binning
            ├── pav/              PAV VCFs — full, small variants, and SVs (reformatted for Emedgene ingestion)
            ├── svafotate/          SVAFotate-annotated SV VCFs
            └── svanna/           SvAnna prioritized SVs and filtered candidates
```
>NOTE: **Emedgene outputs:** VCFs containing `_Emedgene` files are reformatted copies for ingestion into the [Emedgene](https://www.emedgene.com/) clinical interpretation platform. If not using Emedgene, use primary VCFs.

>NOTE: **SvAnna outputs:** Full SvAnna results (csv/html/vcf) are available on PAV from both unfiltered and filtered (AF < 0.01, `*SVAFotate-RARE-UNIQUE.vcf`) SV VCFs. The `SvAnna_PAV_Candidates.csv` file, however, is specifically generated using `*SVAFotate-RARE-UNIQUE.vcf` with a SvAnna psv score >= 4.0. To adjust the psv threshold or change which SvAnna results are used for the candidates file, modify `resources/scripts/filter_svanna_pav.R`. Note that the unfiltered SvAnna run requires more memory than the filtered run; if it fails with OOM, see TROUBLESHOOTING.md section 8.

#### Tool Versions

| Tool | Version | Container | Reference |
|---|---|---|---|
| Yak | 0.1 | `quay.io/biocontainers/yak:0.1--h577a1d6_6` | [github.com/lh3/yak](https://github.com/lh3/yak) |
| HiFiasm | 0.20.0 | `dnalinux/hifiasm:0.20.0` | [github.com/chhylp123/hifiasm](https://github.com/chhylp123/hifiasm) |
| samtools | 1.21 | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` | [github.com/samtools](https://github.com/samtools) |
| PAV | 2.3.4 | `becklab/pav:2.3.4` | [github.com/EichlerLab/pav](https://github.com/EichlerLab/pav) |
| bcftools | 1.19 | `quay.io/biocontainers/bcftools:1.19--h8b25389_1` | [github.com/samtools/bcftools](https://github.com/samtools/bcftools) |
| tabix | 1.11 | `quay.io/biocontainers/tabix:1.11--hdfd78af_0` | [github.com/tabixio/tabix](https://github.com/tabixio/tabix)|
| SVAFotate | 0.2.0 | `jxprismdocker/prism_svafotate:latest` | [github.com/fakedrtom/SVAFotate](https://github.com/fakedrtom/SVAFotate) |
| SvAnna | 1.0.4 | `svanna_v1.0.4.sif`(downloaded via `download_refs.sh` | [github.com/TheJacksonLaboratory/SvAnna](https://github.com/TheJacksonLaboratory/SvAnna) |
| R (tidyverse/dplyr) | 4.3 | `rocker/tidyverse:4.3` | [rocker-project.org](https://rocker-project.org/) |

<br>

---

## Citation

Werren, E.A., Vats, P., Rech, G.E., Peracchio, M., King, C., Charnysh, E.J., Gorham, R.D., Audano, P.A., Robinson, P.N., Kelly, M.A., Matson, A.P., Adams, M.D., & Kalsner, L. (2026). Diagnostic utility of clinical genome reanalysis in rare pediatric disorders using long-read sequencing. *(under review)*

---

## Contact

For questions or issues, please open a GitHub issue.
