#!/usr/bin/env bash
# =============================================================================
#  download_refs.sh — Download large reference files for jax-apml-lrs
# =============================================================================
#
#  Downloads reference files that are too large to bundle in the repository.
#  Small reference files are already bundled under refs/ and do not need
#  to be downloaded.
#
#  Files downloaded:
#      hg38.no_alt.fa                                        ~3.0 GB
#      hg38.no_alt.fa.fai                                    (generated)
#      human_GRCh38_no_alt_analysis_set.trf.bed              ~7.3 MB
#      SVAFotate_core_SV_popAFs.GRCh38.v4.1_LR_v2.0_INS1bp.bed.gz  ~451 MB
#
#  Prerequisites:
#      docker or singularity (used to run samtools faidx — no local
#      samtools installation required)
#
#  Usage:
#      # Download to the default refs/ directory in this repo
#      ./download_refs.sh
#
#      # Download to a custom directory
#      ./download_refs.sh --outdir /path/to/refs/
#
#      # Use Singularity instead of Docker
#      ./download_refs.sh --container-runtime singularity
#
#  Options:
#      --outdir              Output directory for reference files
#                            (default: refs/ in this repo)
#      --container-runtime   Container runtime for samtools faidx:
#                            docker (default) or singularity
#      -h, --help            Show this message
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths and defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTDIR="${SCRIPT_DIR}/refs"
CONTAINER_RUNTIME="apptainer"
SAMTOOLS_IMAGE="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

# Release tag where the SVAFotate BED is hosted as a GitHub Release asset.
# Update this if a new release is published.
RELEASE_TAG="v1.0.0"
GITHUB_REPO="TheJacksonLaboratory/jax-apml-lrs"

FASTA_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz"
TRF_URL="https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed"
SVAFOTATE_URL="https://github.com/${GITHUB_REPO}/releases/download/${RELEASE_TAG}/SVAFotate_core_SV_popAFs.GRCh38.v4.1_LR_v2.0_INS1bp.bed.gz"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'; BOLD='\033[1m'; NC='\033[0m'

err()     { echo -e "${RED}ERROR: $*${NC}" >&2; }
warn()    { echo -e "${YELLOW}WARN:  $*${NC}" >&2; }
info()    { echo -e "${GREEN}INFO:  $*${NC}"; }
section() { echo -e "\n${BOLD}$*${NC}"; }

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)              OUTDIR="$2";            shift 2 ;;
        --container-runtime)   CONTAINER_RUNTIME="$2"; shift 2 ;;
        -h|--help)
            sed -n '/^# ====/,/^# ====/p' "$0" | sed 's/^# \{0,2\}//'
            exit 0
            ;;
        *)
            err "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate container runtime
# ---------------------------------------------------------------------------
if [[ "$CONTAINER_RUNTIME" != "docker" && "$CONTAINER_RUNTIME" != "singularity" && "$CONTAINER_RUNTIME" != "apptainer" ]]; then
    err "--container-runtime must be 'docker', 'singularity', or 'apptainer'"
    exit 1
fi

APPTAINER_BIN=/cm/local/apps/apptainer/current/bin/apptainer
[[ "$CONTAINER_RUNTIME" == "apptainer" && -x "$APPTAINER_BIN" ]] && CONTAINER_RUNTIME="$APPTAINER_BIN"
if ! command -v "$CONTAINER_RUNTIME" &>/dev/null && [[ ! -x "$CONTAINER_RUNTIME" ]]; then
    err "'$CONTAINER_RUNTIME' not found. Please install it or specify --container-runtime."
    exit 1
fi

# ---------------------------------------------------------------------------
# Create output directory
# ---------------------------------------------------------------------------
mkdir -p "$OUTDIR"
info "Output directory: $OUTDIR"

# ---------------------------------------------------------------------------
# Helper: run samtools faidx through container
# ---------------------------------------------------------------------------
run_samtools_faidx() {
    local fasta="$1"
    local fasta_dir
    fasta_dir="$(dirname "$fasta")"
    local fasta_file
    fasta_file="$(basename "$fasta")"

    info "Indexing ${fasta_file} with samtools faidx (via ${CONTAINER_RUNTIME})..."

    if [[ "$CONTAINER_RUNTIME" == "docker" ]]; then
        docker run --rm \
            -v "${fasta_dir}:/data" \
            "$SAMTOOLS_IMAGE" \
            samtools faidx "/data/${fasta_file}"
    elif [[ "$CONTAINER_RUNTIME" == "singularity" ]]; then
        singularity exec
    else
        $CONTAINER_RUNTIME exec \
            --bind "${fasta_dir}:/data" \
            "docker://${SAMTOOLS_IMAGE}" \
            samtools faidx "/data/${fasta_file}"
    fi
}

# ---------------------------------------------------------------------------
# 1. hg38.no_alt.fa
# ---------------------------------------------------------------------------
section "[ 1 / 3 ]  Reference genome (hg38.no_alt.fa)"

FASTA_GZ="${OUTDIR}/hg38.no_alt.fa.gz"
FASTA="${OUTDIR}/hg38.no_alt.fa"
FASTA_FAI="${OUTDIR}/hg38.no_alt.fa.fai"

if [[ -f "$FASTA" ]]; then
    info "hg38.no_alt.fa already exists — skipping download."
else
    info "Downloading hg38.no_alt.fa.gz (~3.0 GB, this may take a while)..."
    wget -c --show-progress -O "$FASTA_GZ" "$FASTA_URL"
    info "Decompressing..."
    gunzip "$FASTA_GZ"
    info "hg38.no_alt.fa download complete."
fi

if [[ -f "$FASTA_FAI" ]]; then
    info "hg38.no_alt.fa.fai already exists — skipping indexing."
else
    run_samtools_faidx "$FASTA"
    info "hg38.no_alt.fa.fai generated."
fi

# ---------------------------------------------------------------------------
# 2. human_GRCh38_no_alt_analysis_set.trf.bed
# ---------------------------------------------------------------------------
section "[ 2 / 3 ]  PBSV tandem repeat annotations (human_GRCh38_no_alt_analysis_set.trf.bed)"

TRF_BED="${OUTDIR}/human_GRCh38_no_alt_analysis_set.trf.bed"

if [[ -f "$TRF_BED" ]]; then
    info "human_GRCh38_no_alt_analysis_set.trf.bed already exists — skipping download."
else
    info "Downloading human_GRCh38_no_alt_analysis_set.trf.bed (~7.3 MB)..."
    wget -c --show-progress -O "$TRF_BED" "$TRF_URL"
    info "Download complete."
fi

# ---------------------------------------------------------------------------
# 3. SVAFotate BED (GitHub Release asset)
# ---------------------------------------------------------------------------
section "[ 3 / 3 ]  SVAFotate population allele frequency database"

SVAFOTATE_BED="${OUTDIR}/SVAFotate_core_SV_popAFs.GRCh38.v4.1_LR_v2.0_INS1bp.bed.gz"

if [[ -f "$SVAFOTATE_BED" ]]; then
    info "SVAFotate BED already exists — skipping download."
else
    info "Downloading SVAFotate BED from GitHub Release ${RELEASE_TAG} (~451 MB)..."
    wget -c --show-progress -O "$SVAFOTATE_BED" "$SVAFOTATE_URL"
    info "Download complete."
fi

# ---------------------------------------------------------------------------
# 4. SvAnna Apptainer image
# ---------------------------------------------------------------------------
section "[ 4 / 5 ]  SvAnna Apptainer image (svanna_v1.0.4.sif)"

SVANNA_SIF="${OUTDIR}/svanna_v1.0.4.sif"
SVANNA_URL="https://github.com/${GITHUB_REPO}/releases/download/${RELEASE_TAG}/svanna_v1.0.4.sif"

if [[ -f "$SVANNA_SIF" ]]; then
    info "svanna_v1.0.4.sif already exists — skipping download."
else
    info "Downloading SvAnna Apptainer image (~1.1 GB, this may take a while)..."
    wget -c --show-progress -O "$SVANNA_SIF" "$SVANNA_URL"
    info "Download complete."
fi

# ---------------------------------------------------------------------------
# 4. Bundled reference files (small files included in the repo)
# ---------------------------------------------------------------------------
section "[ 5 / 5 ]  Bundled reference files"

BUNDLED_REFS="${SCRIPT_DIR}/refs"

if [[ -d "$BUNDLED_REFS" ]]; then
    info "Copying bundled reference files from repo..."
    cp "${BUNDLED_REFS}/"* "$OUTDIR/"
    info "Bundled reference files copied."
else
    warn "Bundled refs directory not found at ${BUNDLED_REFS} — skipping."
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
section "Done"
echo ""
info "Reference files are ready in: $OUTDIR"
echo ""
echo "  Pass this directory to run.sh with:"
echo "      --refs_path $OUTDIR"
echo ""