#!/usr/bin/env bash
# =============================================================================
#  run.sh — Unified launcher for jax-apml-lrs Nextflow workflows
# =============================================================================
#
#  Usage:
#      ./run.sh -w <workflow> -p <profile> --csv_path <samplesheet> \
#               --outputDir <output_dir> [options]
#
#  Required:
#      -w, --workflow   Workflow to run. One of:
#                           lrs_read
#                           lrs_asm_single
#                           lrs_asm_trio
#                           lrs_pbmerge
#
#      -p, --profile    Execution profile. One of:
#                           standard   local execution
#                           hpc        SLURM + Singularity
#                           gcb        Google Cloud Batch
#
#      --csv_path       Path to the input samplesheet CSV.
#
#      --outputDir      Path to the output directory.
#
#  Options:
#      -r, --resume     Resume a previous run (-resume flag passed to Nextflow).
#      -h, --help       Show this help message and exit.
#
#  Any additional arguments are passed directly to Nextflow.
#
#  Examples:
#      ./run.sh -w lrs_read -p hpc \
#               --csv_path /path/to/samples.csv \
#               --outputDir /path/to/output/
#
#      ./run.sh -w lrs_asm_trio -p gcb \
#               --csv_path /path/to/samples.csv \
#               --outputDir gs://my-bucket/results/ \
#               --resume
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW_BIN="${SCRIPT_DIR}/nextflow"
WORKFLOWS_DIR="${SCRIPT_DIR}/workflows"

VALID_WORKFLOWS=(lrs_read lrs_asm_single lrs_asm_trio lrs_pbmerge)
VALID_PROFILES=(standard hpc gcb)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'; NC='\033[0m'

usage() {
cat <<EOF

  Usage: ./run.sh -w <workflow> -p <profile> --csv_path <path> --outputDir <path> [options]

  Required:
    -w, --workflow    ${VALID_WORKFLOWS[*]}
    -p, --profile     ${VALID_PROFILES[*]}
    --csv_path        Path to input samplesheet CSV
    --outputDir       Path to output directory

  Options:
    -r, --resume      Resume a previous Nextflow run
    -h, --help        Show this message

  Any additional arguments are forwarded to Nextflow.

EOF
}

err()  { echo -e "${RED}ERROR: $*${NC}" >&2; }
warn() { echo -e "${YELLOW}WARN:  $*${NC}" >&2; }
info() { echo -e "${GREEN}INFO:  $*${NC}"; }

contains() {
    local target="$1"; shift
    for item in "$@"; do [[ "$item" == "$target" ]] && return 0; done
    return 1
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
WORKFLOW=''
PROFILE=''
CSV_PATH=''
OUTPUT_DIR=''
RESUME=''
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workflow)   WORKFLOW="$2";    shift 2 ;;
        -p|--profile)    PROFILE="$2";     shift 2 ;;
        --csv_path)      CSV_PATH="$2";    shift 2 ;;
        --outputDir)     OUTPUT_DIR="$2";  shift 2 ;;
        -r|--resume)     RESUME="-resume ${2:-}"; [[ -n "${2:-}" && "${2}" != -* ]] && shift; shift   ;;
        -h|--help)       usage; exit 0              ;;
        *)               EXTRA_ARGS+=("$1"); shift  ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
missing=0

if [[ -z "$WORKFLOW" ]]; then
    err "Workflow not specified (-w / --workflow)."
    missing=1
fi

if [[ -z "$PROFILE" ]]; then
    err "Profile not specified (-p / --profile)."
    missing=1
fi

if [[ -z "$CSV_PATH" ]]; then
    err "Samplesheet not specified (--csv_path)."
    missing=1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    err "Output directory not specified (--outputDir)."
    missing=1
fi

[[ $missing -eq 1 ]] && { usage; exit 1; }

# ---------------------------------------------------------------------------
# Validate workflow and profile values
# ---------------------------------------------------------------------------
if ! contains "$WORKFLOW" "${VALID_WORKFLOWS[@]}"; then
    err "Unknown workflow: '$WORKFLOW'. Must be one of: ${VALID_WORKFLOWS[*]}"
    exit 1
fi

if ! contains "$PROFILE" "${VALID_PROFILES[@]}"; then
    err "Unknown profile: '$PROFILE'. Must be one of: ${VALID_PROFILES[*]}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Resolve workflow paths
# ---------------------------------------------------------------------------
WORKFLOW_DIR="${WORKFLOWS_DIR}/${WORKFLOW}"
WORKFLOW_NF="${WORKFLOW_DIR}/${WORKFLOW}.nf"
WORKFLOW_CONFIG="${WORKFLOW_DIR}/nextflow.config"

if [[ ! -f "$WORKFLOW_NF" ]]; then
    err "Workflow file not found: $WORKFLOW_NF"
    exit 1
fi

if [[ ! -f "$WORKFLOW_CONFIG" ]]; then
    err "Config file not found: $WORKFLOW_CONFIG"
    exit 1
fi

# ---------------------------------------------------------------------------
# Validate samplesheet exists (local paths only — skip for gs:// paths)
# ---------------------------------------------------------------------------
if [[ "$CSV_PATH" != gs://* && ! -f "$CSV_PATH" ]]; then
    err "Samplesheet not found: $CSV_PATH"
    exit 1
fi

# ---------------------------------------------------------------------------
# Ensure Nextflow bootstrap is present
# ---------------------------------------------------------------------------
if [[ ! -f "$NEXTFLOW_BIN" ]]; then
    err "Nextflow executable not found at: $NEXTFLOW_BIN"
    echo "Run the following to install it in this directory:" >&2
    echo "    curl -s https://get.nextflow.io | bash" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Print run summary
# ---------------------------------------------------------------------------
info "============================================================"
info "  jax-apml-lrs"
info "------------------------------------------------------------"
info "  Workflow  : $WORKFLOW"
info "  Profile   : $PROFILE"
info "  Samplesheet: $CSV_PATH"
info "  Output dir: $OUTPUT_DIR"
[[ -n "$RESUME" ]] && info "  Resume    : yes"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && info "  Extra args: ${EXTRA_ARGS[*]}"
info "============================================================"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"

# Preserve runDate across resumes using a sentinel file
SENTINEL="${OUTPUT_DIR}/.rundate"
if [[ -n "$RESUME" && -f "$SENTINEL" ]]; then
    RUN_DATE=$(cat "$SENTINEL")
    info "  Resuming run dated: $RUN_DATE"
else
    RUN_DATE=$(date +%Y%m%d)
    echo "$RUN_DATE" > "$SENTINEL"
fi

# Clear stale Nextflow session locks left by abrupt SLURM kills
export NXF_OPTS="-Xms2g -Xmx8g"
"$NEXTFLOW_BIN" -log "${OUTPUT_DIR}.nextflow.${WORKFLOW}.log" run "$WORKFLOW_NF" \
    -w "${OUTPUT_DIR}work/" \
    -c "$WORKFLOW_CONFIG" \
    -profile "$PROFILE" \
    --csv_path "$CSV_PATH" \
    --outputDir "$OUTPUT_DIR" \
    --runDate "$RUN_DATE" \
    ${RESUME} \
    "${EXTRA_ARGS[@]}" || { NXF_EXIT=$?; [[ $NXF_EXIT -eq 2 ]] && true || exit $NXF_EXIT; }

# Post-run publish rescue: copy any unpublished svanna files from work dir
if grep -q "Execution complete" "${OUTPUT_DIR}.nextflow.${WORKFLOW}.log" 2>/dev/null; then
    SAMPLE_ID=$(awk -F',' 'NR==2{print $1}' "$CSV_PATH")
    if [[ "$WORKFLOW" == "lrs_asm_single" || "$WORKFLOW" == "lrs_asm_trio" ]]; then
        OUT_SVANNA="${OUTPUT_DIR}${SAMPLE_ID}/asm/$(cat "${OUTPUT_DIR}.rundate" 2>/dev/null)/svanna"
        mkdir -p "$OUT_SVANNA"
        find "${OUTPUT_DIR}work" \( -name "*_svanna.csv" -o -name "*_svanna.html" -o -name "*_svanna.vcf.gz" \) 2>/dev/null | while read f; do
            bn=$(basename "$f")
            [[ ! -f "${OUT_SVANNA}/${bn}" ]] && cp "$f" "${OUT_SVANNA}/${bn}" && echo "Rescued: ${bn}"
        done
    fi
fi
