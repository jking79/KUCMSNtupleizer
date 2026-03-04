#!/bin/bash
#
# parallelRun.sh - Run cmsRun in parallel for KUCMSNtupleizer miniAOD analysis
#
# Processes ROOT files in parallel using GNU parallel, writing outputs to a
# timestamped EOS directory via xrdcp (similar to CRAB job organisation).
#
# Usage:
#   ./parallelRun.sh <input> [options]
#
# Input (first positional argument):
#   file.root          Single ROOT file (local or xrootd URL)
#   list.txt           Text file with one input path per line (# = comment)
#   '/store/.../*.root' Quoted glob pattern
#
# Examples:
#   ./parallelRun.sh files.txt
#   ./parallelRun.sh files.txt -j 8 -t mySignal
#   ./parallelRun.sh 'root://cmseos.fnal.gov//store/user/me/signal_*.root' -j 4
#   ./parallelRun.sh signal.root -n 100 --no-xrdcp
#   ./parallelRun.sh files.txt --continue          # resume a previous run
#

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
N_JOBS=4
FILES_PER_JOB=1
MAX_EVENTS=-1
EOS_BASE="/store/group/lpcsusylep/anazario"
EOS_SERVER="root://cmseos.fnal.gov"
TAG="kucmsntuple"
FILTER=""           # eventFilter passed to cmsRun (empty = use config default)
LOCAL_DIR=""        # auto-set from /tmp/parallelRun_PID if not given
CONFIG=""           # auto-detected from CMSSW_BASE if not given
DO_MERGE=false
DO_XRDCP=true
DO_CONTINUE=false

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
print_usage() {
    echo "Usage: $0 <input> [options]"
    echo ""
    echo "Input (required, first argument):"
    echo "  file.root              Single ROOT file (local path or xrootd URL)"
    echo "  list.txt               Text file, one input path per line (# = comment)"
    echo "  '/store/.../*.root'    Quoted glob pattern"
    echo ""
    echo "Options:"
    echo "  -c, --config   FILE   cmsRun config (default: auto-detected in CMSSW_BASE)"
    echo "  -f, --filter   NAME   eventFilter passed to cmsRun (default: use config default)"
    echo "  -j, --nJobs    N      Parallel workers (default: 4)"
    echo "  -p, --files-per-job N Input files processed per job (default: 1)"
    echo "  -n, --maxEvents N     Events per job, -1 = all (default: -1)"
    echo "  -o, --outDir   DIR    EOS base directory"
    echo "                        (default: /store/group/lpcsusylep/anazario)"
    echo "  -t, --tag      TAG    Job tag used in output naming (default: kucmsntuple)"
    echo "  -l, --localDir DIR    Local working directory (default: /tmp/parallelRun_PID)"
    echo "      --merge           Run hadd to merge all per-job outputs"
    echo "      --no-xrdcp        Keep outputs local only, skip EOS transfer"
    echo "      --continue        Skip files that succeeded in a previous run"
    echo "  -h, --help            Show this message"
    echo ""
    echo "Output structure on EOS:"
    echo "  <outDir>/<tag>_<identifier>/<YYYYMMDD_HHMMSS>/"
    echo "    <tag>_<identifier>_0.root"
    echo "    <tag>_<identifier>_1.root"
    echo "    ..."
    echo "    <tag>_<identifier>_merged.root   (unless --no-merge)"
}

# ---------------------------------------------------------------------------
# cmsenv check — prompt the user if cmsRun is not in PATH
# ---------------------------------------------------------------------------
check_cmsenv() {
    if command -v cmsRun &>/dev/null; then
        return 0
    fi

    echo -e "${YELLOW}WARNING: cmsRun not found in PATH. Have you run 'cmsenv'?${NC}"
    echo -n "  Attempt to apply cmsenv now? [y/N] "
    read -r _answer
    if [[ "$_answer" =~ ^[Yy]$ ]]; then
        if [[ -z "${CMSSW_BASE:-}" ]]; then
            echo -e "${RED}Error: CMSSW_BASE is not set. Cannot run cmsenv automatically.${NC}"
            echo "  Please cd into your CMSSW release area and run 'cmsenv' manually."
            exit 1
        fi
        # eval the SCRAM runtime into this shell so cmsRun is available here and
        # in all child processes spawned by this script.
        eval "$(cd "$CMSSW_BASE" && scramv1 runtime -sh)"
        if ! command -v cmsRun &>/dev/null; then
            echo -e "${RED}Error: cmsRun still not found after cmsenv attempt.${NC}"
            exit 1
        fi
        echo -e "${GREEN}cmsenv applied successfully.${NC}"
    else
        echo "Aborting. Please run cmsenv and try again."
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Auto-detect the parallel config from CMSSW_BASE
# ---------------------------------------------------------------------------
detect_config() {
    local candidate="${CMSSW_BASE:-}/src/KUCMSNtupleizer/KUCMSNtupleizer/test/llpgana_miniaod_parallel.py"
    if [[ -f "$candidate" ]]; then
        echo "$candidate"
    else
        echo ""
    fi
}

# ---------------------------------------------------------------------------
# Build an input file list from a single .root, a .txt, or a glob.
# Sets globals FILES_LIST and IDENTIFIER directly (no subshell).
# ---------------------------------------------------------------------------
IDENTIFIER=""
FILES_LIST=""
build_file_list() {
    local input="$1"
    FILES_LIST=$(mktemp)

    if [[ "$input" == *.txt && -f "$input" ]]; then
        # Text file list — strip comments and blank lines
        grep -v '^\s*#' "$input" | grep -v '^\s*$' > "$FILES_LIST"
        IDENTIFIER=$(basename "$input" .txt)

    elif [[ "$input" == *.root || "$input" == root://* || "$input" == file:* ]]; then
        # Single ROOT file (local or xrootd)
        echo "$input" > "$FILES_LIST"
        local base
        base=$(basename "$input" .root)
        base="${base#file:}"
        IDENTIFIER="${base:0:50}"   # cap length

    else
        # Treat as a glob pattern — try shell expansion first, then xrdfs
        local found=()
        # shellcheck disable=SC2206
        found=( $input )    # intentional word-splitting for glob expansion
        if [[ ${#found[@]} -eq 0 || ! -e "${found[0]}" ]]; then
            # Try xrdfs glob for EOS patterns
            local server path_pattern
            if [[ "$input" == root://* ]]; then
                server="${input%%//store*}"
                path_pattern="${input#${server}}"
                xrdfs "$server" ls "${path_pattern%/*}" 2>/dev/null \
                    | grep "$(basename "$path_pattern" | tr '*' '.')" \
                    | sed "s|^|${server}|" > "$FILES_LIST" || true
            fi
        else
            printf '%s\n' "${found[@]}" > "$FILES_LIST"
        fi

        if [[ ! -s "$FILES_LIST" ]]; then
            echo -e "${RED}Error: No files matched: $input${NC}" >&2
            rm -f "$FILES_LIST"
            exit 1
        fi
        IDENTIFIER="glob"
    fi
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    print_usage
    exit 1
fi

INPUT="$1"
shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)     CONFIG="$2";           shift 2 ;;
        -f|--filter)     FILTER="$2";           shift 2 ;;
        -j|--nJobs)      N_JOBS="$2";           shift 2 ;;
        -p|--files-per-job) FILES_PER_JOB="$2"; shift 2 ;;
        -n|--maxEvents)  MAX_EVENTS="$2";  shift 2 ;;
        -o|--outDir)     EOS_BASE="$2";    shift 2 ;;
        -t|--tag)        TAG="$2";         shift 2 ;;
        -l|--localDir)   LOCAL_DIR="$2";   shift 2 ;;
        --merge)         DO_MERGE=true;    shift   ;;
        --no-xrdcp)      DO_XRDCP=false;   shift   ;;
        --continue)      DO_CONTINUE=true; shift   ;;
        -h|--help)       print_usage; exit 0 ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Environment and config checks
# ---------------------------------------------------------------------------
check_cmsenv

if [[ -z "$CONFIG" ]]; then
    CONFIG=$(detect_config)
    if [[ -z "$CONFIG" ]]; then
        echo -e "${RED}Error: Could not auto-detect config. Use -c / --config to specify.${NC}"
        exit 1
    fi
fi

if [[ ! -f "$CONFIG" ]]; then
    echo -e "${RED}Error: Config file not found: $CONFIG${NC}"
    exit 1
fi

if ! command -v parallel &>/dev/null; then
    echo -e "${RED}Error: GNU parallel is not installed.${NC}"
    echo "  Install with:  conda install -c conda-forge parallel"
    echo "              or brew install parallel  (macOS)"
    exit 1
fi

# ---------------------------------------------------------------------------
# Build file list
# ---------------------------------------------------------------------------
build_file_list "$INPUT"
N_FILES=$(wc -l < "$FILES_LIST" | tr -d ' ')

if [[ "$N_FILES" -eq 0 ]]; then
    echo -e "${RED}Error: No input files found.${NC}"
    rm -f "$FILES_LIST"
    exit 1
fi

# ---------------------------------------------------------------------------
# Timestamped run directory (CRAB-style)
# ---------------------------------------------------------------------------
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_NAME="${TAG}_${IDENTIFIER}"

if [[ -z "$LOCAL_DIR" ]]; then
    LOCAL_DIR="/tmp/parallelRun_$$/${RUN_NAME}/${TIMESTAMP}"
fi

LOG_DIR="${LOCAL_DIR}/logs"
EOS_RUN_DIR="${EOS_BASE}/${RUN_NAME}/${TIMESTAMP}"

mkdir -p "$LOCAL_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# --continue: filter files that already completed successfully
# ---------------------------------------------------------------------------
N_SKIPPED=0
if [[ "$DO_CONTINUE" == true ]]; then
    FILTERED=$(mktemp)
    JOB_IDX=0
    while IFS= read -r f; do
        LOG="${LOG_DIR}/job_${JOB_IDX}.log"
        if grep -q "CMSRUN_EXIT_SUCCESS" "$LOG" 2>/dev/null; then
            N_SKIPPED=$((N_SKIPPED + 1))
        else
            echo "$f" >> "$FILTERED"
        fi
        JOB_IDX=$((JOB_IDX + 1))
    done < "$FILES_LIST"
    mv "$FILTERED" "$FILES_LIST"
fi

N_TO_RUN=$(wc -l < "$FILES_LIST" | tr -d ' ')

# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  KUCMSNtupleizer Parallel Runner${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "  Input:           $INPUT"
echo "  Total files:     $N_FILES"
[[ "$DO_CONTINUE" == true ]] && echo "  Already done:    $N_SKIPPED"
echo "  To process:      $N_TO_RUN"
echo "  Parallel jobs:   $N_JOBS"
echo "  Files per job:   $FILES_PER_JOB"
echo "  Max events/job:  $MAX_EVENTS"
echo "  Config:          $CONFIG"
echo "  Event filter:    ${FILTER:-<config default>}"
echo "  Run name:        $RUN_NAME"
echo "  Timestamp:       $TIMESTAMP"
echo "  Local work dir:  $LOCAL_DIR"
[[ "$DO_XRDCP" == true ]] && \
echo "  EOS output:      ${EOS_SERVER}/${EOS_RUN_DIR}"
echo ""

if [[ "$N_TO_RUN" -eq 0 ]]; then
    echo -e "${GREEN}Nothing left to process.${NC}"
    rm -f "$FILES_LIST"
    exit 0
fi

# ---------------------------------------------------------------------------
# Create EOS output directory up front
# ---------------------------------------------------------------------------
if [[ "$DO_XRDCP" == true ]]; then
    echo -e "${CYAN}Creating EOS directory...${NC}"
    if ! xrdfs "$EOS_SERVER" mkdir -p "$EOS_RUN_DIR"; then
        echo -e "${RED}Error: Could not create EOS directory: ${EOS_RUN_DIR}${NC}"
        echo "  Is your grid proxy valid? Try: voms-proxy-init --voms cms"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Write a per-job helper script consumed by GNU parallel
# ---------------------------------------------------------------------------
TEMP_JOB=$(mktemp /tmp/kucms_job_XXXXXX.sh)
# Export all variables the job needs via the environment so parallel can
# forward them cleanly without quoting headaches.
export _KP_CONFIG="$CONFIG"
export _KP_LOCAL_DIR="$LOCAL_DIR"
export _KP_TAG="$TAG"
export _KP_IDENTIFIER="$IDENTIFIER"
export _KP_MAX_EVENTS="$MAX_EVENTS"
export _KP_FILTER="$FILTER"
export _KP_EOS_SERVER="$EOS_SERVER"
export _KP_EOS_RUN_DIR="$EOS_RUN_DIR"
export _KP_DO_XRDCP="$DO_XRDCP"

cat > "$TEMP_JOB" << 'JOB_EOF'
#!/bin/bash
# Called by GNU parallel as:  bash TEMP_JOB START_IDX CHUNK_FILE
START_IDX="$1"
CHUNK_FILE="$2"

mkdir -p "${_KP_LOCAL_DIR}/logs"

_FILTER_ARG=()
[[ -n "$_KP_FILTER" ]] && _FILTER_ARG=( "eventFilter=${_KP_FILTER}" )

LOCAL_IDX=0
while IFS= read -r INPUT_FILE; do
    [[ -z "$INPUT_FILE" ]] && continue

    GLOBAL_IDX=$((START_IDX + LOCAL_IDX))
    OUTFILE="${_KP_TAG}_${_KP_IDENTIFIER}_${GLOBAL_IDX}.root"
    LOCAL_OUT="${_KP_LOCAL_DIR}/${OUTFILE}"
    LOG="${_KP_LOCAL_DIR}/logs/job_${GLOBAL_IDX}.log"

    echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: starting $(basename "${INPUT_FILE}")"

    if cmsRun "$_KP_CONFIG" \
          inputFiles="$INPUT_FILE" \
          outputFile="$LOCAL_OUT"  \
          maxEvents="$_KP_MAX_EVENTS" \
          "${_FILTER_ARG[@]}" \
          >> "$LOG" 2>&1; then

        echo "CMSRUN_EXIT_SUCCESS" >> "$LOG"
        echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: done → ${OUTFILE}"

        if [[ "$_KP_DO_XRDCP" == "true" ]]; then
            echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: transferring to EOS..."
            if xrdcp "$LOCAL_OUT" "${_KP_EOS_SERVER}/${_KP_EOS_RUN_DIR}/${OUTFILE}" >> "$LOG" 2>&1; then
                echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: xrdcp done"
                rm -f "$LOCAL_OUT"
            else
                echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: xrdcp FAILED — local file kept at ${LOCAL_OUT}"
            fi
        fi
    else
        EXIT_CODE=$?
        echo "[$(date '+%H:%M:%S')] Job ${GLOBAL_IDX}: FAILED (exit ${EXIT_CODE})"
        if [[ -f "$LOG" ]]; then
            echo "--- error from log (${LOG}) ---"
            grep -A 10 "An exception of category\|Fatal Exception\|%MSG-e\|Throw with message\|cmsRun: error" "$LOG" | head -60
            echo "-------------------------------"
        else
            echo "(log not found: ${LOG})"
        fi
    fi

    LOCAL_IDX=$((LOCAL_IDX + 1))
done < "$CHUNK_FILE"
JOB_EOF
chmod +x "$TEMP_JOB"

# ---------------------------------------------------------------------------
# Split file list into chunks of FILES_PER_JOB, build job list:
# "START_GLOBAL_IDX CHUNK_FILE" one per line.
# Global index starts at N_SKIPPED so output filenames stay consistent when
# --continue is used across multiple invocations.
# ---------------------------------------------------------------------------
CHUNK_DIR=$(mktemp -d)
JOB_LIST=$(mktemp)
chunk_idx=0
file_pos=0
current_chunk=""

while IFS= read -r f; do
    local_idx=$((file_pos % FILES_PER_JOB))
    if [[ $local_idx -eq 0 ]]; then
        current_chunk="${CHUNK_DIR}/chunk_$(printf '%04d' $chunk_idx)"
        start_global=$((N_SKIPPED + file_pos))
        echo "$start_global $current_chunk" >> "$JOB_LIST"
        chunk_idx=$((chunk_idx + 1))
    fi
    echo "$f" >> "$current_chunk"
    file_pos=$((file_pos + 1))
done < "$FILES_LIST"
rm -f "$FILES_LIST"

# ---------------------------------------------------------------------------
# Run with GNU parallel
# ---------------------------------------------------------------------------
N_CHUNKS=$(wc -l < "$JOB_LIST" | tr -d ' ')
echo -e "${YELLOW}Launching ${N_CHUNKS} chunks (${N_TO_RUN} files, ${FILES_PER_JOB} per job) with ${N_JOBS} workers...${NC}"
echo ""

START_TIME=$(date +%s)

parallel --bar -j "$N_JOBS" --colsep ' ' \
    bash "$TEMP_JOB" {1} {2} \
    :::: "$JOB_LIST"

PARALLEL_EXIT=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

rm -f "$TEMP_JOB" "$JOB_LIST"
rm -rf "$CHUNK_DIR"

# ---------------------------------------------------------------------------
# Results summary
# ---------------------------------------------------------------------------
echo ""
N_DONE=$(ls -1 "${LOCAL_DIR}/${TAG}_${IDENTIFIER}_"[0-9]*.root 2>/dev/null | wc -l | tr -d ' ')
N_FAILED=$((N_TO_RUN - N_DONE))

echo -e "${GREEN}Processing finished in ${ELAPSED}s${NC}"
echo "  Succeeded: ${N_DONE} / ${N_TO_RUN}"
[[ $N_FAILED -gt 0 ]] && \
    echo -e "  ${RED}Failed:    ${N_FAILED}${NC}  (logs in ${LOG_DIR})"

# ---------------------------------------------------------------------------
# Merge with hadd
# ---------------------------------------------------------------------------
if [[ "$DO_MERGE" == true && "$N_DONE" -gt 0 ]]; then
    echo ""
    echo -e "${CYAN}Merging with hadd...${NC}"
    MERGED="${LOCAL_DIR}/${TAG}_${IDENTIFIER}_merged.root"
    # shellcheck disable=SC2086
    if hadd -f "$MERGED" "${LOCAL_DIR}/${TAG}_${IDENTIFIER}_"[0-9]*.root; then
        SIZE=$(ls -lh "$MERGED" | awk '{print $5}')
        echo -e "${GREEN}Merged: $(basename "$MERGED")  (${SIZE})${NC}"
        if [[ "$DO_XRDCP" == "true" ]]; then
            echo -e "${CYAN}Transferring merged file to EOS...${NC}"
            if xrdcp "$MERGED" "${EOS_SERVER}/${EOS_RUN_DIR}/$(basename "$MERGED")"; then
                echo -e "${GREEN}EOS: ${EOS_SERVER}/${EOS_RUN_DIR}/$(basename "$MERGED")${NC}"
            else
                echo -e "${RED}xrdcp of merged file failed.${NC}"
            fi
        fi
    else
        echo -e "${RED}hadd failed.${NC}"
    fi
fi

# ---------------------------------------------------------------------------
# Final paths
# ---------------------------------------------------------------------------
echo ""
echo -e "${GREEN}Done!${NC}"
echo "  Run tag:  ${RUN_NAME}/${TIMESTAMP}"
echo "  Local:    ${LOCAL_DIR}/"
[[ "$DO_XRDCP" == "true" ]] && \
    echo "  EOS:      ${EOS_SERVER}/${EOS_RUN_DIR}/"

exit $PARALLEL_EXIT
