#!/bin/bash
#SBATCH --job-name=msolve_drl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_msolve_drl.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
THREADS="${2:-32}"

[[ -f "$INPUT_FILE" ]] || { echo "Input not found: $INPUT_FILE"; exit 2; }

mkdir -p logs results

# Load 'module' if necessary, then try to load msolve if not on PATH
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

if ! command -v msolve >/dev/null 2>&1; then
  module load msolve || {
    echo "msolve not found on PATH and 'module load msolve' failed. Install or load msolve and retry."
    exit 3
  }
fi

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

# Tuning knobs (override at submit: MSOLVE_LINALG=44 MSOLVE_MAX_PAIRS=30000 MSOLVE_HTS=22 sbatch ...)
MSOLVE_LINALG="${MSOLVE_LINALG:-}"     # 1,2,42,44 (default empty = msolve default)
MSOLVE_MAX_PAIRS="${MSOLVE_MAX_PAIRS:-}"   # e.g., 30000
MSOLVE_HTS="${MSOLVE_HTS:-}"               # e.g., 22

# Job-local scratch for intermediates
export TMPDIR="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

STEM="$(basename "$INPUT_FILE")"
STEM="${STEM%.*}"
TS="$(date +%Y%m%d_%H%M%S)"

TMP_IN="$TMPDIR/${STEM}.msin"
TMP_GB_RAW="$TMPDIR/${STEM}_gb_raw.txt"
TMP_SOL="$TMPDIR/${STEM}_sol.out"        # msolve requires -o; we discard solutions
LOG_FILE="logs/${STEM}_msolve_${TS}.log"
OUT_FILE="results/${STEM}_msolve_DRL_${TS}.txt"

# --- Convert your .in to msolve format (comma-separated equations) ---
#   1st line: variables (kept as is)
#   2nd line: p (kept as is)
#   3+: polynomials; we collect non-empty, non-comment lines and add commas between them
awk '
  NR==1 { print; next }
  NR==2 { print; next }
  NR>=3 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/ { arr[++n]=$0 }
  END {
    for (i=1;i<=n;i++) {
      printf "%s%s\n", arr[i], (i<n ? "," : "")
    }
  }
' "$INPUT_FILE" > "$TMP_IN"

# Count #equations for the header (non-empty, non-comment lines from line 3)
NUM_EQS=$(awk 'NR>=3 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/ {n++} END{print n+0}' "$INPUT_FILE")
VARS_LINE=$(sed -n '1p' "$INPUT_FILE")
P_FIELD=$(sed -n '2p' "$INPUT_FILE" | tr -d '[:space:]')

# --- Build msolve CLI ---
args=(-f "$TMP_IN" -o "$TMP_SOL" -t "$THREADS" -g 2 -v 2)
[[ -n "$MSOLVE_LINALG"   ]] && args+=(-l "$MSOLVE_LINALG")
[[ -n "$MSOLVE_MAX_PAIRS" ]] && args+=(-m "$MSOLVE_MAX_PAIRS")
[[ -n "$MSOLVE_HTS"      ]] && args+=(-s "$MSOLVE_HTS")

echo "[$(date)] msolve args: ${args[*]}"
echo "[$(date)] Logging to:  $LOG_FILE"
echo "[$(date)] DRL basis ->  $OUT_FILE"

# --- Run msolve ---
# Per the man page: -g 2 prints the *reduced* grevlex (DRL) basis; -t sets threads; -l/-m/-s tune LA/HTS/pairs.
# We capture the printed DRL basis on STDOUT and the verbose progress on STDERR.
msolve "${args[@]}" > "$TMP_GB_RAW" 2> "$LOG_FILE"

# --- Postprocess and wrap into your pipeline-friendly DRL file ---
{
  echo "# Groebner basis (msolve, DRL) computed for $INPUT_FILE"
  echo "# Variables: $VARS_LINE"
  echo "# Field: GF(${P_FIELD})"
  echo "# Order: degrevlex"
  echo "# Basis: reduced=True"
  echo "# Number of input equations: ${NUM_EQS}"
  echo "# --- Groebner basis ---"
  # msolve prints one polynomial per line, commas separated; strip trailing commas and blank lines
  sed -E 's/,[[:space:]]*$//' "$TMP_GB_RAW" | awk 'NF>0'
} > "$OUT_FILE"

echo "[$(date)] ✅ Done"
echo "Result: $OUT_FILE"
echo "Log:    $LOG_FILE"
