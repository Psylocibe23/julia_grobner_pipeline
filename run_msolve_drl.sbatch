#!/bin/bash
#SBATCH --job-name=msolve_drl
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=980G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_msolve_drl.sbatch <input.in> [threads]"
  exit 1
fi

INPUT_IN="$1"
THREADS="${2:-32}"
SIF="${MSOLVE_SIF:-containers/msolve.sif}"

mkdir -p logs results scratch

# Modules
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load apptainer || module load singularity

export OMP_NUM_THREADS="$THREADS"
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-/scratch/$USER/$SLURM_JOB_ID/apptainer-tmp}
mkdir -p "$APPTAINER_TMPDIR"

[[ -f "$INPUT_IN" ]] || { echo "Input not found: $INPUT_IN"; exit 2; }

stem=$(basename "$INPUT_IN"); stem="${stem%.*}"
msfile="scratch/${stem}.ms"

# Convert your .in -> .ms (vars line, prime line, then polys comma-separated)
# (Overhead is negligible compared to the solve.)
awk '
NR==1 {print; next}
NR==2 {print; next}
{
  s=$0
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
  if (s=="" || s ~ /^#/) next
  polys[++k]=s
}
END {
  for (i=1;i<=k;i++){
    sep=(i<k?",":"")
    print polys[i] sep
  }
}' "$INPUT_IN" > "$msfile"

echo "[$(date)] Input : $INPUT_IN"
echo "[$(date)] Threads: $THREADS"
echo "[$(date)] SIF   : $SIF"
echo "[$(date)] .ms   : $msfile"

# Run msolve: (1) verbose log + GB printed (DRL/grevlex), all to one stream.
# You’ll get detailed matrices/pairs (to derive solving degree, max matrix sizes),
# and at the end the reduced DRL basis (since -g 2).
out_all="logs/${stem}_msolve_full.out"
gb_out="results/${stem}_DRL_msolve.txt"

apptainer exec "$SIF" msolve \
  -f "$msfile" \
  -t "$THREADS" \
  -v 2 \
  -l 44 \
  -m 30000 \
  -s 22 \
  -r 1 \
  -g 2 | tee "$out_all" > "$gb_out"

echo "[$(date)] Done. Full log: $out_all"
echo "[$(date)] DRL basis: $gb_out"
