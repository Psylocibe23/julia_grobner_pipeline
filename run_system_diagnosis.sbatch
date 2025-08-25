#!/bin/bash
#SBATCH --job-name=sysdiag
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
# Write logs in CWD so there’s no missing-dir issue:
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Work in your repo so relative paths are correct:
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline

set -euo pipefail
set -x  # echo commands for debugging

# --- Modules ---
module load sagemath/10.6-singularity

# --- Keep threads sane ---
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

usage() {
  echo "Usage:"
  echo "  Single file:"
  echo "    sbatch run_system_diagnosis.sbatch data/HFE_n32_D64.in [--dim=SECONDS] [--hilbert=SECONDS]"
  echo
  echo "  Directory + array:"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch data [--dim=20] [--hilbert=60]"
  echo
  echo "  Glob + array (quote the glob!):"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch 'data/HFE_n*_D*.in' [--dim=20] [--hilbert=60]"
  echo
  echo "  Manifest + array (one path per line):"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch filelist.txt [--dim=20] [--hilbert=60]"
  echo
  echo "Notes:"
  echo "  - Flags after the file spec are forwarded to the Sage script."
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

INPUT_SPEC="$1"
shift
EXTRA_ARGS=( "$@" )

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CPUs: ${SLURM_CPUS_PER_TASK:-1}  Mem: ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-unknown}}  ArrayID: ${SLURM_ARRAY_TASK_ID:-none}"
echo "[$(date)] CWD: $(pwd)"
which sage || true
sage -v || true

# --- Select file (array-aware) ---
SELECTED_FILE=""
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  declare -a FILES
  if [[ -d "$INPUT_SPEC" ]]; then
    mapfile -t FILES < <(ls -1 "$INPUT_SPEC"/*.in 2>/dev/null | sort || true)
  elif [[ -f "$INPUT_SPEC" && "$INPUT_SPEC" == *.txt ]]; then
    mapfile -t FILES < "$INPUT_SPEC"
  else
    FILES=( $INPUT_SPEC )
    IFS=$'\n' FILES=($(printf "%s\n" "${FILES[@]}" | sort))
  fi
  if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "No input files found for spec: $INPUT_SPEC"
    exit 2
  fi
  IDX=${SLURM_ARRAY_TASK_ID}
  if (( IDX < 0 || IDX >= ${#FILES[@]} )); then
    echo "SLURM_ARRAY_TASK_ID=$IDX out of range (0..$(( ${#FILES[@]}-1 )))"
    exit 3
  fi
  SELECTED_FILE="${FILES[$IDX]}"
else
  SELECTED_FILE="$INPUT_SPEC"
fi

if [[ ! -f "$SELECTED_FILE" ]]; then
  echo "Input file not found: $SELECTED_FILE"
  exit 4
fi

echo "[$(date)] Running diagnosis for: $SELECTED_FILE"
if ((${#EXTRA_ARGS[@]})); then
  echo "[$(date)] Extra args → ${EXTRA_ARGS[*]}"
fi

# Run (script writes its own detailed log alongside printing)
sage scripts/system_diagnosis.sage "$SELECTED_FILE" "${EXTRA_ARGS[@]}"

echo "[$(date)] ✅ Done"
