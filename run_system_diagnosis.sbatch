#!/bin/bash
#SBATCH --job-name=sysdiag
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/sysdiag_%j.out
#SBATCH --error=logs/sysdiag_%j.err

set -euo pipefail

# --- Modules ---
module load sagemath/10.6-singularity

# --- Keep threads sane (Sage will spawn algebra kernels; we don't want BLAS storms) ---
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p logs

usage() {
  echo "Usage:"
  echo "  Single file:"
  echo "    sbatch run_system_diagnosis.sbatch data/hfe_instances/HFE_n32_D64.in [--dim=SECONDS] [--hilbert=SECONDS]"
  echo
  echo "  Directory + array:"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch data/hfe_instances [--dim=20] [--hilbert=60]"
  echo
  echo "  Glob + array (quote the glob!):"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch 'data/hfe_instances/HFE_n*_D*.in' [--dim=20] [--hilbert=60]"
  echo
  echo "  Manifest + array (one path per line):"
  echo "    sbatch --array=0-9 run_system_diagnosis.sbatch filelist.txt [--dim=20] [--hilbert=60]"
  echo
  echo "Notes:"
  echo "  - With a job array, the script picks the Nth file (SLURM_ARRAY_TASK_ID) from:"
  echo "      * all .in files in the directory, OR"
  echo "      * the expanded glob, OR"
  echo "      * the manifest (one path per line)."
  echo "  - Any flags after the file spec are forwarded to the Sage script."
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

INPUT_SPEC="$1"
shift
EXTRA_ARGS=( "$@" )     # e.g. --dim=20 --hilbert=60

SELECTED_FILE=""

# If running as an array job, build the file list and pick the indexed item
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  declare -a FILES
  if [[ -d "$INPUT_SPEC" ]]; then
    # Directory: take all .in files (sorted for stability)
    mapfile -t FILES < <(ls -1 "$INPUT_SPEC"/*.in 2>/dev/null | sort || true)
  elif [[ -f "$INPUT_SPEC" && "$INPUT_SPEC" == *.txt ]]; then
    # Manifest file: one path per line
    mapfile -t FILES < "$INPUT_SPEC"
  else
    # Treat as a glob pattern (quote it at submission time!)
    # shellcheck disable=SC2206
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
  # Non-array: single file path as-is
  SELECTED_FILE="$INPUT_SPEC"
fi

if [[ ! -f "$SELECTED_FILE" ]]; then
  echo "Input file not found: $SELECTED_FILE"
  exit 4
fi

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Running diagnosis for: $SELECTED_FILE"
echo "[$(date)] CPUs: ${SLURM_CPUS_PER_TASK:-1}  Mem: ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-unknown}}  ArrayID: ${SLURM_ARRAY_TASK_ID:-none}"
if ((${#EXTRA_ARGS[@]})); then
  echo "[$(date)] Extra args → ${EXTRA_ARGS[*]}"
fi

# Run the Sage script (prints to stdout and writes logs/*_diagnostics.log)
sage scripts/system_diagnosis.sage "$SELECTED_FILE" "${EXTRA_ARGS[@]}"

echo "[$(date)] ✅ Done"
