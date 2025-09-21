#!/bin/bash
#SBATCH --job-name=gb_F4_long
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=980G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4_long.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
THREADS="${2:-16}"   # fewer F4 workers to lower peak RAM

[[ -f "$INPUT_FILE" ]] || { echo "Input not found: $INPUT_FILE"; exit 2; }

# Avoid oversubscription from math libs (F4 uses its own threads param)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Scratch for temporaries
export TMPDIR=${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}
mkdir -p "$TMPDIR"

# Set F4 tuning knobs (can be overridden at submit time: F4_MAX_PAIRS=4000 sbatch ...)
export F4_THREADS="$THREADS"
export F4_MAX_PAIRS="${F4_MAX_PAIRS:-6000}"
export F4_INITIAL_HTS="${F4_INITIAL_HTS:-19}"
export F4_LA_OPTION="${F4_LA_OPTION:-22}"

# Large stack to avoid recursion crashes in parser or F4
ulimit -s unlimited || true

# Make 'module' available, then load Julia
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load julia

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] F4 threads (nr_thrds): $F4_THREADS"
echo "[$(date)] F4 knobs: la_option=$F4_LA_OPTION initial_hts=$F4_INITIAL_HTS max_nr_pairs=$F4_MAX_PAIRS"

# Use repo as the project; keep depot in HOME by default
export JULIA_PROJECT=$PWD
export JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH:-$HOME/.julia}

# Ensure deps exist
julia --project=$JULIA_PROJECT -e '
using Pkg
Pkg.activate(".")
for pkg in ["AlgebraicSolving","AbstractAlgebra","Nemo"]
    isnothing(Base.find_package(pkg)) && Pkg.add(pkg)
end
Pkg.instantiate(); Pkg.precompile()
Pkg.status() |> println
'

# Run
julia --project=$JULIA_PROJECT scripts/solve_F4_from_file_parsing.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
