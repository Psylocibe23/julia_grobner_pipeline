#!/bin/bash
###############################################################################
# run_grobner_pipeline_F4.sbatch
#
# End-to-end pipeline:
#   1) F4 (Julia)              -> results/<stem>_F4_<timestamp>.txt
#   2) FGLM DRL->LEX (Sage)    -> results/<stem>_F4_<timestamp>_LEX.txt
#   3) Extract solutions (Sage)-> results/<stem>_F4_<timestamp>_LEX_sols.txt
#
# INPUT:
#   sbatch run_grobner_pipeline_F4.sbatch <system>
#   where <system> is EITHER:
#     (a) a stem like "HFE_n6_D8"  -> we use data/<stem>.in
#     (b) a full path to a .in file (e.g., data/HFE_n6_D8.in)
#
# NOTES:
#   - Adjust module names to your cluster environment if needed.
#   - We extract metadata directly from the LEX basis and solutions files
#     (no dependency on a separate "solutions log" file).
###############################################################################

#SBATCH --job-name=grobner_pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --partition=turing-long
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

# -------------------- args & defaults --------------------
SYSTEM_ARG="${1:-HFE_n4_d4_system}"
THREADS="${2:-8}"

mkdir -p logs results

START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
SECONDS=0

echo "[$START_TIME] Starting Gröbner pipeline. Arg: ${SYSTEM_ARG}"

# If arg looks like a .in file and exists, use it; otherwise treat as stem.
if [[ "${SYSTEM_ARG}" == *.in && -f "${SYSTEM_ARG}" ]]; then
  INPUT_FILE="${SYSTEM_ARG}"
  STEM="$(basename "${INPUT_FILE%.*}")"
else
  INPUT_FILE="data/${SYSTEM_ARG}.in"
  STEM="${SYSTEM_ARG}"
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: Input file not found: ${INPUT_FILE}"
  exit 1
fi

echo "Using input: ${INPUT_FILE}"
echo "Threads for F4: ${THREADS}"

# -------------------- STEP 1: F4 (Julia) --------------------
module load julia

echo "[$(date)] Running F4 (Julia)..."
/usr/bin/time -v julia scripts/solve_F4_from_file.jl "${INPUT_FILE}" "${THREADS}" \
  2> "logs/${STEM}_f4_time_mem_${SLURM_JOB_ID}.log"

# Pick newest F4 result matching the input stem
BASIS_FILE=$(ls -1t results/${STEM}_F4_*.txt 2>/dev/null | head -n1 || true)
if [[ -z "${BASIS_FILE}" || ! -f "${BASIS_FILE}" ]]; then
  echo "ERROR: F4 result file not found for stem '${STEM}'."
  exit 2
fi
echo "[$(date)] F4 result: ${BASIS_FILE}"

# -------------------- STEP 2: FGLM (Sage) --------------------
module load sagemath/10.6-singularity

echo "[$(date)] Running FGLM DRL->LEX (Sage)..."
/usr/bin/time -v sage scripts/convert_to_lex_fglm.sage "${BASIS_FILE}" \
  2> "logs/${STEM}_fglm_time_mem_${SLURM_JOB_ID}.log"

LEX_FILE=$(ls -1t results/${STEM}_F4_*_LEX.txt 2>/dev/null | head -n1 || true)
if [[ -z "${LEX_FILE}" || ! -f "${LEX_FILE}" ]]; then
  echo "ERROR: LEX basis file not found for stem '${STEM}'."
  exit 3
fi
echo "[$(date)] LEX result: ${LEX_FILE}"

# -------------------- STEP 3: Extract solutions (Sage) --------------------
echo "[$(date)] Extracting solutions (Sage)..."
/usr/bin/time -v sage scripts/extract_solutions_from_lex.sage "${LEX_FILE}" \
  2> "logs/${STEM}_sol_time_mem_${SLURM_JOB_ID}.log"

SOL_FILE=$(ls -1t results/${STEM}_F4_*_LEX_sols.txt 2>/dev/null | head -n1 || true)
if [[ -z "${SOL_FILE}" || ! -f "${SOL_FILE}" ]]; then
  echo "ERROR: Solutions file not found for stem '${STEM}'."
  exit 4
fi
echo "[$(date)] Solutions: ${SOL_FILE}"

# -------------------- METRICS: parse from files --------------------
# Variables from LEX header
VARNAMES_LINE=$(grep -m1 "^# Variables:" "${LEX_FILE}" || true)
VARS=$(echo "${VARNAMES_LINE#*:}" | tr -d ' ')
NUM_VARS=$(echo "${VARS}" | tr ',' '\n' | sed '/^$/d' | wc -l)

# Field from LEX header (e.g., "# Field: GF(2)")
FIELD_LINE=$(grep -m1 "^# Field:" "${LEX_FILE}" || true)
FIELD=$(echo "${FIELD_LINE#*:}" | tr -d '[:space:]' | sed 's/^GF(//; s/)$//')

# Equations: equal to number of variables in your .in format
NUM_EQS="${NUM_VARS}"

# Num solutions = count lines starting with "{"
NUM_SOL=$(grep -c "^{.*}$" "${SOL_FILE}" || echo 0)

# -------------------- TIME & MEMORY --------------------
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
WALLTIME=$(printf '%02d:%02d:%02d' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))

F4_TIME=$(grep "Elapsed (wall clock) time" "logs/${STEM}_f4_time_mem_${SLURM_JOB_ID}.log" | awk '{print $8}')
F4_RAM=$(grep "Maximum resident set size" "logs/${STEM}_f4_time_mem_${SLURM_JOB_ID}.log" | awk '{print $6}')
FGLM_TIME=$(grep "Elapsed (wall clock) time" "logs/${STEM}_fglm_time_mem_${SLURM_JOB_ID}.log" | awk '{print $8}')
FGLM_RAM=$(grep "Maximum resident set size" "logs/${STEM}_fglm_time_mem_${SLURM_JOB_ID}.log" | awk '{print $6}')
SOL_TIME=$(grep "Elapsed (wall clock) time" "logs/${STEM}_sol_time_mem_${SLURM_JOB_ID}.log" | awk '{print $8}')
SOL_RAM=$(grep "Maximum resident set size" "logs/${STEM}_sol_time_mem_${SLURM_JOB_ID}.log" | awk '{print $6}')

# -------------------- SOFTWARE VERSIONS --------------------
JULIA_VERSION=$(julia -e 'println(VERSION)' 2>/dev/null || echo "unknown")
ALGEBRAIC_VERSION=$(julia -e 'using Pkg; d=Pkg.dependencies(); for (k,v) in d; if string(k)=="AlgebraicSolving"; println(v.version); end; end' 2>/dev/null || echo "unknown")
SAGE_VERSION=$(sage --version 2>/dev/null | head -n1 || echo "unknown")

# -------------------- SUMMARY LOG --------------------
SUMMARY_LOG=logs/experiment_summary.log
{
  echo "# Experiment summary"
  echo "System stem: ${STEM}"
  echo "Input file: ${INPUT_FILE}"
  echo "Variables: ${NUM_VARS}"
  echo "Equations: ${NUM_EQS}"
  echo "Field: GF(${FIELD})"
  echo "Threads: ${THREADS}"
  echo "F4 result: ${BASIS_FILE}"
  echo "LEX file: ${LEX_FILE}"
  echo "Solutions file: ${SOL_FILE}"
  echo "Start time: ${START_TIME}"
  echo "End time: ${END_TIME}"
  echo "Total wall time: ${WALLTIME}"
  echo "F4 wall time: ${F4_TIME}"
  echo "F4 max RAM: ${F4_RAM} kB"
  echo "FGLM wall time: ${FGLM_TIME}"
  echo "FGLM max RAM: ${FGLM_RAM} kB"
  echo "Solution extraction time: ${SOL_TIME}"
  echo "Solution extraction max RAM: ${SOL_RAM} kB"
  echo "Number of solutions: ${NUM_SOL}"
  echo "Software versions: Julia ${JULIA_VERSION}, AlgebraicSolving.jl ${ALGEBRAIC_VERSION}, ${SAGE_VERSION}"
  echo
} >> "${SUMMARY_LOG}"

echo "[$(date)] Pipeline complete for ${STEM}."
