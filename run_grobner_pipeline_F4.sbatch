#!/bin/bash
#SBATCH --job-name=grobner_pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --partition=turing-long
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

SYSTEM_NAME=${1:-HFE_n4_d4_system}  # Set to input system name
THREADS=8

mkdir -p logs

# Start time
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
SECONDS=0

echo "[$START_TIME] Starting Groebner pipeline for system: $SYSTEM_NAME"

# 1. Step 1: Compute Groebner basis with F4 (Julia)
module load julia

echo "[$(date)] Running F4 computation (Julia)..."
/usr/bin/time -v julia scripts/solve_F4_from_file.jl data/${SYSTEM_NAME}.in $THREADS \
    2> logs/f4_time_mem.log

BASIS_FILE=$(ls -1t results/${SYSTEM_NAME}_F4_*.txt 2>/dev/null | head -n1)
if [[ ! -f "$BASIS_FILE" ]]; then
    echo "ERROR: Could not find F4 result file for system '$SYSTEM_NAME'."
    exit 1
fi
echo "[$(date)] F4 result file: $BASIS_FILE"

# 2. Step 2: Convert to LEX using FGLM (SageMath)
module load sagemath/10.6-singularity

echo "[$(date)] Running FGLM conversion (SageMath)..."
/usr/bin/time -v sage -python scripts/convert_to_lex_fglm.sage.py "$BASIS_FILE" \
    2> logs/fglm_time_mem.log


LEX_FILE=$(ls -1t results/${SYSTEM_NAME}_F4_*_LEX.txt 2>/dev/null | head -n1)
if [[ ! -f "$LEX_FILE" ]]; then
    echo "ERROR: Could not find LEX Groebner basis file for system '$SYSTEM_NAME'."
    exit 2
fi
echo "[$(date)] LEX result file: $LEX_FILE"

# 3. Step 3: Extract solutions from LEX basis (SageMath)
echo "[$(date)] Extracting solutions from LEX basis (SageMath)..."
/usr/bin/time -v sage scripts/extract_solutions_from_lex.sage "$LEX_FILE" \
    2> logs/sol_time_mem.log

SOL_FILE=$(ls -1t results/${SYSTEM_NAME}_F4_*_LEX_sols.txt 2>/dev/null | head -n1)
if [[ ! -f "$SOL_FILE" ]]; then
    echo "ERROR: Could not find solution file for system '$SYSTEM_NAME'."
    exit 3
fi
echo "[$(date)] Solution file: $SOL_FILE"

# Locate the most recent solution log file for parsing
SOL_LOG_FILE=$(ls -1t logs/${SYSTEM_NAME}_F4_*_LEX_SOLUTIONS.log 2>/dev/null | head -n1)

# Parse details from the solution log
NUM_VARS=$(grep "# Variables:" "$SOL_LOG_FILE" | awk -F: '{print $2}' | tr -d " []',")
NUM_EQS=$(grep "# Number of equations:" "$SOL_LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')
FIELD=$(grep "# Field characteristic:" "$SOL_LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')
NUM_SOL=$(grep "# Number of solutions:" "$SOL_LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')
VERIF=$(grep "# Verification passed:" "$SOL_LOG_FILE" | awk -F: '{print $2}' | tr -d ' ')

# End time and wall time
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
WALLTIME=$(printf '%02d:%02d:%02d' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))

# Extract timing and memory data
F4_TIME=$(grep "Elapsed (wall clock) time" logs/f4_time_mem.log | awk '{print $8}')
F4_RAM=$(grep "Maximum resident set size" logs/f4_time_mem.log | awk '{print $6}')
FGLM_TIME=$(grep "Elapsed (wall clock) time" logs/fglm_time_mem.log | awk '{print $8}')
FGLM_RAM=$(grep "Maximum resident set size" logs/fglm_time_mem.log | awk '{print $6}')
SOL_TIME=$(grep "Elapsed (wall clock) time" logs/sol_time_mem.log | awk '{print $8}')
SOL_RAM=$(grep "Maximum resident set size" logs/sol_time_mem.log | awk '{print $6}')

# Software versions
JULIA_VERSION=$(julia -e 'println(VERSION)')
ALGEBRAIC_VERSION=$(julia -e 'using Pkg; println(Pkg.dependencies()["AlgebraicSolving"].version)' 2>/dev/null)
SAGE_VERSION=$(sage --version | head -n1)

# Write experiment summary log (append mode)
SUMMARY_LOG=logs/experiment_summary.log
{
echo "# Experiment summary"
echo "System: $SYSTEM_NAME"
echo "Variables: $NUM_VARS"
echo "Equations: $NUM_EQS"
echo "Field: GF($FIELD)"
echo "Threads: $THREADS"
echo "Input file: data/${SYSTEM_NAME}.in"
echo "Start time: $START_TIME"
echo "End time: $END_TIME"
echo "Total wall time: $WALLTIME"
echo "F4 wall time: $F4_TIME"
echo "F4 max RAM: $F4_RAM kB"
echo "FGLM wall time: $FGLM_TIME"
echo "FGLM max RAM: $FGLM_RAM kB"
echo "Solution extraction time: $SOL_TIME"
echo "Solution extraction max RAM: $SOL_RAM kB"
echo "Number of solutions: $NUM_SOL"
echo "Verification: $VERIF"
echo "Software versions: Julia $JULIA_VERSION, AlgebraicSolving.jl $ALGEBRAIC_VERSION, $SAGE_VERSION"
echo
} >> $SUMMARY_LOG

echo "[$(date)] Groebner pipeline complete for system: $SYSTEM_NAME"
echo "Logs and results are available in the current directory and logs/"
