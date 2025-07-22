#!/bin/bash
#SBATCH --job-name=grobner_pipeline_F5
#SBATCH --output=logs/pipeline_F5_%j.log
#SBATCH --partition=turing-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

SYSTEM_NAME="toy_test_prime_field"  # Set to input system name
THREADS=4

mkdir -p logs

# Start time
START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
SECONDS=0

echo "[$START_TIME] Starting Groebner F5 pipeline for system: $SYSTEM_NAME"

# 1. Step 1: Compute Groebner basis with F5 (Julia)
module load julia

echo "[$(date)] Running F5 computation (Julia)..."
/usr/bin/time -v julia scripts/solve_F5_from_file.jl data/${SYSTEM_NAME}.in $THREADS \
    2> logs/f5_time_mem.log

BASIS_FILE=$(ls -1t results/${SYSTEM_NAME}_F5_*.txt 2>/dev/null | head -n1)
if [[ ! -f "$BASIS_FILE" ]]; then
    echo "ERROR: Could not find F5 result file for system '$SYSTEM_NAME'."
    exit 1
fi
echo "[$(date)] F5 result file: $BASIS_FILE"

# 2. Step 2: Convert to LEX using FGLM (SageMath)
module load sagemath/10.6-singularity

echo "[$(date)] Running FGLM conversion (SageMath)..."
/usr/bin/time -v sage scripts/convert_to_lex_fglm.sage "$BASIS_FILE" \
    2> logs/fglm_time_mem.log

LEX_FILE=$(ls -1t results/${SYSTEM_NAME}_F5_*_LEX.txt 2>/dev/null | head -n1)
if [[ ! -f "$LEX_FILE" ]]; then
    echo "ERROR: Could not find LEX Groebner basis file for system '$SYSTEM_NAME'."
    exit 2
fi
echo "[$(date)] LEX result file: $LEX_FILE"

# 3. Step 3: Extract solutions from LEX basis (SageMath)
echo "[$(date)] Extracting solutions from LEX basis (SageMath)..."
/usr/bin/time -v sage scripts/extract_solutions_from_lex.sage "$LEX_FILE" \
    2> logs/sol_time_mem.log

SOL_FILE=$(ls -1t results/${SYSTEM_NAME}_F5_*_LEX_sols.txt 2>/dev/null | head -n1)
if [[ ! -f "$SOL_FILE" ]]; then
    echo "ERROR: Could not find solution file for system '$SYSTEM_NAME'."
    exit 3
fi
echo "[$(date)] Solution file: $SOL_FILE"

# Locate the most recent solution log file for parsing
SOL_LOG_FILE=$(ls -1t logs/${SYSTEM_NAME}_F5_*_LEX_SOLUTIONS.log 2>/dev/null | head -n1)

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
F5_TIME=$(grep "Elapsed (wall clock) time" logs/f5_time_mem.log | awk '{print $8}')
F5_RAM=$(grep "Maximum resident set size" logs/f5_time_mem.log | awk '{print $6}')
FGLM_TIME=$(grep "Elapsed (wall clock) time" logs/fglm_time_mem.log | awk '{print $8}')
FGLM_RAM=$(grep "Maximum resident set size" logs/fglm_time_mem.log | awk '{print $6}')
SOL_TIME=$(grep "Elapsed (wall clock) time" logs/sol_time_mem.log | awk '{print $8}')
SOL_RAM=$(grep "Maximum resident set size" logs/sol_time_mem.log | awk '{print $6}')

# Software versions
JULIA_VERSION=$(julia -e 'println(VERSION)')
ALGEBRAIC_VERSION=$(julia -e 'using Pkg; println(Pkg.dependencies()["AlgebraicSolving"].version)' 2>/dev/null)
SAGE_VERSION=$(sage --version | head -n1)

# Write experiment summary log (append mode)
SUMMARY_LOG=logs/experiment_summary_F5.log
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
echo "F5 wall time: $F5_TIME"
echo "F5 max RAM: $F5_RAM kB"
echo "FGLM wall time: $FGLM_TIME"
echo "FGLM max RAM: $FGLM_RAM kB"
echo "Solution extraction time: $SOL_TIME"
echo "Solution extraction max RAM: $SOL_RAM kB"
echo "Number of solutions: $NUM_SOL"
echo "Verification: $VERIF"
echo "Software versions: Julia $JULIA_VERSION, AlgebraicSolving.jl $ALGEBRAIC_VERSION, $SAGE_VERSION"
echo
} >> $SUMMARY_LOG

echo "[$(date)] Groebner pipeline F5 complete for system: $SYSTEM_NAME"
echo "Logs and results are available in the current directory and logs/"
