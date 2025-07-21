#!/bin/bash
#SBATCH --job-name=grobner_pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --partition=turing-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# Set this to the name of input (without 'data/' prefix)
SYSTEM_NAME="toy_test_prime_field"  # Change this to input system name
THREADS=4

# Ensure a clean logs directory exists
mkdir -p logs

echo "[$(date)] Starting Groebner pipeline for system: $SYSTEM_NAME"

# 1. Step 1: Compute Groebner basis with F4 (Julia)
module load julia

echo "[$(date)] Running F4 computation (Julia)..."
julia scripts/solve_F4_from_file.jl data/${SYSTEM_NAME}.in $THREADS

# Find the latest result file matching your naming pattern
BASIS_FILE=$(ls -1t result_*_${SYSTEM_NAME}_*.out 2>/dev/null | head -n1)
if [[ ! -f "$BASIS_FILE" ]]; then
    echo "ERROR: Could not find F4 result file for system '$SYSTEM_NAME'."
    exit 1
fi
echo "[$(date)] F4 result file: $BASIS_FILE"

# 2. Step 2: Convert to LEX using FGLM (SageMath)
module load sagemath/10.6-singularity

echo "[$(date)] Running FGLM conversion (SageMath)..."
sage scripts/convert_to_lex_fglm.sage "$BASIS_FILE"

# Find the latest LEX result file
LEX_FILE=$(ls -1t lex_*_${SYSTEM_NAME}_*.out 2>/dev/null | head -n1)
if [[ ! -f "$LEX_FILE" ]]; then
    echo "ERROR: Could not find LEX Groebner basis file for system '$SYSTEM_NAME'."
    exit 2
fi
echo "[$(date)] LEX result file: $LEX_FILE"

# 3. Step 3: Extract solutions from LEX basis (SageMath)
echo "[$(date)] Extracting solutions from LEX basis (SageMath)..."
sage scripts/extract_solutions_from_lex.sage "$LEX_FILE"

# Find the latest solution file
SOL_FILE=$(ls -1t solutions_*_${SYSTEM_NAME}_*.txt 2>/dev/null | head -n1)
if [[ ! -f "$SOL_FILE" ]]; then
    echo "ERROR: Could not find solution file for system '$SYSTEM_NAME'."
    exit 3
fi
echo "[$(date)] Solution file: $SOL_FILE"

echo "[$(date)] Groebner pipeline complete for system: $SYSTEM_NAME"
echo "Logs and results are available in the current directory and logs/"
