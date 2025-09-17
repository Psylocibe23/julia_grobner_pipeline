#SBATCH --job-name=ac_extract
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

# Usage:
#   sbatch run_AC_EXTRACT.sbatch <LEX_basis_file.txt> [extract-flags...]
# Examples:
#   # counts only (fast): no enumeration, just AC counts from the eliminant
#   sbatch run_AC_EXTRACT.sbatch results/ANEMOI_p65537_r5_a3_PCICO_..._LEX.txt --skip-enum
#
#   # enumerate roots in minimal extensions up to degree 4 (default):
#   sbatch run_AC_EXTRACT.sbatch results/ANEMOI_p65537_r5_a3_PCICO_..._LEX.txt --enum-maxdeg 4
#
#   # try elimination-ideal if the univariate isn’t present in the basis:
#   sbatch run_AC_EXTRACT.sbatch results/..._LEX.txt --elim --enum-maxdeg 3

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_AC_EXTRACT.sbatch <LEX_basis_file.txt> [extract-flags...]"
  exit 1
fi

INPUT_FILE="$1"; shift
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: file not found: $INPUT_FILE"
  exit 1
fi


module load sagemath/10.6-singularity

mkdir -p logs results

echo "[$(date)] Starting AC extraction on: $INPUT_FILE"
echo "[$(date)] Extra flags: ${*@Q}"

sage scripts/extract_anemoi_solutions.sage "$INPUT_FILE" "$@"

echo "[$(date)] Done."
