#!/usr/bin/env bash
# Submit the 13-wgbs SLURM chain from the repository root:
#   13.1 genome prep -> 13.2 per-sample array -> 13.3 summary
# Usage:  bash code/13-wgbs-submit.sh            # full chain
#         bash code/13-wgbs-submit.sh --test     # first sample only, to check params
# Then run code/13.4-wgbs-methylation.Rmd interactively (R container).

set -euo pipefail
cd "$(dirname "$0")/.."
source code/13.0-wgbs-config.sh
mkdir -p "${LOGS}"   # SLURM will not create the --output directory itself

n=$(( $(wc -l < "${SAMPLES}") - 1 ))
array="1-${n}%${ALIGN_CONCURRENCY}"
[[ "${1:-}" == "--test" ]] && array="1"

where=(--account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}")
prep=$(sbatch --parsable "${where[@]}" code/13.1-wgbs-genome-prep.sh)
align=$(sbatch --parsable "${where[@]}" --dependency=afterok:"${prep}" --array="${array}" code/13.2-wgbs-trim-align.sh)
summ=$(sbatch --parsable "${where[@]}" --dependency=afterok:"${align}" code/13.3-wgbs-qc-summary.sh)

echo "Submitting to ${SLURM_ACCOUNT} / ${SLURM_PARTITION}"
echo "13.1 genome prep   ${prep}"
echo "13.2 align array   ${align} (tasks ${array})"
echo "13.3 summary       ${summ}"
