#!/bin/bash
# Plain driver: one SLURM array task per work unit, one technique's index range at a time (see
# NOTES.md "Per-technique array ranges" -- Mspec costs >10x every other technique, so a single
# uniform resource request is wrong for 3/4 of the work units either way).
#
# Submit once per technique, from the deployed location on bluehive (NOTES.md "Deploying to
# bluehive"), overriding --array/--time/--mem/RESULTS_DIR at submission time, e.g.:
#   sbatch --array=0-379%20   --time=00:15:00 --mem=8G  submit_plain.sbatch
#   sbatch --array=380-759%20 --time=00:30:00 --mem=16G submit_plain.sbatch
# (index ranges from manifest.technique_index_ranges(); the work unit's technique is derived from
# its global index by run_plain.py itself -- there is no separate technique argument, despite what
# an earlier version of this comment claimed; %N caps concurrent tasks, out of consideration for
# other users on the shared partition.)
#
# This cluster's MaxArraySize is 1001, i.e. a raw --array index must stay < 1001 (found live,
# 2026-09-02, trying to resubmit lost Mspec work units at their real global indices 760-1139:
# "Invalid job array specification"). Mspec (760-1139) and MspecBestK (1140-1519) can exceed that
# on their own -- for those, submit with small array indices (e.g. --array=0-379) and set
# IDX_OFFSET to the technique's range start, e.g.:
#   sbatch --array=0-379%15 --time=03:00:00 --mem=120G --export=ALL,IDX_OFFSET=760 submit_plain.sbatch
#
#SBATCH --job-name=dispcurve_plain
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/plain_%A_%a.out
#SBATCH --error=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/plain_%A_%a.err

set -euo pipefail
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch

MANIFEST_CSV="${MANIFEST_CSV:-data/madagascar_stn_conn_ccflist.csv}"
REF_CURVE="${REF_CURVE:-data/SDISPL.ASC}"
RESULTS_DIR="${RESULTS_DIR:-results}"
IDX_OFFSET="${IDX_OFFSET:-0}"
DEPLOY_ROOT="/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"

mkdir -p "$DEPLOY_ROOT/logs" "$DEPLOY_ROOT/$RESULTS_DIR"
cd "$DEPLOY_ROOT/python"  # so `python -m dispcurve_pick_batch...` resolves
python3 -m dispcurve_pick_batch.run_plain \
    "$DEPLOY_ROOT/$MANIFEST_CSV" \
    "$DEPLOY_ROOT/$REF_CURVE" \
    "$DEPLOY_ROOT/$RESULTS_DIR" \
    "$((SLURM_ARRAY_TASK_ID + IDX_OFFSET))"
