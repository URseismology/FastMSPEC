#!/bin/bash
# Plain driver: one SLURM array task per work unit, one technique's index range at a time (see
# NOTES.md "Per-technique array ranges" -- Mspec costs >10x every other technique, so a single
# uniform resource request is wrong for 3/4 of the work units either way).
#
# Submit once per technique, from the deployed location on bluehive (NOTES.md "Deploying to
# bluehive"), overriding --array/--time/--mem/RESULTS_DIR at submission time, e.g.:
#   sbatch --array=0-379%20   --time=00:15:00 --mem=8G  submit_plain.sbatch single-taper
#   sbatch --array=380-759%20 --time=00:30:00 --mem=16G submit_plain.sbatch FastMspec
#   sbatch --array=760-1139%10 --time=02:00:00 --mem=32G submit_plain.sbatch Mspec
#   sbatch --array=1140-1519%20 --time=00:30:00 --mem=16G submit_plain.sbatch MspecBestK
# (index ranges from manifest.technique_index_ranges(); %N caps concurrent tasks, out of
# consideration for other users on the shared partition -- raise once a first small run confirms
# this behaves well; the values above are provisional starting points, not tuned.)
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
DEPLOY_ROOT="/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"

mkdir -p "$DEPLOY_ROOT/logs" "$DEPLOY_ROOT/$RESULTS_DIR"
cd "$DEPLOY_ROOT/python"  # so `python -m dispcurve_pick_batch...` resolves
python3 -m dispcurve_pick_batch.run_plain \
    "$DEPLOY_ROOT/$MANIFEST_CSV" \
    "$DEPLOY_ROOT/$REF_CURVE" \
    "$DEPLOY_ROOT/$RESULTS_DIR" \
    "$SLURM_ARRAY_TASK_ID"
