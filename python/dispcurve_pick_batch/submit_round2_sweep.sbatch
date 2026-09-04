#!/bin/bash
# Round 2 NW-sweep driver: one SLURM array task per (pair, target-NW) point, 60 pairs x 5
# fractions of each pair's own NW_high(R) = 300 total work units. See run_round2_sweep.py and
# docs/notebook5_revamp_progress.md's "Round 2 design" section for the full rationale.
#
# Submit from the deployed location on bluehive:
#   sbatch --array=0-299%15 --time=02:00:00 --mem=120G submit_round2_sweep.sbatch
# (--time/--mem match the FastMspec production budget already proven safe this project's largest
# pairs; %15 concurrency cap out of consideration for other users on the shared partition, same
# as every other job this stage.)
#
#SBATCH --job-name=dispcurve_r2sweep
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/r2sweep_%A_%a.out
#SBATCH --error=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/r2sweep_%A_%a.err

set -euo pipefail
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch

SUBSET_CSV="${SUBSET_CSV:-data/round2_subset.csv}"
MANIFEST_CSV="${MANIFEST_CSV:-data/madagascar_stn_conn_ccflist.csv}"
REF_CURVE="${REF_CURVE:-data/SDISPL.ASC}"
RESULTS_DIR="${RESULTS_DIR:-results_round2_sweep}"
DEPLOY_ROOT="/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"

mkdir -p "$DEPLOY_ROOT/logs" "$DEPLOY_ROOT/$RESULTS_DIR"
cd "$DEPLOY_ROOT/python"
python3 -m dispcurve_pick_batch.run_round2_sweep \
    "$DEPLOY_ROOT/$SUBSET_CSV" \
    "$DEPLOY_ROOT/$MANIFEST_CSV" \
    "$DEPLOY_ROOT/$REF_CURVE" \
    "$DEPLOY_ROOT/$RESULTS_DIR" \
    "$SLURM_ARRAY_TASK_ID"
