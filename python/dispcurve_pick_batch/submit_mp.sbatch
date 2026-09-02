#!/bin/bash
# Multiprocessing driver: one SLURM array task per NODE, each fanning a SLICE of work units out
# across a multiprocessing.Pool sized to that node's --cpus-per-task. See NOTES.md "Plain vs.
# multiprocessing" and "Per-technique array ranges" -- same technique-range logic as
# submit_plain.sbatch, just chunked into node-sized slices instead of one task per work unit.
#
# Submit once per technique with START/END/CHUNK set to that technique's own index range from
# manifest.technique_index_ranges() and a chunk size (work units per array task -- pick so
# n_chunks isn't tiny for a technique with few work units), e.g. for a chunk size of 20:
#   sbatch --array=0-18   --time=00:20:00 --mem=32G  --cpus-per-task=20 submit_mp.sbatch 0 380 20
#   sbatch --array=0-18   --time=00:40:00 --mem=48G  --cpus-per-task=20 submit_mp.sbatch 380 760 20
#   sbatch --array=0-18   --time=02:30:00 --mem=64G  --cpus-per-task=20 submit_mp.sbatch 760 1140 20
#   sbatch --array=0-18   --time=00:40:00 --mem=48G  --cpus-per-task=20 submit_mp.sbatch 1140 1520 20
# (--mem scales with --cpus-per-task since each worker process can independently peak near the
# single-work-unit memory high-water marks measured in Stage 3, e.g. ~18-20GB for Mspec -- do not
# assume workers share that footprint efficiently; size for the worst case x worker count, then
# adjust down once a real run's actual peak is observed.)
#
#SBATCH --job-name=dispcurve_mp
#SBATCH --partition=preempt
#SBATCH --output=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/mp_%A_%a.out
#SBATCH --error=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/mp_%A_%a.err

set -euo pipefail
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch

START="$1"
END="$2"
CHUNK="$3"

TASK_START=$((START + SLURM_ARRAY_TASK_ID * CHUNK))
TASK_END=$((TASK_START + CHUNK))
if [ "$TASK_END" -gt "$END" ]; then TASK_END=$END; fi

MANIFEST_CSV="${MANIFEST_CSV:-data/madagascar_stn_conn_ccflist.csv}"
REF_CURVE="${REF_CURVE:-data/SDISPL.ASC}"
RESULTS_DIR="${RESULTS_DIR:-results}"
DEPLOY_ROOT="/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"

mkdir -p "$DEPLOY_ROOT/logs" "$DEPLOY_ROOT/$RESULTS_DIR"
cd "$DEPLOY_ROOT/python"
python3 -m dispcurve_pick_batch.run_multiprocessing \
    "$DEPLOY_ROOT/$MANIFEST_CSV" \
    "$DEPLOY_ROOT/$REF_CURVE" \
    "$DEPLOY_ROOT/$RESULTS_DIR" \
    "$TASK_START" "$TASK_END" --workers "${SLURM_CPUS_PER_TASK:-4}"
