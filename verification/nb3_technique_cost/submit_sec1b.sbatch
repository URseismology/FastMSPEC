#!/bin/bash
#SBATCH --job-name=nb3_sec1b
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=120G
#SBATCH --output=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/nb3_sec1b_%j.out
#SBATCH --error=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/nb3_sec1b_%j.err

set -euo pipefail
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch
python3 /scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/nb3_precompute/driver_sec1b.py
