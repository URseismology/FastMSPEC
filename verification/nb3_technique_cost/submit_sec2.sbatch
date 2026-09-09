#!/bin/bash
#SBATCH --job-name=nb3_sec2
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --mem=96G
#SBATCH --output=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/nb3_sec2_%j.out
#SBATCH --error=/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/logs/nb3_sec2_%j.err

set -euo pipefail
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch
python3 /scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/nb3_precompute/driver_sec2.py
