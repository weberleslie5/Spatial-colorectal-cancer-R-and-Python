#!/bin/bash
#SBATCH --job-name=joint
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/joint_pipeline.R           # cluster, name, render panels
conda activate pyspatial
python ~/projects/crc_visiumhd/joint_stitch.py             # trim + assemble the final row
