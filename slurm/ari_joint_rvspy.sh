#!/bin/bash
#SBATCH --job-name=ariJ
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/ari_joint_rvspy.py
