#!/bin/bash
#SBATCH --job-name=rclust --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate rspatial
export TMPDIR=~/tmp
Rscript ~/projects/crc_visiumhd/r_cluster_sample.R "$1"
