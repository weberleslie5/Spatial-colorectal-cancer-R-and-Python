#!/bin/bash
#SBATCH --job-name=ari
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/r_export_clusters.R "$1"
conda activate pyspatial
python ~/projects/crc_visiumhd/py_ari_sample.py "$1"
