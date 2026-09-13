#!/bin/bash
#SBATCH --job-name=pyclust --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate pyspatial
export TMPDIR=~/tmp
python ~/projects/crc_visiumhd/py_cluster_sample.py "$1"
