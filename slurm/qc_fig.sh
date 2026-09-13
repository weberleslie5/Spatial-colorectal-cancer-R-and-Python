#!/bin/bash
#SBATCH --job-name=qcfig --cpus-per-task=4 --mem=64G --time=2:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/qc_fig.py
