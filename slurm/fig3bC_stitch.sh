#!/bin/bash
#SBATCH --job-name=stitchC --cpus-per-task=2 --mem=16G --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/fig3bC_legend.R
conda activate pyspatial
python ~/projects/crc_visiumhd/fig3bC_stitch.py
