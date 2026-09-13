#!/bin/bash
#SBATCH --job-name=refC
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/ref_custom.R
