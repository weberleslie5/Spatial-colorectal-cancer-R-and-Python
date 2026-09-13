#!/bin/bash
#SBATCH --job-name=pyjoint
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/py_joint.py
