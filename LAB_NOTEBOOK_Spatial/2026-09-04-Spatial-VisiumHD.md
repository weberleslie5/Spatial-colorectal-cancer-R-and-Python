# 2026-09-04 — Spatial Visium HD

RTCD Level 2

```bash
squeue -u weberl
```

```text
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          11228089       cpu     rctd   weberl  R 1-23:09:40      1 isce003
```

With the official version of spacexr, the rctd analysis is spacing out.
Using now the pull request (PR): proposed changed of the spacexr original code PR#206, re-writing the slow code for Visium HD's bin counts.

Suite

Cleaning for space:

```bash
df -h ~
```

```text
Filesystem          Size  Used Avail Use% Mounted on
issweka1-ib/homefs  100G  100G     0 100% /home
```

```bash
du -sh ~/* 2>/dev/null | sort -h | tail -10
```

```text
268M    /home/weberl/LW-2995_TRI-3
282M    /home/weberl/LW-2995_TRI-6
286M    /home/weberl/LW-2995_TRI-5
796M    /home/weberl/Notebooks
2.9G    /home/weberl/LW-2875_2_LRRK2_dogma_multiome
3.1G    /home/weberl/LW-2875_1_WT_dogma_multiome
3.4G    /home/weberl/LW-2995_TRI-2
15G     /home/weberl/LW-2995_TRI-1
23G     /home/weberl/projects
44G     /home/weberl/data
```

items in data:

```bash
du -sh ~/data/crc_visiumhd/* ~/projects/crc_visiumhd/* 2>/dev/null | sort -h | tail -20
```

```text
12K     /home/weberl/projects/crc_visiumhd/rctd_10872206.log
20K     /home/weberl/projects/crc_visiumhd/rctd_10489792.log
20K     /home/weberl/projects/crc_visiumhd/rctd_11228089.log
28K     /home/weberl/projects/crc_visiumhd/rctd_10196760.log
12M     /home/weberl/projects/crc_visiumhd/r_clusters.csv
16M     /home/weberl/projects/crc_visiumhd/figures
177M    /home/weberl/projects/crc_visiumhd/rctd_reference.rds
402M    /home/weberl/data/crc_visiumhd/HumanColonCancer_VisiumHD
406M    /home/weberl/data/crc_visiumhd/chromium
407M    /home/weberl/projects/crc_visiumhd/rctd_reference_L2.rds
450M    /home/weberl/projects/crc_visiumhd/01_qc.rds
756M    /home/weberl/projects/crc_visiumhd/04_rctd.rds
1.4G    /home/weberl/projects/crc_visiumhd/02_clustered.rds
1.4G    /home/weberl/projects/crc_visiumhd/03_annotated.rds
1.7G    /home/weberl/projects/crc_visiumhd/py_01_qc.h5ad
8.2G    /home/weberl/projects/crc_visiumhd/py_02)clustered.h5ad
8.2G    /home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad
9.8G    /home/weberl/data/crc_visiumhd/binned_outputs
16G     /home/weberl/data/crc_visiumhd/P5CRC
18G     /home/weberl/data/crc_visiumhd/P2CRC
```

```bash
ls ~/data/crc_visiumhd/P2CRC ~/data/crc_visiumhd/P5CRC

rm -f "/home/weberl/projects/crc_visiumhd/py_02)clustered.h5ad"
rm -f /home/weberl/projects/crc_visiumhd/py_01_qc.h5ad
rm -f /home/weberl/projects/crc_visiumhd/01_qc.rds
rm -f /home/weberl/data/crc_visiumhd/binned_outputs/square_008um/cloupe.cloupe
rm -rf /home/weberl/data/crc_visiumhd/binned_outputs/square_008um/raw_feature_bc_matrix*
rm -rf /home/weberl/data/crc_visiumhd/binned_outputs/square_008um/analysis
```

```bash
ls /home/weberl/data/crc_visiumhd/P2CRC/binned_outputs/square_008um/
```

```text
analysis  cloupe.cloupe
```

>> P2CRC extraction is incomplete. Delete partial files.

```bash
rm -rf /home/weberl/data/crc_visiumhd/P2CRC/binned_outputs/

df -h ~
```

```text
Filesystem          Size  Used Avail Use% Mounted on
issweka1-ib/homefs  100G   84G   17G  84% /home
```

extracting again:

```bash
cd ~/data/crc_visiumhd/P2CRC
tar -xzf *.tar.gz --wildcards '*square_008um/filtered_feature_bc_matrix.h5' '*square_008um/spatial/*' '*square_002um/spatial/*.png' '*square_002um/spatial/*.jpg' '*square_002um/spatial/*.tiff'

ls binned_outputs/square_008um/
```

```text
filtered_feature_bc_matrix.h5  spatial
```

```bash
ls binned_outputs/square_008um/spatial/
```

```text
aligned_fiducials.jpg     cytassist_image.tiff       scalefactors_json.json  tissue_lowres_image.png
aligned_tissue_image.jpg  detected_tissue_image.jpg  tissue_hires_image.png  tissue_positions.parquet
```

```bash
rm -f *.tar.gz

cd ~/data/crc_visiumhd/P5CRC/
tar -xzf *.tar.gz --wildcards '*square_008um/filtered_feature_bc_matrix.h5' '*square_008um/spatial/*' '*square_002um/spatial/*.png' '*square_002um/spatial/*.jpg' '*square_002um/spatial/*.tiff'

ls binned_outputs/square_008um/
```

```text
filtered_feature_bc_matrix.h5  spatial
```

```bash
ls binned_outputs/square_008um/spatial/
```

```text
aligned_fiducials.jpg     cytassist_image.tiff       scalefactors_json.json  tissue_lowres_image.png
aligned_tissue_image.jpg  detected_tissue_image.jpg  tissue_hires_image.png  tissue_positions.parquet
```

```bash
cd ~/projects/crc_visiumhd
sbatch r_cluster_sample.sh P2CRC
sbatch py_cluster_sample.sh P2CRC
sbatch r_cluster_sample.sh P5CRC
sbatch py_cluster_sample.sh P5CRC

sleep 30; squeue -u weberl
```

```text
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          11480801       cpu  pyclust   weberl  R       1:05      1 isca208
          11480797       cpu   rclust   weberl  R       1:07      1 isca205
          11480798       cpu  pyclust   weberl  R       1:07      1 isca208
          11480799       cpu   rclust   weberl  R       1:07      1 isca208
```

Installing the spacexr fork (PR #206):

```bash
srun --pty --cpus-per-task=4 --mem=32G --time=1:00:00 bash
conda activate rspatial
Rscript -e 'remotes::install_github("dmcable/spacexr#206")'

cd ~/projects/crc_visiumhd
sbatch run_rctd_L2.sh
squeue -u weberl
```

```text
crc_visiumhd]$ squeue -u weberl
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          11481096       cpu     rctd   weberl  R       0:01      1 isca208
          11480961       cpu     bash   weberl  R       3:41      1 isca072
          11480797       cpu   rclust   weberl  R       8:08      1 isca205
          11480799       cpu   rclust   weberl  R       8:08      1 isca208
```

pyclust jobs died.

```bash
ls -t ~/projects/crc_visiumhd/pyclust_*.log | head -2
```

```text
/home/weberl/projects/crc_visiumhd/pyclust_11480798.log
/home/weberl/projects/crc_visiumhd/pyclust_11480801.log
```

```bash
tail -n 20 ~/projects/crc_visiumhd/pyclust_11480798.log
tail -n 20 ~/projects/crc_visiumhd/pyclust_11480801.log
```

```text
/data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages/scanpy/readwrite.py:248: UserWarning: Variable names are not unique
. To make them unique, call `.var_names_make_unique`.
  adata = adata.copy()
/data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/functools.py:982: UserWarning: zero-centering a sparse array/matrix densifies it.
  return dispatch(args[0].__class__)(*args, **kw)
```

>> no internal error. External ending.

```bash
sacct -j 11480798,11480801 --format=JobID,JobName,State,MaxRSS,Elapsed
```

```text
JobID           JobName      State     MaxRSS    Elapsed
------------ ---------- ---------- ---------- ----------
11480798        pyclust  COMPLETED              00:06:38
11480798.ba+      batch  COMPLETED  36051196K   00:06:38
11480798.ex+     extern  COMPLETED              00:06:38
11480801        pyclust  COMPLETED              00:05:38
11480801.ba+      batch  COMPLETED  30468976K   00:05:38
11480801.ex+     extern  COMPLETED              00:05:38
```

```bash
df -h ~
```

```text
Filesystem          Size  Used Avail Use% Mounted on
issweka1-ib/homefs  100G   70G   31G  70% /home
```

```bash
ls -lh ~/projects/crc_visiumhd/py_02_clustered_P*.h5ad 2>/dev/null
```

```text
-rw-rw-r-- 1 weberl weberl 8.2G Sep  4 17:02 /home/weberl/projects/crc_visiumhd/py_02_clustered_P2CRC.h5ad
-rw-rw-r-- 1 weberl weberl 7.0G Sep  4 17:01 /home/weberl/projects/crc_visiumhd/py_02_clustered_P5CRC.h5ad
```

In fact, completed.

P1/P2/P5 deconvolution with the new fork PR#206:

```bash
cat > ~/projects/crc_visiumhd/rctd_sample.R << "EOF"
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
library(Seurat); library(spacexr)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
ref <- readRDS("~/projects/crc_visiumhd/rctd_reference.rds")
coords <- GetTissueCoordinates(obj)[, c("x","y")]
cts <- GetAssayData(obj, layer = "counts")
puck <- SpatialRNA(coords, cts)
rctd <- create.RCTD(puck, ref, max_cores = 8)
rctd <- run.RCTD(rctd, doublet_mode = "doublet")
saveRDS(rctd, paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
EOF
```

```bash
cat > ~/projects/crc_visiumhd/rctd_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=rctdS --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate rspatial
export TMPDIR=~/tmp
Rscript ~/projects/crc_visiumhd/rctd_sample.R "$1"
EOF
```

Renaming P1:

```bash
cp ~/projects/crc_visiumhd/02_clustered.rds ~/projects/crc_visiumhd/02_clustered_P1CRC.rds

cd ~/projects/crc_visiumhd
sbatch rctd_sample.sh P1CRC
sbatch rctd_sample.sh P2CRC
sbatch rctd_sample.sh P5CRC
sleep 30; squeue -u weberl
```
