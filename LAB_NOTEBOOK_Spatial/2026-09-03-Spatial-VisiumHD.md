# 2026-09-03 — Spatial Visium HD

Checking the progress of RCTD Level 2 job:

```bash
squeue -u weberl
```

## Reproducing Figure 3 from Oliveira et al 2025 Nature Genetics

Downloading samples P2 and P5:

```bash
mkdir -p ~/data/crc_visiumhd/P2CRC ~/data/crc_visiumhd/P5CRC

cd ~/data/crc_visiumhd/P2CRC
nohup wget -c "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P2/Visium_HD_Human_Colon_Cancer_P2_binned_outputs.tar.gz" > wget_P2.log 2>&1 &

cd ~/data/crc_visiumhd/P5CRC
nohup wget -c "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P5/Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz" > wget_P5.log 2>&1 &

jobs
```

## Running samples P2 and P5

Preparing the R and Python scripts for Loading > QC clustering > Cluster map > Saving, before clusters annotation and RCTD (R only) deconvolution:

1. Clustering with R

```bash
cat > ~/projects/crc_visiumhd/r_cluster_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE)    # read sample name from the command line: one script for all samples
sample <- args[1]
library(Seurat); library(ggplot2)
obj <- Load10X_Spatial(paste0("~/data/crc_visiumhd/", sample), bin.size = 8)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, nCount_Spatial.008um >= 100 & nFeature_Spatial.008um >= 50 & percent.mt < 30)
obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj); obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, dims = 1:30); obj <- FindClusters(obj, resolution = 0.6)
df <- GetTissueCoordinates(obj); df$cluster <- obj$seurat_clusters
p <- ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.1) + coord_fixed() + theme_void() + guides(colour = guide_legend(override.aes = list(size = 4)))
ggsave(paste0("~/projects/crc_visiumhd/figures/02_clusters_", sample, ".png"), p, width = 8, height = 6)
saveRDS(obj, paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
EOF
```

R job ticket:

```bash
cat > ~/projects/crc_visiumhd/r_cluster_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=rclust --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate rspatial
export TMPDIR=~/tmp
Rscript ~/projects/crc_visiumhd/r_cluster_sample.R "$1"
EOF
```

2. Clustering with Python

```bash
cat > ~/projects/crc_visiumhd/py_cluster_sample.py << 'EOF'
import sys, scanpy as sc, pandas as pd
sample = sys.argv[1]    # read sample name from the command line: one script for all samples
base = f"/home/weberl/data/crc_visiumhd/{sample}/binned_outputs/square_008um"
adata = sc.read_10x_h5(f"{base}/filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
coords = pd.read_parquet(f"{base}/spatial/tissue_positions.parquet").set_index("barcode").loc[adata.obs_names]
adata.obsm["spatial"] = coords[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
adata = adata[(adata.obs.total_counts >= 100) & (adata.obs.n_genes_by_counts >= 50) & (adata.obs.pct_counts_mt < 30)].copy()
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata; adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10); sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_low", flavor="igraph", n_iterations=2)
adata.write(f"/home/weberl/projects/crc_visiumhd/py_02_clustered_{sample}.h5ad")
EOF
```

Python job ticket:

```bash
cat > ~/projects/crc_visiumhd/py_cluster_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pyclust --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate pyspatial
export TMPDIR=~/tmp
python ~/projects/crc_visiumhd/py_cluster_sample.py "$1"
EOF
```

P2

```bash
tail -2 ~/data/crc_visiumhd/P2CRC/wget_P2.log
cd ~/data/crc_visiumhd/P2CRC
tar -xzf *.tar.gz --wildcards '*square_008um*'
rm *.tar.gz
ls binned_outputs/square_008um/
cd ~/projects/crc_visiumhd
sbatch r_cluster_sample.sh P2CRC
sbatch py_cluster_sample.sh P2CRC
```

P5

```bash
tail -2 ~/data/crc_visiumhd/P5CRC/wget_P5.log
cd ~/data/crc_visiumhd/P5CRC
tar -xzf *.tar.gz --wildcards '*square_008um*'
rm *.tar.gz
ls binned_outputs/square_008um/
cd ~/projects/crc_visiumhd
sbatch r_cluster_sample.sh P5CRC
sbatch py_cluster_sample.sh P5CRC
```
