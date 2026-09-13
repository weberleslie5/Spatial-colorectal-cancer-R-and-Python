# 2026-08-31 — Spatial Visium HD

We are checking on the finer (Level2) deconvolution started on Friday 8/28.

```bash
squeue -u weberl
```

```text
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

```bash
ls -lh ~/projects/crc_visiumhd/04_rctd_L2.rds
```

```text
ls: cannot access '/home/weberl/projects/crc_visiumhd/04_rctd_L2.rds': No such file or directory
```

```bash
ls -lt ~/projects/crc_visiumhd/*.log
```

```text
-rw-rw-r-- 1 weberl weberl 16645 Aug 28 19:01 /home/weberl/projects/crc_visiumhd/rctd_10489792.log
-rw-rw-r-- 1 weberl weberl 25764 Aug 28 00:13 /home/weberl/projects/crc_visiumhd/rctd_10196760.log
```

```bash
tail -30 ~/projects/crc_visiumhd/rctd_10489792.log
```

OOM Killed: Out-of-memory: 30 cell-types models in each of the 8 parallel workers (parallel copies of R that RCTD launches to process bin simultaneously: each hold a copy of the reference data in memory). >> more memory and fewer workers.

```bash
sed -i 's/--mem=128G/--mem=256G/' ~/projects/crc_visiumhd/run_rctd_L2.sh
sed -i 's/max_cores = 8/max_cores = 4/' ~/projects/crc_visiumhd/run_rctd_L2.R

cd ~/projects/crc_visiumhd
sbatch run_rctd_L2.sh
squeue -u weberl
```

```text
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          10872206       cpu     rctd   weberl  R       2:36      1 isca209
```

## In the meantime, we are starting Chapter 2: the analysis done in Python

```bash
srun --pty --cpus-per-task=4 --mem=32G --time=4:00:00 bash
```

Python environment:

```bash
mamba create -n pyspatial -c conda-forge -y scanpy leidenalg python-igraph pyarrow h5py pandas matplotlib
conda activate pyspatial
```

1. Loading and assembling: matrix (sc.read_10x_h5) and parquet coordinates assembled into one AnnData object.

```bash
python
```

COUNTS:

```python
import scanpy as sc #analysis tookit
import pandas as pd #tables

adata = sc.read_10x_h5("/home/weberl/data/crc_visiumhd/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
adata.var_names_make_unique() # de-duplicate gene names (GENE, GENE-1). In R, Load10X_Spatial did it automatically.
adata
```

```text
AnnData object with n_obs × n_vars = 545913 × 18085 # AnnData bins x genes is Seurat (R) convention is genes x bins.
    var: 'gene_ids', 'feature_types', 'genome'
    layers: None (.X)
```

COORDINATES:

```python
coords = pd.read_parquet("/home/weberl/data/crc_visiumhd/binned_outputs/square_008um/spatial/tissue_positions.parquet") # bins positions
coords = coords.set_index("barcode").loc[adata.obs_names] # match coordinate rows to matrix bins by barcode
adata.obsm["spatial"] = coords[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy() # save bin x,y into anndata
adata
```

```text
AnnData object with n_obs × n_vars = 545913 × 18085
    var: 'gene_ids', 'feature_types', 'genome'
    obsm: 'spatial'
    layers: None (.X)
```

2. QC and filter (valid bins if >= 100 counts, >= 50 genes, < 30% mito)

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False) # QC columns computed per bin
adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe() # summary stats
```

```text
total_counts  n_genes_by_counts  pct_counts_mt
count  545913.000000      545913.000000  545850.000000
mean      488.434357         393.164671       4.837324
std       427.640228         332.244261       4.479874
min         0.000000           0.000000       0.000000
25%       124.000000         110.000000       1.941748
50%       380.000000         299.000000       3.921569
75%       770.000000         635.000000       6.340289
max      4145.000000        2594.000000     100.000000
```

50%       380.000000         299.000000       3.921569 >> median matches R's output

Filtering:

```python
n0 = adata.n_obs # bins before
adata = adata[(adata.obs.total_counts >= 100) & (adata.obs.n_genes_by_counts >= 50) & (adata.obs.pct_counts_mt < 30)].copy()
print(n0, "->", adata.n_obs)
```

```text
545913 -> 428405 # Identical to R.
```

```python
adata.write ("/home/weberl/projects/crc_visiumhd/py_01_qc.h5ad")
```

3. Normalize and cluster: log-normalization, PCA, Leiden (Python's standard, equivalent to Louvain)

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata # keep full data for markers
adata = adata[:, adata.var.highly_variable].copy() # only highly variable genes
sc.pp.scale(adata, max_value=10) # center/scale genes
sc.tl.pca(adata, n_comps=30) # compress to 30 dimemsions
```

Building the neighbor graph and Leiden clustering:
(SENT AS A SCRIPT, STARTED AT THE LAST CHECKPOINT: BEFORE NORMALIZATION)

```bash
cat > ~/projects/crc_visiumhd/py_cluster.py << 'EOF'
import scanpy as sc
adata = sc.read_h5ad("/home/weberl/projects/crc_visiumhd/py_01_qc.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata # keep full data for markers
adata = adata[:, adata.var.highly_variable].copy() # only highly variable genes
sc.pp.scale(adata, max_value=10) # center/scale genes
sc.tl.pca(adata, n_comps=30) # compress to 30 dimemsions
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2)
adata.write("/home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad")
EOF
```

```bash
cat > ~/projects/crc_visiumhd/py_cluster.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pyclust --cpus-per-task=8 --mem=128G --time=12:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate pyspatial
export TMPDIR=~/tmp
python ~/projects/crc_visiumhd/py_cluster.py
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch py_cluster.sh
squeue -u weberl
```
