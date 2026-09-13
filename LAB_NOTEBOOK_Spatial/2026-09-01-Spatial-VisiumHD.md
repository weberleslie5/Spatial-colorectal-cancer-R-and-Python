# 2026-09-01 — Spatial Visium HD

Checking on the rctd status:

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
ls -t ~/projects/crc_visiumhd/rctd_*.log | head -2
```

```text
/home/weberl/projects/crc_visiumhd/rctd_10872206.log
/home/weberl/projects/crc_visiumhd/rctd_10489792.log
```

```bash
tail -15 ~/projects/crc_visiumhd/rctd_10872206.log
```

```text
Likelihood value: 1086212.44895894
Sigma value:  0.49
starting worker pid=1167437 on localhost:11608 at 16:04:54.960
starting worker pid=1167460 on localhost:11608 at 16:04:55.148
starting worker pid=1167507 on localhost:11608 at 16:04:55.331
starting worker pid=1167530 on localhost:11608 at 16:04:55.523
Loading required package: spacexr
loaded spacexr and set parent environment
Loading required package: spacexr
loaded spacexr and set parent environment
Loading required package: spacexr
loaded spacexr and set parent environment
Loading required package: spacexr
loaded spacexr and set parent environment
[2026-09-01T15:11:10.101] error: *** JOB 10872206 ON isca209 CANCELLED AT 2026-09-01T15:11:10 DUE TO TIME LIMIT ***
```

Halving the workers doubled the runtime.
New launch: more memory and more time.

```bash
sed -i 's/--time=24:00:00/--time=48:00:00/' ~/projects/crc_visiumhd/run_rctd_L2.sh
sed -i 's/max_cores = 4/max_cores = 8/' ~/projects/crc_visiumhd/run_rctd_L2.R
cd ~/projects/crc_visiumhd && sbatch run_rctd_L2.sh && squeue -u weberl
```

Chapter 2 (python) Sept 1 work did not save.

ARI to compare the clustering of the same bins: R's (Louvain, rds Seurat object) and Python's (Leiden, h5ad object) > extracting both as pdf and joining them by barcodes before running ARI.

```bash
srun --pty --cpus-per-task=4 --mem=96G --time=4:00:00 bash
```

```bash
conda activate rspatial
Rscript -e 'obj <- readRDS("~/projects/crc_visiumhd/02_clustered.rds"); write.csv(data.frame(barcode=colnames(obj), r_cluster=obj$seurat_clusters), "~/projects/crc_visiumhd/r_clusters.csv", row.names=FALSE)' # Louvain labels to csv
```

Python: recluster, save, verify, compare:

```bash
conda activate pyspatial

pip install scikit-learn # ARI function
```

```text
Requirement already satisfied: scikit-learn in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (1.9.0)
Requirement already satisfied: numpy>=1.24.1 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (from scikit-learn) (2.5.2)
Requirement already satisfied: scipy>=1.10.0 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (from scikit-learn) (1.18.0)
Requirement already satisfied: joblib>=1.4.0 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (from scikit-learn) (1.6.0)
Requirement already satisfied: narwhals>=2.0.1 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (from scikit-learn) (2.25.0)
Requirement already satisfied: threadpoolctl>=3.5.0 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (3.6.0)
Requirement already satisfied: cloudpickle>=3.0 in /data1/test03/weberl/miniforge3/envs/pyspatial/lib/python3.14/site-packages (from joblib>=1.4.0->scikit-learn) (3.1.2)
```

```bash
python
```

```python
import scanpy as sc, pandas as pd, os, time
from sklearn.metrics import adjusted_rand_score

adata = sc.read_h5ad("/home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad")
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_low", flavor="igraph", n_iterations=2)
adata.write("/home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad")
print(time.ctime(os.path.getmtime("/home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad")))
```

Comparison:

Note: Shell and R ~ expands to /home/weberl
Python use /home/weberl/...

```python
r = pd.read_csv("/home/weberl/projects/crc_visiumhd/r_clusters.csv").set_index("barcode")   # R labels
both = r.join(pd.DataFrame({"py": adata.obs.leiden_low})).dropna()    # align by barcode
print(adjusted_rand_score(both.r_cluster, both.py))
```

Coarse cluster map:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

xy = adata.obsm["spatial"]
cl = adata.obs.leiden_low.astype(int)
plt.figure(figsize=(9,7))
plt.scatter(xy[:,0], xy[:,1], c=cl % 20, s=0.3, cmap="tab20")
plt.gca().invert_yaxis()
plt.axis("equal"); plt.axis("off")
plt.title("Leiden clusters - Python")
plt.savefig("/home/weberl/projects/crc_visiumhd/figures/py_02_clusters_coarse.png", dpi=150)
```

Immune vs Tumor figure with Level 1 (R):

```r
library(spacexr); library(ggplot2); library(RANN)
rctd <- readRDS("~/projects/crc_visiumhd/04_rctd.rds")
res <- rctd@results$results_df
coords <- as.data.frame(rctd@spatialRNA@coords)[rownames(res), ]
df <- cbind(coords, type = res$first_type)[res$spot_class != "reject", ] # position + cell type

tumor <- df[df$type == "Tumor", c("x", "y")]    # tumor bins positions
df$dist <- nn2(tumor, df[, c("x", "y")], k = 1)$nn.dists[,1]    #each bin's distance to nearest tumor bin
df$band <- cut(df$dist, breaks = c(-1, 0, 50, 100, 200, 400, 800, Inf), labels = c("0", "<50", "50-100", "100-200", "200-400", "400-800", ">800"))  # zones in um

prop <- prop.table(table(df$band, df$type), 1)
pd <- as.data.frame(prop)
names(pd) <- c("band", "type", "fraction")
pd <- pd[pd$type %in% c("T cells", "B cells", "Myeloid", "Fibroblasts"), ]

p7 <- ggplot(pd, aes(band, fraction, colour = type, group = type)) + geom_line(linewidth = 1) + geom_point() + labs(x = "Distance from nearest tumor bin", y = "Fraction of bins", title = "Cell-type composition vs distance from tumor") + theme_minimal()
ggsave("~/projects/crc_visiumhd/figures/06_immune_vs_tumor.png", p7, width = 9, height = 5.5)
```
