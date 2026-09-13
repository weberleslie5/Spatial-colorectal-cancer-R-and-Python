# 2026-09-05 — Spatial Visium HD

Comparing P2 and P5 analyses done in R and Python: ARI comparison (Adjusted Rand Index) of the cluster similarity.

R: export R clusters

```bash
cat > ~/projects/crc_visiumhd/r_export_clusters.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
df <- data.frame(barcode = colnames(obj), cluster = obj$seurat_clusters)
write.csv(df, paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"), row.names = FALSE)
EOF
```

Python: ARI + cluster map

```bash
cat > ~/projects/crc_visiumhd/py_ari_sample.py << 'EOF'
import sys, os, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

sample = sys.argv[1]
home = os.path.expanduser("~")
adata = sc.read_h5ad(f"{home}/projects/crc_visiumhd/py_02_clustered_{sample}.h5ad")
r = pd.read_csv(f"{home}/projects/crc_visiumhd/r_clusters_{sample}.csv").set_index("barcode")

common = adata.obs_names.intersection(r.index)
ari = adjusted_rand_score(r.loc[common, "cluster"], adata.obs.loc[common, "leiden_low"])
print(f"{sample}: {len(common)} shared bins, ARI = {ari:.3f}")

xy = adata.obsm["spatial"]
plt.figure(figsize=(8,8))
plt.scatter(xy[:,0], xy[:,1], s=0.5, c=adata.obs["leiden_low"].astype(int), cmap="tab20")
plt.gca().invert_yaxis(); plt.axis('off'); plt.title(f"{sample} Python clusters")
plt.savefig(f"{home}/projects/crc_visiumhd/figures/py_02_clusters_{sample}.png", dpi = 150)
EOF
```

Job ticket:

```bash
cat > ~/projects/crc_visiumhd/ari_sample.sh << 'EOF'
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
EOF
```

Submitting:

```bash
cd ~/projects/crc_visiumhd
sbatch ari_sample.sh P2CRC
sbatch ari_sample.sh P5CRC

grep "ARI" ~/projects/crc_visiumhd/ari_*.log
```

```text
/home/weberl/projects/crc_visiumhd/ari_11666508.log:P2CRC: 428405 shared bins, ARI = 0.341
/home/weberl/projects/crc_visiumhd/ari_11666509.log:P5CRC: 373156 shared bins, ARI = 0.486
```

```bash
ls -l ~/projects/crc_visiumhd/02_clustered*.rds ~/projects/crc_visiumhd/py_02_clustered*.h5ad
```

```text
-rw-rw-r-- 1 weberl weberl 1487592894 Sep  4 18:01 /home/weberl/projects/crc_visiumhd/02_clustered_P1CRC.rds
-rw-rw-r-- 1 weberl weberl 1487593245 Sep  4 17:13 /home/weberl/projects/crc_visiumhd/02_clustered_P2CRC.rds
-rw-rw-r-- 1 weberl weberl 1233658903 Sep  4 17:08 /home/weberl/projects/crc_visiumhd/02_clustered_P5CRC.rds
-rw-rw-r-- 1 weberl weberl 1487592894 Aug 25 14:09 /home/weberl/projects/crc_visiumhd/02_clustered.rds
-rw-rw-r-- 1 weberl weberl 8771327489 Sep  2 15:42 /home/weberl/projects/crc_visiumhd/py_02_clustered.h5ad
-rw-rw-r-- 1 weberl weberl 8770887644 Sep  4 17:02 /home/weberl/projects/crc_visiumhd/py_02_clustered_P2CRC.h5ad
-rw-rw-r-- 1 weberl weberl 7425664866 Sep  4 17:01 /home/weberl/projects/crc_visiumhd/py_02_clustered_P5CRC.h5ad
```

```bash
wc -l ~/projects/crc_visiumhd/r_clusters*.csv
```

```text
  428406 /home/weberl/projects/crc_visiumhd/r_clusters.csv
  428406 /home/weberl/projects/crc_visiumhd/r_clusters_P2CRC.csv
  373157 /home/weberl/projects/crc_visiumhd/r_clusters_P5CRC.csv
 1229969 total
```

The original first simple analysed was in fact P2:

```bash
for s in P1CRC P2CRC P5CRC; do
  n=$(comm -12 <(cut -d, -f1 r_clusters.csv | tr -d '"' | sort) <(zcat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_${s}.csv.gz | awk -F, '$2 != "NA"' | cut -d, -f1 | sort) | wc -l)
  echo "our_P1 vs authors_$s: $n shared barcodes"
done
```

```text
our_P1 vs authors_P1CRC: 235645 shared barcodes
our_P1 vs authors_P2CRC: 428355 shared barcodes
our_P1 vs authors_P5CRC: 244421 shared barcodes
```

```bash
for s in P1CRC P2CRC P5CRC; do
  n=$(comm -12 <(cut -d, -f1 r_clusters_P5CRC.csv | tr -d '"' | sort) <(zcat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_${s}.csv.gz | awk -F, '$2 != "NA"' | cut -d, -f1 | sort) | wc -l)
  echo "our_P5 vs authors_$s: $n shared barcodes"
done
```

```text
our_P5 vs authors_P1CRC: 162472 shared barcodes
our_P5 vs authors_P2CRC: 233905 shared barcodes
our_P5 vs authors_P5CRC: 373051 shared barcodes
```

Missing P1
Downloading P1:

```bash
mkdir -p ~/data/crc_visiumhd/P1CRC
cd ~/data/crc_visiumhd/P1CRC
nohup wget -c "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer_P1/Visium_HD_Human_Colon_Cancer_P1_binned_outputs.tar.gz" > wget_P1.log 2>&1 &
```

Deleting duplicates:

```bash
cd ~/projects/crc_visiumhd
rm -f 02_clustered_P1CRC.rds 04_rctd_P1CRC.rds 02_clustered.rds py_02_clustered.h5ad r_clusters.csv
rm -rf ~/data/crc_visiumhd/binned_outputs
```

Renaming the files kept:

```bash
mv 03_annotated.rds 03_annotated_P2CRC.rds        # the compartment annotation
mv 04_rctd.rds 04_rctd_stock_P2CRC.rds            # original spacexr run
```

```bash
tail -2 ~/data/crc_visiumhd/P1CRC/wget_P1.log
```
