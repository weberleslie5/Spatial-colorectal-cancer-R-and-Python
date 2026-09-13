# 2026-09-10 — Spatial Visium HD

```bash
cd ~/projects/crc_visiumhd
ls -lh 02_clustered_*.rds py_02_clustered_*.h5ad 04_rctd_*.rds r_clusters_*.csv 2>/dev/null
```

```text
-rw-rw-r-- 1 weberl weberl 578M Sep  6 12:45 02_clustered_P1CRC.rds
-rw-rw-r-- 1 weberl weberl 1.4G Sep  4 17:13 02_clustered_P2CRC.rds
-rw-rw-r-- 1 weberl weberl 1.2G Sep  4 17:08 02_clustered_P5CRC.rds
-rw-rw-r-- 1 weberl weberl 299M Sep  6 16:08 04_rctd_P1CRC.rds
-rw-rw-r-- 1 weberl weberl 746M Sep  5 02:53 04_rctd_P2CRC.rds
-rw-rw-r-- 1 weberl weberl 587M Sep  5 00:30 04_rctd_P5CRC.rds
-rw-rw-r-- 1 weberl weberl 756M Aug 28 00:16 04_rctd_stock_P2CRC.rds
-rw-rw-r-- 1 weberl weberl 5.0G Sep  6 12:36 py_02_clustered_P1CRC.h5ad
-rw-rw-r-- 1 weberl weberl 8.2G Sep  4 17:02 py_02_clustered_P2CRC.h5ad
-rw-rw-r-- 1 weberl weberl 7.0G Sep  4 17:01 py_02_clustered_P5CRC.h5ad
-rw-rw-r-- 1 weberl weberl 7.7M Sep  6 13:09 r_clusters_P1CRC.csv
-rw-rw-r-- 1 weberl weberl  12M Sep  5 15:48 r_clusters_P2CRC.csv
-rw-rw-r-- 1 weberl weberl  11M Sep  5 15:48 r_clusters_P5CRC.csv
```

```bash
grep "ARI" ~/projects/crc_visiumhd/ari_*.log
ls ~/projects/crc_visiumhd/figures
```

```text
01_qc_spatial.png  02_clusters_P1CRC.png  02_clusters_P5CRC.png    03_markers_spatial.png  05_celltypes.png        py_02_clusters_coarse.png  py_02_clusters_P2CRC.png  py_02_clusters.png
01_qc_violins.png  02_clusters_P2CRC.png  02_clusters_spatial.png  04_compartments.png     06_immune_vs_tumor.png  py_02_clusters_P1CRC.png   py_02_clusters_P5CRC.png
```

## Reproducing Fig 3a, using the R pipeline

How many clusters per sample?

```bash
cd ~/projects/crc_visiumhd
for s in P1CRC P2CRC P5CRC; do
    echo -n "$s: "
    cut -d, -f2 r_clusters_${s}.csv | tr -d '"' | grep -v cluster | sort -u | wc -l
done
```

```text
P1CRC: 13
P2CRC: 18
P5CRC: 20
```

Install package for gene enrichment and cluster identity

```bash
srun --cpus-per-task=4 --mem=32G --time=1:00:00 --pty bash
mamba install -n rspatial -y -c conda-forge -c bioconda bioconductor-clusterprofiler bioconductor-org.hs.eg.db
exit
```

Markers per cluster: all significant upregulated genes:

```bash
cat > ~/projects/crc_visiumhd/markers_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
m <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sig <- m[m$p_val_adj < 0.05, ]
write.csv(sig, paste0("~/projects/crc_visiumhd/markers_", sample, ".csv"), row.names = FALSE)
EOF
```

```bash
cat > ~/projects/crc_visiumhd/markers_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=mark
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/markers_sample.R "$1"
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch markers_sample.sh P1CRC
sbatch markers_sample.sh P2CRC
sbatch markers_sample.sh P5CRC
```

On the R clustering: extracting top significantly enriched markers per cluster:

```bash
cat > ~/projects/crc_visiumhd/pathways_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(clusterProfiler); library(org.Hs.eg.db); library(dplyr)
m <- read.csv(paste0("~/projects/crc_visiumhd/markers_", sample, ".csv"))
out <- list()
for (cl in unique(m$cluster)) {
    g <- m$gene[m$cluster == cl]
    e <- enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
    df <- as.data.frame(e)
    if (nrow(df) > 0) out[[as.character(cl)]] <- data.frame(cluster=cl, head(df[, c("Description", "p.adjust", "Count")], 10))
}
write.csv(bind_rows(out), paste0("~/projects/crc_visiumhd/pathways_", sample, ".csv"), row.names = FALSE)
EOF
```

Job ticket:

```bash
cat > ~/projects/crc_visiumhd/pathways_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=path
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/pathways_sample.R $1
EOF
```

```bash
ls markers_*.csv
```

```text
markers_P1CRC.csv  markers_P2CRC.csv  markers_P5CRC.csv
```

```bash
cd ~/projects/crc_visiumhd
sbatch pathways_sample.sh P1CRC
sbatch pathways_sample.sh P2CRC
sbatch pathways_sample.sh P5CRC
```

## Benchmark test

Independently: Benchmark test: following RCTD on each of the patients: comparing the assigned cell identity from my analysis and the paper's (ARI tests clustering: separate values):

```bash
cat > ~/projects/crc_visiumhd/benchmark_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), mine = as.character(res$first_type))[res$spot_class != "reject",]
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
authors <- authors[authors$DeconvolutionClass != "reject" & !is.na(authors$DeconvolutionLabel1), c("barcode", "DeconvolutionLabel1")]
m <- merge(mine, authors, by = "barcode")
cat(sample, "bins compared:", nrow(m), "\n")
cat(sample, "agreement:", mean(m$mine == m$DeconvolutionLabel1), "\n")
print(table(mine = m$mine, authors = m$DeconvolutionLabel1))
EOF
```

```bash
cat > ~/projects/crc_visiumhd/benchmark_sample.sh << 'EOF'
#! /bin/bash
#SBATCH --job-name=bench
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/benchmark_sample.R "$1"
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch benchmark_sample.sh P1CRC
sbatch benchmark_sample.sh P2CRC
sbatch benchmark_sample.sh P5CRC
```

Once run completed:

```bash
grep "agreement" ~/projects/crc_visiumhd/bench_*.log
```

```text
/home/weberl/projects/crc_visiumhd/bench_12019479.log:P1CRC agreement: 0.06541502
/home/weberl/projects/crc_visiumhd/bench_12019480.log:P2CRC agreement: 0.04575582
/home/weberl/projects/crc_visiumhd/bench_12019483.log:P5CRC agreement: 0.1094651
```

Very low... problem with cell labelling concordance?

```bash
tail -30 ~/projects/crc_visiumhd/bench_12019480.log
```

```text
Endothelial                       0      0       7       24       102
  Fibroblast                       15      1       3        7       261
  Intestinal Epithelial             0     26       6        7        98
  Myeloid                           2      1       9       55       321
  Neuronal                          0      0       0        1         0
  Smooth Muscle                   664      0       0        0         1
  T cells                           1      1       1        8        28
  Tumor                             0     23     130      595    211979
                       authors
mine                    Tumor IV Tumor V Unknown III (SM) Vascular Fibroblast
  B cells                     20       1               13                   0
  Endothelial                  2       2                7                   1
  Fibroblast                   5       2              230                  10
  Intestinal Epithelial        6      52                0                   0
  Myeloid                      7       5               31                   1
  Neuronal                     0       0                0                   1
  Smooth Muscle                1       0              470                   9
  T cells                      4       0               15                   2
  Tumor                       63     144                1                   1
                       authors
mine                       vSM
  B cells                  204
  Endothelial               59
  Fibroblast               553
  Intestinal Epithelial      0
  Myeloid                  402
  Neuronal                 133
  Smooth Muscle          11120
  T cells                  181
  Tumor                      3
```

Authors' labelling:

```bash
zcat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_P2CRC.csv.gz | awk -F, 'NR>1 {print $3}' | sort | uniq -c | sort -rn
```

```text
 271240 NA
 213815 Tumor III
  59322 CAF
  30300 Goblet
  18043 Endothelial
  15078 Macrophage
  13225 vSM
  12066 Plasma
  10383 Enterocyte
   9243 Proliferating Macrophages
   9164 Proliferating Fibroblast
   6107 Myofibroblast
   6012 Pericytes
   4647 Neutrophil
   3711 CD4 T cell
   3660 Proliferating Immune II
   3363 Fibroblast
   2447 CD8 T cell
   1935 Lymphatic Endothelial
   1543 Tumor II
   1345 Smooth Muscle
   1112 Unknown III (SM)
    857 Enteric Glial
    726 cDC I
    476 mRegDC
    417 Mast
    343 Neuroendocrine
    276 Adipocyte
    258 Tumor I
    249 Tumor V
    213 Mature B
    171 Epithelial
    141 Tumor IV
    105 pDC
     88 Tuft
     62 Vascular Fibroblast
     60 SM Stress Response
     29 Memory B
     12 NK
```

Finding the Fine (Level2) to coarse (Level1) analysis, grouping cell types in the 9 coarse subtypes the author provide:

```bash
script -e 'meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz"); tab <- table(meta$Level2, meta$Level1); print(data.frame(fine = rownames(tab), coarse = colnames(tab)[max.col(tab)]), row.names = FALSE)'
```

```text
                      fine                coarse
                 Adipocyte            Fibroblast
                       CAF            Fibroblast
                CD4 T cell               T cells
                CD8 T cell               T cells
                     cDC I               Myeloid
               Endothelial           Endothelial
             Enteric Glial              Neuronal
                Enterocyte Intestinal Epithelial
                Epithelial         Smooth Muscle
                Fibroblast            Fibroblast
                    Goblet Intestinal Epithelial
     Lymphatic Endothelial           Endothelial
                Macrophage               Myeloid
                      Mast               Myeloid
                  Mature B               B cells
                  Memory B               B cells
                    mRegDC               Myeloid
             Myofibroblast            Fibroblast
            Neuroendocrine              Neuronal
                Neutrophil               Myeloid
                        NK               T cells
                       pDC               Myeloid
                 Pericytes            Fibroblast
                    Plasma               B cells
  Proliferating Fibroblast            Fibroblast
   Proliferating Immune II               B cells
 Proliferating Macrophages               Myeloid
               QC_Filtered           QC_Filtered
        SM Stress Response         Smooth Muscle
             Smooth Muscle         Smooth Muscle
                      Tuft              Neuronal
                   Tumor I                 Tumor
                  Tumor II                 Tumor
                 Tumor III                 Tumor
                  Tumor IV                 Tumor
                   Tumor V                 Tumor
          Unknown III (SM)         Smooth Muscle
       Vascular Fibroblast            Fibroblast
                       vSM         Smooth Muscle
```

Storing this table:

```bash
Rscript -e 'meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz"); tab <- table(meta$Level2, meta$Level1); write.csv(data.frame(fine = rownames(tab), coarse = colnames(tab)[max.col(tab)]), "~/projects/crc_visiumhd/label_map.csv", row.names = FALSE)'
ls -lh ~/projects/crc_visiumhd/label_map.csv
```

```text
-rw-rw-r-- 1 weberl weberl 1013 Sep 10 15:18 /home/weberl/projects/crc_visiumhd/label_map.csv
```

```bash
ls -lh ~/projects/crc_visiumhd/bench_*.log
```

```text
-rw-rw-r-- 1 weberl weberl 6.8K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019479.log
-rw-rw-r-- 1 weberl weberl 7.2K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019480.log
-rw-rw-r-- 1 weberl weberl 6.8K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019483.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021371.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021372.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021375.log
```

Running the new benchmark call:

```bash
cat > ~/projects/crc_visiumhd/benchmark_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), mine = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
authors <- authors[authors$DeconvolutionClass != "reject" & !is.na(authors$DeconvolutionLabel1), c("barcode", "DeconvolutionLabel1")]

meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
lookup <- setNames(meta$Level1, meta$Level2)               # authors' own subtype-to-class table
authors$L1 <- lookup[authors$DeconvolutionLabel1]          # collapse subtypes to 9 classes
authors$L1[is.na(authors$L1) & grepl("Tumor", authors$DeconvolutionLabel1)] <- "Tumor"          # spatial-only labels
authors$L1[is.na(authors$L1) & grepl("SM", authors$DeconvolutionLabel1)] <- "Smooth Muscle"
authors$L1[is.na(authors$L1) & grepl("Fibro", authors$DeconvolutionLabel1)] <- "Fibroblast"

m <- merge(mine, authors, by = "barcode")
m <- m[!is.na(m$L1), ]                                     # drop anything still unmapped
cat(sample, "bins compared:", nrow(m), "\n")
cat(sample, "agreement:", mean(m$mine == m$L1), "\n")
print(table(mine = m$mine, authors = m$L1))
EOF
cd ~/projects/crc_visiumhd
sbatch benchmark_sample.sh P1CRC
sbatch benchmark_sample.sh P2CRC
sbatch benchmark_sample.sh P5CRC
```

```bash
ls -lh ~/projects/crc_visiumhd/bench_*.log
```

```text
-rw-rw-r-- 1 weberl weberl 6.8K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019479.log
-rw-rw-r-- 1 weberl weberl 7.2K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019480.log
-rw-rw-r-- 1 weberl weberl 6.8K Sep 10 14:43 /home/weberl/projects/crc_visiumhd/bench_12019483.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021371.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021372.log
-rw-rw-r-- 1 weberl weberl 1.6K Sep 10 14:58 /home/weberl/projects/crc_visiumhd/bench_12021375.log
```

```bash
grep agreement ~/projects/crc_visiumhd/bench_1202137*.log
```

```text
/home/weberl/projects/crc_visiumhd/bench_12021371.log:P1CRC agreement: 0.8866033
/home/weberl/projects/crc_visiumhd/bench_12021372.log:P2CRC agreement: 0.8963007
/home/weberl/projects/crc_visiumhd/bench_12021375.log:P5CRC agreement: 0.7598441
```

Back to GO of clusters (R analysis):

```bash
s -lh markers_*.csv
```

```text
-rw-rw-r-- 1 weberl weberl 9.6K Sep 10 13:29 markers_P1CRC.csv
-rw-rw-r-- 1 weberl weberl  14K Sep 10 13:33 markers_P2CRC.csv
-rw-rw-r-- 1 weberl weberl  16K Sep 10 13:33 markers_P5CRC.csv
```

```bash
sbatch pathways_sample.sh P1CRC
sbatch pathways_sample.sh P2CRC
sbatch pathways_sample.sh P5CRC
```
