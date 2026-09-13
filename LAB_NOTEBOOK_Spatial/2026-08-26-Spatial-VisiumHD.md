# 2026-08-26 — Spatial Visium HD

Colorectal cancer Visium HD, Nature Genetics 2025

Compute allocation:

```bash
srun --pty --cpus-per-task=8 --mem=128G --time=8:00:00 bash
```

Environment:

```bash
conda activate rspatial
```

R: restoring last saved object

```r
library(Seurat)
library(ggplot2)
obj <- readRDS("~/projects/crc_visiumhd/02_clustered.rds")
obj #QC, filtered, clustered
```

```text
An object of class Seurat
18085 features across 428405 samples within 1 assay
Active assay: Spatial.008um (18085 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 1 dimensional reduction calculated: pca
 1 spatial field of view present: slice1.008um
```

```r
table(obj$seurat_clusters)
```

```text
   0     1     2     3     4     5     6     7     8     9    10    11    12
71124 61428 47005 40878 32147 26563 25513 23314 19981 19894 15775 10430  7871
   13    14    15    16    17
 7273  6796  5573  4919  1921
```

Identity check: markers panels

```r
genes <- c("EPCAM", "PTPRC", "CD3D", "CD68", "COL1A1")
df2 <- GetTissueCoordinates(obj)                        # collecting each bin's (x,y) position on the slide
df2 <- cbind(df2, FetchData(obj, vars = genes))         # collecting normalized expression fo the 5 genes, per bin.

library(tidyr)
long <- pivot_longer(df2, all_of(genes), names_to = "gene", values_to = "expr") #reshape table: each bin gets 5 row (=long format): format ggplot requires for drawing a panel per gene.
```

Drawing five spatial maps, one per marker; colour = normalized expression

```r
p4 <- ggplot(long, aes(x = x, y = y, colour = expr)) + geom_point(size = 0.05) + scale_colour_viridis_c() + coord_fixed() + theme_void() + facet_wrap(~gene, ncol = 3)

ggsave("~/projects/crc_visiumhd/figures/03_markers_spatial.png", p4, width = 14, height = 9)
```

Mean normalized expression of each marker, per cluster

```r
avg <- AverageExpression(obj, features = genes, group.by = "seurat_clusters", layer = "data")
round(avg$Spatial.008um, 2)
```

```text
PTPRC   0.01  0.01  0.65   0.29  0.22  0.43  0.03  0.02  0.01  0.26  0.02  0.06
EPCAM  21.12 23.27  7.02   1.26  2.86  1.19 28.54 26.91 24.00  0.55 17.52 49.56
CD3D    0.03  0.02  1.10   0.57  0.34  0.92  0.04  0.05  0.04  0.37  0.05  0.10
CD68    0.26  0.25  7.22   1.51  1.35  1.83  0.51  0.77  0.43  1.80  0.54  1.75
COL1A1  0.45  0.82 40.39 272.40 43.38 17.69  1.06  0.55  1.34 12.80  1.98  1.08

PTPRC    0.42  1.29  0.30  0.64  0.46  0.02
EPCAM    0.30  5.97 16.14 15.51  6.79 22.31
CD3D     0.84  5.21  0.43  0.21  0.22  0.05
CD68     2.47  1.44  2.00  1.33 17.82  0.63
COL1A1 140.24 20.17  5.64  1.67 29.85  0.95
```

```r
colnames(avg$Spatial.008um)
```

```text
[1] "g0"  "g1"  "g2"  "g3"  "g4"  "g5"  "g6"  "g7"  "g8"  "g9"  "g10" "g11"
[13] "g12" "g13" "g14" "g15" "g16" "g17"
```

```r
round(AverageExpression(obj, features = c("ACTA2", "DES", "MYH11"), group.by = "seurat_clusters", layer = "data")$Spatial.008um, 2) # assessing if cluster 9 (low in all 5 markers) is smooth-muscle.
```

```text
DES   0.27 0.30 1.94  3.50  4.41 1.96 1.89 0.40 0.32 201.95 0.75 1.24 6.62 3.57
ACTA2 0.08 0.10 4.59 17.62 14.76 2.38 0.44 0.08 0.14  26.38 0.26 0.31 3.76 3.24
MYH11 0.08 0.07 1.23  3.61  4.64 1.28 0.98 0.15 0.11  79.55 0.25 0.52 2.90 1.60

DES   3.22 1.42 0.86 0.24
ACTA2 2.66 0.45 1.49 0.11
MYH11 3.24 0.59 0.42 0.05
```

Annotating clusters:

```r
comp <- c("0"="Epithelial", "1"="Epithelial", "2"="Stroma_immune", "3"="Stroma", "4"="Stroma", "5"="Stroma", "6"="Epithelial", "7"="Epithelial",  "8"="Epithelial", "9"="Stroma", "10"="Epithelial", "11"="Epithelial", "12"="Stroma", "13"="Lymphoid", "14"="Epithelial", "15"="Epithelial", "16"="Stroma_immune", "17"="Epithelial")
obj$compartment <- unname(comp[as.character(obj$seurat_clusters)])
table(obj$compartment)                                              # bin counts per compartment
```

```text
Epithelial      Lymphoid        Stroma Stroma_immune
       241855          7273        127353         51924
```

Compartment map + save checkpoint

```r
df <- GetTissueCoordinates(obj)
df$compartment <- obj$compartment
p5 <- ggplot(df, aes(x = x, y = y, colour = compartment)) + geom_point(size = 0.1) + coord_fixed() + theme_void() + guides(colour = guide_legend(override.aes = list(size = 4)))
ggsave("~/projects/crc_visiumhd/figures/04_compartments.png", p5, width = 8, height = 6)
saveRDS(obj, "~/projects/crc_visiumhd/03_annotated.rds")
```

## scRNAseq referenced-based deconvolution

Cell type composition per bin is estimated using RCTD (spacexr), using the author's scRNAseq reference: provides cell-types proportions independent of cluster boundaries.

Installing spacexr package inside the conda env (previously failed under the system toolchain)
In the terminal (not R):

```bash
conda activate rspatial
Rscript -e 'remotes::install_github("dmcable/spacexr")'
```

Downloading the single-cell reference:

```bash
cd ~/data/crc_visiumhd
git clone https://github.com/10XGenomics/HumanColonCancer_VisiumHD.git
ls HumanColonCancer_VisiumHD/MetaData/
```

```text
DeconvolutionResults_P1CRC.csv.gz  DeconvolutionResults_P5CRC.csv.gz  P2CRC_Metadata.parquet  SingleCell_MetaData.csv.gz
DeconvolutionResults_P2CRC.csv.gz  P1CRC_Metadata.parquet             P5CRC_Metadata.parquet
```

Important files: DeconvolutionResults_P1CRC.csv.gz  SingleCell_MetaData.csv.gz

```bash
zcat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz | head -3
```

```text
Barcode,Patient,BC,QCFilter,Level1,Level2,UMAP1,UMAP2
AAACAAGCAACAGCACACTTTAGG-1,P2CRC,BC1,Remove,QC_Filtered,QC_Filtered,3.98688117424738,4.35304268488498
AAACAAGCAACAGCTAACTTTAGG-1,P2CRC,BC1,Keep,B cells,Plasma,-10.475996281777,-0.267628289586872
```
