# 2026-08-24 — Spatial Visium HD

We need to create the exact folder Seurat expects for bin.size = 8 when using Load10X_Spatial

```bash
cd ~/data/crc_visiumhd
mkdir -p binned_outputs
mv square_008um binned_outputs/
```

We also need more ressources: a compute node with memory:

```bash
srun --pty --cpus-per-task=4 --mem=96G --time=6:00:00 bash
conda activate rspatial
R
```

Then, in R

```r
library(Seurat)
obj <- Load10X_Spatial("~/data/crc_visiumhd", bin.size = 8)
obj
```

```text
An object of class Seurat
18085 features across 545913 samples within 1 assay                 # 18,085 genes and, 545,913 8um bins
Active assay: Spatial.008um (18085 features, 0 variable features)
 1 layer present: counts                                            # raw transcriptional counts
 1 spatial field of view present: slice1.008um                      # tissue image and coordinates
 = Loaded.
```

## QC

1- Mitochondrial % per bin
High mito % = stressed/degraded tissue.
Human mito genes start with "MT-".
They are stored as a new column in the object's per-bin table (meta.data)

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
```

2- Looking at the numbers
nCount = transcripts per bin (depth)
nFeature = genes detected per bin (complexity)

```r
summary(obj$nCount_Spatial.008um)
summary(obj$nFeature_Spatial.008um)
summary(obj$percent.mt)
```

```text
summary(obj$nCount_Spatial.008um)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    0.0   124.0   380.0   488.4   770.0  4145.0
> summary(obj$nFeature_Spatial.008um)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    0.0   110.0   299.0   393.2   635.0  2594.0
> summary(obj$percent.mt)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
  0.000   1.942   3.922   4.837   6.340 100.000      6
```

3- In plots
Violin = distribution shape. Terminal has no screen so saved to pdf.

```r
dir.create("~/projects/crc_visiumhd/figures", recursive = TRUE)

p <- VlnPlot(obj,
            features = c("nCount_Spatial.008um",
                        "nFeature_Spatial.008um",
                        "percent.mt"),
            pt.size = 0)

ggplot2::ggsave("~/projects/crc_visiumhd/figures/01_qc_violins.png", p, width = 12, height = 5)
```

4- Filtering the bins

From these numbers, here are the cutoffs:
Keeping bins with >= 100 transcripts, >= 50 genes and <30% mito

```r
n_before <- ncol(obj)       #ncol = number of bins

obj <- subset(obj, nCount_Spatial.008um >= 100 &
                    nFeature_Spatial.008um >= 50 &
                    percent.mt < 30)

n_after <- ncol(obj)

ncol(obj)
```

```text
[1] 428405
```

5- Tissue image check:

```r
library(ggplot2)
df <- GetTissueCoordinates(obj)
df$depth <- obj$nCount_Spatial.008um

p2 <- ggplot(df, aes(x = x, y = y, colour = depth)) + geom_point(size = 0.1) + scale_colour_viridis_c() + coord_fixed() + theme_void()

ggsave("~/projects/crc_visiumhd/figures/01_qc_spatial.png", p2, width = 7, height = 6)
```

Last, saving the R object on disk: obj filtered, with mito % computed:
Still on R:

```r
saveRDS(obj, "~/projects/crc_visiumhd/01_qc.rds")
```
