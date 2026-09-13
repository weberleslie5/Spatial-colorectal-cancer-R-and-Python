# 2026-08-25 — Spatial Visium HD

Colorectal cancer Visium HD, Nature Genetics 2025

Cluster: Normalize (SCTransform) > PCA > Louvain > plot clusters on tissue.

Compute node with more memory (clustering = heavy)

```bash
srun --pty --cpus-per-task=8 --mem=128G --time=8:00:00 bash
conda activate rspatial
R
```

```r
library(Seurat)
obj <- readRDS("~/projects/crc_visiumhd/01_qc.rds")
```

## A. Normalization

```r
obj <- SCTransform(obj, assay = "Spatial.008um")
```

```text
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 18045 by 428405
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 5000 cells
Found 298 outliers - those will be ignored in fitting/regularization step
Second step: Get residuals using fitted parameters for 18045 genes
Computing corrected count matrix for 18045 genes
Killed
```

SCTransform builds a corrected matrix for the 18045 genes x 428405 bins > too much.
Seurat's HD tutorial uses log-normalization: I switch to that:

```r
library(Seurat)
library(ggplot2)
obj <- readRDS("~/projects/crc_visiumhd/01_qc.rds")

obj <- NormalizeData(obj)                           #log normalize per bin
obj <- FindVariableFeatures(obj, nfeatures = 2000)  #keeping 2000 most informative genes
obj <- ScaleData(obj)    #center: substracts gene average (mean=0). scale: divide by its sd (every gene varies on the same scale)
```

Then: compression and clustering:

```r
obj <- RunPCA(obj, npcs = 30)            # 30 dimensions instead of the 2000 dimensions represented by 2000 highly variable genes
```

```text
PC_ 1
Positive:  COL3A1, VIM, COL1A1, COL1A2, SPARC, COL6A2, MMP2, TAGLN, IGFBP7, LGALS1
           COL5A1, LUM, COL6A3, CALD1, DCN, COL4A1, VCAN, AEBP1, C1S, COL6A1
           TIMP2, C1R, A2M, MYL9, ACTA2, THBS2, FSTL1, THY1, COL5A2, COL4A2
Negative:  PHGR1, NQO1, FABP1, LCN2, TSPAN8, GPX2, EREG, SELENBP1, MCM4, LRATD1
           PTP4A3, TRIM31, FASN, MUC13, KRT23, AREG, HIST1H1B, CXCL2, DHCR24, S100P
           FXYD3, CXCL3, ASS1, FAM83H, VEGFA, MUC1, ID1, TSPOAP1, SPINK1, PI3
PC_ 2
Positive:  PIGR, FCGBP, MUC2, CLCA1, CA2, TSPAN1, TFF3, PLA2G2A, SPINK4, ZG16
           AGR2, OLFM4, KRT20, MUC4, SLC26A3, REG4, DUOX2, B3GALT5, ITLN1, GUCA2A
           CA4, AGR3, SELENOP, PARM1, CDHR5, TSPAN8, CCL28, CA1, ATP8B1, PADI2
Negative:  HIST1H1B, PTP4A3, MCM4, EREG, HIST1H1C, NQO1, MYBL2, FASN, NCOA7, APCDD1
           TPX2, KRT23, TK1, SLC7A5, MCM3, SNRPB, AREG, MMP12, CDC25B, SOD2
           MCM2, TUBA1C, IL32, HIST1H2BH, CXCL2, SCD, AXIN2, HIST1H1D, LRATD1, TSPOAP1
PC_ 3
Positive:  MMP12, CD74, LYZ, SRGN, IFI30, C1QC, CTSB, C1QB, C1QA, FCER1G
           CD68, SPI1, LCP1, RGS1, LAPTM5, CD4, F3, MMP9, MPEG1, PLA2G7
           TYROBP, AIF1, CTSC, PSAP, CD14, SGK1, APOE, IL7R, HMOX1, MS4A6A
Negative:  SFRP4, THBS2, TAGLN, COMP, ELN, AEBP1, GREM1, SFRP2, FN1, COL1A2
           COL11A1, LUM, MMP2, ACTG2, MYL9, TIMP3, COL1A1, DES, CTHRC1, COL3A1
           MGP, EFEMP1, CNN1, MYH11, TPM2, DCN, COL8A1, FBLN1, ITGA11, CCDC80
PC_ 4
Positive:  PI3, LCN2, COL1A1, COL1A2, CEACAM7, DUOX2, TSPAN1, SLC26A3, COL5A1, SPARC
           MUC13, KRT20, COL6A3, COL3A1, ATP1B1, ID1, COL5A2, CLCA4, TRIM31, COL12A1
           CCL20, THBS2, MMP2, LUM, AEBP1, GCNT3, GUCA2A, COL4A1, CD55, CA4
Negative:  JCHAIN, IGKC, IGHG1, MZB1, IGHA1, DES, PIM2, IGHM, DERL3, TENT5C
           IGLC1, CD38, IGHG3, CD79A, MYH11, SRGN, FCRL5, POU2AF1, ADA2, BTG2
           XBP1, LGR5, SEC11C, IRF4, TXNDC5, CYTIP, SYNM, LYZ, TNFRSF17, SYNPO2
PC_ 5
Positive:  VWF, PLVAP, AQP1, PECAM1, CALCRL, EGFL7, RGS5, ADGRL4, ENG, MCAM
           COL4A1, PODXL, CD34, SLCO2A1, ESAM, HSPG2, A2M, COL4A2, CD93, FLT1
           EPAS1, ADAMTS1, CDH5, RAMP3, ADGRF5, NOTCH3, ADAMTS9, SPARCL1, ECSCR, MMRN2
Negative:  LUM, SFRP4, THBS2, DCN, MMP2, C3, AEBP1, SFRP2, COL1A2, COMP
           COL1A1, COL11A1, FBLN1, CTHRC1, COL6A3, GREM1, COL3A1, VCAN, C1S, APOE
           ITGA11, ADAM12, ELN, CTSK, TIMP2, GPNMB, COL8A1, EFEMP1, C1QC, CTSB
```

```r
obj <- FindNeighbors(obj, dims = 1:30)    # dims = 1:30 means using PCA dimensions 1 to 30 # compute nearest neighbor graph
obj <- FindClusters(obj, resolution = 0.6)
```

Plotting and saving:

```r
df <- GetTissueCoordinates(obj)         # building this df from the current obj
df$cluster <- obj$seurat_clusters       # adds the cluster column to df

p3 <- ggplot(df, aes(x = x, y=y, colour = cluster)) + geom_point(size = 0.1) + coord_fixed() + theme_void() + guides(colour = guide_legend(override.aes = list(size = 4)))

ggsave("~/projects/crc_visiumhd/figures/02_clusters_spatial.png", p3, width = 8, height = 6)

saveRDS(obj, "~/projects/crc_visiumhd/02_clustered.rds")
```
