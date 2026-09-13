# Methods

## Data

Visium HD samples P1CRC, P2CRC, P5CRC (Space Ranger binned outputs, 8 um bins) from 10x Genomics public datasets. Single-cell reference: the authors' Chromium Flex matrix (GEO GSE280318) with their `SingleCell_MetaData.csv.gz` (Level1: 9 classes; Level2: 38 subtypes). Authors' per-bin deconvolution and cluster labels: github.com/10XGenomics/HumanColonCancer_VisiumHD. Two conda environments on the IRIS SLURM cluster: `rspatial` (R 4.3, Seurat v5, spacexr) and `pyspatial` (Python, scanpy).

## QC

Bins kept: >= 100 transcripts, >= 50 genes, < 30% mitochondrial. Same thresholds in both pipelines, identical bin sets retained. After QC: 286,973 (P1), 428,405 (P2), 373,156 (P5).

## Clustering

R: log-normalization, 2,000 HVGs, PCA (30 PCs), Louvain (resolution 0.6). Per sample: 13 / 18 / 20 clusters. SCTransform was tried first and dropped (memory).
Python: same QC, log-normalization, 2,000 HVGs, PCA (30), Leiden (resolution 0.4).
Joint: the three samples merged and clustered together (R: 19 clusters; Python: `ad.concat`, same parameters). Used for Figure 3a; tumor clusters are patient-specific, normal cell types shared.

## Annotation

Cluster identity = modal RCTD cell type per cluster (majority vote), supported by FindAllMarkers (adjusted p < 0.05) and GO enrichment (clusterProfiler). Joint tumor clusters assigned to a patient when > 60% of bins come from one sample, then matched to the authors' subtypes by cross-tabulation: P1 = Tumor II, P2 = Tumor III, P5 = Tumor IV. Colors sampled from the published figure.

## Deconvolution

RCTD, doublet mode, per sample, reference from the authors' single cells (QCFilter Keep, classes >= 25 cells). Three references: coarse (9 classes, Level1); Level2 (30 classes, abandoned: out-of-memory and 48 h timeouts); custom (14 classes: Level1 with B cells split into Plasma / Mature B, T cells into CD4 T / CD8 T NK, Fibroblast into CAF / Fibroblast). Figure 3b uses the 14-class reference. All runs use the spacexr PR #206 fork; the stock release stalled at this scale.

## Benchmark

Per-bin agreement of our RCTD `first_type` with the authors' `DeconvolutionLabel1`, joined by barcode, rejects excluded. Two corrections: same patient (the flagship 10x dataset is P2CRC, shown by identical tissue images and 99.99% barcode overlap) and same vocabulary (fine subtypes collapsed to coarse classes by majority-vote lookup, `results/label_map.csv`). Coarse reference: 88.7 / 89.6 / 76.0%. The P5 and P1 losses trace to pooled B-cell profiles absorbing follicular B bins into T cells; the 14-class reference gives 91.7 / 91.8 / 91.4%.

## ARI

Adjusted Rand index on barcode-matched bins. Joint R vs. authors' clusters: 0.816 / 0.605 / 0.730 (P1/P2/P5). Joint R vs. Python: 0.355 overall. Per-sample R vs. Python: 0.456 / 0.341 / 0.486.

## Figures

ggplot2 point maps (size 0.25, coord_fixed), orientation matched to the paper computationally, stitched with PIL. Scale bars from `scalefactors_json.json`: 1 mm main panels, 80 um insets; P1 bars on the right as in the paper. Figure 3c: PIGR, CEACAM6, COL1A1, grey-to-red, limits 0-8, insets picked by maximal signal.

## Deviations from the paper

Log-normalization instead of SCTransform; 9/14-class references instead of the 30-class Level2; Louvain (R) and Leiden (Python) at different resolutions, kept deliberately for cross-pipeline comparison; no sketch subsampling, all bins processed.
