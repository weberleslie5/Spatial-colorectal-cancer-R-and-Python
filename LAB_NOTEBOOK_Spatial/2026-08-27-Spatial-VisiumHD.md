# 2026-08-27 — Spatial Visium HD

Building the RCTD reference: joining the chromium count matrix (GEO GSE280318) to SingleCell_MetaData.csv.gz.
To evaluate the composition per bin, RCTD needs to learn cell-type signatures from annotated single cells.

```bash
srun --pty --cpus-per-task=8 --mem=128G --time=8:00:00 bash
conda activate rspatial
```

Checking spacexr and data:

```bash
Rscript -e 'library(spacexr)'
```

```text
Fatal error: cannot create 'R_TempDir'
```

R cannot write temporary folder. Temp files directed to my home:

```bash
mkdir -p ~/tmp
export TMPDIR=~/tmp
Rscript -e 'library(spacexr)'
```

```bash
zcat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz | head -3
```

```text
Barcode,Patient,BC,QCFilter,Level1,Level2,UMAP1,UMAP2
AAACAAGCAACAGCACACTTTAGG-1,P2CRC,BC1,Remove,QC_Filtered,QC_Filtered,3.98688117424738,4.35304268488498
AAACAAGCAACAGCTAACTTTAGG-1,P2CRC,BC1,Keep,B cells,Plasma,-10.475996281777,-0.267628289586872
```

Barcode = cell ID.
QCFilter: keeping only "Keep" rows
Level1 = coarse cell type (used here for RCTD): B cell
Level2 = finer subtype: Plasma cell

We need the gene expression associated with these barcoded cells: the Chromium count matrix (genes x cells).
We join these two files by barcodes: metadata provides cell identity, matrix supplies cell's gene expression.

```bash
grep -i http ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/README.md | head -20
```

```text
R script used to process the single cell (Flex Gene Expression) data.
https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc
```

From this page, downloaded: Feature barcode matrix HDF5 (filtered) 426 MB:

```bash
mkdir -p ~/data/crc_visiumhd/chromium
cd ~/data/crc_visiumhd/chromium
wget -c https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5
ls -lh
```

406M = 406 MiB (mebibyte) = 426 MB

Building the reference:

```bash
R
```

```r
library(Seurat)
library(spacexr)

counts <- Read10X_h5("~/data/crc_visiumhd/chromium/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5") # single-cell count matrix

meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz") # per-cell annotations
meta <- meta[meta$QCFilter == "Keep", ] # keeping only validated cells

common <- intersect(colnames(counts), meta$Barcode) # barcodes present in both matrix and metadata
length(common)
```

```text
[1] 260506
```

```r
counts <- counts[, common]  #keep only annotated cells

rownames(meta) <- meta$Barcode # index metadata by barcode
celltypes <- factor(meta[common, "Level1"]) # cell-type label per cell
names(celltypes) <- common
table(celltypes) # cells per type
```

```text
celltypes
              B cells           Endothelial            Fibroblast
                33438                  7916                 30681
Intestinal Epithelial               Myeloid              Neuronal
                24699                 25002                  4321
        Smooth Muscle               T cells                 Tumor
                41185                 29364                 63900
```

```r
ref <- Reference(counts, celltypes)
saveRDS(ref, "~/projects/crc_visiumhd/rctd_reference.rds")
```

For the RCTD run, creating two files:
run_rctd.R: the analysis itself: loading the annotated spatial object and the reference, reformating the coordinates and gene counts into RCTD's format (different from Seurat format), running deconvolution, saving the result.
run_rctd.sh: the job ticket for SLURM (cluster's scheduler).

```bash
cat > ~/projects/crc_visiumhd/run_rctd.R << 'EOF'
library(Seurat); library(spacexr)                           # tools
obj <- readRDS("~/projects/crc_visiumhd/03_annotated.rds")  # spatial data
ref <- readRDS("~/projects/crc_visiumhd/rctd_reference.rds") # reference
coords <- GetTissueCoordinates(obj)[, c("x", "y")]          # bin positions
cts <- GetAssayData(obj, layer = "counts")                  # raw counts
puck <- SpatialRNA(coords, cts)                             # spatial query object
rctd <- create.RCTD(puck, ref, max_cores = 8)                # set up
rctd <- run.RCTD(rctd, doublet_mode = "doublet")            # run
saveRDS(rctd, "~/projects/crc_visiumhd/04_rctd.rds")
EOF
```

```bash
cat > ~/projects/crc_visiumhd/run_rctd.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=rctd --cpus-per-task=8 --mem=128G --time=24:00:00
#SBATCH --output=%x_%j.log
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
conda activate rspatial
export TMPDIR=~/tmp
Rscript ~/projects/crc_visiumhd/run_rctd.R
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch run_rctd.sh # submit job to queue
squeue -u weberl # confirm it's queued / running
```

```text
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          10171662       cpu     bash   weberl  R    1:45:07      1 isca062
          10196760       cpu     rctd   weberl  R       0:03      1 isca085
```

To monitor progress:

```bash
tail -20 ~/projects/crc_visiumhd/rctd_10196760.log
```

```text
Loading required package: SeuratObject
Loading required package: sp
Attaching package: 'SeuratObject'
The following objects are masked from 'package:base':
    intersect, t
Begin: process_cell_type_info
process_cell_type_info: number of cells in reference: 82237
process_cell_type_info: number of genes in reference: 18082
              B cells           Endothelial            Fibroblast
                10000                  7916                 10000
Intestinal Epithelial               Myeloid              Neuronal
                10000                 10000                  4321
        Smooth Muscle               T cells                 Tumor
                10000                 10000                 10000
```
