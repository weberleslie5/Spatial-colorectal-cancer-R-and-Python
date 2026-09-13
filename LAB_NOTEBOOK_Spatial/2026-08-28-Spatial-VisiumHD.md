# 2026-08-28 — Spatial Visium HD

```bash
squeue -u weberl
```

```text
  JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

```bash
ls -lh ~/projects/crc_visiumhd/04_rctd.rds
```

```text
-rw-rw-r-- 1 weberl weberl 756M Aug 28 00:16 /home/weberl/projects/crc_visiumhd/04_rctd.rds
```

```bash
tail -20 ~/projects/crc_visiumhd/rctd_10196760.log
```

```text
[1] "gather_results: finished 412000"
[1] "gather_results: finished 413000"
[1] "gather_results: finished 414000"
[1] "gather_results: finished 415000"
[1] "gather_results: finished 416000"
[1] "gather_results: finished 417000"
[1] "gather_results: finished 418000"
[1] "gather_results: finished 419000"
[1] "gather_results: finished 420000"
[1] "gather_results: finished 421000"
[1] "gather_results: finished 422000"
[1] "gather_results: finished 423000"
[1] "gather_results: finished 424000"
[1] "gather_results: finished 425000"
[1] "gather_results: finished 426000"
[1] "gather_results: finished 427000"
[1] "gather_results: finished 428000"
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 7.6 GiB
```

Job completed successfully.

## Verifying the deconvolution output

```bash
srun --pty --cpus-per-task=4 --mem=64GB --time=6:00:00 bash
conda activate rspatial
R
```

```r
library(spacexr)
rctd <- readRDS("~/projects/crc_visiumhd/04_rctd.rds") # loading the results from last night
res <- rctd@results$results_df # from the rctd object, open the results compartment and take the table named results_df, transfered to res: this is a table with one row per 8 µm bin, with RCTD decision: its class (one or more cell), identity of the cells.

dim(res) # prints table's dimensions
```

```text
[1] 428405      9
```

```r
head(res)
```

```text
                     spot_class            first_type second_type first_class
s_008um_00301_00321-1    singlet           Endothelial       Tumor       FALSE
s_008um_00526_00291-1    singlet                 Tumor     B cells       FALSE
s_008um_00128_00278-1    singlet           Endothelial  Fibroblast       FALSE
s_008um_00052_00559-1    singlet Intestinal Epithelial     T cells       FALSE
s_008um_00504_00410-1    singlet                 Tumor     B cells       FALSE
s_008um_00447_00253-1    singlet                 Tumor     B cells       FALSE
                      second_class min_score singlet_score conv_all
s_008um_00301_00321-1        FALSE  227.7695      228.5688     TRUE
s_008um_00526_00291-1        FALSE  765.8615      765.6278     TRUE
s_008um_00128_00278-1        FALSE  331.2901      339.9698     TRUE
s_008um_00052_00559-1        FALSE  336.4178      336.5833     TRUE
s_008um_00504_00410-1        FALSE  653.9680      653.7539     TRUE
s_008um_00447_00253-1        FALSE  504.3937      504.2981     TRUE
                      conv_doublet
s_008um_00301_00321-1         TRUE
s_008um_00526_00291-1         TRUE
s_008um_00128_00278-1         TRUE
s_008um_00052_00559-1         TRUE
s_008um_00504_00410-1         TRUE
s_008um_00447_00253-1         TRUE
```

```r
table(res$spot_class) # bin's classification summary
```

```text
 reject           singlet   doublet_certain doublet_uncertain
             5041            368726             45239              9399
```

```r
table(res$first_type) # primary cell type counts
```

```text
   B cells           Endothelial            Fibroblast
                18317                 24869                 73132
Intestinal Epithelial               Myeloid              Neuronal
                40708                 27396                   987
        Smooth Muscle               T cells                 Tumor
                15300                  8414                219282
```

Creating the deconvolved cell-type map:

```r
library(ggplot2)
coords <- as.data.frame(rctd@spatialRNA@coords) # bin x, y stored into the rctd object (coordinates are not included in the result table > we need to add them)
df <- coords[rownames(res), ] # aligns positions to the results rows
df$type <- res$first_type # attach cell-type labels
df <- df[res$spot_class != "reject", ] # keeps only validated bins

p6 <- ggplot(df, aes(x, y, colour = type)) + geom_point(size = 0.1) + coord_fixed() + theme_void() + guides(colour = guide_legend(override.aes = list(size = 4))) # cell-type map

ggsave("~/projects/crc_visiumhd/figures/05_celltypes.png", p6, width = 9, height = 7)
```

## Comparing to the author's figure

```r
authors <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_P1CRC.csv.gz")
dim(authors)
```

```text
[1] 702244      4
```

```r
head(authors)
```

```text
  barcode DeconvolutionClass DeconvolutionLabel1
1 s_008um_00000_00000-1               <NA>                <NA>
2 s_008um_00000_00001-1               <NA>                <NA>
3 s_008um_00000_00002-1               <NA>                <NA>
4 s_008um_00000_00003-1               <NA>                <NA>
5 s_008um_00000_00004-1               <NA>                <NA>
6 s_008um_00000_00005-1               <NA>                <NA>
  DeconvolutionLabel2
1                <NA>
2                <NA>
3                <NA>
4                <NA>
5                <NA>
6                <NA>
```

```r
mine <- data.frame(barcode = rownames(res), mine = res$first_type, class = res$spot_class)
merged <- merge(mine, authors[, c("barcode", "DeconvolutionLabel1")], by = "barcode") # join on barcode
merged <- merged[merged$class != "reject" & !is.na(merged$DeconvolutionLabel1), ] # bins resolved by both
nrow(merged) # bins compared

mean(merged$mine == merged$DeconvolutionLabel1) # agreement rate
```

```text
[1] 0.009014506
```

```r
table(mine = merged$mine, authors = merged$DeconvolutionLabel1) # confusion matrix
```

```text
      authors
mine                    Adipocyte   CAF CD4 T cell CD8 T cell cDC I Endothelial
  B cells                       9  2346        368        101    39         468
  Endothelial                  22  1878        430        117    46         590
  Fibroblast                   76  3190       2008        524   190        1782
  Intestinal Epithelial        22  6389        447        119    85        1176
  Myeloid                      26  1969        483        121    39         639
  Neuronal                      0    25          8          4     0          17
  Smooth Muscle                 2   420         43         22    10         199
  T cells                       2   559        157         38    12         194
  Tumor                       184 17629       4262       1304   470        5593
                       authors
mine                    Enteric Glial Enterocyte Epithelial Fibroblast Goblet
  B cells                          12        136          3        192    382
  Endothelial                      30        334          4        306   1573
  Fibroblast                       95       2716         49       1508  10597
  Intestinal Epithelial            19         20          0        203     16
  Myeloid                          40        355         11        361   1715
  Neuronal                          0         19          0          2     12
  Smooth Muscle                     0        374          3         21    210
  T cells                          20        148          0         90    428
  Tumor                           420       2342         97       3202  17422
                       authors
mine                    Lymphatic Endothelial Macrophage  Mast Mature B
  B cells                                  56        309    31      266
  Endothelial                              79        420    39      363
  Fibroblast                              213        974   146     1745
  Intestinal Epithelial                   135        808    67       90
  Myeloid                                  83        437    41      413
  Neuronal                                  0         11     0        9
  Smooth Muscle                             6         79    10       23
  T cells                                  23        119     7      149
  Tumor                                   984       3354   389     2825
                       authors
mine                    Memory B mRegDC Myofibroblast Neuroendocrine Neutrophil
  B cells                      9     34           118              2         86
  Endothelial                  4     34           172             10        101
  Fibroblast                  45    154           974             73        246
  Intestinal Epithelial        5     77           351              4        225
  Myeloid                     14     43           178              7        130
  Neuronal                     0      5            12              0          6
  Smooth Muscle                1     12            93              1         46
  T cells                      1     19            64              0         29
  Tumor                      127    451          1488            126       1070
                       authors
mine                       NK   pDC Pericytes Plasma Proliferating Fibroblast
  B cells                   0     8       169    105                      219
  Endothelial               0     4       229    186                      161
  Fibroblast                5    25       595    671                      319
  Intestinal Epithelial     0    12       413    193                      629
  Myeloid                   1     4       214    220                      222
  Neuronal                  0     0         6      2                        5
  Smooth Muscle             0     2        67     48                       48
  T cells                   0     2        81     54                       68
  Tumor                     8    84      2071   1878                     1583
                       authors
mine                    Proliferating Immune II Proliferating Macrophages
  B cells                                    39                       180
  Endothelial                                59                       183
  Fibroblast                                232                       522
  Intestinal Epithelial                      59                       441
  Myeloid                                    70                       202
  Neuronal                                    2                        12
  Smooth Muscle                              17                        82
  T cells                                    18                        57
  Tumor                                     520                      1203
                       authors
mine                    SM Stress Response Smooth Muscle  Tuft Tumor I Tumor II
  B cells                                0             8     4       6     5743
  Endothelial                            5             8     5       7     5943
  Fibroblast                             7            24    32      29    11204
  Intestinal Epithelial                  8             2     6       8    19980
  Myeloid                                3            20     6      10     6315
  Neuronal                               0             0     0       0      302
  Smooth Muscle                          1             0     3       1     3143
  T cells                                0             0     1       1     1812
  Tumor                                 27           181    71      42    34843
                       authors
mine                    Tumor III Tumor IV Tumor V Unknown III (SM)
  B cells                      24        0      78                3
  Endothelial                  32        5      78               21
  Fibroblast                   45       15     100              163
  Intestinal Epithelial        70        3     172                8
  Myeloid                      41        4      72               37
  Neuronal                      0        0       0                0
  Smooth Muscle                 6        2      30                1
  T cells                       7        0      32                3
  Tumor                       456       31     683              202
                       authors
mine                    Vascular Fibroblast   vSM
  B cells                                 0    20
  Endothelial                             0   103
  Fibroblast                              0   672
  Intestinal Epithelial                   0    29
  Myeloid                                 0   167
  Neuronal                                0     2
  Smooth Muscle                           0    12
  T cells                                 0    40
  Tumor                                   5  1252
```

Authors used class II annotation, while I used class I annotations. Re-run RCTD with class II.

## A. Build the level 2 reference

Inspecting the authors's deconvolution script:

```bash
ls ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/
```

```text
ext  Figures  LICENSE.md  MetaData  Methods  README.md  Team
```

```bash
ls ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/Methods
```

```text
AuxFunctions.R  Deconvolution.R  environment_nucleisegmentation.yml  FlexSingleCell.R  NucleiSegmentation.py
```

```bash
cat ~/data/crc_visiumhd/HumanColonCancer_VisiumHD/Methods/Deconvolution.R
```

```r
#### Visium HD Manuscript
## Run SpaceXR (Deconvolution)
# This script can be used as a skeleton to run deconvolution for any sample.
## load libraries
library(arrow)
library(Seurat)
library(spacexr) # Should be the modified version [see PR: ADD LINK TO PR]
FlexOutPath <- "~/AggrOutput/outs" # Path to cellranger aggr output folder
ColonFlex.data <- Read10X_h5(paste0(FlexOutPath,"filtered_feature_bc_matrix.h5"))
## Load Reference Data
FlexRef<-Read10X_h5(paste0(FlexOutPath,"filtered_feature_bc_matrix.h5"))
MetaData<-readRDS('~/Outputs/Flex/FlexSeuratV5_MetaData.rds') # See FlexSingleCell.R if not generated.
# spacexr restriction, Clusters with > 25 cells
KpIdents<-names(which(table(MetaData$Level2)>25))
MetaData<-MetaData[MetaData$Level2%in%KpIdents,]
FlexRef<-FlexRef[,MetaData$Barcode]
## Fix cell type labels as spacexr doesn't allow special characters (i.e. spaces)
CTRef<-MetaData$Level2
CTRef<-gsub("/","_",CTRef)
CTRef<-as.factor(CTRef)
names(CTRef)<-MetaData$Barcode
## Build reference object
reference <- Reference(FlexRef[,names(CTRef)], CTRef , colSums(FlexRef))
# Deconvolve HD Data
counts<-Read10X_h5("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
coords<-read_parquet("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet",as_data_frame = TRUE)
rownames(coords)<-coords$barcode
coords<-coords[colnames(counts),]
coords<-coords[,3:4]
nUMI <- colSums(counts)
puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts)
myRCTD <- create.RCTD(puck, reference, max_cores = 12)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
```

```r
library(Seurat); library(spacexr)
counts <- Read10X_h5("~/data/crc_visiumhd/chromium/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")
meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
meta <- meta[meta$QCFilter == "Keep", ]
keep <- names(which(table(meta$Level2) > 25))
meta <- meta[meta$Level2 %in% keep, ]
common <- intersect(colnames(counts), meta$Barcode)
counts <- counts[, common]
rownames(meta) <- meta$Barcode
ct <-  gsub("/", "_", meta[common, "Level2"])
ct <- factor(ct); names(ct) <- common
ref <- Reference(counts, ct, colSums(counts))
saveRDS(ref, "~/projects/crc_visiumhd/rctd_reference_L2.rds")
q()
```

```bash
sed 's/rctd_reference.rds/rctd_reference_L2.rds/; s/04_rctd.rds/04_rctd_L2.rds/' ~/projects/crc_visiumhd/run_rctd.R > ~/projects/crc_visiumhd/run_rctd_L2.R
sed 's/run_rctd.R/run_rctd_L2.R/' ~/projects/crc_visiumhd/run_rctd.sh > ~/projects/crc_visiumhd/run_rctd_L2.sh

cd ~/projects/crc_visiumhd
sbatch run_rctd_L2.sh
squeue -u weberl
```
