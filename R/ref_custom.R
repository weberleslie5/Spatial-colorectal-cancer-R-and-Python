library(Seurat); library(spacexr)
counts <- Read10X_h5("~/data/crc_visiumhd/chromium/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")
meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
meta <- meta[meta$QCFilter == "Keep", ]
lab <- meta$Level1                                          # start from the coarse classes
lab[meta$Level2 == "Plasma"] <- "Plasma"                    # split B: follicle fix
lab[meta$Level2 %in% c("Mature B", "Memory B")] <- "Mature B"
lab[meta$Level2 == "CD4 T cell"] <- "CD4 T"                 # split T
lab[meta$Level2 %in% c("CD8 T cell", "NK")] <- "CD8 T NK"
lab[meta$Level2 == "CAF"] <- "CAF"                          # split fibroblast: CAF vs normal
names(lab) <- meta$Barcode
common <- intersect(colnames(counts), meta$Barcode)
lab <- lab[common]
keep <- names(table(lab))[table(lab) >= 25]
cells <- names(lab)[lab %in% keep]
ref <- Reference(counts[, cells], as.factor(lab[cells]))
saveRDS(ref, "~/projects/crc_visiumhd/rctd_reference_custom.rds")
cat("classes:\n"); print(table(lab[cells]))
