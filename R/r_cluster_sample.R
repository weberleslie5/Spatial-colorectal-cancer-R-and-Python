args <- commandArgs(trailingOnly = TRUE)    # read sample name from the command line: one script for all samples
sample <- args[1]
library(Seurat); library(ggplot2)
obj <- Load10X_Spatial(paste0("~/data/crc_visiumhd/", sample), bin.size = 8)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, nCount_Spatial.008um >= 100 & nFeature_Spatial.008um >= 50 & percent.mt < 30)
obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj); obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, dims = 1:30); obj <- FindClusters(obj, resolution = 0.6)
df <- GetTissueCoordinates(obj); df$cluster <- obj$seurat_clusters
p <- ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.1) + coord_fixed() + theme_void() + guides(colour = guide_legend(override.aes = list(size = 4)))
ggsave(paste0("~/projects/crc_visiumhd/figures/02_clusters_", sample, ".png"), p, width = 8, height = 6)
saveRDS(obj, paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
