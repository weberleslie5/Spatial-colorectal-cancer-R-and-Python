args <- commandArgs(trailingOnly = TRUE); sample <- args[1] 
library(Seurat)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
df <- data.frame(barcode = colnames(obj), cluster = obj$seurat_clusters)
write.csv(df, paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"), row.names = FALSE)
