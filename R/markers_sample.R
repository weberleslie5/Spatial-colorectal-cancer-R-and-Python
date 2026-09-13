args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat); library(dplyr)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
m <- FindAllMarkers(obj, only.pos = TRUE,
                    max.cells.per.ident = 2000,
                    logfc.threshold = 0.25, min.pct = 0.10)
sig <- m %>% filter(p_val_adj < 0.05)
write.csv(sig, paste0("~/projects/crc_visiumhd/markers_", sample, ".csv"), row.names = FALSE)
