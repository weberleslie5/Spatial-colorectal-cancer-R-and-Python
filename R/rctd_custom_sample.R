args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat); library(spacexr)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
ref <- readRDS("~/projects/crc_visiumhd/rctd_reference_custom.rds")
coords <- GetTissueCoordinates(obj)[, c("x", "y")]
cts <- GetAssayData(obj, layer = "counts")
puck <- SpatialRNA(coords, cts)
rctd <- create.RCTD(puck, ref, max_cores = 8)
rctd <- run.RCTD(rctd, doublet_mode = "doublet")
saveRDS(rctd, paste0("~/projects/crc_visiumhd/04_rctd_custom_", sample, ".rds"))
