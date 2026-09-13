args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
out <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
write.csv(out, paste0("~/projects/crc_visiumhd/rctd_types_", sample, ".csv"), row.names = FALSE)
