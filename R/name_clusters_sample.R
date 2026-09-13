args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
cl <- read.csv(paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"))
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
ct <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
m <- merge(cl, ct, by = "barcode")
tab <- table(m$cluster, m$type)
out <- data.frame(cluster = rownames(tab),
                identity = colnames(tab)[max.col(tab)],
                confidence = round(apply(tab, 1, max) / rowSums(tab), 2))
print(out)
write.csv(out, paste0("~/projects/crc_visiumhd/cluster_identity_", sample, ".csv"), row.names = FALSE)
