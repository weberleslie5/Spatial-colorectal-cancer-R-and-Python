args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(clusterProfiler); library(org.Hs.eg.db); library(dplyr)
m <- read.csv(paste0("~/projects/crc_visiumhd/markers_", sample, ".csv"))
out <- list()
for (cl in unique(m$cluster)) {
    g <- m$gene[m$cluster == cl]
    e <- enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
    df <- as.data.frame(e)
    if (nrow(df) > 0) out[[as.character(cl)]] <- data.frame(cluster=cl, head(df[, c("Description", "p.adjust", "Count")], 10))
}
write.csv(bind_rows(out), paste0("~/projects/crc_visiumhd/pathways_", sample, ".csv"), row.names = FALSE)
