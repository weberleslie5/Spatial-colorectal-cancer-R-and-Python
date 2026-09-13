args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_custom_", sample, ".rds"))
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), mine = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
authors <- authors[authors$DeconvolutionClass != "reject" & !is.na(authors$DeconvolutionLabel1), c("barcode", "DeconvolutionLabel1")]

meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
lab <- meta$Level1
lab[meta$Level2 == "Plasma"] <- "Plasma"
lab[meta$Level2 %in% c("Mature B", "Memory B")] <- "Mature B"
lab[meta$Level2 == "CD4 T cell"] <- "CD4 T"
lab[meta$Level2 %in% c("CD8 T cell", "NK")] <- "CD8 T NK"
lab[meta$Level2 == "CAF"] <- "CAF"
tab <- table(meta$Level2, lab)
lookup <- setNames(colnames(tab)[max.col(tab)], rownames(tab))
authors$cls <- lookup[authors$DeconvolutionLabel1]
authors$cls[is.na(authors$cls) & grepl("Tumor", authors$DeconvolutionLabel1)] <- "Tumor"
authors$cls[is.na(authors$cls) & grepl("SM", authors$DeconvolutionLabel1)] <- "Smooth Muscle"
authors$cls[is.na(authors$cls) & grepl("Fibro", authors$DeconvolutionLabel1)] <- "Fibroblast"

m <- merge(mine, authors, by = "barcode"); m <- m[!is.na(m$cls), ]
cat(sample, "bins compared:", nrow(m), "\n")
cat(sample, "agreement:", mean(m$mine == m$cls), "\n")
print(table(mine = m$mine, authors = m$cls))
