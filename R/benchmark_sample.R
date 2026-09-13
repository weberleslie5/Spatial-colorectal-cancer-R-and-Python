args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), mine = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
authors <- authors[authors$DeconvolutionClass != "reject" & !is.na(authors$DeconvolutionLabel1), c("barcode", "DeconvolutionLabel1")]

meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
lookup <- setNames(meta$Level1, meta$Level2)               # authors' own subtype-to-class table
authors$L1 <- lookup[authors$DeconvolutionLabel1]          # collapse subtypes to 9 classes
authors$L1[is.na(authors$L1) & grepl("Tumor", authors$DeconvolutionLabel1)] <- "Tumor"          # spatial-only labels
authors$L1[is.na(authors$L1) & grepl("SM", authors$DeconvolutionLabel1)] <- "Smooth Muscle"
authors$L1[is.na(authors$L1) & grepl("Fibro", authors$DeconvolutionLabel1)] <- "Fibroblast"

m <- merge(mine, authors, by = "barcode")
m <- m[!is.na(m$L1), ]                                     # drop anything still unmapped
cat(sample, "bins compared:", nrow(m), "\n")
cat(sample, "agreement:", mean(m$mine == m$L1), "\n")
print(table(mine = m$mine, authors = m$L1))
