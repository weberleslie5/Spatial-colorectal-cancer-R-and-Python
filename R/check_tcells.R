library(spacexr)
rctd <- readRDS("~/projects/crc_visiumhd/04_rctd_P1CRC.rds")
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_P1CRC.csv.gz")
m <- merge(mine[mine$type == "T cells", ], authors[, c("barcode", "DeconvolutionLabel1")], by = "barcode")
print(sort(table(m$DeconvolutionLabel1), decreasing = TRUE))   # what the authors called our "T cell" bins
