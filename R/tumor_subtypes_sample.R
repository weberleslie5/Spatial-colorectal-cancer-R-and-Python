args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
cl <- read.csv(paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"))
id <- read.csv(paste0("~/projects/crc_visiumhd/cluster_identity_", sample, ".csv"))
tum <- id$cluster[id$identity == "Tumor"]                  # our tumor clusters only
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
m <- merge(cl[cl$cluster %in% tum, ], authors[, c("barcode", "DeconvolutionLabel1")], by = "barcode")
m <- m[grepl("Tumor", m$DeconvolutionLabel1), ]            # keep bins the authors called a tumor subtype
print(table(ours = m$cluster, authors = m$DeconvolutionLabel1))
