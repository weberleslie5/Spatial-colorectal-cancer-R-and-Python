args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr); library(ggplot2)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_custom_", sample, ".rds"))
res <- rctd@results$results_df
xy <- as.data.frame(rctd@spatialRNA@coords[rownames(res), ])
xy$type <- as.character(res$first_type)
xy <- xy[res$spot_class != "reject", ]
subtype <- c(P1CRC = "Tumor II", P2CRC = "Tumor III", P5CRC = "Tumor IV")[sample]
xy$type[xy$type == "Tumor"] <- subtype
tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp

cols <- c("Tumor II" = "#361E50", "Tumor III" = "#423488", "Tumor IV" = "#C19C41",
          "Intestinal Epithelial" = "#582E76", "CAF" = "#4E6A83", "Fibroblast" = "#D0DFBA",
          "Endothelial" = "#752A63", "Myeloid" = "#7D9A60",
          "Plasma" = "#C1483C", "Mature B" = "#72B0D9",
          "CD4 T" = "#EDE690", "CD8 T NK" = "#CC97A8",
          "Smooth Muscle" = "#B36743", "Neuronal" = "#F6C548")

sf <- paste(readLines(paste0("~/data/crc_visiumhd/", sample, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
mm <- 1000 / mpp
if (sample == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
y0 <- max(xy$y) - 0.04 * diff(range(xy$y))

p <- ggplot(xy, aes(x, y, color = type)) + geom_point(size = 0.25, show.legend = FALSE) +
  scale_color_manual(values = cols) + scale_y_reverse() + coord_fixed() + theme_void() +
  annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "#E8C51A") +
  ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5, size = 28))
ar <- diff(range(xy$y)) / diff(range(xy$x))
ggsave(paste0("~/projects/crc_visiumhd/figures/03bC_", sample, ".png"), p,
       width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")
