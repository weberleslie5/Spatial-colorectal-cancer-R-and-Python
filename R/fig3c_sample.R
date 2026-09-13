args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat); library(ggplot2)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
xy <- GetTissueCoordinates(obj)
tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp

sf <- paste(readLines(paste0("~/data/crc_visiumhd/", sample, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
mm <- 1000 / mpp                                            # 1 mm main-map bar
um80 <- 80 / mpp; ts <- 8 / mpp                             # 80-um inset bar; 8-um tile
if (sample == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
y0 <- max(xy$y) - 0.04 * diff(range(xy$y))

for (g in c("PIGR", "CEACAM6", "COL1A1")) {
  xy$expr <- FetchData(obj, g)[, 1]
  cw <- 250 / mpp
  cellsum <- aggregate(expr ~ I(round(x / cw)) + I(round(y / cw)), data = xy, FUN = sum)
  names(cellsum) <- c("cx", "cy", "s")
  best <- cellsum[which.max(cellsum$s), ]
  ctr <- c(best$cx * cw, best$cy * cw)
  half <- 200 / mpp
  z <- xy[abs(xy$x - ctr[1]) < half & abs(xy$y - ctr[2]) < half, ]

  p <- ggplot(xy[order(xy$expr), ], aes(x, y, color = expr)) + geom_point(size = 0.25) +
    scale_color_gradient(low = "grey85", high = "#C81C0F", limits = c(0, 8), guide = "none") +
    scale_y_reverse() + coord_fixed() + theme_void() +
    annotate("rect", xmin = ctr[1] - half, xmax = ctr[1] + half, ymin = ctr[2] - half, ymax = ctr[2] + half,
             fill = NA, color = "black", linewidth = 0.6) +
    annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "black")
  ar <- diff(range(xy$y)) / diff(range(xy$x))
  ggsave(paste0("~/projects/crc_visiumhd/figures/03c_", g, "_", sample, ".png"), p,
         width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")

  bx <- if (sample == "P1CRC") ctr[1] + half * 0.9 - um80 else ctr[1] - half * 0.9
  pz <- ggplot(z, aes(x, y, fill = expr)) + geom_tile(width = ts, height = ts) +
    scale_fill_gradient(low = "grey85", high = "#C81C0F", limits = c(0, 8), guide = "none") +
    scale_y_reverse() + coord_fixed() + theme_void() +
    annotate("segment", x = bx, xend = bx + um80,
             y = ctr[2] + half * 0.92, yend = ctr[2] + half * 0.92, linewidth = 2.5, color = "#2B5BB5")
  ggsave(paste0("~/projects/crc_visiumhd/figures/03c_zoom_", g, "_", sample, ".png"), pz,
         width = 5, height = 5, dpi = 150, bg = "white")
}
