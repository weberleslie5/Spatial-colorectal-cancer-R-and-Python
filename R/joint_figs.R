library(Seurat); library(ggplot2)
samples <- c("P1CRC", "P2CRC", "P5CRC")
jc <- read.csv("~/projects/crc_visiumhd/joint_clusters.csv")
id <- read.csv("~/projects/crc_visiumhd/cluster_identity_joint.csv")

pal <- list(
  "Tumor II" = "#C17A97", "Tumor III" = "#ECE6A0", "Tumor IV" = "#443E82", "Tumor (shared)" = "#807F3E",
  "Intestinal Epithelial" = c("#7EC4CD", "#DDAC76", "#B5DAC8", "#9BC653", "#000000", "#DC8B7A"),
  "Fibroblast" = c("#812320", "#E48679", "#6684BD"),
  "Endothelial" = c("#C9A785", "#BBD9E0", "#74559E"),
  "Myeloid" = c("#A5594A", "#9CCB85", "#C59E9C"),
  "B cells" = c("#6080B0", "#6594C1", "#6F3182", "#AECA63"),
  "T cells" = c("#374394", "#54816F", "#E79748"),
  "Smooth Muscle" = c("#75AF7F", "#F7E774", "#C4505C", "#8F2F3B"),
  "Neuronal" = c("#6E69AA", "#BDDAA2"))
id <- id[order(id$identity), ]
id$num <- ave(as.numeric(id$cluster), id$identity, FUN = seq_along) - 1
id$label <- paste(id$identity, id$num)
id$color <- unlist(lapply(split(id, id$identity), function(d) rep_len(pal[[d$identity[1]]], nrow(d))))[order(order(id$identity))]
lab_of <- setNames(id$label, as.character(id$cluster)); col_of <- setNames(id$color, id$label)

for (s in samples) {
  obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", s, ".rds"))
  xy <- GetTissueCoordinates(obj); rownames(xy) <- paste0(s, "_", rownames(xy))
  rm(obj); gc()
  sub <- jc[jc$sample == s, ]
  xy <- xy[sub$barcode, ]; xy$label <- lab_of[as.character(sub$cluster)]
  tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp
  sf <- paste(readLines(paste0("~/data/crc_visiumhd/", s, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
  mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
  mm <- 1000 / mpp
  if (s == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
  y0 <- max(xy$y) - 0.04 * diff(range(xy$y))
  p <- ggplot(xy, aes(x, y, color = label)) + geom_point(size = 0.25, show.legend = FALSE) +
    scale_color_manual(values = col_of) + scale_y_reverse() + coord_fixed() + theme_void() +
    annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "#E8C51A") +
    ggtitle(s) + theme(plot.title = element_text(hjust = 0.5))
  ar <- diff(range(xy$y)) / diff(range(xy$x))
  ggsave(paste0("~/projects/crc_visiumhd/figures/joint_panel_", s, ".png"), p,
         width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")
}

# ---- legend only: collapse tumor identities into one numbered "Tumor" block ----
tmap <- c("Tumor III" = 0, "Tumor IV" = 1, "Tumor II" = 2, "Tumor (shared)" = 3)
leg <- id
leg$fam2 <- leg$identity; leg$num2 <- leg$num
istum <- leg$identity %in% names(tmap)
leg$fam2[istum] <- "Tumor"
leg$num2[istum] <- tmap[leg$identity[istum]]
leg <- unique(leg[, c("fam2", "num2", "color")])
leg <- leg[order(leg$fam2, leg$num2), ]
fams <- unique(leg$fam2)
colno <- setNames(rep(1:2, each = ceiling(length(fams) / 2), length.out = length(fams)), fams)
leg$lx <- colno[leg$fam2] * 3
y <- c(0, 0); leg$ly <- NA; hd <- NULL
for (f in fams) {
  k <- sum(leg$fam2 == f); cx <- colno[f]
  leg$ly[leg$fam2 == f] <- y[cx] + seq_len(k)
  hd <- rbind(hd, data.frame(f = f, lx = cx * 3, ly = y[cx] + (1 + k) / 2))
  y[cx] <- y[cx] + k + 1
}
p_leg <- ggplot() +
  geom_tile(data = leg, aes(lx + 0.9, -ly), width = 0.3, height = 0.8, fill = leg$color) +
  geom_text(data = leg, aes(lx + 1.2, -ly, label = num2), hjust = 0, size = 3) +
  geom_text(data = hd, aes(lx + 0.7, -ly, label = f), hjust = 1, size = 3.2) +
  xlim(0.5, 8) + coord_fixed(ratio = 0.35) + theme_void()
ggsave("~/projects/crc_visiumhd/figures/joint_legend.png", p_leg, width = 4, height = 8, dpi = 150, bg = "white")
