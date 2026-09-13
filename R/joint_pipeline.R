library(Seurat); library(spacexr); library(ggplot2)
samples <- c("P1CRC", "P2CRC", "P5CRC")

objs <- list(); coords <- list(); types <- list()
for (s in samples) {
  o <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", s, ".rds"))
  cts <- GetAssayData(o, layer = "counts")
  colnames(cts) <- paste0(s, "_", colnames(cts))
  xy <- GetTissueCoordinates(o); rownames(xy) <- paste0(s, "_", rownames(xy))
  coords[[s]] <- xy
  objs[[s]] <- CreateSeuratObject(cts)
  r <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", s, ".rds"))
  res <- r@results$results_df
  types[[s]] <- data.frame(barcode = paste0(s, "_", rownames(res)), type = as.character(res$first_type))[res$spot_class != "reject", ]
  rm(o, cts, r, res); gc()
}
big <- merge(objs[[1]], objs[2:3]); rm(objs); gc()
big <- JoinLayers(big)
big$sample <- sub("_.*", "", colnames(big))

big <- NormalizeData(big)
big <- FindVariableFeatures(big, nfeatures = 2000)
big <- ScaleData(big)
big <- RunPCA(big, npcs = 30)
big <- FindNeighbors(big, dims = 1:30)
big <- FindClusters(big, resolution = 0.6)
jc <- data.frame(barcode = colnames(big), sample = big$sample, cluster = as.character(big$seurat_clusters))
write.csv(jc, "~/projects/crc_visiumhd/joint_clusters.csv", row.names = FALSE)   # checkpoint
rm(big); gc()

ct <- do.call(rbind, types)
m <- merge(jc, ct, by = "barcode")
tab <- table(m$cluster, m$type)
id <- data.frame(cluster = rownames(tab), identity = colnames(tab)[max.col(tab)],
                 confidence = round(apply(tab, 1, max) / rowSums(tab), 2))
comp <- table(jc$cluster, jc$sample)
subtype <- c(P1CRC = "Tumor II", P2CRC = "Tumor III", P5CRC = "Tumor IV")
for (i in which(id$identity == "Tumor")) {
  dom <- colnames(comp)[which.max(comp[id$cluster[i], ])]
  frac <- max(comp[id$cluster[i], ]) / sum(comp[id$cluster[i], ])
  id$identity[i] <- if (frac > 0.6) subtype[dom] else "Tumor (shared)"
}
write.csv(id, "~/projects/crc_visiumhd/cluster_identity_joint.csv", row.names = FALSE)

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
id$num <- ave(id$cluster, id$identity, FUN = seq_along) - 1
id$label <- paste(id$identity, id$num)
id$color <- unlist(lapply(split(id, id$identity), function(d) rep_len(pal[[d$identity[1]]], nrow(d))))[order(order(id$identity))]
lab_of <- setNames(id$label, id$cluster); col_of <- setNames(id$color, id$label)

for (s in samples) {
  xy <- coords[[s]]
  b <- jc$barcode[jc$sample == s]
  xy <- xy[b, ]; xy$label <- lab_of[jc$cluster[match(b, jc$barcode)]]
  tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp                    # final orientation
  sf <- paste(readLines(paste0("~/data/crc_visiumhd/", s, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
  mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
  mm <- 1000 / mpp
  if (s == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
  y0 <- max(xy$y) - 0.04 * diff(range(xy$y))                # raised scale bar
  p <- ggplot(xy, aes(x, y, color = label)) + geom_point(size = 0.25, show.legend = FALSE) +
    scale_color_manual(values = col_of) + scale_y_reverse() + coord_fixed() + theme_void() +
    annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "#E8C51A") +
    ggtitle(s) + theme(plot.title = element_text(hjust = 0.5))
  ar <- diff(range(xy$y)) / diff(range(xy$x))               # true tissue proportions
  ggsave(paste0("~/projects/crc_visiumhd/figures/joint_panel_", s, ".png"), p,
         width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")
}

fams <- unique(id$identity)                                 # compact two-column legend, numbers on all
colno <- setNames(rep(1:2, each = ceiling(length(fams) / 2), length.out = length(fams)), fams)
id$lx <- colno[id$identity] * 3
y <- c(0, 0); id$ly <- NA; hd <- NULL
for (f in fams) {
  k <- sum(id$identity == f); cx <- colno[f]
  id$ly[id$identity == f] <- y[cx] + seq_len(k)
  hd <- rbind(hd, data.frame(f = f, lx = cx * 3, ly = y[cx] + (1 + k) / 2))
  y[cx] <- y[cx] + k + 1
}
p_leg <- ggplot() +
  geom_tile(data = id, aes(lx + 0.9, -ly), width = 0.3, height = 0.8, fill = id$color) +
  geom_text(data = id, aes(lx + 1.2, -ly, label = num), hjust = 0, size = 3) +
  geom_text(data = hd, aes(lx + 0.7, -ly, label = f), hjust = 1, size = 3.2) +
  xlim(0.5, 8) + coord_fixed(ratio = 0.35) + theme_void()
ggsave("~/projects/crc_visiumhd/figures/joint_legend.png", p_leg, width = 4, height = 8, dpi = 150, bg = "white")
