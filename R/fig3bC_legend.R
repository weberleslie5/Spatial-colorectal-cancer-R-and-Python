library(ggplot2)
cols <- c("Tumor II" = "#361E50", "Tumor III" = "#423488", "Tumor IV" = "#C19C41",
          "Intestinal Epithelial" = "#582E76", "CAF" = "#4E6A83", "Fibroblast" = "#D0DFBA",
          "Endothelial" = "#752A63", "Myeloid" = "#7D9A60",
          "Plasma" = "#C1483C", "Mature B" = "#72B0D9",
          "CD4 T" = "#EDE690", "CD8 T NK" = "#CC97A8",
          "Smooth Muscle" = "#B36743", "Neuronal" = "#F6C548")
df <- data.frame(type = factor(names(cols), levels = names(cols)), y = rev(seq_along(cols)))
p <- ggplot(df, aes(0, y)) + geom_point(color = cols, size = 4.5) +
  geom_text(aes(label = type), hjust = 0, nudge_x = 0.3, size = 4) +
  xlim(-0.3, 4) + theme_void()
ggsave("~/projects/crc_visiumhd/figures/03bC_legend.png", p, width = 2.6, height = 3, dpi = 150, bg = "white")
