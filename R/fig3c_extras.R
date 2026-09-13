library(ggplot2)
for (g in c("PIGR", "CEACAM6", "COL1A1")) {                # vertical italic row labels
  p <- ggplot() + annotate("text", 0, 0, label = g, angle = 90, fontface = "italic", size = 11) + theme_void()
  ggsave(paste0("~/projects/crc_visiumhd/figures/03c_lab_", g, ".png"), p, width = 1.2, height = 6, dpi = 150, bg = "white")
}
d <- data.frame(x = c(1, 1), y = c(1, 1), e = c(0, 8))     # invisible points -> legend only
p <- ggplot(d, aes(x, y, color = e)) + geom_point(alpha = 0) +
  scale_color_gradient(low = "grey85", high = "#C81C0F", limits = c(0, 8), name = "Normalized log\nexpression") +
  theme_void() + theme(legend.key.height = unit(2, "cm"))
ggsave("~/projects/crc_visiumhd/figures/03c_colorbar.png", p, width = 3, height = 7, dpi = 150, bg = "white")
