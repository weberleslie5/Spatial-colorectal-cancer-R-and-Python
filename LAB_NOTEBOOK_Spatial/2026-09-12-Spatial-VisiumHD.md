# 2026-09-12 — Spatial Visium HD

Figure 3b, P2 and P5:

```bash
cd ~/projects/crc_visiumhd
sbatch fig3b_custom_sample.sh P2CRC
sbatch fig3b_custom_sample.sh P5CRC
sbatch benchmark_custom_sample.sh P2CRC
sbatch benchmark_custom_sample.sh P5CRC
```

Edits larger titles:

```bash
sed -i 's/theme(plot.title = element_text(hjust = 0.5))/theme(plot.title = element_text(hjust = 0.5, size = 28))/' ~/projects/crc_visiumhd/fig3b_custom_sample.R
```

Regenerate panels:

```bash
cd ~/projects/crc_visiumhd
sbatch fig3b_custom_sample.sh P1CRC
sbatch fig3b_custom_sample.sh P2CRC
sbatch fig3b_custom_sample.sh P5CRC
```

Building the 14 class legend for 3b:

```bash
cat ~/projects/crc_visiumhd/fig3b_custom_sample.R
```

```r
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
  ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5))
ar <- diff(range(xy$y)) / diff(range(xy$x))
ggsave(paste0("~/projects/crc_visiumhd/figures/03bC_", sample, ".png"), p,
       width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")
```

Legend script:

```bash
cat > ~/projects/crc_visiumhd/fig3bC_legend.R << 'EOF'
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
EOF
```

Stitch script:

```bash
cat > ~/projects/crc_visiumhd/fig3bC_stitch.py << 'EOF'
import os
from PIL import Image, ImageChops
home = os.path.expanduser("~")
fig = f"{home}/projects/crc_visiumhd/figures"

def trim(im):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    return im.crop(ImageChops.difference(im, bg).getbbox())

imgs = [trim(Image.open(f"{fig}/03bC_{s}.png").convert("RGB")) for s in ["P1CRC", "P2CRC", "P5CRC"]]
imgs.append(trim(Image.open(f"{fig}/03bC_legend.png").convert("RGB")))
h = max(im.height for im in imgs[:3])
imgs = [im.resize((round(im.width * h / im.height), h)) for im in imgs]
pad, margin = 40, 60    # pad between panels, margin all around
row = Image.new("RGB", (sum(im.width for im in imgs) + pad * 3 + margin * 2, h + margin * 2), (255, 255, 255))
x = margin
for im in imgs:
    row.paste(im, (x, margin + (h - im.height) // 2)); x += im.width + pad
row.save(f"{fig}/fig3bC_row.png")
print("saved fig3bC_row.png")
EOF
```

Job ticket:

```bash
cat > ~/projects/crc_visiumhd/fig3bC_stitch.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=stitchC --cpus-per-task=2 --mem=16G --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/fig3bC_legend.R
conda activate pyspatial
python ~/projects/crc_visiumhd/fig3bC_stitch.py
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch fig3bC_stitch.sh
```
