# 2026-09-11 — Spatial Visium HD

Assigning cell identity to clusters

```bash
cd ~/projects/crc_visiumhd
tr -d '"' < markers_P1CRC.csv | awk -F, 'NR>1 {a[$6]=a[$6]" "$7} END {for (c in a) print c":"a[c]}' | sort -n | cut -c1-400
```

```text
0: MUC17 PERP PRAP1 HNF4A SDC4 SLC7A1 ALDOB NELFCD PRELID3B CEACAM6
1: COL1A1 COL1A2 MMP2 LUM COL6A3 COL3A1 AEBP1 FBLN1 C3 DCN
2: MUC5B SPINK4 LEFTY1 MUC2 WFDC2 FCGBP PLA2G2A AGR2 PIGR MUC12
3: REG1B REG1A DMBT1 OLFM4 SPINK1 FGGY SLC7A1 ANPEP MUC17 PRAP1
4: PECAM1 AQP1 ENG COL4A1 COL4A2 IGFBP7 A2M VIM ZFP36 SPARC
5: JCHAIN IGKC IGHG1 A2M FOS
6: CXCR4 TRBC2 TRAC TXNIP CD74 PNRC1 B2M PFN1 C3 ACTB
7: C1QC LYZ IFI30 APOE SRGN CD74 PSAP CTSB SAT1 CTSZ
8: F3 CXCL14 VCAN COL5A1 IGFBP7 COL6A1 LGALS1 COL3A1 COL6A2 FTL
9: ATP5F1E PERP POLR1D TSPAN8 CAMK2N1 MISP MT-ND1 FTL UBA52 MT-ND2
10: CEACAM7 GUCA2A CA1 FABP1 SLC26A2 CA2 CKB TSPAN1 SLC26A3 KRT20
11: DES ACTG2 MYH11 MYLK TAGLN MYL9 TPM2 FLNA CSRP1 ACTA2
12: CCL21 AKAP12 VIM ZFP36 IGFBP7 LGALS1 DUSP1 EGR1 FOS JUNB
```

```bash
awk -F, 'NR>1 {print $1": "$2}' pathways_P1CRC.csv
```

```text
1: "extracellular matrix organization"
1: "extracellular structure organization"
1: "external encapsulating structure organization"
1: "collagen fibril organization"
1: "cellular response to amino acid stimulus"
1: "cellular response to acid chemical"
1: "response to amino acid"
1: "response to acid chemical"
1: "transforming growth factor beta production"
1: "bone trabecula formation"
2: "mucus secretion"
2: "secretion by tissue"
2: "body fluid secretion"
3: "antimicrobial humoral immune response mediated by antimicrobial peptide"
3: "antimicrobial humoral response"
3: "humoral immune response"
4: "collagen-activated tyrosine kinase receptor signaling pathway"
4: "regulation of angiogenesis"
4: "regulation of vasculature development"
4: "collagen-activated signaling pathway"
4: "response to glucocorticoid"
4: "response to corticosteroid"
4: "glomerulus vasculature development"
4: "renal system vasculature development"
4: "kidney vasculature development"
4: "branching involved in blood vessel morphogenesis"
5: "humoral immune response"
5: "complement activation"
5: "antibacterial humoral response"
5: "B cell receptor signaling pathway"
5: "antimicrobial humoral response"
5: "immunoglobulin mediated immune response"
5: "B cell mediated immunity"
5: "antigen receptor-mediated signaling pathway"
5: "antibody-dependent cellular cytotoxicity"
5: "response to metal ion starvation"
6: "antigen processing and presentation of exogenous peptide antigen via MHC class II"
6: "antigen processing and presentation of endogenous antigen"
6: "antigen processing and presentation of peptide antigen via MHC class II"
6: "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II"
6: "antigen processing and presentation of exogenous peptide antigen"
6: "leukocyte chemotaxis"
6: "antigen processing and presentation of exogenous antigen"
6: "positive regulation of T cell activation"
6: "positive regulation of immune effector process"
6: "positive regulation of receptor-mediated endocytosis"
7: "antigen processing and presentation of exogenous peptide antigen via MHC class II"
7: "antigen processing and presentation of peptide antigen via MHC class II"
7: "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II"
7: "antigen processing and presentation of exogenous peptide antigen"
7: "antigen processing and presentation of exogenous antigen"
7: "viral life cycle"
8: "collagen fibril organization"
8: "basement membrane organization"
8: "endodermal cell differentiation"
8: "extracellular matrix assembly"
8: "endoderm formation"
8: "connective tissue development"
8: "glial cell migration"
8: "skin development"
8: "extracellular matrix organization"
8: "extracellular structure organization"
9: "proton motive force-driven mitochondrial ATP synthesis"
9: "proton motive force-driven ATP synthesis"
9: "ATP biosynthetic process"
9: "purine ribonucleoside triphosphate biosynthetic process"
9: "purine nucleoside triphosphate biosynthetic process"
9: "ribonucleoside triphosphate biosynthetic process"
9: "nucleoside triphosphate biosynthetic process"
9: "oxidative phosphorylation"
9: "purine ribonucleotide biosynthetic process"
9: "ribonucleotide biosynthetic process"
10: "oxalate transport"
10: "sulfate transmembrane transport"
10: "chloride transport"
10: "inorganic anion transport"
10: "monoatomic anion transport"
10: "sulfur compound transport"
10: "carboxylic acid transport"
10: "organic acid transport"
10: "regulation of intracellular pH"
10: "regulation of cellular pH"
11: "muscle contraction"
11: "muscle system process"
11: "myofibril assembly"
11: "platelet aggregation"
11: "homotypic cell-cell adhesion"
11: "smooth muscle contraction"
11: "cellular component assembly involved in morphogenesis"
11: "cellular anatomical entity morphogenesis"
11: "platelet activation"
11: "membraneless organelle assembly"
12: "response to steroid hormone"
12: "response to glucocorticoid"
12: "response to corticosteroid"
12: "response to ketone"
12: "response to lipopolysaccharide"
12: "cellular response to alcohol"
12: "cellular response to ketone"
12: "response to molecule of bacterial origin"
12: "cellular response to prostaglandin E stimulus"
12: "response to calcium ion"
```

Using the cell type from RCTD, decided by the single-cell reference to label the clusters:

```bash
cat > ~/projects/crc_visiumhd/name_clusters_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
cl <- read.csv(paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"))
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
ct <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
m <- merge(cl, ct, by = "barcode")
tab <- table(m$cluster, m$type)
out <- data.frame(cluster = rownames(tab),
                identity = colnames(tab)[max.col(tab)],
                confidence = round(apply(tab, 1, max) / rowSums(tab), 2))
print(out)
write.csv(out, paste0("~/projects/crc_visiumhd/cluster_identity_", sample, ".csv"), row.names = FALSE)
EOF
```

Ticket:

```bash
sed "s/benchmark_sample.R/name_clusters_sample.R/" ~/projects/crc_visiumhd/benchmark_sample.sh > ~/projects/crc_visiumhd/name_clusters_sample.sh
cd ~/projects/crc_visiumhd
sbatch name_clusters_sample.sh P1CRC
sbatch name_clusters_sample.sh P2CRC
sbatch name_clusters_sample.sh P5CRC
```

```bash
cat ~/projects/crc_visiumhd/cluster_identity_P1CRC.csv
```

```text
"cluster","identity","confidence"
"0","Tumor",1
"1","Fibroblast",0.93
"2","Intestinal Epithelial",0.98
"3","Tumor",1
"4","Endothelial",0.73
"5","B cells",0.31
"6","T cells",0.71
"7","Myeloid",0.65
"8","Tumor",0.59
"9","Tumor",0.99
"10","Intestinal Epithelial",1
"11","Smooth Muscle",0.68
"12","Endothelial",0.57
```

Confidence: ex: 57% of bins are attributed to Endothelial for cluster 12.

```bash
cat ~/projects/crc_visiumhd/cluster_identity_P2CRC.csv
```

```text
"cluster","identity","confidence"
"0","Tumor",1
"1","Tumor",1
"2","Myeloid",0.4
"3","Fibroblast",0.98
"4","Endothelial",0.66
"5","B cells",0.63
"6","Intestinal Epithelial",0.98
"7","Tumor",1
"8","Tumor",1
"9","Smooth Muscle",0.76
"10","Tumor",1
"11","Intestinal Epithelial",1
"12","Fibroblast",0.92
"13","T cells",0.59
"14","Intestinal Epithelial",0.75
"15","Tumor",0.75
"16","Myeloid",0.56
"17","Tumor",1
```

```bash
cat ~/projects/crc_visiumhd/cluster_identity_P5CRC.csv
```

```text
"cluster","identity","confidence"
"0","Tumor",0.98
"1","Smooth Muscle",0.95
"2","Fibroblast",0.97
"3","Intestinal Epithelial",1
"4","Endothelial",0.73
"5","T cells",0.46
"6","B cells",0.61
"7","Intestinal Epithelial",0.99
"8","Myeloid",0.62
"9","Tumor",0.97
"10","Intestinal Epithelial",0.64
"11","Tumor",0.51
"12","Tumor",1
"13","Intestinal Epithelial",0.99
"14","Fibroblast",0.36
"15","Neuronal",0.63
"16","Endothelial",0.73
"17","Tumor",0.39
"18","Smooth Muscle",0.62
"19","Neuronal",0.94
```

Which tumor subtype is each of the tumor cluster identified here:

```bash
cat > ~/projects/crc_visiumhd/tumor_subtypes_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
cl <- read.csv(paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"))
id <- read.csv(paste0("~/projects/crc_visiumhd/cluster_identity_", sample, ".csv"))
tum <- id$cluster[id$identity == "Tumor"]                  # our tumor clusters only
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
m <- merge(cl[cl$cluster %in% tum, ], authors[, c("barcode", "DeconvolutionLabel1")], by = "barcode")
m <- m[grepl("Tumor", m$DeconvolutionLabel1), ]            # keep bins the authors called a tumor subtype
print(table(ours = m$cluster, authors = m$DeconvolutionLabel1))
EOF
sed "s/name_clusters_sample.R/tumor_subtypes_sample.R/" ~/projects/crc_visiumhd/name_clusters_sample.sh > ~/projects/crc_visiumhd/tumor_subtypes_sample.sh
cd ~/projects/crc_visiumhd
sbatch tumor_subtypes_sample.sh P1CRC
sbatch tumor_subtypes_sample.sh P2CRC
sbatch tumor_subtypes_sample.sh P5CRC

grep -A 15 "authors" ~/projects/crc_visiumhd/tumor_*.log
```

```bash
cat /home/weberl/projects/crc_visiumhd/bench_12393128.log
```

```text
    authors
ours Tumor I Tumor II Tumor III Tumor IV Tumor V
  0        4       69      3477    47256    5941
  9      390     2112      1240     8357    2281
  11      23      473        10     6505     110
  12      47       21       558    13111     551
  17       0        7         8      515      96
```

```bash
cat /home/weberl/projects/crc_visiumhd/bench_12393127.log
```

```text
    authors
ours Tumor I Tumor II Tumor III Tumor V
  0        0        1     71113       8
  1        2       13     61404       0
  7        0        1     23219      47
  8        1       11     19946       6
  10       5      384     14935      48
  15       0       18      3542      14
  17       0        0      1920       1
```

```bash
cat /home/weberl/projects/crc_visiumhd/bench_12393126.log
```

```text
    authors
ours Tumor I Tumor II Tumor III Tumor IV Tumor V
   0       6    54370       586        1     998
   3       2    27515        42        0      46
   8       3     7202         1        7      23
   9      13    10272       160        2     210
```

P1: Tumor II, P2: Tumor III, P5: Tumor IV

Assigning names to clusters in figure, croping, adding scale bar to match Fig 3a:

```bash
cat > ~/projects/crc_visiumhd/fig3a_final_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat); library(ggplot2); library(patchwork)

id <- read.csv(paste0("~/projects/crc_visiumhd/cluster_identity_", sample, ".csv"))
cl <- read.csv(paste0("~/projects/crc_visiumhd/r_clusters_", sample, ".csv"))
sizes <- table(cl$cluster)                                  # bins per cluster

subtype <- c(P1CRC = "Tumor II", P2CRC = "Tumor III", P5CRC = "Tumor IV")[sample]  # from the cross-tables
tum <- id$cluster[id$identity == "Tumor"]
tum <- tum[order(-sizes[as.character(tum)])]                # tumor clusters, largest first
id$identity[id$cluster == tum[1]] <- subtype                # largest carries the patient subtype
if (length(tum) > 1) id$identity[id$cluster == tum[2]] <- paste(subtype, "b")
if (length(tum) > 2) id$identity[id$cluster %in% tum[-(1:2)]] <- "Tumor (other)"

id <- id[order(id$identity), ]
id$num <- ave(id$cluster, id$identity, FUN = seq_along) - 1
id$label <- paste(id$identity, id$num)

obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
xy <- GetTissueCoordinates(obj)
xy$label <- setNames(id$label, id$cluster)[as.character(obj$seurat_clusters)]

tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp                      # final orientation, all samples

pal <- list(
  "Tumor II"              = "#C17A97", "Tumor II b"  = "#6A2119",
  "Tumor III"             = "#ECE6A0", "Tumor III b" = "#6A2119",
  "Tumor IV"              = "#443E82", "Tumor IV b"  = "#6A2119",
  "Tumor (other)"         = "#807F3E",
  "Intestinal Epithelial" = c("#7EC4CD", "#DDAC76", "#B5DAC8", "#9BC653", "#000000", "#DC8B7A"),
  "Fibroblast"            = c("#812320", "#E48679", "#6684BD"),
  "Endothelial"           = c("#C9A785", "#BBD9E0", "#74559E"),
  "Myeloid"               = c("#A5594A", "#9CCB85", "#C59E9C"),
  "B cells"               = c("#6080B0", "#6594C1", "#6F3182", "#AECA63"),
  "T cells"               = c("#374394", "#54816F", "#E79748"),
  "Smooth Muscle"         = c("#75AF7F", "#F7E774", "#C4505C", "#8F2F3B"),
  "Neuronal"              = c("#6E69AA", "#BDDAA2"))
id$color <- unlist(lapply(split(id, id$identity), function(d) rep_len(pal[[d$identity[1]]], nrow(d))))[order(order(id$identity))]
cols <- setNames(id$color, id$label)

sf <- paste(readLines(paste0("~/data/crc_visiumhd/", sample, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
mm <- 1000 / mpp
if (sample == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
y0 <- max(xy$y)

p_main <- ggplot(xy, aes(x, y, color = label)) + geom_point(size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = cols) + scale_y_reverse() + coord_fixed() + theme_void() +
  annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "#E8C51A") +
  ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5))

fams <- unique(id$identity)
colno <- setNames(rep(1:2, each = ceiling(length(fams) / 2), length.out = length(fams)), fams)
id$lx <- colno[id$identity] * 3
y <- c(0, 0); id$ly <- NA; hd <- NULL
for (f in fams) {
  k <- sum(id$identity == f); cx <- colno[f]
  id$ly[id$identity == f] <- y[cx] + seq_len(k)
  hd <- rbind(hd, data.frame(f = f, lx = cx * 3, ly = y[cx] + (1 + k) / 2))
  y[cx] <- y[cx] + k + 1
}
id$numlab <- ifelse(grepl("Tumor", id$identity), "", id$num)   # tumor entries: name only, no index
p_leg <- ggplot() +
  geom_tile(data = id, aes(lx + 0.9, -ly), width = 0.3, height = 0.8, fill = id$color) +
  geom_text(data = id, aes(lx + 1.2, -ly, label = numlab), hjust = 0, size = 3) +
  geom_text(data = hd, aes(lx + 0.7, -ly, label = f), hjust = 1, size = 3.2) +
  xlim(0.5, 8) + coord_fixed(ratio = 0.55) + theme_void()

ggsave(paste0("~/projects/crc_visiumhd/figures/03a_final_", sample, ".png"),
       p_main + p_leg + plot_layout(widths = c(3, 1.4)), width = 12, height = 9, dpi = 150, bg = "white")
EOF
cd ~/projects/crc_visiumhd
sbatch fig3a_final_sample.sh P1CRC
sbatch fig3a_final_sample.sh P2CRC
sbatch fig3a_final_sample.sh P5CRC
```

```bash
cd ~/projects/crc_visiumhd
sed -i 's|y0 <- max(xy\$y)$|y0 <- max(xy$y) - 0.04 * diff(range(xy$y))|' fig3a_final_sample.R   # scale bar higher
sed -i 's/ratio = 0.55/ratio = 0.35/' fig3a_final_sample.R                                       # compact legend
sed -i 's/label = numlab/label = num/' fig3a_final_sample.R                                      # numbers on all entries
sbatch fig3a_final_sample.sh P1CRC
sbatch fig3a_final_sample.sh P2CRC
sbatch fig3a_final_sample.sh P5CRC
```

```bash
conda activate pyspatial
python << 'PY'
from PIL import Image
imgs = [Image.open(f"/home/weberl/projects/crc_visiumhd/figures/03a_final_{s}.png") for s in ["P1CRC", "P2CRC", "P5CRC"]]
h = min(i.height for i in imgs)                                   # common height
imgs = [i.resize((int(i.width * h / i.height), h)) for i in imgs]
out = Image.new("RGB", (sum(i.width for i in imgs), h), "white")
x = 0
for i in imgs:
    out.paste(i, (x, 0)); x += i.width
out.save("/home/weberl/projects/crc_visiumhd/figures/fig3a_row.png")
print("saved")
PY
```

Overnight run: redoing the clustering and rctd cell identification by doing batch analysis of the 3 patients:

```bash
cat > ~/projects/crc_visiumhd/joint_pipeline.R << 'EOF'
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
EOF
```

```bash
cat > ~/projects/crc_visiumhd/joint_stitch.py << 'EOF'
from PIL import Image, ImageChops
def trim(im):                                               # crop white borders
    bg = Image.new("RGB", im.size, "white")
    box = ImageChops.difference(im.convert("RGB"), bg).getbbox()
    return im.crop(box)
base = "/home/weberl/projects/crc_visiumhd/figures/"
names = ["joint_panel_P1CRC.png", "joint_panel_P2CRC.png", "joint_panel_P5CRC.png", "joint_legend.png"]
imgs = [trim(Image.open(base + n)) for n in names]
h = min(i.height for i in imgs[:3])                         # common height from the panels
imgs = [i.resize((int(i.width * h / i.height), h)) for i in imgs]
pad = 40
out = Image.new("RGB", (sum(i.width for i in imgs) + pad * 5, h + pad * 2), "white")
x = pad
for i in imgs:
    out.paste(i, (x, pad)); x += i.width + pad
out.save(base + "fig3a_joint.png")
print("stitched")
EOF
```

```bash
cat > ~/projects/crc_visiumhd/joint_pipeline.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=joint
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/joint_pipeline.R           # cluster, name, render panels
conda activate pyspatial
python ~/projects/crc_visiumhd/joint_stitch.py             # trim + assemble the final row
EOF
sbatch ~/projects/crc_visiumhd/joint_pipeline.sh
```

Done early. Good. Simply modifying the legend:

```bash
cat > ~/projects/crc_visiumhd/joint_figs.R << 'EOF'
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
EOF
cd ~/projects/crc_visiumhd
sbatch joint_figs.sh
```

## Figure 3b

```bash
cat > ~/projects/crc_visiumhd/fig3b_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr); library(ggplot2)

rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
xy <- as.data.frame(rctd@spatialRNA@coords[rownames(res), ])
xy$type <- as.character(res$first_type)
xy <- xy[res$spot_class != "reject", ]

subtype <- c(P1CRC = "Tumor II", P2CRC = "Tumor III", P5CRC = "Tumor IV")[sample]
xy$type[xy$type == "Tumor"] <- subtype

tmp <- xy$x; xy$x <- xy$y; xy$y <- tmp                      # final orientation

cols <- c("Tumor II" = "#361E50", "Tumor III" = "#423488", "Tumor IV" = "#C19C41",
          "Intestinal Epithelial" = "#582E76", "Fibroblast" = "#4E6A83",
          "Endothelial" = "#752A63", "Myeloid" = "#7D9A60", "B cells" = "#C1483C",
          "T cells" = "#EDE690", "Smooth Muscle" = "#B36743", "Neuronal" = "#F6C548")

sf <- paste(readLines(paste0("~/data/crc_visiumhd/", sample, "/binned_outputs/square_008um/spatial/scalefactors_json.json")), collapse = "")
mpp <- as.numeric(sub('.*"microns_per_pixel": *([0-9.eE+-]+).*', "\\1", sf))
mm <- 1000 / mpp
if (sample == "P1CRC") { x0 <- max(xy$x) - 0.05 * diff(range(xy$x)) - mm } else { x0 <- min(xy$x) + 0.05 * diff(range(xy$x)) }
y0 <- max(xy$y) - 0.04 * diff(range(xy$y))

p <- ggplot(xy, aes(x, y, color = type)) + geom_point(size = 0.25) +
  scale_color_manual(values = cols) + scale_y_reverse() + coord_fixed() + theme_void() +
  annotate("segment", x = x0, xend = x0 + mm, y = y0, yend = y0, linewidth = 3, color = "#E8C51A") +
  ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ar <- diff(range(xy$y)) / diff(range(xy$x))
ggsave(paste0("~/projects/crc_visiumhd/figures/03b_final_", sample, ".png"), p,
       width = 8 / max(ar, 0.6) + 3, height = 8, dpi = 150, bg = "white")
EOF
cd ~/projects/crc_visiumhd
sbatch fig3b_sample.sh P2CRC
sbatch fig3b_sample.sh P1CRC
sbatch fig3b_sample.sh P5CRC
```

```bash
cd ~/projects/crc_visiumhd
sed -i 's/geom_point(size = 0.25)/geom_point(size = 0.25, show.legend = FALSE)/' fig3b_sample.R

cat > fig3b_legend.R << 'EOF'
library(ggplot2)
cols <- c("Tumor II" = "#361E50", "Tumor III" = "#423488", "Tumor IV" = "#C19C41",
          "Intestinal Epithelial" = "#582E76", "Fibroblast" = "#4E6A83",
          "Endothelial" = "#752A63", "Myeloid" = "#7D9A60", "B cells" = "#C1483C",
          "T cells" = "#EDE690", "Smooth Muscle" = "#B36743", "Neuronal" = "#F6C548")
d <- data.frame(type = names(cols), y = -seq_along(cols))
p <- ggplot(d) + geom_tile(aes(1, y), width = 0.25, height = 0.75, fill = cols) +
  geom_text(aes(1.2, y, label = type), hjust = 0, size = 3.5) +
  xlim(0.8, 4) + coord_fixed(ratio = 0.5) + theme_void()
ggsave("~/projects/crc_visiumhd/figures/fig3b_legend.png", p, width = 3, height = 6, dpi = 150, bg = "white")
EOF
Rscript fig3b_legend.R    # tiny, fine on the login node

sbatch fig3b_sample.sh P1CRC
sbatch fig3b_sample.sh P2CRC
sbatch fig3b_sample.sh P5CRC
```

```bash
conda activate rspatial
Rscript ~/projects/crc_visiumhd/fig3b_legend.R
ls -lh ~/projects/crc_visiumhd/figures/fig3b_legend.png
```

```bash
conda activate pyspatial
python << 'PY'
from PIL import Image, ImageChops
def trim(im):
    bg = Image.new("RGB", im.size, "white")
    return im.crop(ImageChops.difference(im.convert("RGB"), bg).getbbox())
base = "/home/weberl/projects/crc_visiumhd/figures/"
names = ["03b_final_P1CRC.png", "03b_final_P2CRC.png", "03b_final_P5CRC.png", "fig3b_legend.png"]
imgs = [trim(Image.open(base + n)) for n in names]
h = min(i.height for i in imgs[:3])
imgs = [i.resize((int(i.width * h / i.height), h)) for i in imgs[:3]] + [imgs[3]]
pad = 40
out = Image.new("RGB", (sum(i.width for i in imgs) + pad * 5, h + pad * 2), "white")
x = pad
for i in imgs:
    out.paste(i, (x, pad)); x += i.width + pad
out.save(base + "fig3b_row.png")
print("stitched")
PY
```

T cells look like too many

```bash
cat > ~/projects/crc_visiumhd/check_tcells.R << 'EOF'
library(spacexr)
rctd <- readRDS("~/projects/crc_visiumhd/04_rctd_P1CRC.rds")
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_P1CRC.csv.gz")
m <- merge(mine[mine$type == "T cells", ], authors[, c("barcode", "DeconvolutionLabel1")], by = "barcode")
print(sort(table(m$DeconvolutionLabel1), decreasing = TRUE))   # what the authors called our "T cell" bins
EOF
sed "s|benchmark_sample.R \"\$1\"|check_tcells.R|" ~/projects/crc_visiumhd/benchmark_sample.sh > ~/projects/crc_visiumhd/check_tcells.sh
cd ~/projects/crc_visiumhd
sbatch check_tcells.sh

cat ~/projects/crc_visiumhd/bench_12407338.log
```

```text
               CD4 T cell                  Mature B                CD8 T cell
                     7823                      3234                      2153
                      CAF                Fibroblast Proliferating Macrophages
                     1240                       665                       344
                   mRegDC   Proliferating Immune II  Proliferating Fibroblast
                      334                       302                       120
    Lymphatic Endothelial               Endothelial                 Adipocyte
                      109                       104                        93
                Pericytes             Myofibroblast                  Tumor II
                       90                        89                        84
                     Mast                     cDC I                Neutrophil
                       77                        72                        66
               Macrophage                  Memory B                       vSM
                       65                        57                        56
                   Goblet                    Plasma             Enteric Glial
                       46                        40                        27
                      pDC          Unknown III (SM)            Neuroendocrine
                       22                         9                         5
                       NK                  Tumor IV                Enterocyte
                        5                         4                         2
               Epithelial                   Tumor V                      Tuft
                        2                         2                         1
```

Follicular mature B cells were classified under T cells. Re analysing with keeping lymphocytes subsets splits:

```bash
cat > ~/projects/crc_visiumhd/ref_custom.R << 'EOF'
library(Seurat); library(spacexr)
counts <- Read10X_h5("~/data/crc_visiumhd/chromium/HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")
meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
meta <- meta[meta$QCFilter == "Keep", ]
lab <- meta$Level1                                          # start from the coarse classes
lab[meta$Level2 == "Plasma"] <- "Plasma"                    # split B: follicle fix
lab[meta$Level2 %in% c("Mature B", "Memory B")] <- "Mature B"
lab[meta$Level2 == "CD4 T cell"] <- "CD4 T"                 # split T
lab[meta$Level2 %in% c("CD8 T cell", "NK")] <- "CD8 T NK"
lab[meta$Level2 == "CAF"] <- "CAF"                          # split fibroblast: CAF vs normal
names(lab) <- meta$Barcode
common <- intersect(colnames(counts), meta$Barcode)
lab <- lab[common]
keep <- names(table(lab))[table(lab) >= 25]
cells <- names(lab)[lab %in% keep]
ref <- Reference(counts[, cells], as.factor(lab[cells]))
saveRDS(ref, "~/projects/crc_visiumhd/rctd_reference_custom.rds")
cat("classes:\n"); print(table(lab[cells]))
EOF
cat > ~/projects/crc_visiumhd/ref_custom.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=refC
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/ref_custom.R
EOF

cat > ~/projects/crc_visiumhd/rctd_custom_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(Seurat); library(spacexr)
obj <- readRDS(paste0("~/projects/crc_visiumhd/02_clustered_", sample, ".rds"))
ref <- readRDS("~/projects/crc_visiumhd/rctd_reference_custom.rds")
coords <- GetTissueCoordinates(obj)[, c("x", "y")]
cts <- GetAssayData(obj, layer = "counts")
puck <- SpatialRNA(coords, cts)
rctd <- create.RCTD(puck, ref, max_cores = 8)
rctd <- run.RCTD(rctd, doublet_mode = "doublet")
saveRDS(rctd, paste0("~/projects/crc_visiumhd/04_rctd_custom_", sample, ".rds"))
EOF
cat > ~/projects/crc_visiumhd/rctd_custom_sample.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=rctdC
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/rctd_custom_sample.R "$1"
EOF

cd ~/projects/crc_visiumhd
jid=$(sbatch --parsable ref_custom.sh)
sbatch --dependency=afterok:$jid rctd_custom_sample.sh P1CRC
sbatch --dependency=afterok:$jid rctd_custom_sample.sh P2CRC
sbatch --dependency=afterok:$jid rctd_custom_sample.sh P5CRC
```

```bash
cd ~/projects/crc_visiumhd
squeue -u weberl
ls -lh rctd_reference_custom.rds 04_rctd_custom_*.rds 2>/dev/null
tail -15 $(ls -t refC_*.log 2>/dev/null | head -1)
```

```text
Warning message:
In Reference(counts[, cells], as.factor(lab[cells])) :
  Reference: number of cells per cell type is 63229, larger than maximum allowable of 10000. Downsampling number of cells to: 10000
classes:
              B cells                   CAF                 CD4 T
                 2809                 12080                 16830
             CD8 T NK           Endothelial            Fibroblast
                12884                  7884                 18969
Intestinal Epithelial              Mature B               Myeloid
                24641                  8761                 24653
             Neuronal                Plasma         Smooth Muscle
                 4310                 22116                 41059
              T cells                 Tumor
                  281                 63229
```

```bash
ls -lh ~/projects/crc_visiumhd/04_rctd_custom_*.rds
```

```text
-rw-rw-r-- 1 weberl weberl 314M Sep 11 19:51 /home/weberl/projects/crc_visiumhd/04_rctd_custom_P1CRC.rds
```

Starting with P1 files which is ready, with new lymphocytes and fibroblasts annotation:

```bash
cat > ~/projects/crc_visiumhd/fig3b_custom_sample.R << 'EOF'
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
  ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5))
ar <- diff(range(xy$y)) / diff(range(xy$x))
ggsave(paste0("~/projects/crc_visiumhd/figures/03bC_", sample, ".png"), p,
       width = 8 / max(ar, 0.6), height = 8, dpi = 150, bg = "white")
EOF
sed "s|fig3b_sample.R|fig3b_custom_sample.R|; s|--job-name=fig3b|--job-name=fig3bC|" ~/projects/crc_visiumhd/fig3b_sample.sh > ~/projects/crc_visiumhd/fig3b_custom_sample.sh
cd ~/projects/crc_visiumhd
sbatch fig3b_custom_sample.sh P1CRC
```

Comparison on the new 14-class deconvolution, my analysis vs the authors:

```bash
cat > ~/projects/crc_visiumhd/benchmark_custom_sample.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_custom_", sample, ".rds"))
res <- rctd@results$results_df
mine <- data.frame(barcode = rownames(res), mine = as.character(res$first_type))[res$spot_class != "reject", ]
authors <- read.csv(paste0("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/DeconvolutionResults_", sample, ".csv.gz"))
authors <- authors[authors$DeconvolutionClass != "reject" & !is.na(authors$DeconvolutionLabel1), c("barcode", "DeconvolutionLabel1")]

meta <- read.csv("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/SingleCell_MetaData.csv.gz")
lab <- meta$Level1
lab[meta$Level2 == "Plasma"] <- "Plasma"
lab[meta$Level2 %in% c("Mature B", "Memory B")] <- "Mature B"
lab[meta$Level2 == "CD4 T cell"] <- "CD4 T"
lab[meta$Level2 %in% c("CD8 T cell", "NK")] <- "CD8 T NK"
lab[meta$Level2 == "CAF"] <- "CAF"
tab <- table(meta$Level2, lab)
lookup <- setNames(colnames(tab)[max.col(tab)], rownames(tab))
authors$cls <- lookup[authors$DeconvolutionLabel1]
authors$cls[is.na(authors$cls) & grepl("Tumor", authors$DeconvolutionLabel1)] <- "Tumor"
authors$cls[is.na(authors$cls) & grepl("SM", authors$DeconvolutionLabel1)] <- "Smooth Muscle"
authors$cls[is.na(authors$cls) & grepl("Fibro", authors$DeconvolutionLabel1)] <- "Fibroblast"

m <- merge(mine, authors, by = "barcode"); m <- m[!is.na(m$cls), ]
cat(sample, "bins compared:", nrow(m), "\n")
cat(sample, "agreement:", mean(m$mine == m$cls), "\n")
print(table(mine = m$mine, authors = m$cls))
EOF
sed "s|benchmark_sample.R|benchmark_custom_sample.R|" ~/projects/crc_visiumhd/benchmark_sample.sh > ~/projects/crc_visiumhd/benchmark_custom_sample.sh
cd ~/projects/crc_visiumhd
sbatch benchmark_custom_sample.sh P1CRC

grep "agreement" $(ls -t bench_*.log | head -1)
```

```text
P1CRC agreement: 0.9173139
```

## Figure 3c

```bash
cat > ~/projects/crc_visiumhd/fig3c_sample.R << 'EOF'
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
EOF
cd ~/projects/crc_visiumhd
sbatch fig3c_sample.sh P1CRC
sbatch fig3c_sample.sh P2CRC
sbatch fig3c_sample.sh P5CRC
```

Labels and legend:

```bash
cat > ~/projects/crc_visiumhd/fig3c_extras.R << 'EOF'
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
EOF
conda activate rspatial
Rscript ~/projects/crc_visiumhd/fig3c_extras.R
cd ~/projects/crc_visiumhd
sbatch fig3c_sample.sh P1CRC
sbatch fig3c_sample.sh P2CRC
sbatch fig3c_sample.sh P5CRC
```

Stitching:

```bash
conda activate pyspatial
python << 'PY'
from PIL import Image, ImageChops
def trim(im):
    bg = Image.new("RGB", im.size, "white")
    return im.crop(ImageChops.difference(im.convert("RGB"), bg).getbbox())
base = "/home/weberl/projects/crc_visiumhd/figures/"
genes = ["PIGR", "CEACAM6", "COL1A1"]; samples = ["P1CRC", "P2CRC", "P5CRC"]
pad = 30; rows = []
for g in genes:
    imgs = [trim(Image.open(base + f"03c_lab_{g}.png"))]
    for s in samples:
        imgs.append(trim(Image.open(base + f"03c_{g}_{s}.png")))
        imgs.append(trim(Image.open(base + f"03c_zoom_{g}_{s}.png")))
    h = min(i.height for i in imgs[1:])
    imgs = [imgs[0]] + [i.resize((int(i.width * h / i.height), h)) for i in imgs[1:]]
    rows.append(imgs)
bar = trim(Image.open(base + "03c_colorbar.png"))
W = max(sum(i.width for i in r) + pad * 8 for r in rows) + bar.width + pad
H = sum(max(i.height for i in r) for r in rows) + pad * 4
out = Image.new("RGB", (W, H), "white")
ypos = pad
for r in rows:
    rh = max(i.height for i in r); x = pad
    for i in r:
        out.paste(i, (x, ypos + (rh - i.height) // 2)); x += i.width + pad
    ypos += rh + pad
out.paste(bar, (W - bar.width - pad, (H - bar.height) // 2))
out.save(base + "fig3c_grid.png")
print("stitched")
PY
```

Re doing python clustering annotation with rctd single cell labels: batching the three patients like we re-did for R clustering:

```bash
cat > ~/projects/crc_visiumhd/py_joint.py << 'EOF'
import os, pandas as pd, scanpy as sc, anndata as ad
home = os.path.expanduser("~")
samples = ["P1CRC", "P2CRC", "P5CRC"]

parts = []
for s in samples:
    a = sc.read_10x_h5(f"{home}/data/crc_visiumhd/{s}/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()
    pos = pd.read_parquet(f"{home}/data/crc_visiumhd/{s}/binned_outputs/square_008um/spatial/tissue_positions.parquet").set_index("barcode")
    pos = pos.loc[a.obs_names]
    a.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None)
    mt = a.var_names.str.startswith("MT-")
    a.obs["pct_mt"] = (a[:, mt].X.sum(1).A1 / a.obs["total_counts"]) * 100
    a = a[(a.obs["total_counts"] >= 100) & (a.obs["n_genes_by_counts"] >= 50) & (a.obs["pct_mt"] < 30)].copy()  # same QC as before
    a.obs["sample"] = s
    a.obs_names = s + "_" + a.obs_names
    parts.append(a)

adata = ad.concat(parts, join="inner"); del parts           # one object, all patients
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.scale(adata); sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.4, key_added="joint", flavor="igraph", n_iterations=2)
adata.obs[["sample", "joint"]].rename_axis("barcode").to_csv(f"{home}/projects/crc_visiumhd/joint_clusters_py.csv")  # checkpoint
print(adata.obs["joint"].value_counts())
EOF
cat > ~/projects/crc_visiumhd/py_joint.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pyjoint
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/py_joint.py
EOF
cd ~/projects/crc_visiumhd
sbatch py_joint.sh
```

ARI Comparison of R clustering vs Python clustering (by how much are identical bins grouped together with the two pipelines?)

1. Authors's parquet file: what is their clustering like (per bin cluster labels)

```bash
conda activate rspatial
Rscript -e 'library(arrow); d <- read_parquet("~/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/P1CRC_Metadata.parquet"); cat(colnames(d), sep = "\n"); print(as.data.frame(head(d, 3)))'
```

```text
Attaching package: 'arrow'
The following object is masked from 'package:utils':
    timestamp
barcode
tissue
X
Y
DeconvolutionClass
DeconvolutionLabel1
DeconvolutionLabel2
Periphery
UnsupervisedL1
UnsupervisedL2
MacrophageSubtype
GobletSubcluster
                barcode tissue         X         Y DeconvolutionClass
1 s_008um_00000_00000-1      0 -2.042794 -308.7646               <NA>
2 s_008um_00000_00001-1      0 -1.796281 -308.7669               <NA>
3 s_008um_00000_00002-1      0 -1.549768 -308.7692               <NA>
  DeconvolutionLabel1 DeconvolutionLabel2 Periphery UnsupervisedL1
1                <NA>                <NA>      <NA>           <NA>
2                <NA>                <NA>      <NA>           <NA>
3                <NA>                <NA>      <NA>           <NA>
  UnsupervisedL2 MacrophageSubtype GobletSubcluster
1           <NA>              <NA>             <NA>
2           <NA>              <NA>             <NA>
3           <NA>              <NA>             <NA>
```

UnsupervisedL1/ UnsupervisedL2: the author's cluster assignement per bin.

Comparing our R joint clustering vs the paper's, per patient:

```bash
cat > ~/projects/crc_visiumhd/ari_authors.py << 'EOF'
import os, pandas as pd
from sklearn.metrics import adjusted_rand_score
home = os.path.expanduser("~")
joint = pd.read_csv(f"{home}/projects/crc_visiumhd/joint_clusters.csv")

for s in ["P1CRC", "P2CRC", "P5CRC"]:
    a = pd.read_parquet(f"{home}/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/{s}_Metadata.parquet",
                        columns=["barcode", "UnsupervisedL2"]).dropna().set_index("barcode")   # authors' clusters
    j = joint[joint["sample"] == s].copy()                                                     # our joint clusters
    j["bc"] = j["barcode"].str.replace(s + "_", "", regex=False)
    j = j.set_index("bc")
    common = j.index.intersection(a.index)
    print(f"{s} joint vs authors: {len(common)} bins, ARI = "
          f"{adjusted_rand_score(j.loc[common, 'cluster'], a.loc[common, 'UnsupervisedL2']):.3f}")
EOF
cat > ~/projects/crc_visiumhd/ari_authors.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ariA
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/ari_authors.py
EOF
cd ~/projects/crc_visiumhd
sbatch ari_authors.sh
```

```bash
cd ~/projects/crc_visiumhd
cat $(ls -t ariA_*.log | head -1)                      # joint vs authors ARIs
squeue -u weberl                                       # what's still running
ls -lh joint_clusters_py.csv 04_rctd_custom_*.rds 2>/dev/null
ls figures/rvspy_*.png 2>/dev/null
```

```text
P1CRC joint vs authors: 286973 bins, ARI = 0.816
P2CRC joint vs authors: 428405 bins, ARI = 0.605
P5CRC joint vs authors: 373156 bins, ARI = 0.730
```

R vs Python ARI analysis:

ARI job, batch analysis (all patients analysed together):

```bash
cat > ~/projects/crc_visiumhd/ari_joint_rvspy.py << 'EOF'
import os, pandas as pd
from sklearn.metrics import adjusted_rand_score
home = os.path.expanduser("~")
jr = pd.read_csv(f"{home}/projects/crc_visiumhd/joint_clusters.csv").set_index("barcode")      # R joint
jp = pd.read_csv(f"{home}/projects/crc_visiumhd/joint_clusters_py.csv").set_index("barcode")   # Python joint
common = jr.index.intersection(jp.index)
print(f"overall: {len(common)} bins, ARI = {adjusted_rand_score(jr.loc[common, 'cluster'], jp.loc[common, 'joint']):.3f}")
for s in ["P1CRC", "P2CRC", "P5CRC"]:
    c = [b for b in common if b.startswith(s)]
    print(f"{s}: {len(c)} bins, ARI = {adjusted_rand_score(jr.loc[c, 'cluster'], jp.loc[c, 'joint']):.3f}")
EOF
cat > ~/projects/crc_visiumhd/ari_joint_rvspy.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ariJ
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/ari_joint_rvspy.py
EOF
cd ~/projects/crc_visiumhd
sbatch ari_joint_rvspy.sh

cat $(ls -t ~/projects/crc_visiumhd/ariJ_*.log | head -1)
```

```text
overall: 1088534 bins, ARI = 0.355
P1CRC: 286973 bins, ARI = 0.454
P2CRC: 428405 bins, ARI = 0.418
P5CRC: 373156 bins, ARI = 0.498
```

Cell types per bins exports (rctd deconvolution):

```bash
cat > ~/projects/crc_visiumhd/export_types.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE); sample <- args[1]
library(spacexr)
rctd <- readRDS(paste0("~/projects/crc_visiumhd/04_rctd_", sample, ".rds"))
res <- rctd@results$results_df
out <- data.frame(barcode = rownames(res), type = as.character(res$first_type))[res$spot_class != "reject", ]
write.csv(out, paste0("~/projects/crc_visiumhd/rctd_types_", sample, ".csv"), row.names = FALSE)
EOF
cat > ~/projects/crc_visiumhd/export_types.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=types
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate rspatial
Rscript ~/projects/crc_visiumhd/export_types.R "$1"
EOF
cd ~/projects/crc_visiumhd
sbatch export_types.sh P1CRC
sbatch export_types.sh P2CRC
sbatch export_types.sh P5CRC
```

Plotting script and ticket:

```bash
cat > ~/projects/crc_visiumhd/py_rvspy_sample.py << 'EOF'
import sys, os, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

sample = sys.argv[1]
home = os.path.expanduser("~")
tumor = {"P1CRC": ("Tumor II", "#C17A97"), "P2CRC": ("Tumor III", "#ECE6A0"), "P5CRC": ("Tumor IV", "#443E82")}[sample]
cols = {"Tumor": tumor[1], "Intestinal Epithelial": "#7EC4CD", "Fibroblast": "#812320",
        "Endothelial": "#C9A785", "Myeloid": "#A5594A", "B cells": "#6080B0",
        "T cells": "#374394", "Smooth Muscle": "#75AF7F", "Neuronal": "#6E69AA"}

adata = sc.read_h5ad(f"{home}/projects/crc_visiumhd/py_02_clustered_{sample}.h5ad")
r = pd.read_csv(f"{home}/projects/crc_visiumhd/r_clusters_{sample}.csv").set_index("barcode")
rid = pd.read_csv(f"{home}/projects/crc_visiumhd/cluster_identity_{sample}.csv").set_index("cluster")
types = pd.read_csv(f"{home}/projects/crc_visiumhd/rctd_types_{sample}.csv").set_index("barcode")

jr = pd.read_csv(f"{home}/projects/crc_visiumhd/joint_clusters.csv").set_index("barcode")
jp = pd.read_csv(f"{home}/projects/crc_visiumhd/joint_clusters_py.csv").set_index("barcode")
jc = [b for b in jr.index.intersection(jp.index) if b.startswith(sample)]
ari_rp = adjusted_rand_score(jr.loc[jc, "cluster"], jp.loc[jc, "joint"])            # R vs Python

authors = pd.read_parquet(f"{home}/data/crc_visiumhd/HumanColonCancer_VisiumHD/MetaData/{sample}_Metadata.parquet",
                          columns=["barcode", "UnsupervisedL2"]).dropna().set_index("barcode")
js = jr[jr.index.str.startswith(sample)].copy()
js.index = js.index.str.replace(sample + "_", "", regex=False)
ca = js.index.intersection(authors.index)
ari_pub = adjusted_rand_score(js.loc[ca, "cluster"], authors.loc[ca, "UnsupervisedL2"])  # R vs published

common = adata.obs_names.intersection(r.index)
py = pd.DataFrame({"cl": adata.obs["leiden_low"].astype(str)}, index=adata.obs_names).join(types)
py_id = py.dropna().groupby("cl")["type"].agg(lambda s: s.mode()[0])

sub = adata.obs_names.isin(common)
xy = adata.obsm["spatial"][sub][:, [1, 0]]
r_col = [cols.get(rid.loc[c, "identity"], "grey70") for c in r.loc[adata.obs_names[sub], "cluster"]]
p_col = [cols.get(py_id.get(c, ""), "#BBBBBB") for c in adata.obs.loc[sub, "leiden_low"].astype(str)]

fig, axes = plt.subplots(1, 2, figsize=(16, 8.6))
axes[0].scatter(xy[:, 0], xy[:, 1], s=0.5, c=r_col); axes[0].set_title(f"{sample} — R (Seurat)", fontsize=20)
axes[1].scatter(xy[:, 0], xy[:, 1], s=0.5, c=p_col); axes[1].set_title(f"{sample} — Python (scanpy)", fontsize=20)
for ax in axes: ax.invert_yaxis(); ax.axis("off"); ax.set_aspect("equal")
plt.figtext(0.5, 0.02,
            f"Joint clustering concordance — R vs. Oliveira et al. (2025): ARI = {ari_pub:.3f}   |   R vs. Python: ARI = {ari_rp:.3f}",
            ha="center", fontsize=15)
plt.subplots_adjust(bottom=0.09)
plt.savefig(f"{home}/projects/crc_visiumhd/figures/rvspy_{sample}.png", dpi=150, bbox_inches="tight", facecolor="white")
EOF
cd ~/projects/crc_visiumhd
sbatch py_rvspy_sample.sh P1CRC
sbatch py_rvspy_sample.sh P2CRC
sbatch py_rvspy_sample.sh P5CRC
```

Stitching:

```bash
conda activate pyspatial
python << 'PY'
from PIL import Image, ImageChops
def trim(im):
    bg = Image.new("RGB", im.size, "white")
    return im.crop(ImageChops.difference(im.convert("RGB"), bg).getbbox())
base = "/home/weberl/projects/crc_visiumhd/figures/"
imgs = [trim(Image.open(base + f"rvspy_{s}.png")) for s in ["P1CRC", "P2CRC", "P5CRC"]]
w = min(i.width for i in imgs)
imgs = [i.resize((w, int(i.height * w / i.width))) for i in imgs]
pad = 40
out = Image.new("RGB", (w + pad * 2, sum(i.height for i in imgs) + pad * 4), "white")
y = pad
for i in imgs:
    out.paste(i, (pad, y)); y += i.height + pad
out.save(base + "fig_rvspy.png")
print("stitched")
PY
```

Edits:

```bash
cd ~/projects/crc_visiumhd
sed -i 's/plt.figtext(0.5, 0.02,/plt.figtext(0.5, 0.045,/' py_rvspy_sample.py
sed -i 's/plt.subplots_adjust(bottom=0.09)/plt.subplots_adjust(bottom=0.07)/' py_rvspy_sample.py

cat > rvspy_stitch.py << 'EOF'
from PIL import Image, ImageChops
def trim(im):
    bg = Image.new("RGB", im.size, "white")
    return im.crop(ImageChops.difference(im.convert("RGB"), bg).getbbox())
base = "/home/weberl/projects/crc_visiumhd/figures/"
imgs = [trim(Image.open(base + f"rvspy_{s}.png")) for s in ["P1CRC", "P2CRC", "P5CRC"]]
w = min(i.width for i in imgs)
imgs = [i.resize((w, int(i.height * w / i.width))) for i in imgs]
pad = 40
out = Image.new("RGB", (w + pad * 2, sum(i.height for i in imgs) + pad * 4), "white")
y = pad
for i in imgs:
    out.paste(i, (pad, y)); y += i.height + pad
out.save(base + "fig_rvspy.png")
print("stitched")
EOF
cat > rvspy_stitch.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=stitch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/rvspy_stitch.py
EOF

j1=$(sbatch --parsable py_rvspy_sample.sh P1CRC)
j2=$(sbatch --parsable py_rvspy_sample.sh P2CRC)
j3=$(sbatch --parsable py_rvspy_sample.sh P5CRC)
sbatch --dependency=afterok:$j1:$j2:$j3 rvspy_stitch.sh
```

Edits 2:

```bash
cd ~/projects/crc_visiumhd
sed -i 's/Joint clustering concordance/Clustering concordance/' py_rvspy_sample.py
sed -i -E 's/plt\.figtext\(0\.5, 0\.[0-9]+,/plt.figtext(0.5, 0.05,/' py_rvspy_sample.py
sed -i -E 's/plt\.subplots_adjust\(bottom=0\.[0-9]+\)/plt.subplots_adjust(bottom=0.06)/' py_rvspy_sample.py
grep -n "Clustering concordance\|figtext\|subplots_adjust" py_rvspy_sample.py   # verify all three lines

j1=$(sbatch --parsable py_rvspy_sample.sh P1CRC)
j2=$(sbatch --parsable py_rvspy_sample.sh P2CRC)
j3=$(sbatch --parsable py_rvspy_sample.sh P5CRC)
sbatch --dependency=afterok:$j1:$j2:$j3 rvspy_stitch.sh
```
