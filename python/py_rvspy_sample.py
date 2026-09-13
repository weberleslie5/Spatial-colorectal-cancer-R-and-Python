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
plt.figtext(0.5, 0.05,
            f"Clustering concordance — R vs. Oliveira et al. (2025): ARI = {ari_pub:.3f}   |   R vs. Python: ARI = {ari_rp:.3f}",
            ha="center", fontsize=15)
plt.subplots_adjust(bottom=0.06)
plt.savefig(f"{home}/projects/crc_visiumhd/figures/rvspy_{sample}.png", dpi=150, bbox_inches="tight", facecolor="white")
