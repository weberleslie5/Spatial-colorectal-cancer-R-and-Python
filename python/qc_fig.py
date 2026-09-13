import numpy as np, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

samples = ["P1CRC", "P2CRC", "P5CRC"]
metrics = {"total_counts": ("Transcripts per bin", 100),
           "n_genes_by_counts": ("Genes per bin", 50),
           "pct_counts_mt": ("Mitochondrial fraction per bin (%)", 30)}
data = {}
for s in samples:
    a = sc.read_10x_h5(f"/home/weberl/data/crc_visiumhd/{s}/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    keep = (a.obs.total_counts >= 100) & (a.obs.n_genes_by_counts >= 50) & (a.obs.pct_counts_mt < 30)
    print(f"{s}: {a.n_obs} bins -> {keep.sum()} kept")
    data[s] = a.obs[list(metrics)].copy()

fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, (m, (label, thr)) in zip(axes, metrics.items()):
    vals = [np.maximum(data[s][m].to_numpy(), 0.1) for s in samples]    # avoid log(0)
    ax.violinplot(vals, showextrema=False)
    ax.axhline(thr, color="red", linestyle="--", linewidth=1)
    if m != "pct_counts_mt":
        ax.set_yscale("log")
        ax.set_ylabel(f"{label} (log scale)", fontsize=11)
    else:
        ax.set_ylabel(label, fontsize=11)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(samples)
    ax.set_xlabel("Sample", fontsize=11)
    ax.set_title(label, fontsize=12)
fig.suptitle("QC per 8 um bin — violin width: number of bins; dashed line: filter threshold (>= 100 transcripts, >= 50 genes, < 30% mito)", fontsize=10)
fig.tight_layout()
fig.savefig("/home/weberl/projects/crc_visiumhd/figures/fig_qc.png", dpi=150, bbox_inches="tight")
