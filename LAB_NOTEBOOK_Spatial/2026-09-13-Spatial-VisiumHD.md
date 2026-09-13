# 2026-09-13 — Spatial Visium HD

Comparing the cell type attribution between us vs authors: for each 8 µm bin, comparison of the RCTD generated cell labelling, matched by bin with the barcodes.

```bash
grep "agreement" ~/projects/crc_visiumhd/bench*.log | tail -3
```

```text
/home/weberl/projects/crc_visiumhd/bench_12488504.log:P1CRC agreement: 0.9173139
/home/weberl/projects/crc_visiumhd/bench_12527436.log:P2CRC agreement: 0.9183618
/home/weberl/projects/crc_visiumhd/bench_12527445.log:P5CRC agreement: 0.9137777
```

QC figure:

```bash
cat > ~/projects/crc_visiumhd/qc_fig.py << 'EOF'
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
EOF
```

```bash
cat > ~/projects/crc_visiumhd/qc_fig.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=qcfig --cpus-per-task=4 --mem=64G --time=2:00:00
#SBATCH --output=%x_%j.log
source /data1/test03/weberl/miniforge3/etc/profile.d/conda.sh
export TMPDIR=~/tmp
conda activate pyspatial
python ~/projects/crc_visiumhd/qc_fig.py
EOF
```

```bash
cd ~/projects/crc_visiumhd
sbatch qc_fig.sh
```

## Preparing the github repository

organizing files from the cluster:

```bash
cd ~/projects/crc_visiumhd
mkdir -p ~/crc-visiumhd-repro/{R,python,slurm,figures,results,notes}

cp r_cluster_sample.R r_export_clusters.R rctd_sample.R ref_custom.R rctd_custom_sample.R \
   markers_sample.R pathways_sample.R benchmark_sample.R benchmark_custom_sample.R \
   name_clusters_sample.R tumor_subtypes_sample.R check_tcells.R joint_pipeline.R joint_figs.R \
   fig3b_custom_sample.R fig3bC_legend.R fig3c_sample.R fig3c_extras.R export_types.R \
   ~/crc-visiumhd-repro/R/ 2>/dev/null

cp py_cluster_sample.py py_ari_sample.py py_joint.py ari_authors.py ari_joint_rvspy.py \
   py_rvspy_sample.py rvspy_stitch.py joint_stitch.py fig3bC_stitch.py qc_fig.py \
   ~/crc-visiumhd-repro/python/ 2>/dev/null

cp r_cluster_sample.sh py_cluster_sample.sh rctd_sample.sh ref_custom.sh rctd_custom_sample.sh \
   markers_sample.sh pathways_sample.sh benchmark_sample.sh benchmark_custom_sample.sh \
   name_clusters_sample.sh tumor_subtypes_sample.sh check_tcells.sh joint_pipeline.sh \
   fig3b_custom_sample.sh fig3bC_stitch.sh fig3c_sample.sh ari_sample.sh ari_authors.sh \
   ari_joint_rvspy.sh py_joint.sh py_rvspy_sample.sh rvspy_stitch.sh export_types.sh qc_fig.sh \
   ~/crc-visiumhd-repro/slurm/ 2>/dev/null

cp figures/fig3a_joint.png figures/fig3bC_row.png figures/fig3c_grid.png \
   figures/fig_rvspy.png figures/fig_qc.png ~/crc-visiumhd-repro/figures/ 2>/dev/null

cp cluster_identity_*.csv label_map.csv ~/crc-visiumhd-repro/results/ 2>/dev/null

ls -R ~/crc-visiumhd-repro
du -sh ~/crc-visiumhd-repro
```

```text
/home/weberl/crc-visiumhd-repro:
figures  notes  python  R  results  slurm
/home/weberl/crc-visiumhd-repro/figures:
fig3a_joint.png  fig3bC_row.png  fig3c_grid.png  fig_qc.png  fig_rvspy.png
/home/weberl/crc-visiumhd-repro/notes:
/home/weberl/crc-visiumhd-repro/python:
ari_authors.py      fig3bC_stitch.py  py_ari_sample.py      py_joint.py         qc_fig.py
ari_joint_rvspy.py  joint_stitch.py   py_cluster_sample.py  py_rvspy_sample.py  rvspy_stitch.py
/home/weberl/crc-visiumhd-repro/R:
benchmark_custom_sample.R  export_types.R         fig3c_extras.R  joint_pipeline.R        pathways_sample.R     rctd_sample.R        tumor_subtypes_sample.R
benchmark_sample.R         fig3bC_legend.R        fig3c_sample.R  markers_sample.R        r_cluster_sample.R    ref_custom.R
check_tcells.R             fig3b_custom_sample.R  joint_figs.R    name_clusters_sample.R  rctd_custom_sample.R  r_export_clusters.R
/home/weberl/crc-visiumhd-repro/results:
cluster_identity_joint.csv  cluster_identity_P1CRC.csv  cluster_identity_P2CRC.csv  cluster_identity_P5CRC.csv  label_map.csv
/home/weberl/crc-visiumhd-repro/slurm:
ari_authors.sh              benchmark_sample.sh  fig3b_custom_sample.sh  name_clusters_sample.sh  py_rvspy_sample.sh     rctd_sample.sh
ari_joint_rvspy.sh          check_tcells.sh      fig3c_sample.sh         pathways_sample.sh       qc_fig.sh              ref_custom.sh
ari_sample.sh               export_types.sh      joint_pipeline.sh       py_cluster_sample.sh     r_cluster_sample.sh    rvspy_stitch.sh
benchmark_custom_sample.sh  fig3bC_stitch.sh     markers_sample.sh       py_joint.sh              rctd_custom_sample.sh  tumor_subtypes_sample.sh
```

```text
18M     /home/weberl/crc-visiumhd-repro
```

Upload from the cluster to Github:

```bash
cd ~/crc-visiumhd-repro
git init -b main
git add .
git commit -m "Reproduction of Oliveira et al. 2025 Figure 3 - R and Python pipelines"
git remote add origin https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git
git push -u origin main
```

```text
Username for 'https://github.com': weberleslie5
Password for 'https://weberleslie5@github.com':
remote: Write access to repository not granted.
fatal: unable to access 'https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git/': The requested URL returned error: 403
```

Token permissions fixed (Contents: read and write), token embedded in the remote URL:

```bash
git remote set-url origin https://weberleslie5:<TOKEN>@github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git
```

```bash
git pull origin main --allow-unrelated-histories --no-edit
git push -u origin main
```

```text
remote: Enumerating objects: 3, done.
remote: Counting objects: 100% (3/3), done.
remote: Compressing objects: 100% (2/2), done.
remote: Total 3 (delta 0), reused 0 (delta 0), pack-reused 0 (from 0)
Unpacking objects:  33% (1/3)
Unpacking objects: 100% (3/3), 921 bytes | 34.00 KiB/s, done.
From https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python
 * branch            main       -> FETCH_HEAD
 * [new branch]      main       -> origin/main
hint: You have divergent branches and need to specify how to reconcile them.
fatal: Need to specify how to reconcile divergent branches.
```

```text
git push -u origin main
To https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git
 ! [rejected]        main -> main (non-fast-forward)
error: failed to push some refs to 'https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git'
hint: Updates were rejected because the tip of your current branch is behind
hint: its remote counterpart. If you want to integrate the remote changes,
hint: use 'git pull' before pushing again.
```

Importing to git:

```bash
git pull origin main --allow-unrelated-histories --no-rebase --no-edit
git push -u origin main
```

```text
From https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python
 * branch            main       -> FETCH_HEAD
Merge made by the 'ort' strategy.
 README.md | 2 ++
 1 file changed, 2 insertions(+)
 create mode 100644 README.md
```

```text
Enumerating objects: 73, done.
Counting objects: 100% (73/73), done.
Delta compression using up to 56 threads
Compressing objects: 100% (72/72), done.
Writing objects: 100% (72/72), 16.99 MiB | 19.48 MiB/s, done.
Total 72 (delta 25), reused 0 (delta 0), pack-reused 0
remote: Resolving deltas: 100% (25/25), done.
To https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git
   1ef2f49..199e344  main -> main
branch 'main' set up to track 'origin/main'.
```

Reset URL with no token:

```bash
git remote set-url origin https://github.com/weberleslie5/Spatial-colorectal-cancer-R-and-Python.git
```
