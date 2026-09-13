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
