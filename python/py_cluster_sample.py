import sys, scanpy as sc, pandas as pd
sample = sys.argv[1]    # read sample name from the command line: one script for all samples
base = f"/home/weberl/data/crc_visiumhd/{sample}/binned_outputs/square_008um"
adata = sc.read_10x_h5(f"{base}/filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
coords = pd.read_parquet(f"{base}/spatial/tissue_positions.parquet").set_index("barcode").loc[adata.obs_names]
adata.obsm["spatial"] = coords[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
adata = adata[(adata.obs.total_counts >= 100) & (adata.obs.n_genes_by_counts >= 50) & (adata.obs.pct_counts_mt < 30)].copy()
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata; adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10); sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_low", flavor="igraph", n_iterations=2)
adata.write(f"/home/weberl/projects/crc_visiumhd/py_02_clustered_{sample}.h5ad")
