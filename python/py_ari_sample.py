import sys, os, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

sample = sys.argv[1]
home = os.path.expanduser("~")
adata = sc.read_h5ad(f"{home}/projects/crc_visiumhd/py_02_clustered_{sample}.h5ad")
r = pd.read_csv(f"{home}/projects/crc_visiumhd/r_clusters_{sample}.csv").set_index("barcode")

common = adata.obs_names.intersection(r.index)
ari = adjusted_rand_score(r.loc[common, "cluster"], adata.obs.loc[common, "leiden_low"])
print(f"{sample}: {len(common)} shared bins, ARI = {ari:.3f}")

xy = adata.obsm["spatial"]
plt.figure(figsize=(8,8))
plt.scatter(xy[:,0], xy[:,1], s=0.5, c=adata.obs["leiden_low"].astype(int), cmap="tab20")
plt.gca().invert_yaxis(); plt.axis('off'); plt.title(f"{sample} Python clusters")
plt.savefig(f"{home}/projects/crc_visiumhd/figures/py_02_clusters_{sample}.png", dpi = 150)
