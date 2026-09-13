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
