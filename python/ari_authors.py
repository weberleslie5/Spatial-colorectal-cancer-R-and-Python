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
