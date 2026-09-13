# 2026-09-06 — Spatial Visium HD

Extracting P1CRC

```bash
tail -2 ~/data/crc_visiumhd/P1CRC/wget_P1.log
```

```text
2026-09-05 22:17:30 (3.53 MB/s) - 'Visium_HD_Human_Colon_Cancer_P1_binned_outputs.tar.gz' saved [12319027011/12319027011]
```

```bash
cd ~/data/crc_visiumhd/P1CRC

ls binned_outputs/square_008um/ binned_outputs/square_008um/spatial/
```

```text
binned_outputs/square_008um/:
filtered_feature_bc_matrix.h5  spatial
binned_outputs/square_008um/spatial/:
aligned_fiducials.jpg     cytassist_image.tiff       scalefactors_json.json  tissue_lowres_image.png
aligned_tissue_image.jpg  detected_tissue_image.jpg  tissue_hires_image.png  tissue_positions.parquet
```

```bash
conda activate pyspatial
python -c "import h5py,os; f=h5py.File(os.path.expanduser('~/data/crc_visiumhd/P1CRC/binned_outputs/square_008um/filtered_feature_bc_matrix.h5')); print(f['matrix/barcodes'].shape)"
```

```text
(507684,)
```

```bash
rm -f ~/data/crc_visiumhd/P1CRC/*.tar.gz
cd ~/projects/crc_visiumhd
sbatch r_cluster_sample.sh P1CRC
sbatch py_cluster_sample.sh P1CRC
```

deconvolution:

```bash
sbatch rctd_sample.sh P1CRC
```

R vs Python concordance and cluster map

```bash
sbatch ari_sample.sh P1CRC
```
