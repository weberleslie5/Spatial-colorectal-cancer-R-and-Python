# 2026-08-12 — Spatial Visium HD

Step 1: connect to iris through terminal app:

```bash
ssh <MSKusername>@islogin01.mskcc.org
```

Step 2: Iris on Demand
In the Iris file browser: creating two folders: data/crc_visiumhd and projects.

```bash
mkdir -p ~/data/crc_visiumhd ~/projects
```

From the crc dataset: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc
In the Output and supplemental files folder, downloaded 'Binned outputs (all bin levels)'. Opened and compressed square_008um.

Using the Upload button on Iris on Demand: uploaded square_008um.
Upload did not work using the Iris Upload button.
On the Mac Terminal:

```bash
scp ~/Downloads/square_008um.zip weberl@islogin01.mskcc.org:~/data/crc_visiumhd/
```

Once scp done, in Iris Terminal:

```bash
cd ~/data/crc_visiumhd
ls -lh
```

```text
total 9.5G
-rw-r--r-- 1 weberl weberl 9.5G Aug 13 13:40 square_008um.zip
```

```bash
unzip *.zip
```

```text
error: invalid zip file with overlapped components (possible zip bomb)
To unzip the file anyway, rerun the command with UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE environmnent variable
```

Used instead:

```bash
python3 -m zipfile -e *.zip .
ls
```

```text
__MACOSX  square_008um  square_008um.zip
```

```bash
ls square_008um/spatial
```

```text
aligned_fiducials.jpg     cytassist_image.tiff       scalefactors_json.json  tissue_lowres_image.png
aligned_tissue_image.jpg  detected_tissue_image.jpg  tissue_hires_image.png  tissue_positions.parquet
```

Meanwhile on Iris terminal:

```bash
module avail R
```

```text
R/4.3.0
```

Request a compute node and make R and conda available:

```bash
srun --pty --cpus-per-task=8 --mem=32G --time=3:00:00 bash
module load R/4.3.0
module load miniforge3
R --version
```

```text
R version 4.3.0 (2023-04-21) -- "Already Tomorrow"
```

Install the R packages:

```bash
mkdir -p ~/R/library
echo 'R_LIBS_USER=~/R/library' >> ~/.Renviron # R to always use ~/R/library as my personal package folder.

Rscript -e 'install.packages(c("Seurat", "hdf5r", "arrow", "tidyverse", "patchwork", "RANN", "remotes"), repos="https://cloud.r-project.org"); remotes::install_github("dmcable/spacexr")'
```

Two errors:
1. Rfast and spacexr: compiler issue
2. igraph/SeuratObject/Seurat compiled but could not load

Cleanup of the failed attempt:

```bash
rm -rf ~/R/library
rm -f ~/.Renviron
```

I built a conda environment (named rspatial) with the following packages (prebuilt binaries from the conda-forge repository):

```bash
mamba create -n rspatial -c conda-forge -y r-base=4.3 r-seurat r-arrow r-hdf5r r-tidyverse r-patchwork r-rann r-remotes c-compiler cxx-compiler gfortran make

conda activate rspatial
hash -r # Wipes the shell memory of the last R path used, so that it finds the conda one.
```

Seurat check:

```bash
R
```

Then, in R:

```r
library(Seurat)
```

= OK.
