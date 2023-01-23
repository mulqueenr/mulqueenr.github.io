---
title: Environment Setup
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /mda_environment/
category: mdanderson
---

Environment Setup for Navin Lab compute cluster. Some installs may be redundant with the cluster...

Login
```bash
ssh mulqueen@10.132.80.157
pwd
#/volumes/USR2/Ryan
```
Make directories for src (usually my own code) and tools (other peoples analysis pipelines). 

```bash
mkdir tools
mkdir src
```
Download and install conda
```bash
cd tools
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_22.11.1-1-Linux-x86_64.sh
bash Miniconda3-py310_22.11.1-1-Linux-x86_64.sh
```
Conda updates the .bashrc, so next time you log in it is properly in your PATH.

Logout and login to activate conda
```
conda create -n r4.2 
conda activate r4.2
conda install -c conda-forge r-base scipy numpy
#add conda environment to .bashrc
echo "conda activate r4.2" >> ~/.bashrc
```
Log out and login to automatic load in of conda environment.

### General RNA and ATAC processing can be done with Seurat and Signac

Now install all the R packages.
```R
install.packages("Seurat","Signac","ComplexHeatmap","devtools","BiocManager")

BiocManager::install("cicero")
devtools::install_github('cole-trapnell-lab/monocle3')

```

### HiC Data Analysis can be done with the DipC group's released hickit
https://github.com/lh3/hickit

```bash
cd ~/tools

# Download precompiled binaries for Linux
curl -L https://github.com/lh3/hickit/releases/download/v0.1/hickit-0.1_x64-linux.tar.bz2 | tar -jxf -
cd hickit-0.1_x64-linux

# Map Dip-C reads and extract contacts (skip if you use your own pipeline)
./seqtk mergepe read1.fq.gz read2.fq.gz | ./pre-dip-c - | bwa mem -5SP -p hs37d5.fa - | gzip > aln.sam.gz
./k8 hickit.js vcf2tsv phased.vcf > phased_SNP.tsv   # extract phased SNPs from VCF
./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ./k8 hickit.js chronly - | ./k8 hickit.js bedflt par.bed - | gzip > contacts.seg.gz # for male
#./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ./k8 hickit.js chronly -y - | gzip > contacts.seg.gz # for female
./hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz  # optional

# Impute phases (-i also works with contacts.seg.gz)
./hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz
./hickit -i contacts.pairs.gz --out-val=impute.val     # estimate imputation accuracy by holdout
# Infer 3D structure
./hickit -i impute.pairs.gz -Sr1m -c1 -r10m -c5 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O imput.3dg

# 2D contact map in PNG (bin size determined by the image width)
./hickit -i impute.pairs.gz --out-png impute.png
# Compute CpG density (optional)
./hickit.js gfeat -r hs37d5.fa.gz imput.3dg | gzip > imput.cpg.3dg.gz
# Visualize 3D structure (requiring a graphical card)
./hickit-gl -I imput.cpg.3dg.gz --view
```

### Methylation Analysis with scbs
https://github.com/LKremer/scbs
https://www.bioconductor.org/packages/release/bioc/vignettes/Melissa/inst/doc/process_files.html

```bash
python3 -m pip install scbs
```

```R
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("Melissa")

```
## File Location