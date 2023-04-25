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

## General tools
### FastQ Generation:
bcl2fastq2

### Alignment:
bwa
bowtie2
hisat2
bismark (methylation)

### Fastq/sam/bam manipulation
samtools

### Bed file manipulation
bedtools

### ATAC peak calling
macs3

```bash
#bcl2fastq2
conda install -c bih-cubi bcl2fastq2

#bwa
cd ~/tools
git clone https://github.com/lh3/bwa.git
cd bwa; make

#bowtie2
cd ~/tools
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip

#hisat2
conda install -c bioconda hisat2

#bismark
cd ~/tools
wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.24.0.tar.gz

#samtools
conda install -c bioconda samtools

#bedtools
conda install -c bioconda bedtools

#macs3
pip install macs3

#Add to PATH (for those not conda installed)
echo "PATH=$PATH:~/tools/bowtie2-2.5.1-linux-x86_64" >> ~/.bashrc
echo "PATH=$PATH:~/tools/bwa" >> ~/.bashrc
echo "PATH=$PATH:~/tools/Bismark-0.24.0/"  >> ~/.bashrc

#multiqc
pip install multiqc
```

### General RNA and ATAC processing can be done with Seurat and Signac in R

Now install all the R packages.

Install package dependencies in the environment
```bash
conda install -c conda-forge r-xml
conda install -c conda-forge r-gert
conda install -c conda-forge r-ragg
conda install -c conda-forge r-spdep
conda install -c conda-forge r-terra
conda install -c conda-forge zlib
```

I'm sure there are a ton that I'm missing. But I'm getting started with these. 

```R
install.packages(c("BiocManager","devtools")) 
install.packages(c("Seurat","Signac","harmony")) #biocmanager install is necessary for signac
BiocManager::install(c("cicero","ComplexHeatmap"))
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("navinlabcode/copykit")
```

### HiC Data Analysis can be done with the DipC group's released hickit
https://github.com/4dn-dcic/pairix#installation-for-pairix
Pairix download and installation
```bash
cd ~/tools
git clone https://github.com/4dn-dcic/pairix
cd pairix
make
#and then add to path
```
Installing cooler on alternative conda env

```bash
conda deactivate #get out of r3.4 env
conda create -n "cooler_env" python=3.9.15
conda activate cooler_env
pip install cooler
pip install Cython
pip install cooltools #need to downgrate python for this
pip install pypairix
pip install plotly seaborn
```


# Testing of SV detectors

Installing EagleC
https://github.com/XiaoTaoWang/EagleC
```bash
conda install mamba -n base -c conda-forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba create -n EagleC scikit-learn statsmodels matplotlib cooler pyBigWig pyensembl python=3.8 joblib=1.0.1 tensorflow=2 cython=0.29.24

conda activate EagleC
pip install eaglec
download-pretrained-models
```


Installing HiSV in cooler_env
https://github.com/GaoLabXDU/HiSV
```bash
conda activate cooler_env
conda install gcc
pip install numpy pysam pandas prox_tv
#prox_tv was annoying but installed from git source
git clone https://github.com/gaolabXDU/HiSV
```


Installing HiC_breakfinder in cooler_env
https://github.com/dixonlab/hic_breakfinder
```bash
cd ~/tools
git clone https://gitlab.com/libeigen/eigen.git #https://eigen.tuxfamily.org/dox/GettingStarted.html
git clone https://github.com/pezmaster31/bamtools.git #https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
conda install -c conda-forge jsoncpp

git clone https://github.com/dixonlab/hic_breakfinder.git

./configure --prefix=/volumes/USR2/Ryan/tools/hic_breakfinder/ CPPFLAGS="-I /volumes/USR2/Ryan/tools/bamtools/include -I /volumes/USR2/Ryan/tools/eigen" LDFLAGS="-L/volumes/USR2/Ryan/tools/bamtools/lib/"

git clone https://github.com/dixonlab/hic_breakfinder

 conda install -c bioconda bamtools
 conda install -c conda-forge eigen
./configure CPPFLAGS="-I /volumes/USR2/Ryan/miniconda3/envs/cooler_env/include/bamtools -I /volumes/USR2/Ryan/miniconda3/envs/cooler_env/include/eigen3/" LDFLAGS="-L/volumes/USR2/Ryan/miniconda3/envs/cooler_env/lib/bamtools"
make

```
<!--
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
-->

### Samblaster for eccDNA
Installation of samblaster
https://github.com/GregoryFaust/samblaster
```bash
cd ~/tools
wget https://github.com/GregoryFaust/samblaster/archive/refs/heads/master.zip
unzip master.zip
cd samblaster-master
make
```

### Circle-finder for eccDNA
https://github.com/pk7zuva/Circle_finder
```bash
cd ~/tools
mkdir circle-finder
wget https://github.com/pk7zuva/Circle_finder/blob/3eb333db2ea6277dde36cbf640be9afeb710c717/circle_finder-pipeline-bwa-mem-samblaster.sh
```

# Separate Conda Env: scbs_env
Need a separate environment for scbs analysis because python version has to be lower.

### Methylation Analysis with scbs
https://github.com/LKremer/scbs
https://www.bioconductor.org/packages/release/bioc/vignettes/Melissa/inst/doc/process_files.html

Naming a scbs conda environment as scbs_env
```bash
conda deactivate #get out of r3.4 env
conda create -n "scbs_env" python=3.9.15
conda activate scbs_env
python3 -m pip install scbs
conda install -c conda-forge r-xml r-gert r-ragg r-spdep r-terra r-stringi
scbs --version                                                                                    
#scbs, version 0.5.4            
```

Install R packages for downstream analysis.
```R
install.packages(c("BiocManager","devtools")) 
BiocManager::install("Melissa")
```

# References Download
Downloading mm10 and hg38 reference genomes. Using the cellranger ones for consistency with commercial products.

```bash
mkdir ~/ref
cd ~/ref
wget https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

tar -xvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -xvf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

```

## Prepare Reference Genomes For Aligners
Skipping mouse genome for now, but its similar
```bash
bwa index ~/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa &

bismark_genome_preparation ~/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta &
```