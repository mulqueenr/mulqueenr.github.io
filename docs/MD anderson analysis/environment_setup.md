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

#seqkit
conda install -c bioconda seqkit

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
install.packages(c("Seurat","Signac","harmony","vcfR")) #biocmanager install is necessary for signac
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

Installing schicluster
https://github.com/zhoujt1994/scHiCluster
```bash
conda create -n schicluster python==3.6.8
conda activate schicluster
pip install git+https://github.com/zhoujt1994/scHiCluster.git
https://zhoujt1994.github.io/scHiCluster/intro.html
```
# Testing of SV detectors

Installing EagleC
https://github.com/XiaoTaoWang/EagleC

Additionally installing NeoLoopFinder in the same EagleC environment.
https://github.com/XiaoTaoWang/NeoLoopFinder

```bash
conda install mamba -n base -c conda-forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba create -n EagleC scikit-learn statsmodels matplotlib cooler pyBigWig pyensembl python=3.8 joblib=1.0.1 tensorflow=2 cython=0.29.24

conda activate EagleC
mamba install cooler matplotlib pyensembl pybigwig intervaltree scikit-learn=1.1.2 joblib=1.1.0 rpy2 r-mgcv

mamba install -c bioconda bedtools

pip install eaglec
pip install numpy==1.21
download-pretrained-models

#mamba install -c anaconda pomegranate=0.14.4
pip install pomegranate==0.14.4
pip install -U neoloop TADLib
pip install cooltools
conda install sniffles=2.2 #for ONT data

#quick start testing of eaglec
cd ~/ref
wget -O SKNAS-MboI-allReps-filtered.mcool -L https://www.dropbox.com/s/f80bgn11d7wfgq8/SKNAS-MboI-allReps-filtered.mcool?dl=0
predictSV --hic-5k SKNAS-MboI-allReps-filtered.mcool::/resolutions/5000 \
            --hic-10k SKNAS-MboI-allReps-filtered.mcool::/resolutions/10000 \
            --hic-50k SKNAS-MboI-allReps-filtered.mcool::/resolutions/50000 \
            -O SK-N-AS -g hg38 --balance-type CNV --output-format full \
            --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

#testing install of neoloopfinder
cd ~/ref
wget -O SKNMC-MboI-allReps-filtered.mcool -L https://www.dropbox.com/s/tuhhrecipkp1u8k/SKNMC-MboI-allReps-filtered.mcool?dl=0
calculate-cnv -H SKNMC-MboI-allReps-filtered.mcool::resolutions/25000 -g hg38 \
                -e MboI --output SKNMC_25k.CNV-profile.bedGraph

segment-cnv --cnv-file ~/ref/SKNMC_25k.CNV-profile.bedGraph --binsize 25000 \
              --ploidy 2 --output ~/ref/SKNMC_25k.CNV-seg.bedGraph --nproc 4
```

Installing bedGraphtoBigWig in EagleC environment (for plotting)

```bash
cd ~/tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedGraphToBigWig
chmod a+x bedGraphToBigWig
#add to path
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

Installation of scDNA Replication tools
https://github.com/shahcompbio/scdna_replication_tools/tree/main
```bash
git clone git@github.com:shahcompbio/scdna_replication_tools.git
cd scdna_replication_tools-main
conda create -n scdna_replication_tools python==3.7.4
conda activate scdna_replication_tools
python -m venv venv/
source venv/bin/activate
pip install numpy cython
pip install -r requirements3.txt
python setup.py develop

https://luminousmen.com/post/resolve-cython-and-numpy-dependencies

import numpy
from Cython.Build import cythonize
from setuptools import setup, Extension

setup(
    ...
    setup_requires=[
        'setuptools>=18.0',
        'cython',
    ],
    ext_modules=[
        Extension('package.cython_code1', sources=['package/cython_code1.pyx']),
        Extension('package.cython_code2', sources=['package/cython_code2.pyx']),
    ],
    include_dirs=[numpy.get_include()],
)
```

## Installing Singularity, Nextflow and ONT Analysis Pipeline
Necessary tools for ONT validation of structural variant calls

```bash
#install singularity (conda)
conda install -c conda-forge singularity



#initialize epi2me
./nextflow run epi2me-labs/wf-human-variation --help

#test
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --snp --sv \
    --bam demo_data/demo.bam \
    --bed demo_data/demo.bed \
    --ref demo_data/demo.fasta \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
    --sample_name MY_SAMPLE \
    --out_dir ${OUTPUT} \
    -with-singularity \
    -without-docker

#success!
```

# run on local terminal but use 10.132.80.157 files and directories

```zsh
bash

#install dorado (prebuilt binary)
"""
mkdir ~/src
cd ~/src
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz
tar -xvf dorado-0.3.4-linux-x64.tar.gz
#add to path in ~/.bashrc
#PATH=$PATH:/volumes/USR2/Ryan/tools/dorado-0.3.4-linux-x64/bin
#download dorado models
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 #5khz
"""

#install nextflow
"""
cd ~/
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.6-amzn
curl -s https://get.nextflow.io | bash
./nextflow self-update
mv ~/nextflow ~/tools #moving to in PATH
"""


#install docker from https://www.docker.com/products/docker-desktop/ using apple chip
"""
#run docker and it should be in the path
docker version
"""
#connect to server both USR2 and seq
ref="/Volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
wd_out="/Volumes/seq/projects/gccACT/230808_mdamb231_ONT"
output_name="20230726_1239_2D_PAO38369_output" #change to each flowcell
pod5_dir="/Volumes/seq/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231/20230726_1239_2D_PAO38369_dde6ac95" #change to each flowcell

#connect through Finder to seq and USR2 (for reference genome.fa)

nextflow run epi2me-labs/wf-basecalling \
    --input $pod5_dir \
    --dorado_ext "pod5" \
    --output_bam \
    --cuda-device cuda:all \
    --ref $ref \
    --verbose \
    --out_dir ${wd_out}/${output_name} \
    --basecaller_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.2.0" \
    --remora_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2"


```
trying Higashi for scHiC
```bash
cd ~/tools 
git clone https://github.com/ma-compbio/Higashi/
cd Higashi
python setup.py install

ssh seadragon
bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=1:gmem=4 -M 16 -R rusage[mem=16] /bin/bash #get interactive gpu node


bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
module load miniconda3/39_23.5.0; eval "$(/risapps/rhel7/miniconda3/py39_4.12.0/bin/conda shell.bash hook)"
module load cuda11.5/toolkit/11.5.1


conda create --name hic python=3.11 #installing to conda base
#install higashi from 4dn
conda activate hic
conda install -c conda-forge mamba
mamba install pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 cudatoolkit=11.3 -c pytorch
mamba install -c bioconda pybedtools cooler
conda install -c bioconda pybedtools

mamba install -c conda-forge zlib

mamba install -c ruochiz fasthigashi
# CUDA 10.2 (version of cuda in modules)
conda install pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 cudatoolkit=10.2 -c pytorch

mkdir ~/tools
cd ~/tools 
git clone https://github.com/ma-compbio/Higashi/
cd Higashi
conda install python==3.9 numpy==1.16.5
python setup.py install
```


conda create --name hic 
conda activate hic
mamba install -c ruochiz fasthigashi
    q

wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz


Download data via transfer node
```bash
bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

#ONT data
rsync  \
-LPr mulqueen@10.132.80.157:/volumes/seq/projects/gccACT/230808_mdamb231_ONT ~/projects/gccACT

#dorado prebuilt
rsync \
-LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/tools/dorado-0.3.4-linux-x64.tar.gz ~/tools
#dorado reference genome
#download model for base calling
#wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 #5khz #cpg??
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0
rsync \
-LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ~/
rsync \
-LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 ~/

#nextflow epi2me download
rsync \
-LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/wf-human-variation-master ~/

#my 10.132.80.157 references
rsync \
-LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/ref ~/


```

Use CUDA node
```bash
bsub -Is -W 6:00 -q gpu-medium -n 1 -gpu num=2:gmem=4 -M 16 -R rusage[mem=16] /bin/bash #get interactive gpu node

module load singularity/3.7.0
module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1
module load samtools/1.15 
#module load dorado
#module load wf-human-variation

OUTPUT=output
nextflow run ./wf-human-variation-master/main.nf \
    -w ${OUTPUT}/workspace \
    --snp --sv --bam demo_data/demo.bam \
    --bed demo_data/demo.bed \
    --ref demo_data/demo.fasta \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0' \
    --sample_name MY_SAMPLE \
    --out_dir ${OUTPUT} \
    -with-singularity -without-docker


# If you want to start with a FASTQ, you can generate an unaligned BAM from FASTQ with samtools import http://www.htslib.org/doc/samtools-import.html and allow the workflow to take care of mapping.

#CPU Node Run
bsub -Is -W 6:00 -q medium -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive cpu node

#module load singularity/3.7.0
#module load nextflow/23.04.3
#module load cuda11.5/toolkit/11.5.1
module load samtools/1.15 
#module load dorado

ref="/rsrch4/home/genetics/rmulqueen/ref/genome.fa"
wd_out="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT"
output_name="20230726_1239_2D_PAO38369_output" #change to each flowcell
pod5_dir="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231/20230726_1239_2D_PAO38369_dde6ac95" #change to each flowcell

dorado basecaller \
    --verbose \
    --device cpu \
    --reference ${ref} \
    --emit-sam \
    --max-reads 100 \
    --batchsize 64 \
    'dna_r10.4.1_e8.2_400bps_hac@v4.2.0' \
    ${pod5_dir}/pod5_pass/ \

#untested for this part, but use bam from wf-basecalling output as input
nextflow run ~/wf-human-variation-master/main.nf \
    -w ${wd_out}/${output_name}/workspace \
    -profile singularity \
    --snp --sv --cnv --methyl \
    --ref ${ref} \
    --bam ${wd_out}/${output_name}.bam \
    --dorado_ext pod5 \
    --basecaller_basemod_threads 40 \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2' \
    --sample_name ${output_name} \
    --out_dir ${wd_out}/${output_name}/ \
    -with-singularity \
    -without-docker

```



```bash
ssh r1prpsciapp13 
module purge
module load nextflow/23.04.3
module load singularity/3.7.0

#manually pull singularity image?
nextflow pull epi2me-labs/wf-human-variation
singularity pull --name ontresearch-wf-human-variation-shac4db03c19b6ff1277a24ec28a19e564d628d478f.img.pulling.1669977561040 docker://ontresearch/wf-human-variation:shac4db03c19b6ff1277a24ec28a19e564d628d478f
 > /dev/null

sif_in="ontresearch-wf-human-variation-shac4db03c19b6ff1277a24ec28a19e564d628d478f.img.pulling.1669977561040"
nextflow run ~/wf-human-variation-master/main.nf 
OUTPUT=output
nextflow run ~/wf-human-variation-master/main.nf \
     -w ${OUTPUT}/workspace \
     -profile standard \
     --snp --sv --bam demo_data/demo.bam \
     --bed demo_data/demo.bed \
     --ref demo_data/demo.fasta \
     --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0' \
     --sample_name MY_SAMPLE --out_dir ${OUTPUT} \
     -without-docker -with-singularity $sif_in

tar -xzf /risapps/tps_source/tps_source/wf-human-variation/demo_data.tar.gz
OUTPUT=output
nextflow run wf-human-variation-master -w ${OUTPUT}/workspace -profile standard --snp --sv --bam demo_data/demo.bam --bed demo_data/demo.bed --ref demo_data/demo.fasta --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0' --sample_name MY_SAMPLE --out_dir ${OUTPUT} -with-singularity -without-docker
```


bsub -q "gpu-medium" dorado_test.bsub