---
title: gccACT ONT Analysis
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /gccACT_ONT/
category: mdanderson
---

Generated multiple flowcell lanes of data on MDA-MB-231 cell line to check SVs detected in HiC data.

nextflow analysis pipeline for nanopore
https://github.com/epi2me-labs/wf-human-variation

Using ONT software for analysis of data.

First testing pipeline on test data

```bash
ssh r1prpsciapp13
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz
tar -xzvf demo_data.tar.gz
```

Basecalling: https://github.com/nanoporetech/dorado
Methylation calling: https://github.com/nanoporetech/remora
small variants: https://www.github.com/HKU-BAL/Clair3
QDNA seq CNV caller: https://bioconductor.org/packages/release/bioc/html/QDNAseq.html

Now our data:

## Set up seadragon for data processing
Download data via transfer node
```bash
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

#ONT data
rsync  -LPr mulqueen@10.132.80.157:/volumes/seq/projects/gccACT/230808_mdamb231_ONT ~/projects/gccACT

#dorado prebuilt
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/tools/dorado-0.3.4-linux-x64.tar.gz ~/tools
#dorado reference genome
#download model for base calling
#wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 #5khz #cpg??
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ~/
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 ~/

#nextflow epi2me download
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/wf-human-variation-master ~/

#my 10.132.80.157 references
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/ref ~/


```

```bash
## Dorado basecalling and alignment

This one works! 
20230731_1851_3E_PAO38479_822d79b2.dorado.bsub
```bash
#BSUB -J 20230731_1851_3E_PAO38479_822d79b2
#BSUB -W 6:00
#BSUB -o /rsrch4/home/genetics/rmulqueen/
#BSUB -e /rsrch4/home/genetics/rmulqueen/
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/
#BSUB -q gpu-medium
#BSUB -gpu num=2:gmem=4 
#BSUB -M 160
#BSUB -R "rusage[mem=160]"
#BSUB -B
#BSUB -N
#BSUB -u rmulqueen@mdanderson.org

#bsub -Is -W 4:00 -q short -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

pwd
module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1
module load samtools/1.15 
echo $(hostname)

ref="/rsrch4/home/genetics/rmulqueen/ref/genome.fa"
wd_out="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT"
output_name="20230731_1851_3E_PAO38479_822d79b2_output" #change to each flowcell
pod5_dir="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231/20230731_1851_3E_PAO38479_822d79b2" #change to each flowcell

~/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
    --verbose \
    --reference ${ref} \
    --emit-sam \
    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
    ${pod5_dir}/pod5_pass/ | samtools view -b - > ${wd_out}/${output_name}.bam

```

20230802_1920_2D_PAO38925_a09c109d.dorado.bsub
```bash
#BSUB -J 20230802_1920_2D_PAO38925_a09c109d
#BSUB -W 6:00
#BSUB -o /rsrch4/home/genetics/rmulqueen/
#BSUB -e /rsrch4/home/genetics/rmulqueen/
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/
#BSUB -q gpu-medium
#BSUB -gpu num=2:gmem=4 
#BSUB -M 160
#BSUB -R "rusage[mem=160]"
#BSUB -B
#BSUB -N
#BSUB -u rmulqueen@mdanderson.org

#bsub -Is -W 4:00 -q short -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

pwd
module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1
module load samtools/1.15 
echo $(hostname)

ref="/rsrch4/home/genetics/rmulqueen/ref/genome.fa"
wd_out="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT"
output_name="20230802_1920_2D_PAO38925_a09c109d_output" #change to each flowcell
pod5_dir="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231_2/MDA_MB_231/20230802_1920_2D_PAO38925_a09c109d" #change to each flowcell

~/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
    --verbose \
    --reference ${ref} \
    --emit-sam \
    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
    ${pod5_dir}/pod5_pass/ | samtools view -b - > ${wd_out}/${output_name}.bam

```

20230726_1239_2D_PAO38369_dde6ac95.dorado.bsub
```bash
#BSUB -J 20230726_1239_2D_PAO38369_dde6ac95
#BSUB -W 6:00
#BSUB -o /rsrch4/home/genetics/rmulqueen/
#BSUB -e /rsrch4/home/genetics/rmulqueen/
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/
#BSUB -q gpu-medium
#BSUB -gpu num=2:gmem=4 
#BSUB -M 160
#BSUB -R "rusage[mem=160]"
#BSUB -B
#BSUB -N
#BSUB -u rmulqueen@mdanderson.org

#bsub -Is -W 4:00 -q short -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

pwd
module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1
module load samtools/1.15 
echo $(hostname)

ref="/rsrch4/home/genetics/rmulqueen/ref/genome.fa"
wd_out="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT"
output_name="20230726_1239_2D_PAO38369_output" #change to each flowcell
pod5_dir="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231/20230726_1239_2D_PAO38369_dde6ac95" #change to each flowcell

~/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
    --verbose \
    --reference ${ref} \
    --emit-sam \
    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
    ${pod5_dir}/pod5_pass/ | samtools view -b - > ${wd_out}/${output_name}.bam

```






Local install
```bash
#bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=2:gmem=4 -R rusage[mem=160] /bin/bash #get interactive gpu node

module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1

ref="/Volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
wd_out="/Volumes/seq/projects/gccACT/230808_mdamb231_ONT"
output_name="20230726_1239_2D_PAO38369_output" #change to each flowcell
pod5_dir="/Volumes/seq/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231/20230726_1239_2D_PAO38369_dde6ac95" #change to each flowcell

#download model for base calling
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 #5khz #cpg??
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0

#untested for this part, but use bam from wf-basecalling output as input
nextflow run epi2me-labs/wf-human-variation \
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

<!--

#connect through Finder to seq and USR2 (for reference genome.fa)

# cd ~
# nextflow run epi2me-labs/wf-basecalling \
#     -w ${wd_out}/${output_name}/workspace \
#     --input $pod5_dir \
#     --dorado_ext pod5 \
#     --ref $ref \
#     --out_dir ${wd_out}/${output_name} \
#     --basecaller_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.2.0" \
#     --basecaller_basemod_threads 80 \
#     --remora_cfg "dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2"


```

```bash
ssh seadragon
sftp mulqueen@10.132.80.157
lcd ~/projects/gccACT #seadragon directory
get -R /volumes/seq/projects/gccACT/230808_mdamb231_ONT/ #download ONT data
get -R ~/wf-human-variation #downloaded 

#transfer
bqueues
bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=1:gmem=4 -M 16 -R rusage[mem=16] /bin/bash #get interactive gpu node

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive gpu node


#test demo data
wget -O demo_data.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz #download
tar -xzvf demo_data.tar.gz #extract


module load singularity/3.5.2 #load singularity
#module load cuda10.1/toolkit/10.1.243 #load cuda
module load nextflow/22.10.6 #load nextflow
#install nextflow instead of using module
bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive gpu node
cd ~/

#on 10.132.80.157  server do 
curl -s "https://get.sdkman.io" | bash
curl -s https://get.nextflow.io | bash
#transfer to seadragon

sftp mulqueen@10.132.80.157 #transfer sdk to home directory
get -R .sdkman
get -R .nextflow
get nextflow
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.6-amzn
#transfer nextflow to home directory
./nextflow self-update
mv ~/nextflow ~/tools #moving to in PATH


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
```






Our servers don't have GPUs so I'm just going to use my macbook
https://labs.epi2me.io/downloads/
and install docker desktop
https://www.docker.com/products/docker-desktop/

installing nextflow
curl -s "https://get.sdkman.io" | bash
source "/Users/rmulqueen/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.6-tem
wget -qO- https://get.nextflow.io | bash

-->