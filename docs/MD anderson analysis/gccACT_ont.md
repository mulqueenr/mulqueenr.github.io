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

First testing pipeline on test data. This is all written for seadragon.

Basecalling: https://github.com/nanoporetech/dorado
Methylation calling: https://github.com/nanoporetech/remora
small variants: https://www.github.com/HKU-BAL/Clair3
QDNA seq CNV caller: https://bioconductor.org/packages/release/bioc/html/QDNAseq.html

Now our data:

## Set up seadragon for data processing
Download data via transfer node
```bash
ssh seadragon

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

#ONT data
#rsync  -LPr mulqueen@10.132.80.157:/volumes/seq/projects/gccACT/230808_mdamb231_ONT ~/projects/gccACT
#just focusing on the big pod5 files
rsync -LPr mulqueen@10.132.80.157:/volumes/seq/projects/gccACT/230808_mdamb231_ONT/MDA_MB_231_2/MDA_MB_231/20230802_1920_2D_PAO38925_a09c109d/pod5_pass  ~/projects/gccACT

#dorado prebuilt
#rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/tools/dorado-0.3.4-linux-x64.tar.gz ~/tools
#or
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.4-linux-x64.tar.gz

#epi2me test data
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz
tar -xzvf demo_data.tar.gz

#dorado reference genome
#download model for base calling
#ran on navin 10.132.80.157 cluster first and pulled to seadragon, can also be done on transfer node
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ~/
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 ~/

#nextflow epi2me download (pulled from github originally)
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/wf-human-variation-master ~/

#my 10.132.80.157 references (used just for genome.fa of hg38 which i pulled from the 10x website for consistency)
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/ref ~/

cd ~/singularity
module load singularity/3.7.0
export NXF_SINGULARITY_CACHEDIR="~/singularity/"

#manual pull of singularity containers so i can run on gpu nodes (taking these from output log of test data ran on seadragon transfer node to see what docker containers it was pulling.) I'm not sure if this step is necessary anymore, since i set env variables to tell it not to pull these in the job submissions

singularity pull docker://ontresearch/wf-human-variation:sha0337567399c09ef14d1ab9cc114f77de86398e12 

#cnv
singularity pull docker://ontresearch/wf-cnv:sha428cb19e51370020ccf29ec2af4eead44c6a17c2 

#sv
singularity pull docker://ontresearch/wf-human-variation-sv:shabc3ac908a14705f248cdf49f218956ec33e93ef9 

#snp
singularity pull docker://ontresearch/wf-human-variation-snp:sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e 
#methyl
singularity pull docker://ontresearch/wf-human-variation-methyl:sha44a13bcf48db332b2277bb9f95b56d64e393a1d5 > /dev/null

singularity pull  --name ontresearch-wf-human-variation-methyl-sha44a13bcf48db332b2277bb9f95b56d64e393a1d5.img.pulling.1695478159368 docker://ontresearch/wf-human-variation-methyl:sha44a13bcf48db332b2277bb9f95b56d64e393a1d5 > /dev/null

#all together at once
module load nextflow/23.04.3

nextflow download /home/rmulqueen/wf-human-variation-master/main.nf --container singularity
```

## Running ONT nextflow pipeline.

Local install and run using GPUs on seadragon

Note for some, the pod5 files are so big (600gb) that they need to be uploaded and cleared from seadragon since the disk quota is 1TB.

Written as interactive node, but can be formatted for bsub job submisison as well. Run in a screen so you don't have to keep an active terminal through it.

20230726_1239_2D_PAO38369_dde6ac95.nextflow.bsub

```bash
#BSUB -J 20230726_1239_2D_PAO38369_dde6ac95_nextflow
#BSUB -W 12:00
#BSUB -o /rsrch4/home/genetics/rmulqueen/
#BSUB -e /rsrch4/home/genetics/rmulqueen/
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/
#BSUB -q gpu
#BSUB -gpu num=4:gmem=16 
#BSUB -M 160
#BSUB -R "rusage[mem=160]"
#BSUB -B
#BSUB -N
#BSUB -u rmulqueen@mdanderson.org


#bsub -Is -W 4:00 -q gpu -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

module load nextflow/23.04.3
module load cuda11.5/toolkit/11.5.1
module load singularity/3.7.0
module load samtools/1.15 

ref="/rsrch4/home/genetics/rmulqueen/ref/genome.fa"
wd_out="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT"
output_name="20230726_1239_2D_PAO38369_output" #change to each flowcell

pod5_dir="/rsrch4/home/genetics/rmulqueen/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_dde6ac95" #change to each flowcell
#download model for base calling
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 #5khz #cpg??
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0

#these make it run locally, and use the singularity containers we pulled manually above
export NXF_SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity/"
export SINGULARITY_CACHEDIR="/rsrch4/home/genetics/rmulqueen/singularity/"
#mkdir $SINGULARITY_CACHEDIR #make sure these directories are made
#mkdir $SINGULARITY_CACHEDIR/tmp
#mkdir $SINGULARITY_CACHEDIR/pull
export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
export SINGULARITY_PULLDIR=$SINGULARITY_CACHEDIR/pull
export CWL_SINGULARITY_CACHE=$SINGULARITY_PULLDIR
export NXF_OFFLINE='TRUE' #https://nf-co.re/docs/usage/offline

#dorado run (succeeded previously so commenting out here)
#output bam file from dorado caller has to be sorted before it can be used in the pipeline.
#~/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
#    --verbose \
#    --reference ${ref} \
#    --emit-sam \
#    --modified-bases-models dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 \
#    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
#    ${pod5_dir}/pod5_pass/ | samtools sort -@ 10 -T $HOME | samtools view -b - > ${wd_out}/${output_name}.sorted.bam

nextflow run /home/rmulqueen/wf-human-variation-master/main.nf \
    -w ${wd_out}/${output_name}/workspace \
    -profile singularity \
    --snp --sv --cnv \
    --ref ${ref} \
    --bam ${wd_out}/${output_name}.sorted.bam \
    --dorado_ext pod5 \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2' \
    --sample_name ${output_name} \
    --out_dir ${wd_out}/${output_name}/ \
    -with-singularity \
    -without-docker \
    -offline

#add methyl in future
#note sif files may need to be manually pulled (as above) for updates to wf-human-variation-master in the future
```

```bash
bsub < 20230726_1239_2D_PAO38369_dde6ac95.nextflow.bsub
```

Transfer back to navin cluster with results.

```bash
ssh seadragon
screen

bsub -Is -W 4:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up

sftp mulqueen@10.132.80.157
cd /volumes/seq/projects/gccACT/230808_mdamb231_ONT 
put -R ~/projects/gccACT/230808_mdamb231_ONT
```