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
#ran on navin 10.132.80.157 cluster first and pulled to seadragon
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 
#dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.2.0
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ~/
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 ~/

#nextflow epi2me download
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/wf-human-variation-master ~/

#my 10.132.80.157 references
rsync -LPr mulqueen@10.132.80.157:/volumes/USR2/Ryan/ref ~/

mkdir ~/singularity
export NXF_SINGULARITY_CACHEDIR="~/singularity/"

cd ~/singularity
module load singularity/3.7.0

#manual pull of singularity containers so i can run on gpu nodes (taking these from output log of test data ran on seadragon transfer node to see what docker containers it was pulling.)
singularity pull docker://ontresearch/wf-human-variation-sv:shabc3ac908a14705f248cdf49f218956ec33e93ef9 
singularity pull docker://ontresearch/wf-human-variation:sha0337567399c09ef14d1ab9cc114f77de86398e12 
singularity pull docker://ontresearch/wf-cnv:sha428cb19e51370020ccf29ec2af4eead44c6a17c2 
singularity pull docker://ontresearch/wf-human-variation-snp:sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e 

```

## Dorado basecalling and alignment
Run these scripts with command:
```bash
bsub < 20230731_1851_3E_PAO38479_822d79b2.dorado.bsub 
```
or use the commented out call for an interative gpu node and run line by line.

20230731_1851_3E_PAO38479_822d79b2.dorado.bsub
```bash
#BSUB -J 20230731_1851_3E_PAO38479_822d79b2
#BSUB -W 12:00
#BSUB -o /rsrch4/home/genetics/rmulqueen/
#BSUB -e /rsrch4/home/genetics/rmulqueen/
#BSUB -cwd /rsrch4/home/genetics/rmulqueen/
#BSUB -q gpu-medium
#BSUB -gpu num=4:gmem=16 
#BSUB -M 160
#BSUB -R "rusage[mem=160]"
#BSUB -B
#BSUB -N
#BSUB -u rmulqueen@mdanderson.org

#bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

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
<!--
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

#bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

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

#bsub -Is -W 4:00 -q gpu-medium -n 1 -gpu num=2:gmem=4 -M 160 -R rusage[mem=160] /bin/bash #get interactive gpu node

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
    --modified-bases-models dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 \
    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
    ${pod5_dir}/pod5_pass/ | samtools view -b - > ${wd_out}/${output_name}.bam

```
-->
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
#BSUB -q gpu-medium
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

#dorado run
#output bam file from dorado caller has to be sorted before it can be used in the pipeline.
~/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
    --verbose \
    --reference ${ref} \
    --emit-sam \
    --modified-bases-models dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2 \
    dna_r10.4.1_e8.2_400bps_hac@v4.2.0 \
    ${pod5_dir}/pod5_pass/ | samtools sort -@ 10 -T $HOME | samtools view -b - > ${wd_out}/${output_name}.sorted.bam

nextflow run /home/rmulqueen/wf-human-variation-master/main.nf \
    -w ${wd_out}/${output_name}/workspace \
    -profile singularity \
    --snp --sv --cnv --methyl \
    --ref ${ref} \
    --bam ${wd_out}/${output_name}.sorted.bam \
    --dorado_ext pod5 \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2' \
    --sample_name ${output_name} \
    --out_dir ${wd_out}/${output_name}/ \
    -with-singularity \
    -without-docker


```

```bash
bsub < 20230726_1239_2D_PAO38369_dde6ac95.nextflow.bsub
```