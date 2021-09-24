---
title: HMMcopy on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /hmmcopy_nmt/
category: CEDAR
---

# Script to test Shah lab CNV calling using HMMcopy on SE bismark aligned reads. Starting at split fastq files for alignment to hg38.

## Source data
```bash
/home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scNMT_NOMeWorkFlow/samples/raw/*fastq.gz
/home/groups/CEDAR/doe/projects/my_NMT/new_MCF7/scNMT_NOMeWorkFlow/samples/raw/*fastq.gz
```

### Setting up fastq files, and run QC

Some MCF7 cells were run twice, so I'm going to merge them at the fastq level before alignment. 

```bash
mkdir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams
cd /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams
cp /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scNMT_NOMeWorkFlow/samples/raw/*fastq.gz . & #copy files over


for i in /home/groups/CEDAR/doe/projects/my_NMT/new_MCF7/scNMT_NOMeWorkFlow/samples/raw/*fastq.gz; do outname=`basename $i`; cat $i $outname > merged.${outname}; done & #merge files that were resequenced
for i in merged.*fastq.gz; do unmerged_file=${i:7}; rm $unmerged_file; done & #remove duplicate files and keep merged ones
rename "merged." "" merged*gz #rename to make it consistent
ls *fastq.gz | parallel -j15 fastqc {} & #run 15 parallel instances of fastqc
multiqc -f . &
```

Setting up bismark alignment as single-end since NMT is known to generate chimeric reads. Setting trimming settings for R1 and R2 separately based on MultiQC output.

I decided to trim the first 25 base pairs for both read 1 and read 2.

## Bismark genome preparation

```bash
bismark_genome_preparation \
--path_to_aligner /opt/installed/bowtie2/bin/bowtie2 \
 --verbose \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/
```

## Trimming fastq reads with trim_galore
Using trim_galore, which uses cutadapt for trimming reads.
Running 15 parallel instances.

```bash
ls *fastq.gz | parallel -j15 trim_galore --clip_R1 25 {} & #run 15 parallel instances of fastqc
multiqc -f . & #remake multiqc
```
## Batch script for alignment of trimmed fastq files
Using bismark single-end mode.

```bash
cd /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams
bismark \
--genome /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/ \
--non_directional \
--gzip \
--parallel 8 \
--single_end \
`ls -m *trimmed.fq.gz | tr -d " " | tr -d "\n" ` &

#ls command generates a comma separated list of files
```
## Merge Read 1 and Read 2 files per cell

Merged single end bam files

```bash
for i in `ls *_R1_*bam`;
do R1=$i;
R2=`echo $R1 | sed 's/R1/R2/g'`;
outname=`echo $R1 | sed 's/R1/merged/g'`;
samtools cat -@ 20 -o $outname $R1 $R2; done &
```

## Batch script for deduplication
Using the bismark deduplication script.

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-573
#SBATCH --tasks-per-node=30 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams/*merged_trimmed_bismark_bt2.bam | wc -l`

bam_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams/*merged_trimmed_bismark_bt2.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun deduplicate_bismark \
--single \
--output_dir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams/dedup_bams \
${bam_in}

```

## Perform modal copy number bin normalization similar to Shah sc-wgs work
This takes the input of files that are generated in the hmm_generator function in the R script (files starting with "counts_")

```python
'''
Created on Feb 21, 2018
@author: dgrewal
https://github.com/shahcompbio/single_cell_pipeline/blob/master/single_cell/workflows/hmmcopy/scripts/correct_read_count.py
Modified by mulqueenr to fit HMMcopy correction output for our own pipeline.
'''

from __future__ import division
import argparse
import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats.mstats import mquantiles

import glob
import warnings

'''
Compute quantile regression curves and select the modal quantile.
'''
# 2nd order polynomial quantile regression

def modal_quantile_regression(df_regression, lowess_frac=0.2):
    q_range = range(10, 91, 1)
    quantiles = np.array(q_range) / 100
    quantile_names = [str(x) for x in q_range]
    if len(df_regression) < 10:     #need at least 3 values to compute the quantiles
        return df_regression
    poly2_quantile_model = smf.quantreg('reads ~ gc + I(gc ** 2.0)', data=df_regression)
    poly2_quantile_fit = [poly2_quantile_model.fit(q=q) for q in quantiles]
    poly2_quantile_predict = [poly2_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))]
    poly2_quantile_params = pd.DataFrame()
    for i in range(len(quantiles)):
        df_regression[quantile_names[i]] = poly2_quantile_predict[i]
        poly2_quantile_params[quantile_names[i]] = poly2_quantile_fit[i].params
    gc_min = df_regression['gc'].quantile(q=0.10)     # integration and mode selection
    gc_max = df_regression['gc'].quantile(q=0.90)
    poly2_quantile_integration = np.zeros(len(quantiles)+1)
    for i in range(len(quantiles)):
        params = poly2_quantile_params[quantile_names[i]].tolist()
        params.reverse()
        poly2 = np.poly1d(params)
        integ = poly2.integ()
        integrand = integ(gc_max) - integ(gc_min)
        poly2_quantile_integration[i+1] = integrand
    distances = poly2_quantile_integration[1:] - poly2_quantile_integration[:-1]     # find the modal quantile
    df_dist = pd.DataFrame({'quantiles': quantiles, 'quantile_names': quantile_names, 'distances': distances})
    dist_max = df_dist['distances'].quantile(q=0.95)
    df_dist_filter = df_dist[df_dist['distances']<dist_max]
    df_dist_filter['lowess'] = lowess(df_dist_filter['distances'], df_dist_filter['quantiles'], frac=lowess_frac, return_sorted=False)
    modal_quantile = quantile_names[np.argmin(df_dist_filter['lowess'])]
    df_regression['modal_quantile'] = modal_quantile     # add values to table
    df_regression['modal_curve'] = df_regression[modal_quantile]
    df_regression['modal_corrected'] = df_regression['reads'] / df_regression[modal_quantile]
    return df_regression

warnings.filterwarnings("ignore")

working_dir="/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test"

files=[name for name in glob.glob(working_dir+'/counts*tsv')]

for infile in files:
    df=pd.read_csv(infile,sep="\t")
    df['modal_quantile'] = 'NaN'
    df['modal_curve'] = 'NaN'
    df['modal_corrected'] = 'NaN'
    df_ideal = df[df['ideal']==True] # filtering and sorting
    df_regression = pd.DataFrame.copy(df_ideal)
    df_regression.sort_values(by='gc', inplace=True)
    df_regression = modal_quantile_regression(df_regression=df_regression, lowess_frac=0.2) # modal quantile regression
    df.loc[df_regression.index, 'modal_quantile'] = df_regression['modal_quantile'] # map results back to full data frame
    df.loc[df_regression.index, 'modal_curve'] = df_regression['modal_curve']
    df.loc[df_regression.index, 'modal_corrected'] = df_regression['modal_corrected']
    df['copy'] = df['modal_corrected']
    df = df.rename(columns=({ "modal_corrected" : "cor_gc"}))
    df.to_csv(infile[0:-3]+"modal_corrected.tsv", index=False, sep=',', na_rep="NA")
    print("Completed:",infile)

```

## Preparing bulk WGS data for cell lines to check CNVs
Files are from a previous publication by Hisham (here.)[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68355]

```bash
#conda install sra-tools #have to set up SRA tools to get fastq files
#conda install bwa-mem

mkdir /home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/
cd /home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/

#T47D 3 replicates of input DNA
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000691/SRR2000691.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000692/SRR2000692.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000693/SRR2000693.1

#MCF7 3 replicates of input DNA
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000750/SRR2000750.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000751/SRR2000751.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2000752/SRR2000752.1


for i in SRR*1; do 
    fastq-dump --split-3 $i & done & #using sra to fastq-dump data

for i in *fastq; do
    gzip $i & done & #gzip those fastqs because I'm not a criminal

mv SRR2000691.1.fastq.gz t74d_input.1.fq.gz
mv SRR2000692.1.fastq.gz t74d_input.2.fq.gz
mv SRR2000693.1.fastq.gz t74d_input.3.fq.gz
mv SRR2000750.1.fastq.gz mcf7_input.1.fq.gz
mv SRR2000751.1.fastq.gz mcf7_input.2.fq.gz
mv SRR2000752.1.fastq.gz mcf7_input.3.fq.gz


#using bowtie2 for alignment
#preparing human genome
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz

#alignment of reads using either bwa mem or bowtie2
for i in *fq.gz; do
bwa mem -P -S -t 10 -T 20 /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa $i 2>>align.log | samtools sort -@ 5 -T . -m 5G - | samtools rmdup -s - $outname; done &

for i in *fq.gz; do
outname=${i::-6}.bam
bowtie2 --threads 20 -x /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome -U $i | samtools sort -@ 5 -T . -m 5G - | samtools rmdup -s - $outname; done &

```

## Run HMMcopy using single-end deduplicated bam files as input
Some functions taken from SCOPE for convenient bam file read in


```R

library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)
library(SCOPE)
library(HMMcopy)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(reshape2)
library(philentropy)
library(dendextend)
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams/dedup_bams")


#set up single cell bam files from given directory
prepare_bam_bed_obj<-function(resol.){
    #Initalization
    bamfolder <- "/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/redo_bams/dedup_bams"
    bamFile <- list.files(bamfolder, pattern = 'deduplicated.bam$')
    #add reference bam files
    ref_bamfolder<-"/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata"
    ref_bamFile <- list.files(ref_bamfolder, pattern = '.bam$')
    bamdir <- file.path(bamfolder, bamFile)
    ref_bamdir<-file.path(ref_bamfolder,ref_bamFile)
    sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
    refname_raw<-substr(ref_bamFile,1, nchar(ref_bamFile)-4)
    
    bamdir<-append(bamdir,ref_bamdir)
    sampname_raw<-append(sampname_raw,refname_raw)

    #set genomic window size to 500kb
    bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38",resolution=resol,sex=T)#resolution of 100 = 100kbp
    #Compute GC content and mappability for each reference bin.
    mapp <- get_mapp(bambedObj$ref, hgref = "hg38")
    gc <- get_gc(bambedObj$ref, hgref = "hg38")
    values(bambedObj$ref) <- cbind(values(bambedObj$ref), DataFrame(gc, mapp))
    return(bambedObj)
}

#modified scope get_coverage_scDNA
read_in_reads<-function(i,seq="single-end"){
    bamurl <- bamdir[i]
    what <- c("rname", "pos", "mapq", "qwidth")
    if (seq == "paired-end") {
        flag <- scanBamFlag(isPaired = TRUE,
                isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE,
                isFirstMateRead = TRUE)
        param <- ScanBamParam(what = what, flag = flag)
        bam <- scanBam(bamurl, param = param)[[1]]
        }
    else if (seq == "single-end") {
            flag <- scanBamFlag(isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE) #removed checking for unpaired 
            param <- ScanBamParam(what = what, flag = flag)
            bam <- scanBam(bamurl, param = param)[[1]]
        }
        message("Getting coverage for sample ", i, ": ", sampname[i],"...", sep = "")
            bam.ref <- GRanges(seqnames = bam$rname, ranges = IRanges(start = bam[["pos"]],width = bam[["qwidth"]]))
            bam.ref <- suppressWarnings(bam.ref[countOverlaps(bam.ref,mask.ref) == 0])
            Y<- countOverlaps(ref, bam.ref)
            ref$counts<-Y
    return(Y)
}

#standard read correction
hmm_generator<-function(i){
    cellname<-colnames(Y)[i]
    dat<-cbind(ref,Y[,i])
    colnames(dat)[ncol(dat)]<-"reads"
    dat_copy <- suppressWarnings(correctReadcount(dat))
    write.table(dat_copy,
    file=paste0("counts_",cellname,".tsv"),
    sep="\t",quote=F,row.names=T,col.names=T)
    return(dat_copy)
}

#standard segmentation
copy_seg_func<-function(i){
    dat<-copy_estimate[[i]]
    dat<-HMMsegment(dat) #run HMMsegment
    segs<-makeGRangesFromDataFrame(dat$segs,keep.extra.columns=T)   #make granges object
    out_segs_idx<-as.data.frame(findOverlaps(segs,ref_granges)) #overlap GRanges for plotting consistent matrix bins
    out_segs<-segs[out_segs_idx$queryHits,]$state #return list of states length of reference bins
    return(out_segs)
}

#parameter initials
strength<-rep(1e+30,7)
#e<-rep(0.995,7) #shah suggested for scwgs
e<-rep(0.995,7)
mu<-c(-0.1, 1, 2, 1.5, 2.0, 2.5, 3.0)
#mu <- c(0, 1,2,3,4,5,6) #https://www.nature.com/articles/nmeth.4140.pdf?origin=ppub
lambda<-rep(20,7)
#nu<-rep(2.1,7) #shah suggested for scwgs
nu<-rep(2.1,7)
#kappa<-c(25, 50, 800, 50, 25, 25, 25)
kappa<-c(25, 50, 800, 50, 25, 25, 25)
m<-c(0, 1, 2, 3, 4, 5, 6)
eta<-rep(5e+04,7)
gamma<-rep(3,7)
S<-rep(35,7)
#S<-rep(0.01858295,7)
param<-data.frame(strength=strength,e=e,mu=mu,lambda=lambda,nu=nu,kappa=kappa,m=m,eta=eta,gamma=gamma,S=S)

#Seven copy number states were used; the m values were set to 0, 0.5, 1.0, 1.5, 2.0, 2.5, and 3.0 for copy number states 0, 1, 2, 3, 4, 5, and 6, respectively;
#the μ values were set to 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 for copy number states 0, 1, 2, 3, 4, 5, and 6, respectively;
#the κ values were set to  25, 50, 800, 50, 25, 25, 25 for copy number states 0, 1, 2, 3, 4, 5, and 6, respectively; 
#the e value was set to 0.995; and the S value was set to 35.
#https://www.biorxiv.org/content/10.1101/696179v1.full.pdf suggests e>=0.999999, nu=4 and strength=1e+07

#based on https://www.nature.com/articles/nmeth.4140#Sec8
param$m <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
param$mu <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
param$kappa <- c(25, 50, 800, 50, 25, 25, 25)
param$e <- 0.995
param$S <- 35

#Set up genome masking of duplicated regions and gaps
# Get segmental duplication regions
seg.dup <- read.table(system.file("extdata", "GRCh38GenomicSuperDup.tab", package = "WGSmapp"))
# Get hg38 gaps
gaps <- read.table(system.file("extdata", "hg38gaps.txt", package = "WGSmapp"))
#Set up genome masking from SCOPE get_masked_ref function
seg.dup <- seg.dup[!is.na(match(seg.dup[,1], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
seg.dup <- GRanges(seqnames = seg.dup[,1], ranges = IRanges(start = seg.dup[,2], end = seg.dup[,3]))
gaps <- gaps[!is.na(match(gaps[,2], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
gaps <- GRanges(seqnames = gaps[,2], ranges = IRanges(start = gaps[,3], end = gaps[,4]))
mask.ref <- sort(c(seg.dup, gaps)) # Generate mask region

#generate bam bed object
    resol=500
    bambedObj<-prepare_bam_bed_obj(resol.=resol) #500 is 500kbp
    saveRDS(bambedObj,file="bambedObj.500kbp.rds")
    bambedObj<-readRDS(file="bambedObj.500kbp.rds")
    ref <- bambedObj$ref
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname

#Get read count per bin per cell
    Y_out<-mclapply(1:length(sampname),read_in_reads,mc.cores=20)
    Y<-do.call("cbind",Y_out)
    colnames(Y)<-sampname
    ref<-as.data.frame(ref)
    row.names(Y)<-paste0(ref$seqnames,":",ref$start,"_",ref$end)
    saveRDS(Y,file="rawcount.500kbp.rds")
    Y<-readRDS(file="rawcount.500kbp.rds")
    
    #filter Y to cells with >500000 reads
    Y<-Y[,colSums(Y)>50000] #with fewer reads i was getting errors, was at 500k, but i want to try to go lower

#Prepare windows with gc and mappability data
    ref <- as.data.frame(bambedObj$ref)
    ref$gc<-ref$gc/100 #fitting to HMMcopy analysis
    ref$chr<-ref$seqname
    #ref$chr<-substr(ref$seqname,4,6)
    ref<-ref[c("chr","start","end","gc","mapp")]
    colnames(ref)<-c("chr","start","end","gc","map")
    ref_granges<-makeGRangesFromDataFrame(ref,keep.extra.columns=T) #used in the next few functions for overlapping consistent windows

#Normalize bins
    copy_estimate<-lapply(1:ncol(Y),hmm_generator) #correct bins by gc content, 
    names(copy_estimate)<-colnames(Y)
    saveRDS(copy_estimate,file="copyestimate.500kb.rds")
    copy_estimate<-readRDS(file="copyestimate.500kb.rds")

#segmentation for plotting
    copy_segmentation<-lapply(1:length(copy_estimate),copy_seg_func) #segment genome by copy estimate (log2 of cor.map)
    copy_segmentation<-lapply(copy_segmentation,as.numeric)
    names(copy_segmentation)<-colnames(Y)
    copy_segmentation<-as.data.frame(do.call("cbind",copy_segmentation))
    row.names(copy_segmentation)<-paste0(ref$chr,":",ref$start,"_",ref$end)
    copy_segmentation<-as.data.frame(t(copy_segmentation))
    saveRDS(copy_segmentation,file="copysegmentation.500kb.rds")
    copy_segmentation<-readRDS(file="copysegmentation.500kb.rds")

    copy_segmentation<-copy_segmentation[,!grepl("Y",colnames(copy_segmentation))] #remove Y

#filter to qc cells and set up annotations
    Y_plot<-as.data.frame(t(Y))
    Y_plot<-Y_plot[,!grepl("Y",colnames(Y_plot))]

    Y_plot<-as.data.frame(Y_plot %>% mutate_all(., ~ ifelse(.!=0, log10(.+0.00000000001), log10(1))))

    #add sample cell line names as metadata
    cell_annot<-data.frame(basename=unlist(lapply(strsplit(row.names(Y_plot),"_"),"[",1)))

    #set up treatment conditions
    cell_annot$treatment<-substr(unlist(lapply(strsplit(row.names(Y_plot),"_"),"[",2)),1,1) #set up T47D cell treatment first
    cell_annot[startsWith(cell_annot$basename,"M7"),]$treatment<-substr(cell_annot[startsWith(cell_annot$basename,"M7"),]$basename,3,3)
    cell_annot[startsWith(cell_annot$basename,"BSM7"),]$treatment<-substr(cell_annot[startsWith(cell_annot$basename,"BSM7"),]$basename,5,5)

    cell_annot$cellLine_treatment<-paste(cell_annot$basename,cell_annot$treatment,sep="_")

    #clean up metadata a bit more
    cell_annot$cellLine<-"T47D"
    cell_annot[cell_annot$basename %in% c("BSM7E6","M7C1A","M7C2B","M7E4C","mcf7bulk"),]$cellLine<-"MCF7"
    cell_annot[cell_annot$treatment %in% c("C"),]$treatment<-"control"
    cell_annot[cell_annot$treatment %in% c("E"),]$treatment<-"estrogen"
    cell_annot$bulk<-"single-cell"
    cell_annot[cell_annot$basename=="mcf7bulk",]$bulk<-"MCF7"
    cell_annot[cell_annot$basename=="t47dbulk",]$bulk<-"T47D"

    col_fun_cn=structure(c("#2166ac", "#67a9cf", "#f7f7f7","#fddbc7","#ef8a62","#b2182b","#630410"), names = c("1", "2", "3", "4", "5", "6","7"))

    ha = HeatmapAnnotation(which="row",cellline=cell_annot$cellLine,treatment=cell_annot$treatment,assay=cell_annot$cellLine_treatment,readcount=log10(colSums(Y)),
        col = list(assay = c("BSM7E6_E"="red","M7C1A_C"="blue","M7C2B_C"="purple","M7E4C_E"="brown","T_C"="green","T_E"="gray","mcf7bulk_i"="yellow","t47dbulk_i"="orange"),
            cellline=c("T47D"="red","MCF7"="blue"),
            treatment=c("control"="white","estrogen"="black","i"="green")))

    cellcount=data.frame(sampname=row.names(Y_plot),cellno=rep(1,nrow(Y_plot)))
    cellcount[grepl("input",cellcount$sampname),]$cellno<-"bulk"

    col_fun_reads=colorRamp2(quantile(unlist(Y_plot),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
    c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))


#optimize segmentation 
#https://github.com/shahcompbio/single_cell_pipeline/blob/master/single_cell/workflows/hmmcopy/scripts/hmmcopy.R
#https://advances.sciencemag.org/content/6/50/eabd6454
#Set up clustering of cells
    dist_method="euclidean"
    dist_x<-philentropy::distance(copy_segmentation,method=dist_method,as.dist.obj=T,use.row.names=T)
    dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
    k_search<-find_k(dend,krange=4:10) #search for optimal K from 5-10
    k_clus_number<-k_search$nc
    k_clus_id<-k_search$pamobject$clustering
    dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters

#Read count raw
plt1<-Heatmap(Y_plot,
    show_row_names=F,
    row_split=cell_annot$bulk,
    show_column_names=F,
    column_order=1:ncol(Y_plot),
    col=col_fun_reads,
    #cluster_rows=dend,
    row_title="Raw read count",
    name="Log10 Reads",#,
    left_annotation=ha,
    column_split=factor(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1))))
    )

#HMM segmentation
plt2<-Heatmap(copy_segmentation,
    show_row_names=F,
    show_column_names=F,
    row_split=cell_annot$bulk,
    column_order=1:ncol(copy_segmentation),
    col=col_fun_cn,
    #cluster_rows=dend,
    row_title="Segmentation",
    name="copy state",
    left_annotation=ha,
    column_split=factor(unlist(lapply(strsplit(colnames(copy_segmentation),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(copy_segmentation),":"),"[",1))))
    )


pdf("HMMcopy_test.pdf",width=10,height=3)
par(mfrow=c(2,1))
plt1
plt2
dev.off()

system("slack -F HMMcopy_test.pdf ryan_todo")

#custom mapd function
mapd<-function(cell){
  d<-unlist(lapply(1:nrow(Y)-1, function(x) (Y[x,cell]-Y[x+1,cell])/mean(Y[,cell])))
  mad<-median(abs(d-median(d)))
  return(mad)
}


mapd_list<-unlist(mclapply(colnames(Y),mapd,mc.cores=20)) #25 cores

dist_df<-data.frame("cellID"=colnames(Y),"mapd"=mapd_list) #making a MAPD matrix
dist_df$read_count<-colSums(Y)
dist_df$sample<-unlist(lapply(strsplit(dist_df$cellID,"_"),"[",1))
library(dplyr)
dist_df %>% group_by(sample) %>% summarize(mean_MAD=mean(mapd),median_MAD=median(mapd),sd_MAD=sd(mapd),mean_readcount=mean(read_count),median_readcount=median(read_count))


## A tibble: 7 x 6
#  sample   mean_MAD median_MAD  sd_MAD mead_readcount median_readcount
#  <chr>       <dbl>      <dbl>   <dbl>          <dbl>            <dbl>
#1 BSM7E6     0.196      0.190  0.0293        3187439.          2226462
#2 M7C1A      0.288      0.290  0.0301        1128324.          1069844
#3 M7C2B      0.333      0.344  0.0411         718630.           592200
#4 M7E4C      0.301      0.299  0.0284         304512.           327448
#5 mcf7bulk   0.0905     0.0921 0.0140       22368933.         20262739
#6 T          0.272      0.262  0.0473         752268.           570749
#7 t47dbulk   0.0901     0.0917 0.00733      30191789.         32174336



plt1<-ggplot(dist_df,aes(x=paste(sample),y=mapd,color=paste(sample)))+geom_jitter()+geom_boxplot(aes(fill=NULL))+theme_bw()+ylim(c(0,1))
plt2<-ggplot(dist_df,aes(x=paste(sample),y=log10(read_count),color=paste(sample)))+geom_jitter()+geom_boxplot(aes(fill=NULL))+theme_bw()+ylim(c(0,9))
ggsave(plt1/plt2,file="mapd_scores.pdf")
system("slack -F mapd_scores.pdf ryan_todo")


#didnt run below yet



### Merge single-cell data by clades for read depth visualization
#Restructure bins to finer resolution

    bambedObj<-readRDS(file="bambedObj.500kbp.rds")
    ref <- as.data.frame(bambedObj$ref)
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname
    ref$gc<-ref$gc/100 #fitting to HMMcopy analysis
    ref$chr<-ref$seqname
    #ref$chr<-substr(ref$seqname,4,6)
    ref<-ref[c("chr","start","end","gc","mapp")]
    colnames(ref)<-c("chr","start","end","gc","map")
    ref_granges<-makeGRangesFromDataFrame(ref,keep.extra.columns=T) #used in the next few functions for overlapping consistent windows

    Y<-readRDS(file="rawcount.500kbp.rds")

    #Merge data by row means
    #arbitrarily set bulk data as separate clades
    k_clus_id[which(grepl("mcf7bulk",names(k_clus_id)))]<-"mcf7bulk"
    k_clus_id[which(grepl("t47dbulk",names(k_clus_id)))]<-"t47dbulk"
    clade_dat<-list()
    for (i in unique(k_clus_id)){
        clade<-names(k_clus_id[k_clus_id==i])
        clade_combined<-rowSums(Y[,colnames(Y) %in% clade],na.rm=T)
        clade_dat[[paste0("clade_",i)]]<-clade_combined
    }

    Ybulk<-do.call("cbind",clade_dat)
    colnames(Ybulk)<-names(clade_dat)


#standard read correction
hmm_generator_bulk<-function(i){
    cellname<-colnames(Ybulk)[i]
    dat<-cbind(ref,Ybulk[,i])
    colnames(dat)[ncol(dat)]<-"reads"
    dat_copy <- suppressWarnings(correctReadcount(dat))
    write.table(dat_copy,
    file=paste0("bulkcounts_",cellname,".tsv"),
    sep="\t",quote=F,row.names=T,col.names=T)
    return(dat_copy)
}

#Normalize bins
    copy_estimate<-lapply(1:ncol(Ybulk),hmm_generator_bulk) #correct bins by gc content, 
    names(copy_estimate)<-colnames(Ybulk)
    saveRDS(copy_estimate,file="copyestimate.bulk.500kb.rds")

#segmentation for plotting
    copy_segmentation<-lapply(1:length(copy_estimate),copy_seg_func) #segment genome by copy estimate (log2 of cor.map)
    copy_segmentation<-lapply(copy_segmentation,as.numeric)
    names(copy_segmentation)<-colnames(Ybulk)
    copy_segmentation<-as.data.frame(do.call("cbind",copy_segmentation))
    row.names(copy_segmentation)<-paste0(ref$chr,":",ref$start,"_",ref$end)
    copy_segmentation<-as.data.frame(t(copy_segmentation))
    saveRDS(copy_segmentation,file="copysegmentation.bulk.500kb.rds")


    copy_segmentation<-readRDS(file="copysegmentation.bulk.500kb.rds")
    copy_estimate<-readRDS(file="copyestimate.bulk.500kb.rds")

#prepare data frames for plotting
    #remove chrY and set up data frames for plotting
    copy_segmentation<-copy_segmentation[,!grepl("Y",colnames(copy_segmentation))]
    copy_segmentation<-as.data.frame(t(copy_segmentation))
    copy_estimate_plot<-as.data.frame(do.call("cbind",lapply(1:length(copy_estimate),function(x) copy_estimate[[x]]$cor.map)))
    colnames(copy_estimate_plot)<-colnames(Ybulk)
    row.names(copy_estimate_plot)<-row.names(Ybulk)
    copy_estimate_plot<-copy_estimate_plot[row.names(copy_estimate_plot)%in%row.names(copy_segmentation),]

    #convert copy estimates to integer copy number state via https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4140/MediaObjects/41592_2017_BFnmeth4140_MOESM175_ESM.pdf
    copy_estimate_plot<-lapply(1:ncol(copy_estimate_plot), function(x) {
        int_val<-median(as.numeric(copy_estimate_plot[which(copy_segmentation[,x]==2),x]),na.rm=T)/2
        return(copy_estimate_plot[,x]/int_val)
        })
    copy_estimate_plot<-as.data.frame(do.call("cbind",copy_estimate_plot))
    colnames(copy_estimate_plot)<-colnames(copy_segmentation)
    row.names(copy_estimate_plot)<-row.names(copy_segmentation)

    est_plt<-melt(as.matrix(copy_estimate_plot))
    seg_plt<-melt(as.matrix(copy_segmentation))
    colnames(seg_plt)<-c("gloc","clade","state")
    plt_melt<-cbind(seg_plt,copy_est=est_plt$value)

    cols<-c("1"="#2166ac", "2"="#d0e7f5", "3"="#f5f5f5","4"="#fddbc7","5"="#ef8a62","6"="#b2182b")
    #plot bins for a cell across a chromosome

    range_gc<-quantile(plt_melt$copy_est, na.rm = TRUE, prob = c(0.05,0.95))
    plt_melt$contig<-unlist(lapply(strsplit(as.character(plt_melt$gloc),":"),"[",1))
    plt_melt$contig<-factor(plt_melt$contig,level=unique(plt_melt$contig))
    plt_melt$start<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(plt_melt$gloc),"_"),"[",1)),":"),"[",2)))
    plt_melt$end<-as.numeric(unlist(lapply(strsplit(as.character(plt_melt$gloc),"_"),"[",2)))

    plt_melt$row_order<-1:(nrow(plt_melt)/length(unique(plt_melt$clade)))

    plt<-ggplot(plt_melt,aes(x=row_order,y=copy_est))+
    geom_rect(aes(fill=as.character(state),xmin=row_order,xmax=row_order+1,ymin=0,ymax=6,alpha=0.01))+
    scale_fill_manual(values=cols)+
    geom_point(color="black",size=1,alpha=0.2)+
    ylab("")+
    xlab("")+
    facet_grid(plt_melt$clade~plt_melt$contig,space="free",scales="free_x")+
    theme_minimal()+
    scale_y_continuous(breaks=c(0,1,2,3,4,5,6),limits=c(0,6))+
    theme(axis.text.y = element_text(size=30),
        axis.text.x = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        panel.spacing.x=unit(0.1,"lines"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position="none",
        panel.border = element_rect(colour = "black", fill = NA,size=3))


    ggsave(plt,file="HMMcopy.merged.500kb.png",width=2000,height=1000,units="mm",limitsize=F)

system("slack HMMcopy.merged.500kb.png")

```

