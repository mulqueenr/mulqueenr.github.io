---
title: HMMcopy on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /hmmcopy_nmt/
category: CEDAR
---

# Script to test Shah lab CNV calling using HMMcopy on SE bismark aligned reads

## Preparing bam files from Aaron Doe's preprocessing
Need to merge resequenced bam files and combine read 1 and read 2 into a single-end bam

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
mkdir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams
cd /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams
ls /home/groups/CEDAR/doe/projects/my_NMT/MCF7_reseq/merge_all_MCF7/scNMT_NOMeWorkFlow/bismarkSE/dedup/*merged.bam > files.txt
ls /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scNMT_NOMeWorkFlow/bismarkSE/*bam | grep -v "M7" >> files.txt
for i in `cat files.txt`; do outname=`basename $i`; ln -s $i $outname; done
for i in *_R1.*bt2.bam; do 
r1=$i;
r2=`echo $i | sed "s/R1/R2/g"`;
r2=`echo $r2 | sed "s/val_1_/val_2_/g"`;
outname=`echo $i | awk 'OFS="_" {split($1,a,"_");print a[1],a[2],a[3]}'`;
samtools cat -o ${outname}_merged.bam $r1 $r2 ;
done &
```
{% endcapture %} {% include details.html %} 

## Batch script for deduplication
Using the bismark deduplication script.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-374
#SBATCH --tasks-per-node=30 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams/*merged.bam | wc -l`

bam_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams/*merged.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun deduplicate_bismark \
--single \
--output_dir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/dedup_bams \
${bam_in}

```
{% endcapture %} {% include details.html %} 

## Perform modal copy number bin normalization similar to Shah sc-wgs work
This takes the input of files that are generated in the hmm_generator function in the R script (files starting with "counts_")
{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 

## Preparing bulk WGS data for cell lines to check CNVs
Files are from a previous publication by Hisham (here.)[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68355]
{% capture summary %} Code {% endcapture %} {% capture details %}  
```bash
conda install sra-tools #have to set up SRA tools to get fastq files
conda install bwa-mem

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
{% endcapture %} {% include details.html %} 

## Run HMMcopy using single-end deduplicated bam files as input
Some functions taken from SCOPE for convenient bam file read in

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
library(philentropy)
library(dendextend)
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test")


#set up single cell bam files from given directory
prepare_bam_bed_obj<-function(resol.){
    #Initalization
    bamfolder <- "/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/dedup_bams"
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
copy_segmentation<-function(i){
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
    #filter Y to cells with >100000 reads
    Y<-Y[,colSums(Y)>100000] #with fewer reads i was getting errors

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
    copy_segmentation<-lapply(1:length(copy_estimate),copy_segmentation) #segment genome by copy estimate (log2 of cor.map)
    copy_segmentation<-lapply(copy_segmentation,as.numeric)
    names(copy_segmentation)<-colnames(Y)
    copy_segmentation<-as.data.frame(do.call("cbind",copy_segmentation))
    row.names(copy_segmentation)<-paste0(ref$chr,":",ref$start,"_",ref$end)
    copy_segmentation<-as.data.frame(t(copy_segmentation))
    saveRDS(copy_segmentation,file="copysegmentation.500kb.rds")
    copy_segmentation<-readRDS(file="copysegmentation.500kb.rds")

#filter to qc cells and set up annotations
    Y_plot<-as.data.frame(t(Y))
    Y_plot<-as.data.frame(Y_plot %>% mutate_all(., ~ ifelse(.!=0, log10(.+0.00000000001), log10(1))))

    col_fun_cn=structure(c("#2166ac", "#67a9cf", "#f7f7f7","#fddbc7","#ef8a62","#b2182b","#630410"), names = c("1", "2", "3", "4", "5", "6","7"))

    ha = HeatmapAnnotation(which="row",assay=unlist(lapply(strsplit(row.names(Y_plot),"_"),"[",1)),
        col = list(assay = c("BSM7E6"="red","M7C1A"="blue","M7C2B"="purple","M7E4C"="brown","T"="green","mcf7bulk"="yellow","t47dbulk"="orange")))

    cellcount=data.frame(sampname=row.names(Y_plot),cellno=rep(1,nrow(Y_plot)))
    cellcount[grepl("input",cellcount$sampname),]$cellno<-"bulk"

    col_fun_reads=colorRamp2(quantile(unlist(Y_plot),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
    c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))


#optimize segmentation 
#https://github.com/shahcompbio/single_cell_pipeline/blob/master/single_cell/workflows/hmmcopy/scripts/hmmcopy.R
#https://advances.sciencemag.org/content/6/50/eabd6454

#Set up clustering of cells
    dist_method="jaccard"
    dist_x<-philentropy::distance(copy_segmentation,method=dist_method,as.dist.obj=T,use.row.names=T)
    dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
    k_search<-find_k(dend,krange=3:10) #search for optimal K from 2-10
    k_clus_number<-k_search$nc
    k_clus_id<-k_search$pamobject$clustering
    dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters

#Read count raw
plt1<-Heatmap(Y_plot,
    show_row_names=F,
    show_column_names=F,
    column_order=1:ncol(Y_plot),
    col=col_fun_reads,
    cluster_rows=dend,
    row_title="Raw read count",
    name="Log10 Reads",#,
    left_annotation=ha,
    column_split=factor(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1))))
    )

#HMM segmentation
plt2<-Heatmap(copy_segmentation,
    show_row_names=F,
    show_column_names=F,
    column_order=1:ncol(copy_segmentation),
    col=col_fun_cn,
    cluster_rows=dend,
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

#custom mapd function
mapd<-function(cell){
  d<-unlist(lapply(1:nrow(Y)-1, function(x) (Y[x,cell]-Y[x+1,cell])/mean(Y[,cell])))
  mad<-median(abs(d-median(d)))
  return(mad)
}

mapd_corrected<-function(cell){
  d<-unlist(lapply(1:ncol(shah_copy)-1, function(x) 
    if(is.finite(shah_copy[cell,x]) & is.finite(shah_copy[cell,x+1])){
        (shah_copy[cell,x]-shah_copy[cell,x+1])/mean(shah_copy[cell,])}
  ))
  mad<-median(abs(d-median(d)))
  return(mad)
}

mapd_list<-unlist(mclapply(colnames(Y),mapd,mc.cores=20)) #25 cores
mapd_corrected_list<-unlist(mclapply(row.names(shah_copy),mapd_corrected,mc.cores=20)) #25 cores

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
```

{% endcapture %} {% include details.html %} 
