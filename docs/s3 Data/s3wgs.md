---
title: s3WGS Analysis
layout: s3data
author: Ryan Mulqueen
permalink: /s3wgs/
category: s3processing
---

# Processing for s3WGS portion of s3 paper.
This notebook is a continuation of ["s3 Preprocessing"](https://mulqueenr.github.io/s3preprocess/) and details the processing of s3WGS libraries. This notebook starts with a merged, barcode-based removal of duplicate bam file.

{% capture summary %} Initial Files and Directory Structure {% endcapture %} {% capture details %}  
```bash
#Initial directory structure (filtered to relevant files):
#Note GCC type libraries will also be processes as WGS files, so are included here.
/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis
├── filtered_bam
│   ├── CRC-4442CRC-4671_GCC_H.bbrd.q10.filt.cellIDs.list
│   ├── CRC4442_WGS_E.bbrd.q10.filt.cellIDs.list
│   ├── CRC4442_WGS_F.bbrd.q10.filt.cellIDs.list
│   ├── CRC4671_WGS_G.bbrd.q10.filt.cellIDs.list
│   ├── GM12878_WGS_B.bbrd.q10.filt.cellIDs.list
│   ├── GM12878_WGS_D.bbrd.q10.filt.cellIDs.list
├── raw_alignment
│   ├── CRC-4442CRC-4671_GCC_H.nsrt.bam
│   ├── CRC4442_WGS_E.nsrt.bam
│   ├── CRC4442_WGS_F.nsrt.bam
│   ├── CRC4671_WGS_G.nsrt.bam
│   ├── GM12878_WGS_B.nsrt.bam
│   ├── GM12878_WGS_D.nsrt.bam
├── s3wgs_data
├── s3gcc_data

```
{% endcapture %} {% include details.html %} 


# Filter bam files to just cells passing original QC
Performing this on pre-barcode based remove duplicate bams for future complexity plotting and projection modeling.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
#First going to subset un-rmdup bam files to just those passing QC
DIR="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis"
OUTPUT_DIR=$DIR"/s3wgs_data/singlecell_bam"
mkdir $OUTPUT_DIR

#concatenating lists of cellIDs passing QC per experiment on initial QC pass
cat $DIR/filtered_bam/GM12878_WGS*filt.cellIDs.list > $DIR/s3wgs_data/s3wgs_gm12878.bbrd.q10.filt.cellIDs.list 
cat $DIR/filtered_bam/CRC*WGS*filt.cellIDs.list > $DIR/s3wgs_data/s3wgs_crc.bbrd.q10.filt.cellIDs.list 
cat $DIR/filtered_bam/CRC*GCC*filt.cellIDs.list > $DIR/s3gcc_data/s3gcc_crc.bbrd.q10.filt.cellIDs.list 

#Filter pre-rmdup bam files to just cellIDs passing filter first
scitools bam-filter -L $DIR/s3wgs_data/s3wgs_gm12878.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3wgs_gm12878_D $DIR/raw_alignment/GM12878_WGS_D.nsrt.bam &
scitools bam-filter -L $DIR/s3wgs_data/s3wgs_gm12878.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3wgs_gm12878_B $DIR/raw_alignment/GM12878_WGS_B.nsrt.bam &
scitools bam-filter -L $DIR/s3wgs_data/s3wgs_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3wgs_crc4671_G $DIR/raw_alignment/CRC4671_WGS_G.nsrt.bam &
scitools bam-filter -L $DIR/s3gcc_data/s3gcc_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3gcc_crc_H $DIR/raw_alignment/CRC-4442CRC-4671_GCC_H.nsrt.bam &
scitools bam-filter -L $DIR/s3wgs_data/s3wgs_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3wgs_crc4442_E $DIR/raw_alignment/CRC4442_WGS_E.nsrt.bam &
scitools bam-filter -L $DIR/s3wgs_data/s3wgs_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3wgs_data/s3wgs_crc4442_F $DIR/raw_alignment/CRC4442_WGS_F.nsrt.bam &

#Add read group for cellID specific splitting
for i in *filt.bam ; do scitools bam-addrg $i & done &

#Split out pre-deduplicate bam files for single_cell projections
mkdir $DIR/single_cell_projections

#cellID is few enough that there are no IO errors, split by RG
#for i in s3gcc_crc_H.filt.RG.bam s3wgs_crc4442_E.filt.RG.bam  s3wgs_crc4442_F.filt.RG.bam  s3wgs_crc4671_G.filt.RG.bam  s3wgs_gm12878_B.filt.RG.bam; do samtools split -@20 -f './single_cell_projections/%*_%!.%.' $i ; done & #perform samtools split on RG (cellIDs)
for i in s3wgs_gm12878_D.filt.RG.bam; do samtools split -@20 -f './single_cell_projections/%*_%!.%.' $i ; done & #perform samtools split on RG (cellIDs)

#perform parallelized duplicate removal
find . -type f -name '*.bam' | parallel -j 20 scitools bam-rmdup {} 

#perform parallelized read projection
find . -type f -name '*.bam' | parallel -j 20 scitools bam-project -r 1000 -X -e {} &

for i in s3wgs_gm12878_D.filt.RG_CTCTGGAGTTGGTTCTAATGCCTC.bam s3wgs_gm12878_D.filt.RG_CTCTGGAGTTGGTTCTCAGTAGGC.bam s3wgs_gm12878_D.filt.RG_CTCTGGAGTTGGTTCTTTGCCTAG.bam s3wgs_gm12878_D.filt.RG_GACTAGCACAATTCGTAGTTCAGG.bam s3wgs_gm12878_D.filt.RG_GACTAGCACAATTCGTGTCGGAGC.bam s3wgs_gm12878_D.filt.RG_GACTAGCATTGGTTCTCAGTAGGC.bam; do scitools bam-project -r 1000 -X -e $i & done

find . -type f -name 's3wgs_gm12878_D*.bam' | parallel -j 5 scitools bam-project -r 1000 -X -e {} &

#Then move folders into s3wgs_data directory 
mkdir complexity
mv *complexity* ./complexity/ #hold complexity output
mkdir projections
mv *projection* *rand* ./projections/ #hold single cell projections output
mv *q10.bam ../singlecell_bam/ #hold deduplicated single cell bams for SCOPE

#perform parallelized insert size analysis in deduplicated single cell bams
cd ../singlecell_bam
find . -type f -name '*q10.bam' | parallel -j 15 scitools isize -M 20000000 {} & #Setting size limit high for s3gcc

#Sort single cell bam files
find . -type f -name '*.bam' | parallel -j 20 samtools sort -o {}.sorted.bam {} &

mkdir isize
mv *isize* ./isize


#Annotation file containing apriori knowledge from experimental setup (cellID, assay and sample)
/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/s3wgs_gcc.annot

#head s3wgs_gcc.annot
#cellID  experiment      assay   sample  cellID_idx
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTAGTTCAGG     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTAGTTCAGG
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTCATAGAGT     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTCATAGAGT
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTCTAGCGCT     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTCTAGCGCT
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTCTCTCGTC     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTCTCTCGTC
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTGGCATTCT     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTGGCATTCT
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTTACCGAGG     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTTACCGAGG
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTTACTCATA     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTTACTCATA
#s3wgs_crc4442_E.RG_ACGCGACGAGAGGACTTATTAGCT     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGAGAGGACTTATTAGCT
#s3wgs_crc4442_E.RG_ACGCGACGCAATGAGAAGTTCAGG     s3wgs_crc4442_E s3wgs   crc4442 ACGCGACGCAATGAGAAGTTCAGG
```
{% endcapture %} {% include details.html %} 

## QC Directory Structure 

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data
├── singlecell_bam #contains deduplicated single cell bams for SCOPE analysis
│   └── isize #contains isize output for insert size distribution
└── single_cell_projections #contains pre-dedupliacted single cell bams
    ├── complexity #contains deduplication complexity logs per cell
    ├── projections #contains projections of complexity per cell
    ├── undedup_bams #contains pre-bbrd bams
    
```
{% endcapture %} {% include details.html %} 

## Collate data into single files for plotting

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash

cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data

#read projections
for i in ./single_cell_projections/projections/*.read_projections; 
do cellid=${i:38:-17};
awk -v cellid=$cellid 'OFS="\t" {print $1,$2,$3,cellid}' $i/cell_summaries.txt;
done > ./s3wgs_projected_reads.txt 

#Complexity
for i in ./single_cell_projections/complexity/*complexity.txt;
do cellid=${i:37:-15};
awk -v cellid=$cellid 'OFS="\t" {print $2,$3,$4,$5,cellid}' $i;
done > ./s3wgs_complexity.txt
```
Combine data for complexity files.

```R
dat<-read.table("s3wgs_complexity.txt",header=F)
annot<-read.table("s3wgs_gcc.annot",header=T)
colnames(dat)<-c("cellID_idx","total_reads","uniq_reads","perc_uniq","cellID")
dat<-merge(dat,annot,by="cellID_idx")
write.table(dat,file="s3wgs_complexity.txt",sep="\t",col.names=T,quote=F)
```
Collate insert size distribution.
```bash
#isize
for i in ./singlecell_bam/isize/*values;
do cellid=${i:23:-22};
awk -v cellid=$cellid 'OFS="\t" {print $1,cellid}' $i;
done > ./s3wgs_isize.txt 

```
{% endcapture %} {% include details.html %} 

## Generate complexity projections per cell

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(dplyr)


setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
#Updating gcc cells with cell line sample name (a priori from pcr indexes)
annot<-read.table("s3wgs_gcc.annot",header=T)

proj<-read.table("s3wgs_projected_reads.txt",sep="\t")
colnames(proj)<-c("proj_perc_uniq","proj_reads","proj_unique_reads","cellID")
proj$cellID_idx<-unlist(lapply(strsplit(proj$cellID,"RG_"),"[",2))

dat<-merge(annot,proj,by="cellID_idx")
dat$plate<-unlist(lapply(strsplit(dat$experiment,"_"),"[",3))


dat<- as.data.frame(dat %>% 
                    group_by(sample,assay,proj_perc_uniq) %>% 
                    summarize(mean=mean(log10(proj_unique_reads)),
                              sd=sd(log10(proj_unique_reads)),
                              median=median(log10(proj_unique_reads))))


dat[dat$proj_perc_uniq==0.05,] #looking at read depth at estimated 95% duplication rate.
#     sample assay proj_perc_uniq     mean        sd   median
#1   crc4442 s3gcc           0.05 5.845659 0.3314061 5.751168
#95  crc4442 s3wgs           0.05 6.509959 0.3144060 6.489320
#191 crc4671 s3gcc           0.05 6.050565 0.3503318 6.041202
#284 crc4671 s3wgs           0.05 6.393986 0.3490392 6.337743
#380 gm12878 s3wgs           0.05 7.087167 0.3415564 7.044962

ggplot(dat,aes(x=as.numeric(proj_perc_uniq),fill = paste(sample,assay),color=paste(sample,assay)))+
geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),alpha=0.1) +
geom_line(aes(y=mean,linetype="solid"))+
geom_line(aes(y=median,linetype="dashed")) +
theme_bw() + scale_x_reverse() + facet_grid(assay ~ .)

ggsave(file="projected_readcount.png")
ggsave(file="projected_readcount.pdf")

system("slack -F projected_readcount.pdf ryan_todo")
system("slack -F projected_readcount.png ryan_todo")

#Plot Boxplots
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
annot<-read.table("s3wgs_gcc.cellsummary.txt",header=T)

ggplot(annot,aes(x=paste(assay,sample),y=log10(mapped),color=paste(assay,sample)))+geom_jitter()+geom_boxplot()+ylim(c(0,8))+theme_bw()
ggsave(file="readcount.png")
ggsave(file="readcount.pdf")

system("slack -F readcount.pdf ryan_todo")
```
{% endcapture %} {% include details.html %} 

##Using SCOPE to analyze single-cell data on single cell bam directory
Library used for scWGS data analysis is [SCOPE](https://github.com/rujinwang/SCOPE) available as a [preprint](https://www.biorxiv.org/content/10.1101/594267v1.full). SCOPE works on pre-aligned deduplicated bam files. So I split files post-deduplication into a subdirectory to load in (above).

Using R 4.0.0
Note a lot of this code and even the comments and explanation is directly taken from the SCOPE example given.


## Read in files from directory
First reading in the split bam files and setting up the reference genome.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(parallel)

#Initalization
bamfolder <- "/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/singlecell_bam/sorted_bams"
bamFile <- list.files(bamfolder, pattern = 'sorted.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- paste(sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1),sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 3),sep=".")

#set genomic window size to 500kb
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38",resolution=500,sex=T)#resolution of 100 = 100kbp

#Compute GC content and mappability for each reference bin.
mapp <- get_mapp(bambedObj$ref, hgref = "hg38")
gc <- get_gc(bambedObj$ref, hgref = "hg38")
values(bambedObj$ref) <- cbind(values(bambedObj$ref), DataFrame(gc, mapp))

#For 500kb bins
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 10, seq = 'paired-end', hgref = "hg38") #using a Q10 score for filtering

saveRDS(coverageObj,"scope_covobj.500kb.rds")
saveRDS(bambedObj,"scope_bambedObj.500kb.rds")
```
{% endcapture %} {% include details.html %} 


##Quality control of bins and cells via SCOPE

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(parallel)

coverageObj<-readRDS("scope_covobj.500kb.rds")
bambedObj<-readRDS("scope_bambedObj.500kb.rds")
#Perform QC and Ploidy Estimate
#For 500kb bins
#Perform QC
#Quality control to remove samples/cells with low proportion of mapped reads, bins that have extreme GC content (less than 20% and greater than 80%) and low mappability (less than 0.9) to reduce artifacts.
QCmetric_raw <- get_samp_QC(bambedObj)
saveRDS(QCmetric_raw,"scope_qcmetric.500kb.rds")
QCmetric_raw<-readRDS("scope_qcmetric.500kb.rds")

qcObj <- perform_qc(Y_raw = coverageObj$Y,sampname_raw = row.names(QCmetric_raw), ref_raw = bambedObj$ref, QCmetric_raw = QCmetric_raw)
"""Removed 0 samples due to failed library preparation.
Removed 133 samples due to failure to meet min coverage requirement.
Removed 0 samples due to low proportion of mapped reads.
Excluded 276 bins due to extreme GC content.
Excluded 450 bins due to low mappability.
Removed 0 samples due to excessive zero read counts in
            library size calculation.
There are 1268 samples and 5449 bins after QC step."""
saveRDS(qcObj,"scope_qcObj.500kb.rds")

#Generate gini coeffcients for 500kb binned data
Gini<-get_gini(qcObj$Y) #SCOPE function
saveRDS(Gini,"scope_gini.500kb.rds")

#500 kb first pass normalization
#First pass at estimating ploidy using no latent factors
norm_index_gini<-which(grepl("gm12878",colnames(qcObj$Y)))#grabbing 20 cells lowest on gini index

# first-pass CODEX2 run with no latent factors
normObj_gini<- normalize_codex2_ns_noK(Y_qc = qcObj$Y, gc_qc = qcObj$ref$gc, norm_index = norm_index_gini)
saveRDS(normObj_gini,"scope_noKnorm.gini.500kb.rds")

normObj_gini<-readRDS("scope_noKnorm.gini.500kb.rds")
##Ploidy initialization 500kb
#Generate ploidy estimate per cell from gini index
ploidy <- initialize_ploidy(Y = qcObj$Y, Yhat = normObj_gini$Yhat, ref = qcObj$ref, SoS.plot = F)
saveRDS(ploidy,"scope_ploidy.gini.500kb.rds")
```
{% endcapture %} {% include details.html %} 


## Using ploidy estimates to assess copy number change

SCOPE uses gini index to determine normal ploidy cells. I run here, but a priori set GM12878 cells at normal ploidy later in script.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(SCOPE)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
#Performing normalization via gini euploid estimates
qcObj<-readRDS("scope_qcObj.500kb.rds")
ploidy<-readRDS("scope_ploidy.gini.500kb.rds")
norm_index<-which(grepl("gm12878",colnames(qcObj$Y)))
normObj<-readRDS("scope_noKnorm.gini.500kb.rds")

#Normalize with no Latent Factors
normObj.scope.gini <- normalize_scope_foreach(
	Y_qc = qcObj$Y,
	gc_qc = qcObj$ref$gc,
	K = 5,
	ploidyInt = ploidy,
	norm_index = norm_index,
	T = 1:6,
	beta0 = normObj$beta.hat)
    
saveRDS(normObj.scope.gini,"scope_normforeach.gini.500kb.rds")
```
{% endcapture %} {% include details.html %} 

## Generate final segmentation plots.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R

#Clustering and plotting
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(philentropy)
library(ape)
library(ggdendro)
library(dendextend)
library(rpart)
library(cluster)
library(parallel)
library(SCOPE)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")


#Apriori Ploidy Estimate
annot<-read.table("s3wgs_gcc.cellsummary.txt")
qcObj<-readRDS("scope_qcObj.500kb.rds")

#Set up apriori ploidy estimate
apriori_ploidy<-data.frame(sample=annot[match(qcObj$sampname, annot$cellID),]$sample, ploidy=c(2))
apriori_ploidy[apriori_ploidy$sample=="crc4442",]$ploidy<-2.6 #4442 ploidy of 2.6 from karyotyping
apriori_ploidy[apriori_ploidy$sample=="crc4671",]$ploidy<-2.6 #4671 ploidy of 2.6 from karyotyping
saveRDS(apriori_ploidy,"scope_ploidy.apriori.500kb.rds")
apriori_ploidy<-readRDS("scope_ploidy.apriori.500kb.rds")

norm_index<-which(apriori_ploidy$sample=="gm12878") #effectively gm12878 cells

# first-pass CODEX2 run with no latent factors
normObj<- normalize_codex2_ns_noK(Y_qc = qcObj$Y, gc_qc = qcObj$ref$gc, norm_index = norm_index)
saveRDS(normObj,"scope_noKnorm.500kb.rds")
normObj<-readRDS("scope_noKnorm.500kb.rds")

#Performing normalization via a priori euploid states
#Normalize with 5 Latent Factors
normObj.scope <- normalize_scope_foreach(Y_qc = qcObj$Y,gc_qc = qcObj$ref$gc,
    K = 5,ploidyInt = apriori_ploidy$ploidy,
    norm_index = norm_index,T = 1:6,beta0 = normObj$beta.hat,minCountQC=20,nCores=30)
saveRDS(normObj.scope,"scope_normforeach.500kb_k5.rds") #K=5,ploidyInt = apriori_ploidy, norm_index=gm12878, 
normObj.scope<-readRDS("scope_normforeach.500kb_k5.rds")

#Perform CBS
chr_cbs <- function(x,scope=normObj.scope) {
    print(paste("Running for ",x))
    chr_seg<-segment_CBScs(Y = qcObj$Y,
    Yhat = as.data.frame(scope$Yhat[[which.max(scope$BIC)]]),
    sampname = colnames(qcObj$Y),
    ref = qcObj$ref,
    chr = x,
    mode = "integer", max.ns = 2)
    return(chr_seg)
}

#Running CBS segmentation on 500kb windows (gini ploidy estimates)
chrs <- unique(as.character(seqnames(qcObj$ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs

segment_cs<-list()
segment_cs<-mclapply(chrs,FUN=chr_cbs,mc.cores=length(chrs))
names(segment_cs) <- chrs #mclapply returns jobs in same order as specified
saveRDS(segment_cs,"scope_segmentcs.500kb_k5.ns2.rds") #ns 2
segment_cs<-readRDS("scope_segmentcs.500kb_k5.ns2.rds") #ns 2

#Plot raw reads
Y <- as.data.frame(qcObj$Y)
Y<-as.data.frame(Y %>% mutate_all(., ~ ifelse(.!=0, log10(.+0.00000000001), log10(1))))
row.names(Y)<-row.names(qcObj$Y)
Y<-as.data.frame(t(Y))

#Plot normalized reads
Yhat = as.data.frame(normObj.scope$Yhat[[which.max(normObj.scope$BIC)]])
Yhat<-as.data.frame(Yhat %>% mutate_all(., ~ ifelse(.!=0, log10(.+0.00000000001), log10(1))))
colnames(Yhat)<-qcObj$sampname
row.names(Yhat)<-row.names(qcObj$Y)
Yhat<-as.data.frame(t(Yhat))

col_fun_reads_raw=colorRamp2(quantile(unlist(Y),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))
col_fun_reads_normalized=colorRamp2(quantile(unlist(Yhat),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))

#Set up bins for plotting
bin_loci<-data.frame(qcObj$ref)
bin_loci$bin_name<-row.names(bin_loci)
bin_loci_subject<-makeGRangesFromDataFrame(bin_loci,keep.extra.columns=T)
iCN <- do.call(rbind, lapply(segment_cs, function(z){names(z)}))
iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
row.names(iCN)<-paste(bin_loci$seqnames,bin_loci$start,bin_loci$end,sep="_")
iCN[which(iCN>5,arr.ind=T)]<-5 #limit copy number amplification to 5+
iCN<-as.data.frame(t(iCN))
col_fun_iCN = structure(c("#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b"), names = c("0","1", "2", "3", "4", "5"))

#Set up clustering of cells
x<-iCN
x[which(x>2,arr.ind=T)]<-3
x[which(x<2,arr.ind=T)]<-1 #categorize copy numer data to account for ploidy changes
dist_method="jaccard"
dist_x<-philentropy::distance(x,method=dist_method,as.dist.obj=T,use.row.names=T)
#dist_x<-as.dist(cor(t(x),method=dist_method))
dend <- dist_x %>%  hclust(method="ward.D") %>% as.dendrogram(edge.root=F,h=3) 
k_search<-find_k(dend,krange=6:12)
k_clus_number<-k_search$nc
k_clus_id<-k_search$pamobject$clustering
dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters

#Perform k=3 clustering as well
k_search_3<-find_k(dend,krange=3)
k_clus_number_3<-k_search_3$nc
k_clus_id_3<-k_search_3$pamobject$clustering

saveRDS(dend,"scope_segmentcs.500kb_k2.maxns2.dendrogram.rds")
saveRDS(k_clus_id,"scope_segmentcs.500kb_k2.maxns2.clusterID.rds")
saveRDS(k_clus_id_3,"scope_segmentcs.500kb_k2.maxns2.clus3.clusterID.rds")

annot<-read.table("s3wgs_gcc.cellsummary.txt")
assay_annot<-annot$assay
sample_annot<-annot$sample
cluster_id<-as.factor(as.character(k_clus_id[match(annot$cellID,names(k_clus_id))]))
cluster_id_3<-as.factor(as.character(k_clus_id_3[match(annot$cellID,names(k_clus_id_3))]))

    #Read count raw
    plt1<-Heatmap(Y,
        show_row_names=F,
        show_column_names=F,
        column_order=1:ncol(Y),
        cluster_rows=dend,
        column_split=factor(unlist(lapply(strsplit(colnames(Y),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(Y),":"),"[",1)))),
        col=col_fun_reads_raw,
        row_title="Raw read count",
        name="Log10 Reads",
        left_annotation=rowAnnotation(assay=assay_annot,sample=sample_annot))

    #Read count corrected
    plt2<-Heatmap(Yhat,
        show_row_names=F,
        show_column_names=F,
        column_order=1:ncol(Yhat),
        cluster_rows=dend,
        column_split=factor(unlist(lapply(strsplit(colnames(Yhat),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(Yhat),":"),"[",1)))),
        col=col_fun_reads_normalized,
        row_title="Normalized read count",
        name="Log10 Reads",
        left_annotation=rowAnnotation(assay=assay_annot,sample=sample_annot))


    plt3<-Heatmap(iCN,
        show_row_names=F,
        show_column_names=F,
        column_order=1:ncol(iCN),
        cluster_rows=dend,
        column_split=factor(unlist(lapply(strsplit(colnames(iCN),"_"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(iCN),"_"),"[",1)))),
        col=col_fun_iCN,
        row_title="SCOPE Segmentation",  
        name="copy_number",
        left_annotation=rowAnnotation(assay=assay_annot,sample=sample_annot,cluster_3=cluster_id_3,cluster=cluster_id))

pdf("scope_final_scplot.pdf",width=10,height=3)
par(mfrow=c(3,1))
#plt1
#plt2
plt3
dev.off()

system("slack -F scope_final_scplot.pdf ryan_todo")

write.table(iCN,file="SourceData_Fig5b.tsv",sep="\t",row.names=T,quote=F)
system("slack -F SourceData_Fig5b.tsv ryan_todo")

```
{% endcapture %} {% include details.html %} 

## Coverage uniformity
Next generate MAPD and DIMAPD scores for all cells, using a custom script based on description from this website
from https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/interpret/metrics

MAPD, or the Median Absolute deviation of Pairwise Differences, is defined for a vector v[i] as follows:

    -construct the vector of consecutive pairwise differences scaled by the mean d[i] = ( v[i] - v[i+1] )/Mean(v),
    -compute the median absolute deviation of d[i] defined as Median(|d[i] - Median(d[i])|).
    -Performing genome coverage distribution on post GC corrected and mappability limited matrix (Y)


{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(SCOPE)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

#qcObj (raw reads, filtered bins)
qcObj<-readRDS("scope_qcObj.500kb.rds")
#normObj (normalization and bin-specific factors)
normObj.scope.gini<-readRDS("scope_normforeach.gini.500kb.rds")

#Sample information
cell_info<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/s3wgs_gcc.annot",header=T)

#CHANGE Read count matrix TO NORMALIZED Y VALUES
#Mad calculated on bins prenormalized with library size and bin-specific factor based on SCOPE paper
#Y(ij)/[N(j)B(j)]
Bj <- unlist(normObj.scope.gini$beta.hat)
Nj <- qcObj$QCmetric$mapq20
Y<-qcObj$Y
Y<-sweep(Y, MARGIN=1, Bj, `*`) #beta bin normalization
Y<-sweep(Y, MARGIN=2, Nj, `*`) #library size normalization

sampname <- qcObj$sampname #get cell name
ref <- qcObj$ref 

#custom mapd function
mapd<-function(cell){
  d<-unlist(lapply(1:nrow(Y)-1, function(x) (Y[x,cell]-Y[x+1,cell])/mean(Y[,cell])))
  mad<-median(abs(d-median(d)))
  return(mad)
}

library(parallel) #parallellize function per cell
mapd_list<-unlist(mclapply(colnames(Y),mapd,mc.cores=10)) #25 cores

dist_df<-data.frame("cellID"=colnames(Y),"mapd"=mapd_list) #making a MAPD matrix
dat<-merge(dist_df,cell_info,by="cellID")
qcmetric<-qcObj$QCmetric
qcmetric$cellID<-row.names(qcmetric)
dat<-merge(dat,qcmetric,by="cellID")

write.table(dat,file="s3wgs_gcc.cellsummary.txt",col.names=T,sep="\t",quote=F)

library(ggplot2)
plt<-ggplot(dat,aes(x=paste(assay,sample),y=mapd,color=paste(assay,sample)))+geom_jitter()+geom_boxplot(aes(fill=NULL))+theme_bw()+ylim(c(0,1))
ggsave(plt,file="mapd_scores.pdf")

#system("slack -F mapd_scores.pdf ryan_todo")
```
{% endcapture %} {% include details.html %} 

Generation of a cell summary file.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
# Group-wise ploidy initialization
#Updating annotation file
library(ggplot2)
library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(patchwork)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
bamfolder <- "/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/singlecell_bam"

annot<-read.table("s3wgsgcc_cellsummary.500kb.tsv",header=T)
annot$cellID_idx<-unlist(lapply(strsplit(annot$cellID,"_"),"[",3))

wgs_cell_summary<-read.table("s3wgs_cell_summary.tsv",header=T) #contains cell line information for wgs
#location of file to be changed later
gcc_cell_summary<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/for_gurkan/s3gcc_cellsummary.tsv",header=T)

wgs_cell_summary<-wgs_cell_summary[c("cellID","i5_idx_seq","i5_idx_name",
                                     "i7_idx_seq","tn5_idx_seq","total_reads_q10",
                                     "uniq_reads_q10","perc_uniq","plate",
                                     "sample.y","assay")]
gcc_cell_summary<-gcc_cell_summary[c("cellID","i5_idx_seq","i5_idx_name",
                                     "i7_idx_seq","tn5_idx_seq","total_reads",
                                     "uniq_reads","perc_uniq","plate",
                                     "sample")]
gcc_cell_summary$assay<-"GCC"
colnames(wgs_cell_summary)<-c("cellID_idx","i5_idx_seq","i5_idx_name","i7_idx_seq","tn5_idx_seq","total_reads_q10","uniq_reads_q10","perc_uniq","plate","sample","assay")
colnames(gcc_cell_summary)<-c("cellID_idx","i5_idx_seq","i5_idx_name","i7_idx_seq","tn5_idx_seq","total_reads_q10","uniq_reads_q10","perc_uniq","plate","sample","assay")

cell_summary<-rbind(wgs_cell_summary,gcc_cell_summary)
cell_summary<-cell_summary[!duplicated(cell_summary$cellID),]
annot<-merge(annot,cell_summary,by="cellID_idx")
write.table(annot,"s3wgsgcc_cellsummary.500kb.tsv",col.names=T)


```
{% endcapture %} {% include details.html %} 

## Additional complexity plots.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
#annotation for cell line processing
#annotation files stored on google drive
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

compl_J<-read.table("../complexity_data/CellLine_WGS_J.complexity.txt",header=F)
compl_K<-read.table("../complexity_data/CellLine_WGS_K.complexity.txt",header=F)
compl_E<-read.table("../complexity_data/CRC4442_WGS_E.complexity.txt",header=F)
compl_F<-read.table("../complexity_data/CRC4442_WGS_F.complexity.txt",header=F)
compl_G<-read.table("../complexity_data/CRC4671_WGS_G.complexity.txt",header=F)
compl_B<-read.table("../complexity_data/GM12878_WGS_B.complexity.txt",header=F)
compl_D<-read.table("../complexity_data/GM12878_WGS_D.complexity.txt",header=F)

colnames(compl_J)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_K)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_E)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_F)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_G)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_B)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
colnames(compl_D)<-c("row.number","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")

tn5_annot<-read.table("s3WGSandGCC_cellline_tn5_annotation.txt",header=T)
pcr_annot<-read.table("s3WGS_pcr_annotation.txt",header=T)

dat<-rbind(compl_J,compl_K,compl_E,compl_F,compl_G,compl_B,compl_D)

dat$i7_idx_seq<-substr(dat$cellID,0,8)
dat$i5_idx_seq<-substr(dat$cellID,9,16)
dat$tn5_idx_seq<-substr(dat$cellID,17,24)

dat<-merge(dat,pcr_annot,by="i5_idx_seq")

dat_notcellline<-dat[dat$sample!="cellline",]
dat_cellline<-dat[dat$sample=="cellline",]
dat_cellline<-merge(dat_cellline,tn5_annot,by="tn5_idx_seq")

dat_cellline<-dat_cellline[c("i5_idx_seq","cellID","total_reads_q10","uniq_reads_q10","perc_uniq","i7_idx_seq","tn5_idx_seq","i5_idx_name","i5_idx_cycle","sample.x","assay","plate","sample.y")]
dat_notcellline<-dat_notcellline[c("i5_idx_seq","cellID","total_reads_q10","uniq_reads_q10","perc_uniq","i7_idx_seq","tn5_idx_seq","i5_idx_name","i5_idx_cycle","sample","assay","plate")]
dat_notcellline$sample.y<-dat_notcellline$sample
colnames(dat_notcellline)<-c("i5_idx_seq","cellID","total_reads_q10","uniq_reads_q10","perc_uniq","i7_idx_seq","tn5_idx_seq","i5_idx_name","i5_idx_cycle","sample.x","assay","plate","sample.y")

dat<-rbind(dat_cellline,dat_notcellline)

#filtering data frame by listed cellIDs
celline_cellID_filt<-read.table("s3wgs_cellline.bbrd.q10.filt.cellIDs.list",header=F)
crc_cellID_filt<-read.table("s3wgs_crc.bbrd.q10.filt.cellIDs.list",header=F)
gm12878_cellID_filt<-read.table("s3wgs_gm12878.bbrd.q10.filt.cellIDs.list",header=F)

cellID_filt<-rbind(celline_cellID_filt,crc_cellID_filt,gm12878_cellID_filt)
colnames(cellID_filt)<-c("cellID")

dat<-dat[dat$cellID %in% cellID_filt$cellID,]

write.table(dat,file="s3wgs_cell_summary.tsv",sep="\t",quote=F,col.names=T,row.names=F)

dat<-read.table("s3wgs_cell_summary.tsv",header=T)
library(dplyr)
dat %>% group_by(plate,sample.y) %>% summarize(mean=mean(uniq_reads_q10),median=median(uniq_reads_q10))

library(ggplot2)

plt<-ggplot(dat,aes(y=log10(uniq_reads_q10),x=paste(dat$plate,dat$sample.y),color=perc_uniq))+
geom_jitter()+
geom_boxplot(outlier.shape=NA)+
theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="s3wgs_complexity_saturation_boxplot.pdf")

plt<-ggplot(dat,aes(y=log10(uniq_reads_q10),x=paste(dat$plate,dat$sample.y),color=paste(dat$plate,dat$sample.y)))+
geom_jitter()+
geom_boxplot(outlier.shape=NA)+
theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="s3wgs_complexity_platecellline_boxplot.pdf")
```
{% endcapture %} {% include details.html %} 

## Now plotting those QC metrics out with ggplot and R
This is a big old copy and paste script and can probably be parsed down by user defined functions.
{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(reshape2)
###This is all kind of a mess since there is a lot of copy paste code and reassigning the same variable names

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200219_s3WGS_CRC_reprocessing/preprocessing")

s3_s3GCC_4442<-read.table("hg38.s3GCC_4442.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4442)<-c("row_carryover","tot_reads","uniq_reads","perc_uniq")
s3_s3GCC_4442$sample<-c("4442")
s3_s3GCC_4442$assay<-c("GCC")

s3_s3GCC_4671<-read.table("hg38.s3GCC_4671.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4671)<-c("row_carryover","tot_reads","uniq_reads","perc_uniq")
s3_s3GCC_4671$sample<-c("4671")
s3_s3GCC_4671$assay<-c("GCC")

s3_s3WGS_4442<-read.table("hg38.s3WGS_4442.complexity.txt",header=F,row.names=2)
colnames(s3_s3WGS_4442)<-c("row_carryover","tot_reads","uniq_reads","perc_uniq")
s3_s3WGS_4442$sample<-c("4442")
s3_s3WGS_4442$assay<-c("WGS")

s3_s3WGS_4671<-read.table("hg38.s3WGS_4671.complexity.txt",header=F,row.names=2)
colnames(s3_s3WGS_4671)<-c("row_carryover","tot_reads","uniq_reads","perc_uniq")
s3_s3WGS_4671$sample<-c("4671")
s3_s3WGS_4671$assay<-c("WGS")

proj_uniq_reads <- function(x,df,proj_perc_num){
df_temp<-df[df$cellID==x,]
val<-df_temp[which.min(abs(proj_perc_num-df_temp$proj_perc)),]$proj_uniq
return(val)
}

s3_s3GCC_4442_proj<-read.table("./hg38_prefilt.s3GCC.s3GCC_4442.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4442_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4442_proj_est<-data.frame(cellID=unique(s3_s3GCC_4442_proj$cellID),
uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4442_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_proj,proj_perc_num=0.5)),
uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4442_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_proj,proj_perc_num=0.05)))

s3_s3GCC_4671_proj<-read.table("./hg38_prefilt.s3GCC.s3GCC_4671.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4671_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4671_proj_est<-data.frame(cellID=unique(s3_s3GCC_4671_proj$cellID),
uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4671_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_proj,proj_perc_num=0.5)),
uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4671_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_proj,proj_perc_num=0.05)))


s3_s3WGS_4442_proj<-read.table("./hg38.s3WGS_4442.read_projections/cell_projections.txt",header=F)
colnames(s3_s3WGS_4442_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3WGS_4442_proj_est<-data.frame(cellID=unique(s3_s3WGS_4442_proj$cellID),
uniq_reads_50=unlist(lapply(X=unique(s3_s3WGS_4442_proj$cellID),proj_uniq_reads,df=s3_s3WGS_4442_proj,proj_perc_num=0.5)),
uniq_reads_95=unlist(lapply(X=unique(s3_s3WGS_4442_proj$cellID),proj_uniq_reads,df=s3_s3WGS_4442_proj,proj_perc_num=0.05)))


s3_s3WGS_4671_proj<-read.table("./hg38_prefilt.s3WGS_4671.read_projections/cell_projections.txt",header=F)
colnames(s3_s3WGS_4671_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3WGS_4671_proj_est<-data.frame(cellID=unique(s3_s3WGS_4671_proj$cellID),
uniq_reads_50=unlist(lapply(X=unique(s3_s3WGS_4671_proj$cellID),proj_uniq_reads,df=s3_s3WGS_4671_proj,proj_perc_num=0.5)),
uniq_reads_95=unlist(lapply(X=unique(s3_s3WGS_4671_proj$cellID),proj_uniq_reads,df=s3_s3WGS_4671_proj,proj_perc_num=0.05)))

#Also adding GM12878 to this


s3_s3WGS_GM12878<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/191118_sciWG_96plex/191118_sciWG_96plex.complexity.txt",header=F,row.names=2)
colnames(s3_s3WGS_GM12878)<-c("row_carryover","tot_reads","uniq_reads","perc_uniq")
s3_s3WGS_GM12878$sample<-c("GM12878")
s3_s3WGS_GM12878$assay<-c("WGS")


s3_s3WGS_GM12878_proj<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/191118_sciWG_96plex/191118_sciWG_96plex.maxPopulation.filt.read_projections/cell_projections.txt",header=F)
colnames(s3_s3WGS_GM12878_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3WGS_GM12878_proj_est<-data.frame(cellID=unique(s3_s3WGS_GM12878_proj$cellID),
uniq_reads_50=unlist(lapply(X=unique(s3_s3WGS_GM12878_proj$cellID),proj_uniq_reads,df=s3_s3WGS_GM12878_proj,proj_perc_num=0.5)),
uniq_reads_95=unlist(lapply(X=unique(s3_s3WGS_GM12878_proj$cellID),proj_uniq_reads,df=s3_s3WGS_GM12878_proj,proj_perc_num=0.05)))

#And finally the distal GCC reads

s3_s3GCC_4671_distal<-read.table("hg38.s3GCC_4671.distal.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4671_distal)<-c("row_carryover","dist_tot_reads","dist_uniq_reads","dist_perc_uniq")
s3_s3GCC_4671_distal$sample<-c("4671")
s3_s3GCC_4671_distal$assay<-c("GCC")

s3_s3GCC_4442_distal<-read.table("hg38.s3GCC_4442.distal.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4442_distal)<-c("row_carryover","dist_tot_reads","dist_uniq_reads","dist_perc_uniq")
s3_s3GCC_4442_distal$sample<-c("4442")
s3_s3GCC_4442_distal$assay<-c("GCC")

s3_s3GCC_4671_transchr<-read.table("hg38.s3GCC_4671.transchr.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4671_transchr)<-c("row_carryover","transchr_tot_reads","transchr_uniq_reads","transchr_perc_uniq")
s3_s3GCC_4671_transchr$sample<-c("4671")
s3_s3GCC_4671_transchr$assay<-c("GCC")

s3_s3GCC_4442_transchr<-read.table("hg38.s3GCC_4442.transchr.complexity.txt",header=F,row.names=2)
colnames(s3_s3GCC_4442_transchr)<-c("row_carryover","transchr_tot_reads","transchr_uniq_reads","transchr_perc_uniq")
s3_s3GCC_4442_transchr$sample<-c("4442")
s3_s3GCC_4442_transchr$assay<-c("GCC")

s3_s3GCC_4442_distal_proj<-read.table("./hg38.s3GCC_4442.distal.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4442_distal_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4442_distal_proj_est<-data.frame(cellID=unique(s3_s3GCC_4442_distal_proj$cellID),type=c("distal"),
dist_uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4442_distal_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_distal_proj,proj_perc_num=0.5)),
dist_uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4442_distal_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_distal_proj,proj_perc_num=0.05)))

s3_s3GCC_4671_distal_proj<-read.table("./hg38.s3GCC_4671.distal.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4671_distal_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4671_distal_proj_est<-data.frame(cellID=unique(s3_s3GCC_4671_distal_proj$cellID),type=c("distal"),
dist_uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4671_distal_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_distal_proj,proj_perc_num=0.5)),
dist_uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4671_distal_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_distal_proj,proj_perc_num=0.05)))

s3_s3GCC_4442_transchr_proj<-read.table("./hg38.s3GCC_4442.transchr.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4442_transchr_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4442_transchr_proj_est<-data.frame(cellID=unique(s3_s3GCC_4442_transchr_proj$cellID),type=c("transchr"),
trans_uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4442_transchr_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_transchr_proj,proj_perc_num=0.5)),
trans_uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4442_transchr_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4442_transchr_proj,proj_perc_num=0.05)))

s3_s3GCC_4671_transchr_proj<-read.table("./hg38.s3GCC_4671.transchr.read_projections/cell_projections.txt",header=F)
colnames(s3_s3GCC_4671_transchr_proj)<-c("cellID","proj_depth","proj_tot","proj_uniq","proj_perc")
s3_s3GCC_4671_transchr_proj_est<-data.frame(cellID=unique(s3_s3GCC_4671_transchr_proj$cellID),type=c("transchr"),
trans_uniq_reads_50=unlist(lapply(X=unique(s3_s3GCC_4671_transchr_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_transchr_proj,proj_perc_num=0.5)),
trans_uniq_reads_95=unlist(lapply(X=unique(s3_s3GCC_4671_transchr_proj$cellID),proj_uniq_reads,df=s3_s3GCC_4671_transchr_proj,proj_perc_num=0.05)))


dat<-rbind(s3_s3GCC_4442,s3_s3GCC_4671,s3_s3WGS_4442,s3_s3WGS_4671,s3_s3WGS_GM12878)
dat_proj<-rbind(s3_s3GCC_4442_proj_est,s3_s3GCC_4671_proj_est,s3_s3WGS_4442_proj_est,s3_s3WGS_4671_proj_est,s3_s3WGS_GM12878_proj_est)

dat$cellID<-row.names(dat)
dat_m<-merge(dat,dat_proj,by="cellID")
dat_m<-dat_m[,c("uniq_reads","uniq_reads_50","uniq_reads_95","sample","assay","cellID")]
dat_m<-dat_m[dat_m$uniq_reads>=100000,]


dat_m<-melt(dat_m)

dat_m_wgs<-dat_m[dat_m$assay=="WGS",]

ggplot(dat_m_wgs,aes(x=as.factor(paste(variable,assay,sample)),y=log10(value),color=paste(assay,sample)))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3WGS_projected_complexity.svg")

dat_m_gcc<-dat_m[dat_m$assay=="GCC",]

ggplot(dat_m_gcc,aes(x=as.factor(paste(variable,assay,sample)),y=log10(value),color=as.factor(paste(assay,sample))))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3GCC_projected_complexity.svg")

#And now distal and trans read projections


dat<-rbind(s3_s3GCC_4442,s3_s3GCC_4671)
dat_proj<-rbind(s3_s3GCC_4442_proj_est,s3_s3GCC_4671_proj_est)

dat$cellID<-row.names(dat)
dat_m<-merge(dat,dat_proj,by="cellID")
dat_m<-dat_m[,c("uniq_reads","uniq_reads_50","uniq_reads_95","sample","assay","cellID")]
dat_m<-dat_m[dat_m$uniq_reads>=100000,]
dat_m_gcc<-dat_m[dat_m$assay=="GCC",]



dat<-rbind(s3_s3GCC_4671_distal,s3_s3GCC_4442_distal)
dat_proj<-rbind(s3_s3GCC_4671_distal_proj_est,s3_s3GCC_4442_distal_proj_est)

dat$cellID<-row.names(dat)
dat_m_dist<-merge(dat,dat_proj,by="cellID") #merge to get cis

dat<-rbind(s3_s3GCC_4442_transchr,s3_s3GCC_4671_transchr)
dat_proj<-rbind(s3_s3GCC_4671_transchr_proj_est,s3_s3GCC_4442_transchr_proj_est)
dat$cellID<-row.names(dat)
dat_m_trans<-merge(dat,dat_proj,by="cellID") #merge to get trans

cell_accepted<-c(row.names(s3_s3GCC_4442[s3_s3GCC_4442$uniq_reads>=100000,]),row.names(s3_s3GCC_4671[s3_s3GCC_4671$uniq_reads>=100000,]))
dat_m_dist<-dat_m_dist[dat_m_dist$cellID %in% cell_accepted,] #filter distal reads to same cellIDs are unique defined by total reads
dat_m_trans<-dat_m_trans[dat_m_trans$cellID %in% cell_accepted,] #filter trans reads to same cellIDs are unique defined by total reads

dat_gcc_final<-merge(dat_m_gcc,dat_m_dist,by="cellID")
dat_gcc_final<-merge(dat_gcc_final,dat_m_trans,by="cellID")


dat_gcc_final<-dat_gcc_final[,c("cellID","sample","assay","uniq_reads","uniq_reads_50","uniq_reads_95","dist_tot_reads","dist_uniq_reads","dist_perc_uniq","dist_uniq_reads_50","dist_uniq_reads_95","transchr_tot_reads", "transchr_uniq_reads","transchr_perc_uniq","trans_uniq_reads_50", "trans_uniq_reads_95"  )]
write.table(dat_gcc_final,file="s3GCC_cellQC.txt",sep="\t",col.names=T,row.names=F,quote=F)


dat_gcc_trans<-dat_gcc_final[,c("cellID","sample","assay","uniq_reads","transchr_uniq_reads" )]
dat_gcc_dist<-dat_gcc_final[,c("cellID","sample","assay","uniq_reads","dist_uniq_reads")]

ggplot(dat_gcc_dist,aes(x=as.factor(paste(assay,sample)),y=dist_uniq_reads/uniq_reads,color=as.factor(paste(assay,sample))))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3GCC_dist_percentreads_complexity.svg")

ggplot(dat_gcc_trans,aes(x=as.factor(paste(assay,sample)),y=transchr_uniq_reads/uniq_reads,color=as.factor(paste(assay,sample))))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3GCC_trans_percentreads_complexity.svg")

dat_gcc_trans<-dat_gcc_final[,c("cellID","sample","assay","transchr_uniq_reads","trans_uniq_reads_50", "trans_uniq_reads_95"  )]
dat_gcc_dist<-dat_gcc_final[,c("cellID","sample","assay","dist_uniq_reads","dist_perc_uniq","dist_uniq_reads_50","dist_uniq_reads_95"  )]

dat_gcc_dist<-melt(dat_gcc_dist)
ggplot(dat_gcc_dist,aes(x=as.factor(paste(variable,assay,sample)),y=log10(value),color=as.factor(paste(assay,sample))))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3GCC_dist_projected_complexity.svg")

dat_gcc_trans<-melt(dat_gcc_trans)
ggplot(dat_gcc_trans,aes(x=as.factor(paste(variable,assay,sample)),y=log10(value),color=as.factor(paste(assay,sample))))+geom_jitter()+geom_boxplot()+theme_bw()
ggsave(file="s3GCC_trans_projected_complexity.svg")

for (i in unique(dat_m$sample)){
	print(i)
	print(nrow(dat_m[dat_m$sample==i,]))
	print(summary(dat_m[dat_m$sample==i,]$uniq_reads_95))
}

```
{% endcapture %} {% include details.html %} 


## Merging Cells for Clade Analysis and Plotting
Combing clades of cell lines for high resolution annotation using 50kb windows

First reading in libraries and setting up functions.

{% capture summary %} Code {% endcapture %} {% capture details %}  
```R
library(SCOPE)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(philentropy)
library(ape)
library(ggdendro)
library(dendextend)
library(rpart)
library(cluster)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)
library(ggrepel)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")


#set up single cell bam files from given directory
prepare_bam_bed_obj<-function(resol=resol){
    #Initalization
    bamfolder <- "/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/singlecell_bam/sorted_bams"
    bamFile <- list.files(bamfolder, pattern = 'sorted.bam$')
    bamdir <- file.path(bamfolder, bamFile)
    sampname_raw <- paste(sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1),sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 3),sep=".")
    #set genomic window size to 500kb
    bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38",resolution=resol,sex=T)#resolution of 100 = 100kbp

    #Compute GC content and mappability for each reference bin.
    mapp <- get_mapp(bambedObj$ref, hgref = "hg38")
    gc <- get_gc(bambedObj$ref, hgref = "hg38")
    values(bambedObj$ref) <- cbind(values(bambedObj$ref), DataFrame(gc, mapp))
    return(bambedObj)
}

#modified scope get_coverage_scDNA
read_in_reads<-function(i,seq="paired-end"){
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

```

Now collating reads to 50kbp regions, normalizing and segmenting genomes.

```R
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

#Read in annotation files
annot<-read.table("s3wgs_gcc.cellsummary.txt")
resol=50
bambedObj<-prepare_bam_bed_obj(resol=50) #500 is 500kbp
saveRDS(bambedObj,file="bambedObj.50kbp.rds")
bambedObj<-readRDS(file="bambedObj.50kbp.rds")
ref <- bambedObj$ref
bamdir <- bambedObj$bamdir
sampname <- bambedObj$sampname

#Get read count per bin per cell
Y_out<-mclapply(1:length(sampname),read_in_reads,mc.cores=20) #can be set to "single-end" for single-end reads
Y<-do.call("cbind",Y_out)
colnames(Y)<-sampname
ref<-as.data.frame(ref)
row.names(Y)<-paste0(ref$seqnames,":",ref$start,"_",ref$end)
saveRDS(Y,file="rawcount.50kbp.rds")
Y<-readRDS(file="rawcount.50kbp.rds")

#Apriori Ploidy Estimate
annot<-read.table("s3wgs_gcc.cellsummary.txt")
#reorder Y to annot
Y<-Y[,which(sampname %in% annot$cellID)]
#qcObj<-readRDS("scope_qcObj.500kb.rds")
dend<-readRDS("scope_segmentcs.500kb_k2.maxns2.dendrogram.rds")
k_clus_id<-readRDS("scope_segmentcs.500kb_k2.maxns2.clusterID.rds")
ploidy<-readRDS("scope_ploidy.apriori.500kb.rds")
ploidy<-ploidy$ploidy

#determine which is gm12878 clade
for (i in unique(k_clus_id)){
    clade<-names(k_clus_id[k_clus_id==i])
    print(paste(i,length(clade),sum(startsWith(clade,"s3wgs_gm12878"))/length(clade),mean(ploidy[colnames(Y) %in% clade],na.rm=T)))
}
#clade 1

clade_dat<-list()
for (i in unique(k_clus_id)){
    clade<-names(k_clus_id[k_clus_id==i])
    clade_combined<-rowMeans(Y[,colnames(Y) %in% clade],na.rm=T)
    clade_dat[[paste0("clade_",i)]]<-clade_combined
}

dat<-do.call("cbind",clade_dat)
colnames(dat)<-names(clade_dat)


#clade ploidy
clade_ploidy<-list()
for (i in unique(k_clus_id)){
clade<-names(k_clus_id[k_clus_id==i])
clade_ploidy[[paste0("clade_",i)]]<-median(ploidy[colnames(Y) %in% clade],na.rm=T)
}
clade_ploidy

#manually filter bins for mappability and gc content (using the scope defaults)
ref<-as.data.frame(ref)
ref<-ref[ref$mapp>0.9 & ref$gc>=20 & ref$gc<=80,]
dat<-dat[paste0(ref$seqnames,":",ref$start,"_",ref$end),]

#normalize clade data, clade 6 contains most gm12878 cells and is considered the normal index
normObj<- normalize_codex2_ns_noK(Y_qc = dat, gc_qc = ref$gc, norm_index = 6)
saveRDS(normObj,"scope.50kb.normalizedobj.mean.rds")
normObj<-readRDS("scope.50kb.normalizedobj.mean.rds")


#Perform CBS originally run with ns=2
chr_cbs <- function(x,scope=normObj,ns=1) {
    print(paste("Running for ",x))
    chr_seg<-segment_CBScs(Y = dat,
    Yhat = as.data.frame(scope$Yhat),
    sampname = colnames(dat),
    ref = makeGRangesFromDataFrame(ref),
    chr = x,
    mode = "integer", max.ns = ns)
    return(chr_seg)
}

chrs <- unique(as.character(seqnames(bambedObj$ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs

segment_cs<-list()
segment_cs<-mclapply(chrs,FUN=chr_cbs,mc.cores=25)
names(segment_cs) <- chrs #mclapply returns jobs in same order as specified

#saveRDS(segment_cs,"scope_segmentcs_clade_50kb.rds")
#saveRDS(segment_cs,"scope_segmentcs_clade_100kb.ns1.rds")
saveRDS(segment_cs,"scope_segmentcs_clade_50kb.ns1.rds")
segment_cs<-readRDS("scope_segmentcs_clade_50kb.ns1.rds")

qcObj<-readRDS("scope_qcObj.1mb.rds")
ploidy<-readRDS("scope_ploidy.gini.1mb.rds")
norm_index<-which(grepl("gm12878",colnames(qcObj$Y)))
normObj<-readRDS("scope_noKnorm.gini.1mb.rds")

# #Normalize with Latent Factors
# normObj.scope <- normalize_scope_foreach(
#     Y_qc = dat,
#     gc_qc = as.data.frame(ref)$gc,
#     K = 1,
#     ploidyInt = unlist(clade_ploidy),
#     norm_index = 6,
#     T = 1:6,
#     beta0 = normObj$beta.hat)
    
# saveRDS(normObj.scope,"scope_normforeach.50kb.rds")

#saveRDS(segment_cs,"scope_segmentcs_clade_50kb.ns6.rds")
segment_cs<-readRDS("scope_segmentcs_clade_50kb.rds")

#remove chrY, it failed because the read count is low
segment_cs<-segment_cs[1:23]
bin_loci<-data.frame(ref)
bin_loci<-bin_loci[bin_loci$seqnames!="chrY",]
bin_loci$bin_name<-row.names(bin_loci)
bin_loci_subject<-makeGRangesFromDataFrame(bin_loci,keep.extra.columns=T)
iCN <- do.call(rbind, lapply(segment_cs, function(z){names(z)}))
iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
row.names(iCN)<-paste(bin_loci$seqnames,bin_loci$start,bin_loci$end,sep="_")
iCN[which(iCN>5,arr.ind=T)]<-5
iCN<-as.data.frame(t(iCN))


iCN<-iCN[1:6,]


library(ggplot2)
library(patchwork)
library(reshape2)

image.orig<-NULL
image.orig <- do.call(rbind, lapply(segment_cs, function(z){names(z)}))
image.orig <- do.call(rbind, lapply(segment_cs, function(z){z[["image.orig"]]}))
row.names(image.orig)<-paste(bin_loci$seqnames,bin_loci$start,bin_loci$end,sep="_")
colnames(image.orig)<-colnames(dat)
plt_melt<-melt(image.orig)
plt_melt$gloc<-row.names(image.orig)
plt_melt<-cbind(plt_melt,melt(t(iCN)))
plt_melt<-plt_melt[,c(1,2,3,6)]
colnames(plt_melt)<-c("clade","copy","gloc","state")
plt_melt$row_order<-1:nrow(bin_loci)

cols<-c("0"="#2166ac", "1"="#d0e7f5", "2"="#f5f5f5","3"="#fddbc7","4"="#ef8a62","5"="#b2182b")
#plot bins for a cell across a chromosome

range_gc<-quantile(plt_melt$copy, na.rm = TRUE, prob = c(0.05,0.95))
plt_melt$contig<-unlist(lapply(strsplit(plt_melt$gloc,"_"),"[",1))
plt_melt$contig<-factor(plt_melt$contig,level=unique(plt_melt$contig))

plt_melt$start<-as.numeric(unlist(lapply(strsplit(plt_melt$gloc,"_"),"[",2)))
plt_melt$end<-as.numeric(unlist(lapply(strsplit(plt_melt$gloc,"_"),"[",3)))
plt_melt$gene<-c("")
plt_melt_6<-plt_melt[plt_melt$clade=="clade_6",] #just adding annotation track to clade 6
plt_melt<-plt_melt[plt_melt$clade!="clade_6",]

#Making another track to plot PDAC linked genes
pdac_genes<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/pdac_gene.bed",header=F)
colnames(pdac_genes)<-c("chr","start","end","gene")
pdac_genes<-makeGRangesFromDataFrame(pdac_genes,keep.extra.columns=T)
pdac_genes_idx<-as.data.frame(findOverlaps(pdac_genes,makeGRangesFromDataFrame(ref))) #overlap GRanges for plotting consistent matrix bins
pdac_genes_idx<-pdac_genes_idx[!duplicated(pdac_genes_idx$queryHits),]
pdac_genes_idx<-cbind(pdac_genes_idx,pdac_genes)
pdac_genes_idx<-merge(bin_loci,pdac_genes_idx,by.x="bin_name",by.y="subjectHits")
pdac_genes_idx$gloc<-paste(pdac_genes_idx$seqnames.x,pdac_genes_idx$start.x,pdac_genes_idx$end.x,sep="_")
plt_melt_6[match(pdac_genes_idx$gloc,plt_melt_6$gloc,nomatch=0),]$gene<-pdac_genes_idx$gene

plt_melt<-rbind(plt_melt,plt_melt_6)

library(ggrepel)

cols_clade<-c("na"="grey","clade_1"="#cccc99","clade_2"="#ffcc33","clade_3"="#ff99cc","clade_4"="#ff9966","clade_5"="#66cc99","clade_6"="#cccccc")

plt<-ggplot(plt_melt,aes(x=row_order,y=(2^(copy)),label=gene))+
geom_rect(aes(fill=as.character(state),xmin=row_order,xmax=row_order+1,ymin=0,ymax=6,alpha=0.01))+
scale_fill_manual(values=cols)+
#geom_hline(yintercept=c(0,1,2,3,4,5),linetype="dashed")+
geom_point(color="black",size=1,alpha=0.2)+
scale_color_manual(values=cols_clade)+
ylab("")+
xlab("")+
geom_text_repel(size=20,y=1,nudge_y= 2, direction="x",angle= 90, hjust= 0,segment.size = 2,max.iter = 1e4,min.segment.length=0)+
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



ggsave(plt,file="scope.merged.50kb.final.final.png",width=2000,height=1000,units="mm",limitsize=F)
system("slack -F scope.merged.50kb.final.final.png ryan_todo")

```
