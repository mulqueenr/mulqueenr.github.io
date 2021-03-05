---
title: s3WGS Analysis
layout: s3data
author: Ryan Mulqueen
permalink: /s3wgs/
category: s3processing
---

# Processing for s3WGS portion of s3 paper.
This notebook is a continuation of "s3 Preprocessing to Bam files" and details the processing of s3WGS libraries. This notebook starts with a merged, barcode-based removal of duplicate, >Q10 filtered bam file which was subset to barcodes >=10K unique reads


```python
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

# Filter bam files to just cells passing original QC
Performing this on pre-barcode based remove duplicate bams for future complexity plotting and projection modeling.


```python
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

### QC Directory Structure 


```python
/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data
├── singlecell_bam #contains deduplicated single cell bams for SCOPE analysis
│   └── isize #contains isize output for insert size distribution
└── single_cell_projections #contains pre-dedupliacted single cell bams
    ├── complexity #contains deduplication complexity logs per cell
    ├── projections #contains projections of complexity per cell
    ├── undedup_bams #contains pre-bbrd bams
    
```

### Collate data into single files for plotting
 


```python

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

R
dat<-read.table("s3wgs_complexity.txt",header=F)
annot<-read.table("s3wgs_gcc.annot",header=T)
colnames(dat)<-c("cellID_idx","total_reads","uniq_reads","perc_uniq","cellID")
dat<-merge(dat,annot,by="cellID_idx")
write.table(dat,file="s3wgs_complexity.txt",sep="\t",col.names=T,quote=F)

#isize
for i in ./singlecell_bam/isize/*values;
do cellid=${i:23:-22};
awk -v cellid=$cellid 'OFS="\t" {print $1,cellid}' $i;
done > ./s3wgs_isize.txt 

```

Generate complexity projections per cell


```python
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


dat[dat$proj_perc_uniq==0.05,]
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


```python
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
annot<-read.table("s3wgs_gcc.cellsummary.txt",header=T)

ggplot(annot,aes(x=paste(assay,sample),y=mapd,color=paste(assay,sample)))+geom_jitter()+geom_boxplot()+ylim(c(0,1))+theme_bw()
ggsave(file="mapd.png")
ggsave(file="mapd.pdf")

system("slack -F mapd.pdf ryan_todo")
```

### Generate Insert Size Distributions 


```python
R
library(ggplot2)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

dat<-read.table("s3wgs_isize.txt",header=F)
colnames(dat)<-c("size","cellID")
dat$cellID_idx<-unlist(lapply(strsplit(dat$cellID,"RG_"),"[",2))
saveRDS(dat,"s3wgs_isize.rds")

annot<-read.table("s3wgs_gcc.cellsummary.txt",header=T)
annot<-annot[c("cellID_idx","assay","sample")]
dat<-merge(dat,annot,by="cellID_idx")
saveRDS(dat,"s3wgs_isize.annot.rds")


```

## Using SCOPE to analyze single-cell data on single cell bam directory
Library used for scWGS data analysis is [SCOPE](https://github.com/rujinwang/SCOPE) available as a [preprint](https://www.biorxiv.org/content/10.1101/594267v1.full). SCOPE works on pre-aligned deduplicated bam files. So I split files post-deduplication into a subdirectory to load in (above).

Using R 4.0.0
Note a lot of this code and even the comments and explanation is directly taken from the SCOPE example given.


### Read in files from directory
First reading in the split bam files and setting up the reference genome.


```python
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
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, hgref = "hg38",resolution=250,sex=T)#resolution of 100 = 100kbp

#Compute GC content and mappability for each reference bin.
mapp <- get_mapp(bambedObj$ref, hgref = "hg38")
gc <- get_gc(bambedObj$ref, hgref = "hg38")
values(bambedObj$ref) <- cbind(values(bambedObj$ref), DataFrame(gc, mapp))

#For 500kb bins
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 10, seq = 'paired-end', hgref = "hg38") #using a Q10 score for filtering

saveRDS(coverageObj,"scope_covobj.250kb.rds")
saveRDS(bambedObj,"scope_bambedObj.250kb.rds")
```


```python
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

coverageObj<-readRDS("scope_covobj.1mb.rds")
bambedObj<-readRDS("scope_bambedObj.1mb.rds")
#Perform QC and Ploidy Estimate
#For 500kb bins
#Perform QC
#Quality control to remove samples/cells with low proportion of mapped reads, bins that have extreme GC content (less than 20% and greater than 80%) and low mappability (less than 0.9) to reduce artifacts.
QCmetric_raw <- get_samp_QC(bambedObj)
saveRDS(QCmetric_raw,"scope_qcmetric.250kb.rds")
QCmetric_raw<-readRDS("scope_qcmetric.1mb.rds")

qcObj <- perform_qc(Y_raw = coverageObj$Y,sampname_raw = row.names(QCmetric_raw), ref_raw = bambedObj$ref, QCmetric_raw = QCmetric_raw)
"""Removed 0 samples due to failed library preparation.
Removed 133 samples due to failure to meet min coverage requirement.
Removed 0 samples due to low proportion of mapped reads.
Excluded 276 bins due to extreme GC content.
Excluded 450 bins due to low mappability.
Removed 0 samples due to excessive zero read counts in
            library size calculation.
There are 1268 samples and 5449 bins after QC step."""
saveRDS(qcObj,"scope_qcObj.250kb.rds")

#Generate gini coeffcients for 500kb binned data
Gini<-get_gini(qcObj$Y) #SCOPE function
saveRDS(Gini,"scope_gini.1mb.rds")

#500 kb first pass normalization
#First pass at estimating ploidy using no latent factors
norm_index_gini<-which(grepl("gm12878",colnames(qcObj$Y)))#grabbing 20 cells lowest on gini index

# first-pass CODEX2 run with no latent factors
normObj_gini<- normalize_codex2_ns_noK(Y_qc = qcObj$Y, gc_qc = qcObj$ref$gc, norm_index = norm_index_gini)
saveRDS(normObj_gini,"scope_noKnorm.gini.1mb.rds")

normObj_gini<-readRDS("scope_noKnorm.gini.1mb.rds")
###Ploidy initialization 500kb
#Generate ploidy estimate per cell from gini index
ploidy <- initialize_ploidy(Y = qcObj$Y, Yhat = normObj_gini$Yhat, ref = qcObj$ref, SoS.plot = F)
saveRDS(ploidy,"scope_ploidy.gini.1mb.rds")
```

### Using ploidy estimates to assess copy number change


```python
library(SCOPE)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")
#Performing normalization via gini euploid estimates
qcObj<-readRDS("scope_qcObj.1mb.rds")
ploidy<-readRDS("scope_ploidy.gini.1mb.rds")
norm_index<-which(grepl("gm12878",colnames(qcObj$Y)))
normObj<-readRDS("scope_noKnorm.gini.1mb.rds")

#Normalize with no Latent Factors
normObj.scope.gini <- normalize_scope_foreach(
	Y_qc = qcObj$Y,
	gc_qc = qcObj$ref$gc,
	K = 1,
	ploidyInt = ploidy,
	norm_index = norm_index,
	T = 1:6,
	beta0 = normObj$beta.hat)
    
saveRDS(normObj.scope.gini,"scope_normforeach.gini.1mb.rds")
```


```python
library(SCOPE)


#Perform CBS
chr_cbs <- function(x,scope=normObj.scope.gini) {
    print(paste("Running for ",x))
    chr_seg<-segment_CBScs(Y = qcObj$Y,
    Yhat = as.data.frame(scope$Yhat[[which.max(scope$BIC)]]),
    sampname = colnames(qcObj$Y),
    ref = qcObj$ref,
    chr = x,
    mode = "integer", max.ns = 1)
    return(chr_seg)
}

###Processing gini index strategy
    #Read in files generated above
    qcObj<-readRDS("scope_qcObj.1mb.rds")
    ploidy<-readRDS("scope_ploidy.gini.1mb.rds")
    norm_index<-head(order(Gini),n=20) #grabbing 20 cells lowest on gini index
    normObj.scope.gini<-readRDS("scope_normforeach.gini.1mb.rds")

    #Running CBS segmentation on 500kb windows (gini ploidy estimates)
    chrs <- unique(as.character(seqnames(qcObj$ref)))
    segment_cs <- vector('list',length = length(chrs))
    names(segment_cs) <- chrs

    segment_cs<-list()
    segment_cs<-mclapply(chrs,FUN=chr_cbs,mc.cores=length(chrs))
    names(segment_cs) <- chrs #mclapply returns jobs in same order as specified
    saveRDS(segment_cs,"scope_segmentcs.gini.1mb.rds")
```


```python
library(SCOPE)
library(ComplexHeatmap)
library(dplyr)
library(ape)
library(RColorBrewer)
library(reshape2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data")

#PDAC common mutated genes 
pdac_genes<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200422_SearsCRCData/TCGA_PDAC_genes.txt")
colnames(pdac_genes)<-c("gene","chr","start","end")
pdac_genes<-makeGRangesFromDataFrame(pdac_genes,keep.extra.columns=T)

#Sample information
cell_info<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/s3wgs_gcc.annot",header=T)

###Processing gini index strategy
    #Read in files generated above
    qcObj<-readRDS("scope_qcObj.500kb.rds")
    segment_cs<-readRDS("scope_segmentcs.gini.500kb.rds")

#Inferred Copy Number 500kb
    bin_loci<-data.frame(qcObj$ref)
    bin_loci$bin_name<-row.names(bin_loci)
    bin_loci_subject<-makeGRangesFromDataFrame(bin_loci,keep.extra.columns=T)

    chrs <- unique(as.character(seqnames(qcObj$ref)))
    iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
    row.names(iCN)<-paste(bin_loci$seqnames,bin_loci$start,bin_loci$end,sep="_")


#Breakpoint data from 500kb bins
   bin_loci<-data.frame(qcObj$ref)
    bin_loci$bin_name<-row.names(bin_loci)

    chrs <- unique(as.character(seqnames(qcObj$ref)))
    breakpoint <-lapply(segment_cs, function(z){z[["finalcall"]]})
    breakpoint<-do.call("rbind",breakpoint)    
    breakpoint$chr<-unlist(lapply(strsplit(row.names(breakpoint),"[.]"),"[",1))
    breakpoint$bin<-paste0(breakpoint$chr,":",
                           breakpoint$st_bin,"-",
                           breakpoint$ed_bin)
    breakpoint<-breakpoint[c("sample_name","chr","bin","cnv_no")]
    breakpoint<-dcast(breakpoint,sample_name~bin,value.var="cnv_no")
    #uwot::umap(breakpoint)
            
    row.names(breakpoint)<-breakpoint$sample_name
    breakpoint<-breakpoint[,2:ncol(breakpoint)]
    breakpoint[which(is.na(breakpoint),arr.ind=T)]<-2
    breakpoint_nj<-nj(dist(breakpoint,method="canberra"))

#Set up annotations
    assay_annot<-unlist(lapply(strsplit(breakpoint_nj$tip.label,"_"),"[",1))
    sample_annot<-unlist(lapply(strsplit(breakpoint_nj$tip.label,"_"),"[",2))
    annot<-paste(assay_annot,sample_annot)
    levels(annot)<-unique(annot)
    annot_colors = as.data.frame(cbind(levels(annot),I(brewer.pal(nlevels(annot),name="Set1"))))
    
#    pdf("scope.gini.breakpoint_phylo.pdf",height=20,width=20)
#    plot.phylo(breakpoint_nj,type="fan",show.tip.label=T,tip.color=annot_colors[match(annot,annot_colors$V1),]$V2)
#    dev.off()
#    system("slack -F scope.gini.breakpoint_phylo.pdf ryan_todo")

#Plotting 500kb bins
    row.names(iCN)<-row.names(qcObj$Y)
    chr_list<-unlist(lapply(strsplit(row.names(iCN),":"),"[",1))

    cellid_order<-colnames(qcObj$Y)
    assay_annot<-setNames(cell_info[match(cellid_order,cell_info$cellID),]$assay,colnames(qcObj$Y))
    sample_annot<-setNames(cell_info[match(cellid_order,cell_info$cellID),]$sample,colnames(qcObj$Y))
    
    gene_loc<-as.data.frame(findOverlaps(pdac_genes,bin_loci_subject))
    gene_names<-as.list(pdac_genes[order(gene_loc$subjectHits),]$gene)
    gene_loc<-as.list(gene_loc[order(gene_loc$subjectHits),]$subjectHits)
    gene_annot = columnAnnotation(PDAC_Genes = anno_mark(at = as.numeric(gene_loc), labels = gene_names))
    column_order<-1:nrow(iCN)
    library(circlize)
    col_fun = colorRamp2(c(0,1,2,3,4,5), c("#2166ac", "#67a9cf", "#f7f7f7","#fddbc7","#ef8a62","#b2182b"))

    
plt<-Heatmap(as.data.frame(t(iCN)),
    cluster_columns=F,
    name="copy_number",
    column_split=factor(chr_list,levels=unique(chr_list)),
    column_order=column_order,
    show_row_names=F,
    left_annotation=rowAnnotation(assay=assay_annot,sample=sample_annot),
    show_column_names=F,
    col=col_fun,
    bottom_annotation=gene_annot,
    #row_km=3,
    show_row_den=T,row_dend_side="right",use_raster=T
)

pdf("scope.gini.CopyNumberHeatmap.500kb.pdf",width=30)
plt
dev.off()
system("slack -F scope.gini.CopyNumberHeatmap.500kb.pdf ryan_todo")


    
#Plotting 500kb bins
    #gcc only
    iCN_gcc<-iCN[,startsWith(colnames(iCN),"s3gcc")]
    chr_list<-unlist(lapply(strsplit(row.names(iCN_gcc),":"),"[",1))

    qcObj$Y_gcc<-qcObj$Y[,startsWith(colnames(qcObj$Y),"s3gcc_")]
    cellid_order<-colnames(qcObj$Y_gcc)
    assay_annot<-setNames(cell_info[match(cellid_order,cell_info$cellID),]$assay,colnames(qcObj$Y_gcc))
    sample_annot<-setNames(cell_info[match(cellid_order,cell_info$cellID),]$sample,colnames(qcObj$Y_gcc))
    
    gene_loc<-as.data.frame(findOverlaps(pdac_genes,bin_loci_subject))
    gene_names<-as.list(pdac_genes[order(gene_loc$subjectHits),]$gene)
    gene_loc<-as.list(gene_loc[order(gene_loc$subjectHits),]$subjectHits)
    gene_annot = columnAnnotation(PDAC_Genes = anno_mark(at = as.numeric(gene_loc), labels = gene_names))
    column_order<-1:nrow(iCN_gcc)
    library(circlize)
    col_fun = colorRamp2(c(0,1,2,3,4,5), c("#2166ac", "#67a9cf", "#f7f7f7","#fddbc7","#ef8a62","#b2182b"))

    
plt<-Heatmap(as.data.frame(t(iCN_gcc)),
    cluster_columns=F,
    name="copy_number",
    column_split=factor(chr_list,levels=unique(chr_list)),
    column_order=column_order,
    show_row_names=F,
    left_annotation=rowAnnotation(assay=assay_annot,sample=sample_annot),
    show_column_names=F,
    col=col_fun,
    bottom_annotation=gene_annot,
    #row_km=3,
    show_row_den=T,row_dend_side="right",use_raster=T
)


pdf("scope.gini.CopyNumberHeatmap.gcconly.500kb.pdf",width=30)
plt
dev.off()
system("slack -F scope.gini.CopyNumberHeatmap.gcconly.500kb.pdf ryan_todo")

#Plot gene location CNV calls per cell line
library(ggplot2)
library(patchwork)

# Barplot change to percentage rather than pie chart
plot_piecharts<-function(x){
    row_idx<-unlist(gene_loc[x])
    gene_name<-gene_names[x]
    tmp<-as.data.frame(cbind("cnv"=iCN_gcc[row_idx,],"sample"=sample_annot))
    #tmp<-melt(table(tmp))
    tmp$cnv<-as.integer(tmp$cnv)
    plt_tmp<- ggplot(dat=tmp,aes(x=cnv,color=sample))+
        geom_density(adjust=1/10,alpha=0.2)+
        ggtitle(paste(gene_name,row.names(iCN_gcc)[row_idx])) + 
        scale_x_discrete(name ="Copy Number", 
                    limits=c(seq(0,6,1))) + xlim(c(0,6))
        theme_bw()
    return(plt_tmp)
    }
plt_list<-lapply(1:length(gene_loc),plot_piecharts)
plt<-wrap_plots(plt_list,ncol=2)
ggsave(plt,file="test.png",width=5,height=15)
system("slack -F test.png ryan_todo")

```


```python
## Single cell VCF Calling
```


```python
cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data
#VarScan for WES
module load varscan/2.3.7
#http://varscan.sourceforge.net/copy-number-calling.html
ref="/home/groups/oroaklab/refs/hg38/hg38.fa"

#Load and run GATK
module load GATK/3.2.2
ref="/home/groups/oroaklab/refs/hg38/hg38.fa"

for i in *RG.bam; do
iname=${i::-4} ;
samtools sort -m 5G -@ 30 -n $i > ${iname}.sorted.bam ;
samtools fixmate -r -@ 10 ${iname}.sorted.bam ${iname}.sorted.fixedmatededup.bam ;
samtools sort -@ 30 ${iname}.sorted.fixedmatededup.bam > ${iname}.sorted.fixedmatededup.sort.bam; &

for i in *fixedmatededup.sort.bam; do
samtools reheader -c 'grep -v ^@PG' $i >  ${i::-4}.reheader.bam;
samtools index -@ 10 ${i::-4}.reheader.bam ${i::-4}.reheader.bai ; done & 

#Load and run GATK
gatk="java -jar -Xmx2g gatk-package-4.1.9.0-local.jar"
ref="/home/groups/oroaklab/refs/hg38/hg38.fa"

for i in *sort.reheader.bam; do
gatk -T HaplotypeCaller \
-R $ref \
-I $i \
-bamout ${i::-4}.gatk_hc.bam > ${i::-4}.gatk_hc.vcf ; done &

```

### Coverage uniformity
Next generate MAPD and DIMAPD scores for all cells, using a custom script based on description from this website
from https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/interpret/metrics

MAPD, or the Median Absolute deviation of Pairwise Differences, is defined for a vector v[i] as follows:

    -construct the vector of consecutive pairwise differences scaled by the mean d[i] = ( v[i] - v[i+1] )/Mean(v),
    -compute the median absolute deviation of d[i] defined as Median(|d[i] - Median(d[i])|).
    -Performing genome coverage distribution on post GC corrected and mappability limited matrix (Y)



```python
R
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

Upon completion of normalization and segmentation at a first pass, SCOPE includes the option to cluster cells based on the matrix of normalized z-scores, estimated copy numbers, or estimated changepoints.
Given the inferred subclones, SCOPE can opt to perform a second round of group-wise ploidy initialization and normalization
I haven't doen this yet on the samples but is a good follow up once we have more sequencing.



```python
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

SCOPE provides the cross-sample segmentation, which outputs shared breakpoints across cells from the same clone. This step processes the entire genome chromosome by chromosome. Shared breakpoints and integer copy-number profiles will be returned.
Using circular binary segmentation (CBS) for breakpoint analysis


```python


normObj<-readRDS("SCOPE_normObj.noK.500kb.rds")
qcObj<-readRDS("qcObj_coverage_500kb.rds")
normObj.scope<-readRDS("SCOPE_normObj.scope.k1.500kb.rds")

#cluster by absolute copy number or corrected bins
Y <- normObj$Y
Yhat <- as.data.frame(normObj.scope$Yhat[[which.max(normObj.scope$BIC)]])
alpha_hat<-as.data.frame(normObj.scope$alpha.hat[[which.max(normObj.scope$BIC)]])
ref <- as.data.frame(qcObj$ref)

row.names(Yhat)<-paste(ref$seqnames,ref$start,ref$end,sep="_")
colnames(Yhat)<-colnames(Y)
row.names(alpha_hat)<-paste(ref$seqnames,ref$start,ref$end,sep="_")
colnames(alpha_hat)<-colnames(Y)

library(dbscan)
library(patchwork)
#Clustering cells by normalized read count per bin
Y_dims<-prcomp(alpha_hat,scale.=F) #PCA dim reduction
dims<-as.data.frame(uwot::umap(t(Y_dims)[1:10],n_components=2))#just umap the whole thing?
clus<-dbscan(t(alpha_hat),eps=1)
dims$clus<-clus$cluster
dims$cellID<-colnames(alpha_hat)
dims$cellID_idx<-unlist(lapply(strsplit(dims$cellID,"_"),"[",3))
dims<-merge(dims,annot,by="cellID_idx")
plt1<-ggplot(dat=dims,aes(x=V1,y=V2,color=as.factor(paste(assay,sample))))+geom_point()+theme_bw()+ggtitle("PC1 and 2, Experiment")
plt2<-ggplot(dat=dims,aes(x=V1,y=V2,color=as.factor(paste(clus))))+geom_point()+theme_bw()+ggtitle("PC1 and 2, cluster")
plt<-plt1+plt2
ggsave(plt,file="normalized_copynumber_PC.png",width=20)
system("slack -F normalized_copynumber_PC.png ryan_todo")

var_explained<-data.frame(summary(Y_dims)$importance) #variance explained per PC
var_explained<-var_explained[3,]
var_explained<-melt(var_explained)
var_explained$pc<-1:nrow(var_explained)

plt<-ggplot(dat=var_explained,aes(x=pc,y=value))+geom_point()+theme_bw()+ylim(c(0,1))+xlim(c(1,50))
ggsave(plt,file="genomic_segment_pca_varexplained.500kb.svg") #plot variance explained
system("slack -F genomic_segment_pca_varexplained.500kb.svg ryan_todo")


#using multiple PCs
#setting up a umap clustering function to be parallelized

Y_dims<-prcomp(alpha_hat,scale.=F) #PCA dim reduction


umap_clus<-function(z){
Y_pc<-as.data.frame(Y_dims$rotation)[,1:z]
clus<-dbscan(Y_pc,eps=0.05)
dims<-as.data.frame(uwot::umap(Y_pc,n_components=2))
row.names(dims)<-row.names(Y_pc)
colnames(dims)<-c("x","y")    
dims$clus<-clus$cluster
dims$cellID_idx<-unlist(lapply(strsplit(row.names(Y_pc),"_"),"[",3))
dims<-merge(dims,annot,by="cellID_idx")
plt1<-ggplot(dat=dims,aes(x=x,y=y,color=as.factor(paste(assay,sample))))+geom_point()+theme_bw()+ggtitle(paste("sample.500kb.PC.",str(z)))
plt2<-ggplot(dat=dims,aes(x=x,y=y,color=as.factor(paste(clus))))+geom_point()+theme_bw()+ggtitle(paste0("cluster.500kb.PC.",str(z)))
plt<-plt1+plt2
return(plt)
}

library(parallel)
plt_list<-mclapply(c(5,10,15,20,25,30),FUN=umap_clus,mc.cores=5)
plt_list<-wrap_plots(plt_list)
ggsave(plt_list,file="umap_muliplePC_clustering.pdf",height=20,width=30)
system("slack -F umap_muliplePC_clustering.pdf ryan_todo")

Y_dims<-as.data.frame(uwot::umap(t(Y)))
row.names(Y_dims)<-colnames(Y)
Y_dims$cellID<-row.names(Y_dims)

Y_dims<-merge(Y_dims,annot,by="cellID")
plt<-ggplot(dat=Y_dims,aes(x=V1,y=V2,color=as.factor(experiment)))+geom_point()+theme_bw()
ggsave(plt,file="genome_segment.umap.500kb.svg") #plot umap projection, honestly this isn't the most rigorous way to do this but I'm curious


###Perform group-wise normalization as alternative to for each normalization###

qcObj<-readRDS("qcObj_coverage_500kb.rds")
sampname <- qcObj$sampname
ref <- qcObj$ref
gc_qc<-ref$gc

normObj<-readRDS("SCOPE_normObj.noK.500kb.rds")
Y <- normObj$Y
beta.hat.noK <- normObj$beta.hat

#updating values with corrected amounts
normObj.scope<-readRDS("SCOPE_normObj.scope.k1.500kb.rds")
Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]

annot<-read.table("s3wgsgcc_cellsummary.500kb.tsv",header=T)
groups<-annot[match(annot$cellID,sampname),]$sample
norm_index<-which(annot[match(annot$cellID,sampname),]$gini<=0.3)

ploidy.group <- initialize_ploidy_group(Y = Y, Yhat = Yhat,
                                ref = ref, groups = groups)
saveRDS(ploidy.group,file="SCOPE_ploidygroup.rds")

# Group-wise normalization
normObj.scope.group <- normalize_scope_group(Y_qc = Y,
                                    gc_qc = gc,
                                    K = 1, ploidyInt = ploidy.group,
                                    norm_index = norm_index,
                                    groups = groups,
                                    T = 1:5,
                                    beta0 = beta.hat.noK)
saveRDS(normObj.scope.group,file="SCOPE_normalizationgroup.rds")

Yhat.group <- normObj.scope.group$Yhat[[which.max(
                                    normObj.scope.group$BIC)]]
fGC.hat.group <- normObj.scope.group$fGC.hat[[which.max(
                                    normObj.scope.group$BIC)]]

```


```python
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


## Now plotting those QC metrics out with ggplot and R
This is a big old copy and paste script and can probably be parsed down by user defined functions.



```python
R
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

## Using bcftools to generate vcf files.

Going to look at mutations in 4442 (Sample 1) and 4671 (Sample 2) mutations. Both have a KRAS driver mutation at chr12:25245350. Sample 1 is G12D (C>T transistion), Sample 2 is G12C.


```python

module load bcftools/1.1

#vcf version
#pileup
pdac_genes="/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/pdac_gene.bed"
ref="/home/groups/oroaklab/refs/hg38/hg38.fa"

#Add RG to bam file to output on sc level
#To add RG (cell IDs) to a bam file

#generate list of cellIDs from input bam as temporary annotation
#then make a RG bam from input
N=30 #N cores for parallel sort

for bam in hg38.s3WGS_GM12878.bbrd.q10.filt.bam hg38.s3GCC_4442.bbrd.q10.filt.bam hg38.s3GCC_4671.bbrd.q10.filt.bam hg38.s3WGS_4442.bbrd.q10.filt.bam hg38.s3WGS_4671.bbrd.q10.filt.bam;
do samtools view $bam | awk 'OFS="\t" {split($1,a,":"); print a[1],"temp"}' | sort -k1,1 -T . --parallel=$N -S 2G | uniq > cellID_list.tmp.annot;
{(samtools view -H $bam) ; 
(awk 'OFS="\t" {print "@RG","ID:"$1,"SM:"$1,"LB:"$1,"PL:SCI"}' cellID_list.tmp.annot) ; 
(echo -e "@PG\tID:scitools_bam-addrg\tVN:dev");
(samtools view $bam | awk 'OFS="\t" {split($1,a,":"); print $0,"RG:Z:"a[1]}');
 } | samtools view -bS - > ${bam%.bam}.RG.bam;
 done &

#Index output bam file
for bam in hg38.s3WGS_GM12878.bbrd.q10.filt.bam hg38.s3GCC_4442.bbrd.q10.filt.bam hg38.s3GCC_4671.bbrd.q10.filt.bam hg38.s3WGS_4442.bbrd.q10.filt.bam hg38.s3WGS_4671.bbrd.q10.filt.bam; do
samtools index -b -@ 10 ${bam%.bam}.RG.bam ; done &


#Generate pileup and vcf
for bam in hg38.s3WGS_GM12878.bbrd.q10.filt.bam hg38.s3GCC_4442.bbrd.q10.filt.bam hg38.s3GCC_4671.bbrd.q10.filt.bam hg38.s3WGS_4442.bbrd.q10.filt.bam hg38.s3WGS_4671.bbrd.q10.filt.bam; do
/home/groups/oroaklab/src/bcftools/bcftools/bcftools mpileup -f $ref -R $pdac_genes --threads 10 -O z ${bam%.bam}.RG.bam > ${bam%.bam}.vcf.gz ; done

#sort
/home/groups/oroaklab/src/bcftools/bcftools/bcftools sort -T . -O z -o test_4442.pdac.mpileup.sort.vcf test_4442.pdac.mpileup.vcf

#index for visualiztion in IGV
~/tools/IGVTools/igvtools index test_4442.pdac.mpileup.sort.vcf 

#Look at by hand
awk '$1=="chr12" && $2=25235382 {print $0}' ${bam%.bam}.vcf | less-S

```
