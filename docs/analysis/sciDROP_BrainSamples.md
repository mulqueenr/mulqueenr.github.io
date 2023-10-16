---
title: txci-ATAC
layout: analysis
permalink: /scidrop/
category: alternative
---


## Processing for txci-ATAC Project

![sciDROP Overview](/assets/images/sciDROP.png){:width="80%"}


Note: sciDROP was renamed to txci-ATAC for final manuscript submission. 


This notebook details the processing of the "20K" and "70K" loaded mouse brain and human cortex samples. It begins with scitools wrapper functions for intial alignment to a concatenated mouse and human genome, following with splitting of reads and realignment to separate human and mouse genomes. It then follows the established scitools formation of a counts matrix and Signac processing.
https://github.com/adeylab/scitools

```bash
#libraries were generated as two separate lanes of a NovaSeq S4 flowcell.
#bcl2fastq was run prior to the transfer
FASTQ_DIR="/home/groups/oroaklab/fastq/201103_NovaSeq_sciDropATAC"
#fastq files downloaded to $FASTQ_DIR

#Using a perl written scitools function for demultiplexing reads
cd $FASTQ_DIR
scitools fastq-dump-10x -V \
-1 70K_S4_L004_R1_001.fastq.gz \
-2 70K_S4_L004_R2_001.fastq.gz \
-i 70K_S4_L004_I1_001.fastq.gz \
-j 70K_S4_L004_I2_001.fastq.gz \
-o sciDROP_70k -R $FASTQ_DIR &

scitools fastq-dump-10x -V \
-1 20K_S3_L003_R1_001.fastq.gz \
-2 20K_S3_L003_R2_001.fastq.gz \
-i 20K_S3_L003_I1_001.fastq.gz \
-j 20K_S3_L003_I2_001.fastq.gz \
-o sciDROP_20K -R $FASTQ_DIR &

#Also including the 10% of 20k library loading which was done in house
FASTQ_DIR="/home/groups/oroaklab/fastq/201007_NS500556_0428_AHGFMMAFX2"
cd $FASTQ_DIR
 #no -V option since seq chem is different on nextseq
scitools fastq-dump-10x \
-1 Undetermined_S0_R1_001.fastq.gz \
-2 Undetermined_S0_R2_001.fastq.gz \
-i Undetermined_S0_I1_001.fastq.gz \
-j Undetermined_S0_I2_001.fastq.gz \
-o sciDROP_20k_10perc -R $FASTQ_DIR &

#Files automatically output to a directory with the prefix as the directory name
sciDROP_20K_demux="/home/groups/oroaklab/demultiplex/sciDROP_20K"
sciDROP_70k_demux="/home/groups/oroaklab/demultiplex/sciDROP_70k"
sciDROP_20k_10perc_demux="/home/groups/oroaklab/demultiplex/sciDROP_20k_10perc"

#Aligning fastq reads to a concatenated human and mouse genome with bwa-mem scitools wrapper
#Using a setting in scitools to sort by cellID (-n) this will decrease the memory footprint for duplicate removal later
scitools fastq-align -m 5G -n -t 20 -r 20 hg38mm10 \
sciDROP_20k \
$sciDROP_20K_demux/sciDROP_20K.1.fq.gz \
$sciDROP_20K_demux/sciDROP_20K.2.fq.gz &

scitools fastq-align -m 5G -n -t 20 -r 20 hg38mm10 \
sciDROP_70k \
$sciDROP_70k_demux/sciDROP_70k.1.fq.gz \
$sciDROP_70k_demux/sciDROP_70k.2.fq.gz &

scitools fastq-align -m 5G -n -t 20 -r 20 hg38mm10 \
sciDROP_20k_10perc \
$sciDROP_20k_10perc_demux/sciDROP_20k_10perc.1.fq.gz \
$sciDROP_20k_10perc_demux/sciDROP_20k_10perc.2.fq.gz &

#Remove duplicates based on cellID, chromosome and start sites per read
#using the name sorted (-n) barcode based removal of duplicates
scitools bam-rmdup -n -t 12 sciDROP_70k.nsrt.bam 
scitools bam-rmdup -n -t 12 sciDROP_20k.nsrt.bam 
scitools bam-rmdup -n -t 12 sciDROP_20k_10perc.nsrt.bam
#combine 20k_10perc with novaseq data

#Run a barnyard comparison for the duplicate removal bams
#Counts read alignments based on chromosomes and plots
scitools barnyard-compare sciDROP_70k.bbrd.q10.bam
scitools barnyard-compare sciDROP_20k.bbrd.q10.bam

#Move it all to a new project directory
wd="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard"
tree $wd
.
├── sciDROP_20K
│   ├── sciDROP_20K.1.fq.gz
│   ├── sciDROP_20K.2.fq.gz
│   ├── sciDROP_20k.align.log
│   ├── sciDROP_20k.bbrd.q10.bam
│   ├── sciDROP_20k.bbrd.q10.barnyard_cells.plot.pdf
│   ├── sciDROP_20k.bbrd.q10.barnyard_cells.plot.png
│   ├── sciDROP_20k.bbrd.q10.barnyard_cells.plot.r
│   ├── sciDROP_20k.bbrd.q10.barnyard_cells.txt
│   ├── sciDROP_20k.bbrd.q10.barnyard_stats.txt
│   ├── sciDROP_20k.complexity.hist.pdf
│   ├── sciDROP_20k.complexity.hist.png
│   ├── sciDROP_20k.complexity.log
│   ├── sciDROP_20k.complexity.pdf
│   ├── sciDROP_20k.complexity.plot.r
│   ├── sciDROP_20k.complexity.plot.txt
│   ├── sciDROP_20k.complexity.png
│   ├── sciDROP_20k.complexity.txt
│   ├── sciDROP_20K.fail.1.fq.gz
│   ├── sciDROP_20K.fail.2.fq.gz
│   ├── sciDROP_20k.nsrt.bam
└── sciDROP_70k
    ├── sciDROP_70k.1.fq.gz
    ├── sciDROP_70k.2.fq.gz
    ├── sciDROP_70k.align.log
    ├── sciDROP_70k.bbrd.q10.bam
    ├── sciDROP_70k.bbrd.q10.barnyard_cells.plot.pdf
    ├── sciDROP_70k.bbrd.q10.barnyard_cells.plot.png
    ├── sciDROP_70k.bbrd.q10.barnyard_cells.plot.r
    ├── sciDROP_70k.bbrd.q10.barnyard_cells.txt
    ├── sciDROP_70k.bbrd.q10.barnyard_stats.txt
    ├── sciDROP_70k.complexity.hist.pdf
    ├── sciDROP_70k.complexity.hist.png
    ├── sciDROP_70k.complexity.log
    ├── sciDROP_70k.complexity.pdf
    ├── sciDROP_70k.complexity.plot.r
    ├── sciDROP_70k.complexity.plot.txt
    ├── sciDROP_70k.complexity.png
    ├── sciDROP_70k.complexity.txt
    ├── sciDROP_70k.fail.1.fq.gz
    ├── sciDROP_70k.fail.2.fq.gz
    ├── sciDROP_70k.nsrt.bam

2 directories, 76 files

```

## Combine 10% Sampling of 20K Experiment with current run data


```bash
dir_20k_10perc="/home/groups/oroaklab/adey_lab/projects/sciDROP/201007_BrainBarnyard_Test"
sciDROP_20K_demux="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_20K"

#Complexity file
cat $sciDROP_20K_demux/sciDROP_20k.complexity.txt \
$dir_20k_10perc/scidrop_barnyard.complexity.txt > $sciDROP_20K_demux/sciDROP_20k.complexity.merge.txt

#Barnyard cells file
cat $sciDROP_20K_demux/sciDROP_20k.bbrd.q10.barnyard_cells.txt \
$dir_20k_10perc/scidrop_barnyard.bbrd.q10.barnyard_cells.txt > $sciDROP_20K_demux/sciDROP_20k.bbrd.q10.barnyard_cells.merge.txt

#Fastq files
dir_20k_10perc_demux="/home/groups/oroaklab/demultiplex/201007_NS500556_0428_AHGFMMAFX2"
cat $sciDROP_20K_demux/sciDROP_20K.1.fq.gz \
$dir_20k_10perc_demux/201007_NS500556_0428_AHGFMMAFX2.1.fq.gz > $sciDROP_20K_demux/sciDROP_20K.1.merge.fq.gz &

cat $sciDROP_20K_demux/sciDROP_20K.2.fq.gz \
$dir_20k_10perc_demux/201007_NS500556_0428_AHGFMMAFX2.2.fq.gz > $sciDROP_20K_demux/sciDROP_20K.2.merge.fq.gz &


```

## Calculate collision rate from barnyard experiment

```R
#Processing of barnyard comparisons
library(ggplot2)
library(Biostrings)
#Read in index file to assign well position to indexes
index_file<-read.table("/home/groups/oroaklab/src/scitools/scitools-dev/SCI_stdchem_Indexes.txt")

idx_pcr<-index_file[index_file$V2==1,]
colnames(idx_pcr)<-c("i7_idx_name","i7_idx_cycle","i7_idx_sequence")
idx_tn5<-index_file[index_file$V2==3,]
colnames(idx_tn5)<-c("tn5_idx_name","tn5_idx_cycle","tn5_idx_sequence")
idx_tn5$tn5_column<-c(1:12)
idx_tn5$tn5_row<-rep(c("A","B","C","D","E","F","G","H"),c(rep(12,8)))

#Processing of 70k samples
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_70k")
dat<-read.table("sciDROP_70k.bbrd.q10.barnyard_cells.txt",header=F)
colnames(dat)<-c("cellID","total_reads_q20","hg38_count","mm10_count","percent_human","species_call")

compl<-read.table("sciDROP_70k.complexity.txt",header=F)
colnames(compl)<-c("row_carryover","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
compl<-compl[compl$uniq_reads_q10>=10000 & compl$perc_uniq <= 90,] #filter cells based on minimum read count and library complexity

dat<-merge(dat,compl,by="cellID")
dat$pcr_idx<-substr(dat$cellID,1,8)
dat$gem_idx<-substr(dat$cellID,9,24)
dat$tn5_idx<-substr(dat$cellID,25,32)

dat<-merge(dat,idx_pcr,by.x="pcr_idx",by.y="i7_idx_sequence") #add pcr i7 index, defining 10% or 90% pool

table(dat$i7_idx_name)
#PCR_i7_P7.S707 PCR_i7_P7.S708
#  6460          54928

dat<-merge(dat,idx_tn5,by.x="tn5_idx",by.y="tn5_idx_sequence")

library(dplyr)
summary(as.data.frame(dat %>% group_by(tn5_column,tn5_row) %>% summarize(count=n()))$count) #cell count distribution per tn5 well
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  104.0   512.0   615.0   646.2   767.0  1171.0

nrow(dat)/length(unique(dat$gem_idx)) #count of cells within unique GEMs
#[1] 1.665256
gem_count<-unlist(lapply(unique(dat$gem_idx),function(x) nrow(dat[dat$gem_idx==x,])))
    plt<-ggplot()+geom_histogram(aes(x=gem_count),binwidth=1)+theme_bw()+xlim(c(0,10))
    ggsave(plt,file="70k.gem_count.pdf")
    #system("slack -F 70k.gem_count.pdf ryan_todo")

dat$condition<-"Mix"
dat[dat$tn5_row %in% c("A","B"),]$condition<-"Human"
dat[dat$tn5_row %in% c("C","D"),]$condition<-"Mouse"

ggplot(dat,aes(x=perc_uniq,y=log10(uniq_reads_q10),color=species_call))+geom_point()+theme_bw()+ylim(c(0,7))+xlim(c(0,100))
ggsave("complexity.pdf")
#system("slack -F complexity.pdf ryan_todo")

library(patchwork)
plt_barnyard<-ggplot(dat[dat$condition=="Mix",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Barnyard")
plt_human<-ggplot(dat[dat$condition=="Human",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Human")
plt_mouse<-ggplot(dat[dat$condition=="Mouse",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Mouse")
plt<-plt_barnyard+plt_human+plt_mouse
ggsave(plt,file="barnyard.pdf",width=20)
#system("slack -F barnyard.pdf ryan_todo")

#count of barnyard species calls
table(dat[dat$condition=="Mix",]$species_call)
#Human Mixed Mouse
#10040   170 18985
(table(dat[dat$condition=="Mix",]$species_call)[["Mixed"]]/nrow(dat))*2 #get estimated collision rate
#[1] 0.005538542 so 0.55%
write.table(dat,"70k_barnyard.summary.txt",quote=F,sep="\t",col.names=T,row.names=F)


#annot species
dat<-dat[(dat$condition=="Mouse" & dat$species_call=="Mouse") | (dat$condition=="Human" & dat$species_call=="Human") | dat$condition=="Mix", ]
table(paste(dat$condition,dat$species_call))
#Human Human   Mix Human   Mix Mixed   Mix Mouse Mouse Mouse
#12525       10040         170       18985       19526
annot<-dat[c("cellID","species_call")]
write.table(annot,"species.annot",quote=F,sep="\t",col.names=F,row.names=F)

#Processing of 20k samples
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_20K")
dat<-read.table("sciDROP_20k.bbrd.q10.barnyard_cells.merge.txt",header=F)
colnames(dat)<-c("cellID","total_reads_q20","hg38_count","mm10_count","percent_human","species_call")

compl<-read.table("sciDROP_20k.complexity.merge.txt",header=F)
colnames(compl)<-c("row_carryover","cellID","total_reads_q10","uniq_reads_q10","perc_uniq")
compl<-compl[compl$uniq_reads_q10>=10000 & compl$perc_uniq <= 90,]

dat<-merge(dat,compl,by="cellID")
dat$pcr_idx<-substr(dat$cellID,1,8)
dat$gem_idx<-substr(dat$cellID,9,24)
dat$tn5_idx<-substr(dat$cellID,25,32)

dat<-merge(dat,idx_pcr,by.x="pcr_idx",by.y="i7_idx_sequence") #add pcr i7 index, defining 10% or 90% pool

table(dat$i7_idx_name)
#PCR_i7_P7.S701 PCR_i7_P7.S702
#  1848          17293

dat<-merge(dat,idx_tn5,by.x="tn5_idx",by.y="tn5_idx_sequence")

library(dplyr)
summary(as.data.frame(dat %>% group_by(tn5_column,tn5_row) %>% summarize(count=n()))$count) #cell count distribution per tn5 well
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#63.0   161.0   195.0   201.5   246.5   331.0   
nrow(dat)/length(unique(dat$gem_idx)) #count of cells within unique GEMs
#1.257456
gem_count<-unlist(lapply(unique(dat$gem_idx),function(x) nrow(dat[dat$gem_idx==x,])))
plt<-ggplot()+geom_histogram(aes(x=gem_count),binwidth=1)+theme_bw()+xlim(c(0,10))
ggsave(plt,file="20k.gem_count.pdf")
#system("slack -F 20k.gem_count.pdf ryan_todo")

dat$condition<-"Mix"
dat[dat$tn5_row %in% c("A","B"),]$condition<-"Human"
dat[dat$tn5_row %in% c("C","D"),]$condition<-"Mouse"

ggplot(dat,aes(x=perc_uniq,y=log10(uniq_reads_q10),color=species_call))+geom_point()+theme_bw()+ylim(c(0,7))+xlim(c(0,100))
ggsave("complexity.pdf")
#system("slack -F complexity.pdf ryan_todo")

library(patchwork)
plt_barnyard<-ggplot(dat[dat$condition=="Mix",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Barnyard")
plt_human<-ggplot(dat[dat$condition=="Human",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Human")
plt_mouse<-ggplot(dat[dat$condition=="Mouse",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Mouse")
plt<-plt_barnyard+plt_human+plt_mouse
ggsave(plt,file="barnyard.pdf",width=20)
#system("slack -F barnyard.pdf ryan_todo")

write.table(dat,"20k_barnyard.summary.txt",quote=F,sep="\t",col.names=T,row.names=F)

#count of barnyard species calls
table(dat[dat$condition=="Mix",]$species_call)
#Human Mixed Mouse
#10040   170 18985
(table(dat[dat$condition=="Mix",]$species_call)[["Mixed"]]/nrow(dat))*2 #get estimated collision rate
#[1] 0.002660036 so 0.2% collision rate

#Filter to data that matches apriori assumption or single cell
#annot species
dat<-dat[(dat$condition=="Mouse" & dat$species_call=="Mouse") | (dat$condition=="Human" & dat$species_call=="Human") | dat$condition=="Mix", ]

#Human Human   Mix Human   Mix Mixed   Mix Mouse Mouse Mouse
# 3687        3001          24        6572        5841
annot<-dat[c("cellID","species_call")]
write.table(annot,"species.annot",quote=F,sep="\t",col.names=F,row.names=F)


```

## Split out species from barnyard experiments


```bash
sciDROP_20k_dir="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_20K"
sciDROP_70k_dir="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_70k"

#Split fastq files based on barnyard analysis
scitools split-fastq -A $sciDROP_20k_dir/species.annot $sciDROP_20k_dir/sciDROP_20K.1.merge.fq.gz $sciDROP_20k_dir/sciDROP_20K.2.merge.fq.gz &
scitools split-fastq -A $sciDROP_70k_dir/species.annot $sciDROP_70k_dir/sciDROP_70k.1.fq.gz $sciDROP_70k_dir/sciDROP_70k.2.fq.gz &

#Realign fastq files to proper genome
scitools fastq-align -t 20 -r 20 -n hg38 hg38 $sciDROP_20k_dir/species.Human.1.fq.gz $sciDROP_20k_dir/species.Human.2.fq.gz &
scitools fastq-align -t 20 -r 20 -n mm10 mm10 $sciDROP_20k_dir/species.Mouse.1.fq.gz $sciDROP_20k_dir/species.Mouse.2.fq.gz &

scitools fastq-align -t 20 -r 20 -n hg38 hg38 $sciDROP_70k_dir/species.Human.1.fq.gz $sciDROP_70k_dir/species.Human.2.fq.gz &
scitools fastq-align -t 20 -r 20 -n mm10 mm10 $sciDROP_70k_dir/species.Mouse.1.fq.gz $sciDROP_70k_dir/species.Mouse.2.fq.gz &

#Remove duplicates
scitools bam-rmdup -n $sciDROP_20k_dir/mm10.bam &
scitools bam-rmdup -n $sciDROP_20k_dir/hg38.bam &
scitools bam-rmdup -n $sciDROP_70k_dir/mm10.bam &
scitools bam-rmdup -n $sciDROP_70k_dir/hg38.bam &

#move merged bam files up a directory
out_dir="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard"
samtools merge -f -@ 20 -h $sciDROP_70k_dir/hg38.bbrd.q10.bam $out_dir/hg38.merged.bbrd.q10.bam $sciDROP_70k_dir/hg38.bbrd.q10.bam $sciDROP_20k_dir/hg38.bbrd.q10.bam &  
samtools merge -f -@ 20 -h $sciDROP_70k_dir/mm10.bbrd.q10.bam $out_dir/mm10.merged.bbrd.q10.bam $sciDROP_70k_dir/mm10.bbrd.q10.bam $sciDROP_20k_dir/mm10.bbrd.q10.bam &  

#Call peaks by read pileups
cd $out_dir
scitools callpeaks hg38.merged.bbrd.q10.bam &
scitools callpeaks mm10.merged.bbrd.q10.bam &

wc -l hg38.merged.bbrd.q10.500.bed
#332864 hg38.merged.bbrd.q10.500.bed

wc -l mm10.merged.bbrd.q10.500.bed
#256912 mm10.merged.bbrd.q10.500.bed


scitools atac-counts -O hg38 hg38.merged.bbrd.q10.bam \
hg38.merged.bbrd.q10.500.bed &
scitools atac-counts -O mm10 mm10.merged.bbrd.q10.bam \
mm10.merged.bbrd.q10.500.bed &

#scitools wrapper for samtools isize
scitools isize hg38.merged.bbrd.q10.bam &
scitools isize mm10.merged.bbrd.q10.bam &

#scitools wrapper for tss enrichment
scitools bam-tssenrich mm10.merged.bbrd.q10.bam mm10 &
scitools bam-tssenrich hg38.merged.bbrd.q10.bam hg38 &

```

### Tabix fragment file generation

Tabix file format is a tab separated multicolumn data structure.

| Column Number | Name | Description |
|:--------|:-------:|:--------|
|1 |chrom |  Reference genome chromosome of fragment |
|2 |chromStart | Adjusted start position of fragment on chromosome. |
|3 |chromEnd   | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
|4 |barcode | The 10x (or sci) cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment. |
|5 |duplicateCount |The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeq instrument are excluded from this count. |

```bash
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
#human processing
input_bam="hg38.merged.bbrd.q10.bam"
output_name=${input_bam::-13}
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=10 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz; wait ;
$tabix -p bed $output_name.fragments.tsv.gz &
#mouse processing
input_bam="mm10.merged.bbrd.q10.bam"
output_name=${input_bam::-13}
samtools view --threads 20 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=20 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
$tabix -p bed $output_name.fragments.tsv.gz &
```

# sciDROP Full Processing

### Generating Seurat Objects

Using R v4.0.0 and Signac v1.0

```R
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(Matrix)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

#function to read in sparse matrix format from atac-count
read_in_sparse<-function(x){ #x is character file prefix followed by .bbrd.q10.500.counts.sparseMatrix.values.gz
IN<-as.matrix(read.table(paste0(x,".counts.sparseMatrix.values.gz")))
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table(paste0(x,".counts.sparseMatrix.cols.gz"))
colnames(IN)<-COLS$V1
ROWS<-read.table(paste0(x,".counts.sparseMatrix.rows.gz"))
row.names(IN)<-ROWS$V1
writeMM(IN,file=paste0(x,".counts.mtx")) #this is to generate counts matrices in scrublet friendly format
return(IN)
}

hg38_counts<-read_in_sparse("hg38") # make hg38 counts matrix from sparse matrix
mm10_counts<-read_in_sparse("mm10") # make mm10 counts matrix from sparse matrix

#write out as MM format
#Read in fragment path for coverage plots
mm10_fragment.path="./mm10.merged.fragments.tsv.gz"
hg38_fragment.path="./hg38.merged.fragments.tsv.gz"

# extract gene annotations from EnsDb
hg38_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
mm10_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style
seqlevelsStyle(hg38_annotations) <- 'UCSC'
seqlevelsStyle(mm10_annotations) <- 'UCSC'

genome(hg38_annotations) <- "hg38"
genome(mm10_annotations) <- "mm10"

#Generate ChromatinAssay Objects
hg38_chromatinassay <- CreateChromatinAssay(
  counts = hg38_counts,
  genome="hg38",
  min.cells = 1,
  annotation=hg38_annotations,
  sep=c("_","_"),
  fragments=hg38_fragment.path
)

mm10_chromatinassay <- CreateChromatinAssay(
  counts = mm10_counts,
  genome="mm10",
  min.cells = 1,
  annotation=mm10_annotations,
  sep=c("_","_"),
  fragments=mm10_fragment.path
)


#Create Seurat Objects
hg38_atac <- CreateSeuratObject(
  counts = hg38_chromatinassay,
  assay = "peaks"
)
mm10_atac <- CreateSeuratObject(
  counts = mm10_chromatinassay,
  assay = "peaks"
)

#Meta.data to be updated after clustering

#saving unprocessed SeuratObjects
saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")
saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")
```

### Perform Scrublet on Data to Ensure Single-cells

Code from tutorial here.[https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb]

```python
#using a conda environment set up by ARSN
#source /home/groups/oroaklab/nishida/scitools_env/bin/activate
#Installing scrublet
#pip install scrublet
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.sparse import coo_matrix
import gzip
import pandas as pd

#Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.

input_dir = '/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/'

#Perform scrublet on mm10 cells
counts_matrix = scipy.io.mmread(input_dir + 'mm10.counts.mtx').T.tocsc() #generated during the initialization of the Seurat Object

peaks= np.array(gzip.open(input_dir+'mm10.counts.sparseMatrix.rows.gz', 'rt').read().split()) #This is read in to check that our data frame is in the correct orientation
cellid= gzip.open(input_dir+'mm10.counts.sparseMatrix.cols.gz', 'rt').read().split() #This is read in to check that our data frame is in the correct orientation
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(peaks)))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
#Run the default pipeline, which includes:
#Doublet simulation
#Normalization, gene filtering, rescaling, PCA
#Doublet score calculation
#Doublet score threshold detection and doublet calling
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

df = pd.DataFrame({'cellid':cellid, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv('mm10.scrublet.tsv', index=False, sep="\t")

#Perform on hg38 cells
counts_matrix = scipy.io.mmread(input_dir + 'hg38.counts.mtx').T.tocsc() #generated during the initialization of the Seurat Object

peaks= np.array(gzip.open(input_dir+'hg38.counts.sparseMatrix.rows.gz', 'rt').read().split()) #This is read in to check that our data frame is in the correct orientation
cellid= gzip.open(input_dir+'hg38.counts.sparseMatrix.cols.gz', 'rt').read().split() #This is read in to check that our data frame is in the correct orientation
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(peaks)))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
#Run the default pipeline, which includes:
#Doublet simulation
#Normalization, gene filtering, rescaling, PCA
#Doublet score calculation
#Doublet score threshold detection and doublet calling
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

df = pd.DataFrame({'cellid':cellid, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv('hg38.scrublet.tsv', index=False, sep="\t")
```

### Add library complexity data to RDS files.


```bash
 cat ./sciDROP_20K/hg38.complexity.txt ./sciDROP_20K/mm10.complexity.txt ./sciDROP_70k/hg38.complexity.txt ./sciDROP_70k/mm10.complexity.txt > complexity.txt
```

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
hg38_atac$cellid<-row.names(hg38_atac@meta.data)
mm10_atac$cellid<-row.names(mm10_atac@meta.data)

hg38_doublet<-read.table("hg38.scrublet.tsv",head=T)
hg38_atac<-AddMetaData(object=hg38_atac,col="doublet_scores",metadata=setNames(hg38_doublet$doublet_scores,hg38_doublet$cellid))
hg38_atac<-AddMetaData(object=hg38_atac,col="predicted_doublets",metadata=setNames(hg38_doublet$predicted_doublets,hg38_doublet$cellid))

mm10_doublet<-read.table("mm10.scrublet.tsv",head=T)
mm10_atac<-AddMetaData(object=mm10_atac,col="doublet_scores",metadata=setNames(mm10_doublet$doublet_scores,mm10_doublet$cellid))
mm10_atac<-AddMetaData(object=mm10_atac,col="predicted_doublets",metadata=setNames(mm10_doublet$predicted_doublets,mm10_doublet$cellid))

compl<-read.table(file="complexity.txt",header=F)
colnames(compl)<-c("row_carryover","cellid","total_reads","unique_reads","percent_uniq")
row.names(compl)<-compl$cellid

hg38_atac$cellid<-row.names(hg38_atac@meta.data)
hg38_atac <- AddMetaData(object = hg38_atac, col="total_reads",metadata = setNames(compl$total_reads,row.names(compl)))
hg38_atac <- AddMetaData(object = hg38_atac, col="unique_reads",metadata = setNames(compl$unique_reads,row.names(compl)))
hg38_atac <- AddMetaData(object = hg38_atac, col="percent_uniq",metadata = setNames(compl$percent_uniq,row.names(compl)))
hg38_atac$pcr_idx<-substr(hg38_atac$cellid,1,8)
hg38_atac$gem_idx<-substr(hg38_atac$cellid,9,24)
hg38_atac$tn5_idx<-substr(hg38_atac$cellid,25,32)

mm10_atac <- AddMetaData(object = mm10_atac, col="total_reads",metadata = setNames(compl$total_reads,row.names(compl)))
mm10_atac <- AddMetaData(object = mm10_atac, col="unique_reads",metadata = setNames(compl$unique_reads,row.names(compl)))
mm10_atac <- AddMetaData(object = mm10_atac, col="percent_uniq",metadata = setNames(compl$percent_uniq,row.names(compl)))
mm10_atac$pcr_idx<-substr(mm10_atac$cellid,1,8)
mm10_atac$gem_idx<-substr(mm10_atac$cellid,9,24)
mm10_atac$tn5_idx<-substr(mm10_atac$cellid,25,32)


saveRDS(mm10_atac,"mm10_SeuratObject.Rds")
saveRDS(hg38_atac,"hg38_SeuratObject.Rds")

plt1<-ggplot(mm10_atac@meta.data[mm10_atac@meta.data$predicted_doublets=="False",], aes(x=as.numeric(percent_uniq), y=log10(as.numeric(unique_reads))) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,100)) +
  scale_y_continuous(expand = c(0, 0),limits=c(3,6)) +
  theme(legend.position='none')
plt2<-ggplot(hg38_atac@meta.data[hg38_atac@meta.data$predicted_doublets=="False",], aes(x=as.numeric(percent_uniq), y=log10(as.numeric(unique_reads))) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,100)) +
  scale_y_continuous(expand = c(0, 0),limits=c(3,6)) +
  theme(legend.position='none')
ggsave(plt1+plt2,file="mm10.hg38.complexity.2d.pdf")
#system("slack -F mm10.hg38.complexity.2d.pdf ryan_todo")

#hard coded these numbers just because i had them for a meeting
#cell_count_75k<-data.frame(count=c(10040,18985,170,12663-(141+19),19530-(484+50),141+19,484+50),names=c("by_h","by_m","by_mix","hum","mus","scrub_h","scrub_m"),loading=c("75k"))
#cell_count_20k<-data.frame(count=c(3001,6572,24,3703-(27),5841-(449),27,449),names=c("by_h","by_m","by_mix","hum","mus","scrub_h","scrub_m"),loading=c("20k"))
#cell_count<-rbind(cell_count_75k,cell_count_20k)
#plt<-ggplot(cell_count,aes(x=loading,y=count,fill=factor(names,levels=rev(c("mus","hum","by_h","by_m","by_mix","scrub_h","scrub_m")))))+geom_bar(position="stack",stat="identity")
#ggsave(plt,file="mm10.hg38.cellcount.pdf")
##system("slack -F mm10.hg38.cellcount.pdf ryan_todo")

```

## Running Initial Clustering of Cells
Using CisTopic for Dimensionality reduction and UMAP for projection.

```R
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(harmony,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/lib_backup_210125")

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")

hg38_cistopic_counts_frmt<-hg38_atac$peaks@counts

#renaming row names to fit granges expectation of format
row.names(hg38_cistopic_counts_frmt)<-sub("-", ":", row.names(hg38_cistopic_counts_frmt))

#set up CisTopicObjects
hg38_atac_cistopic<-cisTopic::createcisTopicObject(hg38_cistopic_counts_frmt)

#Run warp LDA on objects
hg38_atac_cistopic_models<-cisTopic::runWarpLDAModels(hg38_atac_cistopic,topic=c(22,24,26,28,30),nCores=5,addModels=FALSE)

#Saving all models for posterity
saveRDS(hg38_atac_cistopic_models,file="hg38_CisTopicObject.Rds")

#Performing same operation for mm10
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
mm10_cistopic_counts_frmt<-mm10_atac$peaks@counts
row.names(mm10_cistopic_counts_frmt)<-sub("-", ":", row.names(mm10_cistopic_counts_frmt))
mm10_atac_cistopic<-cisTopic::createcisTopicObject(mm10_cistopic_counts_frmt)
mm10_atac_cistopic_models<-cisTopic::runWarpLDAModels(mm10_atac_cistopic,topic=c(24,26,28,30),nCores=4,addModels=FALSE)
saveRDS(mm10_atac_cistopic_models,file="mm10_CisTopicObject.Rds")

mm10_atac_cistopic_models<-readRDS(file="mm10_CisTopicObject.Rds")
hg38_atac_cistopic_models<-readRDS(file="hg38_CisTopicObject.Rds")

#Setting up topic count selection
pdf("hg38_atac_model_selection.pdf")
par(mfrow=c(3,3))
hg38_atac_cistopic_models <- selectModel(hg38_atac_cistopic_models, type='derivative')
dev.off()

pdf("mm10_atac_model_selection.pdf")
par(mfrow=c(3,3))
mm10_atac_cistopic_models <- selectModel(mm10_atac_cistopic_models, type='derivative')
dev.off()

#system("slack -F hg38_atac_model_selection.pdf ryan_todo")
#system("slack -F mm10_atac_model_selection.pdf ryan_todo")

#set topics based on derivative
#selected topics subject to change
mm10_selected_topic=28
hg38_selected_topic=30
mm10_cisTopicObject<-cisTopic::selectModel(mm10_atac_cistopic_models,select=mm10_selected_topic,keepModels=F)
hg38_cisTopicObject<-cisTopic::selectModel(hg38_atac_cistopic_models,select=hg38_selected_topic,keepModels=F)

#saving model selected RDS
saveRDS(hg38_cisTopicObject,file="hg38_CisTopicObject.Rds")
saveRDS(mm10_cisTopicObject,file="mm10_CisTopicObject.Rds")

#Read in cisTopic objects
hg38_cisTopicObject<-readRDS("hg38_CisTopicObject.Rds")
mm10_cisTopicObject<-readRDS("mm10_CisTopicObject.Rds")
#read in seurat format object
hg38_atac<-readRDS("hg38_SeuratObject.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.Rds")



#run UMAP on topics
hg38_topic_df<-as.data.frame(hg38_cisTopicObject@selected.model$document_expects)
row.names(hg38_topic_df)<-paste0("Topic_",row.names(hg38_topic_df))
hg38_dims<-as.data.frame(uwot::umap(t(hg38_topic_df),n_components=2))
row.names(hg38_dims)<-colnames(hg38_topic_df)
colnames(hg38_dims)<-c("x","y")
hg38_dims$cellID<-row.names(hg38_dims)
hg38_dims<-merge(hg38_dims,hg38_atac@meta.data,by.x="cellID",by.y="row.names")

mm10_topic_df<-as.data.frame(mm10_cisTopicObject@selected.model$document_expects)
row.names(mm10_topic_df)<-paste0("Topic_",row.names(mm10_topic_df))
mm10_dims<-as.data.frame(uwot::umap(t(mm10_topic_df),n_components=2))
row.names(mm10_dims)<-colnames(mm10_topic_df)
colnames(mm10_dims)<-c("x","y")
mm10_dims$cellID<-row.names(mm10_dims)
mm10_dims<-merge(mm10_dims,mm10_atac@meta.data,by.x="cellID",by.y="row.names")

#plot heatmaps of topics
pdf("mm10_atac_cistopic_heatmap.pdf")
cellTopicHeatmap(mm10_cisTopicObject, method='Probability')
dev.off()

pdf("hg38_atac_cistopic_heatmap.pdf")
cellTopicHeatmap(hg38_cisTopicObject, method='Probability')
dev.off()

#Add cell embeddings into seurat
hg38_cell_embeddings<-as.data.frame(hg38_cisTopicObject@selected.model$document_expects)
colnames(hg38_cell_embeddings)<-hg38_cisTopicObject@cell.names
hg38_n_topics<-nrow(hg38_cell_embeddings)
row.names(hg38_cell_embeddings)<-paste0("topic_",1:hg38_n_topics)
hg38_cell_embeddings<-as.data.frame(t(hg38_cell_embeddings))

mm10_cell_embeddings<-as.data.frame(mm10_cisTopicObject@selected.model$document_expects)
colnames(mm10_cell_embeddings)<-mm10_cisTopicObject@cell.names
mm10_n_topics<-nrow(mm10_cell_embeddings)
row.names(mm10_cell_embeddings)<-paste0("topic_",1:mm10_n_topics)
mm10_cell_embeddings<-as.data.frame(t(mm10_cell_embeddings))

#Add feature loadings into seurat
hg38_feature_loadings<-as.data.frame(hg38_cisTopicObject@selected.model$topics)
row.names(hg38_feature_loadings)<-paste0("topic_",1:hg38_n_topics)
hg38_feature_loadings<-as.data.frame(t(hg38_feature_loadings))

mm10_feature_loadings<-as.data.frame(mm10_cisTopicObject@selected.model$topics)
row.names(mm10_feature_loadings)<-paste0("topic_",1:mm10_n_topics)
mm10_feature_loadings<-as.data.frame(t(mm10_feature_loadings))

#combined cistopic results (cistopic loadings and umap with seurat object)
hg38_cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(hg38_cell_embeddings),loadings=as.matrix(hg38_feature_loadings),assay="peaks",key="topic_")
hg38_umap_dims<-as.data.frame(as.matrix(hg38_dims[2:3]))
colnames(hg38_umap_dims)<-c("UMAP_1","UMAP_2")
row.names(hg38_umap_dims)<-hg38_dims$cellID
hg38_cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(hg38_umap_dims),assay="peaks",key="UMAP_")
hg38_atac@reductions$cistopic<-hg38_cistopic_obj
hg38_atac@reductions$umap<-hg38_cistopic_umap

mm10_cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(mm10_cell_embeddings),loadings=as.matrix(mm10_feature_loadings),assay="peaks",key="topic_")
mm10_umap_dims<-as.data.frame(as.matrix(mm10_dims[2:3]))
colnames(mm10_umap_dims)<-c("UMAP_1","UMAP_2")
row.names(mm10_umap_dims)<-mm10_dims$cellID
mm10_cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(mm10_umap_dims),assay="peaks",key="UMAP_")
mm10_atac@reductions$cistopic<-mm10_cistopic_obj
mm10_atac@reductions$umap<-mm10_cistopic_umap

hg38_n_topics<-ncol(Embeddings(hg38_atac,reduction="cistopic"))

hg38_atac <- FindNeighbors(
  object = hg38_atac,
  reduction = 'cistopic',
  dims = 1:hg38_n_topics
)
hg38_atac <- FindClusters(
  object = hg38_atac,
  verbose = TRUE,
  resolution=0.01 #8 communities
)

mm10_n_topics<-ncol(Embeddings(mm10_atac,reduction="cistopic"))

mm10_atac <- FindNeighbors(
  object = mm10_atac,
  reduction = 'cistopic',
  dims = 1:mm10_n_topics
)
mm10_atac <- FindClusters(
  object = mm10_atac,
  verbose = TRUE,
  resolution=0.02 #targetting roughly 10 communities for gross cell clustering
)

###save Seurat files
saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")
saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

#Plotting 2d projection and clusters

plt<-DimPlot(hg38_atac,group.by=c('seurat_clusters','predicted_doublets',"pcr_idx"))
ggsave(plt,file="hg38.umap.i7idx.pdf",width=10)
#system("slack -F hg38.umap.i7idx.pdf ryan_todo")

plt<-DimPlot(mm10_atac,group.by=c('seurat_clusters','predicted_doublets',"pcr_idx"))
ggsave(plt,file="mm10.umap.i7idx.pdf",width=10)
#system("slack -F mm10.umap.i7idx.pdf ryan_todo")

plt<-FeaturePlot(hg38_atac,feature=c('doublet_scores'))
ggsave(plt,file="hg38.umap.scrub.pdf")
#system("slack -F hg38.umap.scrub.pdf ryan_todo")

plt<-FeaturePlot(mm10_atac,feature=c('doublet_scores'))
ggsave(plt,file="mm10.umap.scrub.pdf")
#system("slack -F mm10.umap.scrub.pdf ryan_todo")


```

## Correcting for systematic bias with harmony


```R
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(harmony,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/lib_backup_210125")

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

#Correcting bias with harmony
pdf("hg38.harmony.convergence.pdf")
harm_mat<-HarmonyMatrix(hg38_atac@reductions$cistopic@cell.embeddings, hg38_atac@meta.data$pcr_idx,do_pca=FALSE,nclust=14,plot_convergence=T)
dev.off()
#system("slack -F hg38.harmony.convergence.pdf ryan_todo")
hg38_atac@reductions$harmony<-CreateDimReducObject(embeddings=as.matrix(harm_mat),assay="peaks",key="topic_")
hg38_atac<-RunUMAP(hg38_atac, reduction = "harmony",dims=1:ncol(hg38_atac@reductions$harmony))
hg38_atac <- FindNeighbors(object = hg38_atac,reduction = 'harmony')
hg38_atac <- FindClusters(object = hg38_atac,verbose = TRUE,resolution=0.05)

plt<-DimPlot(hg38_atac,group.by=c("pcr_idx","seurat_clusters"))
ggsave(plt,file="hg38.umap.i7idx.harm.pdf",width=15)
#system("slack -F hg38.umap.i7idx.harm.pdf ryan_todo")

saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

#Correcting bias with harmony
pdf("mm10.harmony.convergence.pdf")
harm_mat<-HarmonyMatrix(mm10_atac@reductions$cistopic@cell.embeddings, mm10_atac@meta.data$pcr_idx,do_pca=FALSE,nclust=15)
dev.off()

mm10_atac@reductions$harmony<-CreateDimReducObject(embeddings=as.matrix(harm_mat),assay="peaks",key="topic_")
mm10_atac<-RunUMAP(mm10_atac, reduction = "harmony",dims=2:ncol(mm10_atac@reductions$harmony))
mm10_atac <- FindNeighbors(object = mm10_atac,reduction = 'harmony')
mm10_atac <- FindClusters(object = mm10_atac,verbose = TRUE,resolution=0.075 )

plt<-DimPlot(mm10_atac,group.by=c("pcr_idx","seurat_clusters"))
ggsave(plt,file="mm10.umap.i7idx.harm.pdf",width=10)
#system("slack -F mm10.umap.i7idx.harm.pdf ryan_todo")

plt1<-DimPlot(hg38_atac,group.by="seurat_clusters")
plt2<-DimPlot(mm10_atac,group.by="seurat_clusters")
plt<-plt1+plt2
ggsave(plt,file="hg38_mm10_seurat.clusters.pdf")
#system("slack -F hg38_mm10_seurat.clusters.pdf ryan_todo")

saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")


```

### Subclustering 

Going to exclude cells in subclustering that are identified by scrublet as potential doublets.

```R
library(Signac)
library(Seurat)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

#Read in data and modify to monocle CDS file
#read in RDS file.
hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
dir.create("subcluster")

subset_seurat<-function(x,i,prefix){
#Perform cistopic on subclusters of data
    atac_sub<-subset(x,subset=seurat_clusters==i)
    atac_sub<-subset(atac_sub,subset=predicted_doublets=="False") #remove doublets from subclustering
    outname<-paste0(prefix,"_",i,".subset.SeuratObject.Rds")
    saveRDS(atac_sub,paste0("./subcluster/",outname))
    atac_sub<-subset(atac_sub,subset=pcr_idx %in% c("CAGAGAGG","CTCTCTAC")) #limit to just 75k loading
    outname<-paste0(prefix,"_",i,".75k.subset.SeuratObject.Rds")
    saveRDS(atac_sub,paste0("./subcluster/",outname))
}

for (i. in unique(hg38_atac$seurat_clusters)){subset_seurat(hg38_atac,i=i.,prefix="hg38")}
for (i. in unique(mm10_atac$seurat_clusters)){subset_seurat(mm10_atac,i=i.,prefix="mm10")}

```

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
set.seed(1234)

wd<-"/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/subcluster"
setwd(wd)
args = commandArgs(trailingOnly=TRUE)

cistopic_generation<-function(x){
    #Perform cistopic on subclusters of data 
    outname<-paste0(strsplit(x,"[.]")[[1]][1],".",strsplit(x,"[.]")[[1]][2])
    atac_sub<-readRDS(x)
    cistopic_counts_frmt<-atac_sub$peaks@counts
    row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
    sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
    print("made cistopic object")
    sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(22:28),nCores=7,addModels=FALSE)
    saveRDS(sub_cistopic_models,
        file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
    print("finshed running cistopic")

    pdf(paste0(wd,"/",outname,"_model_selection.pdf"))
    par(mfrow=c(3,3))
    sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
    dev.off()
    rm(sub_cistopic_models)
    rm(sub_cistopic)
    rm(atac_sub)
    }

cistopic_generation(x=args[1])

```

```bash
for i in *75k.subset.SeuratObject.Rds;
do Rscript subcluster.R $i ; done &
```

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(cisTopic)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(dplyr)
library(ggrepel)
library(clustree)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/subcluster")

#UMAP Projection and clustering on selected cistopic model
clustering_loop<-function(topicmodel_list.=topicmodel_list,sample,topiccount_list.=topic_count_list){
    #topicmodel_list.=topicmodel_list;topiccount_list.=topic_count_list
    #set up outname
    topicmodel_list.<-topicmodel_list.[sample]

    outname<-strsplit(topicmodel_list.,split="[.]")[[1]][1]
    object_input<-readRDS(paste0(outname,".75k.subset.SeuratObject.Rds"))
    #select_topic
    models_input<-readRDS(topicmodel_list.)
    cisTopicObject<-cisTopic::selectModel(models_input,select=topiccount_list.[topicmodel_list.],keepModels=F)
    
    #perform UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))

    #get cell embeddings
    cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
    colnames(cell_embeddings)<-cisTopicObject@cell.names
    n_topics<-nrow(cell_embeddings)
    row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
    cell_embeddings<-as.data.frame(t(cell_embeddings))
    
    #get feature loadings
    feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
    row.names(feature_loadings)<-paste0("topic_",1:n_topics)
    feature_loadings<-as.data.frame(t(feature_loadings))
    
    #combine with seurat object for celltype seuratobject  
    cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
    object_input@reductions$cistopic<-cistopic_obj
    n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))
    object_input@assays$peaks@key<-"peaks_"
    object_input<-RunUMAP(object_input,reduction="cistopic",dims=1:n_topics)    #finally recluster
    #Clustering with multiple resolutions to account for different celltype complexities
    object_input <- FindNeighbors(object = object_input, reduction = 'cistopic', dims = 1:n_topics)
    object_input <- FindClusters(object = object_input,resolution=0.01)
    object_input <- FindClusters(object = object_input,resolution=0.025)
    object_input <- FindClusters(object = object_input,resolution=0.05)
    object_input <- FindClusters(object = object_input,resolution=0.1)
    object_input <- FindClusters(object = object_input,resolution=0.2)
    object_input <- FindClusters(object = object_input,resolution=0.5)
    object_input <- FindClusters(object = object_input,resolution=0.9)
    
    saveRDS(object_input,paste0(outname,".75k.subset.SeuratObject.Rds"))
    plt<-DimPlot(object_input,group.by=c('peaks_snn_res.0.01','peaks_snn_res.0.025','peaks_snn_res.0.05','peaks_snn_res.0.1','peaks_snn_res.0.2','peaks_snn_res.0.5','peaks_snn_res.0.9'))
    ggsave(plt,file=paste(outname,"clustering.pdf",sep="."))
    system(paste0("slack -F ",paste(outname,"clustering.pdf",sep=".")," ryan_todo"))
    plt<-clustree(object_input, prefix = "peaks_snn_res.")
    ggsave(plt,file=paste(outname,"clustree.pdf",sep="."))
    system(paste0("slack -F ",paste(outname,"clustree.pdf",sep=".")," ryan_todo"))
    object_input<-AddMetaData(object_input,object_input@reductions$umap@cell.embeddings,col.name=c("UMAP_1","UMAP_2"))
    plt<-clustree_overlay(object_input, prefix = "peaks_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2",red_dim="umap")
    ggsave(plt,file=paste(outname,"clustree.overlay.pdf",sep="."))
    system(paste0("slack -F ",paste(outname,"clustree.overlay.pdf",sep=".")," ryan_todo"))

}

####################################
### Processing ###
    topicmodel_list<-list.files(pattern="75k.CisTopicObject.Rds$")

    #determine model count to use for each cell type
    for (i in list.files(pattern="model_selection.pdf")){system(paste0("slack -F ",i," ryan_todo"))}

    #selecting topics based on derivative, making a named vector, note some of these seem to suggest a topic count over 28, but were artificially capped
    topic_count_list<-setNames(c(25,26,27,27,28,25,26,26,25,24,25,27,25,25,27,28,27),topicmodel_list)

    #Running clustering loop
    for (i in 1:length(topic_count_list)){clustering_loop(sample=i)}

    #selecting resolution by plots
    for (i in list.files(pattern="clustering.pdf")){system(paste0("slack -F ",i," ryan_todo"))}
    
    setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

    #Read in data and modify to monocle CDS file
    #read in RDS file.
    hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")

    #adding all subclustering info back into main RDS object
    hg38_atac$seurat_subcluster<-"NA"
    hg38_atac$subcluster_x<-"NA"
    hg38_atac$subcluster_y<-"NA"

    #Assign clustering resolution based on clustering.pdf output
    celltype_list<-list.files(path="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/subcluster",pattern="clustering.pdf")
    hg38_celltype_list<-celltype_list[startsWith(prefix="hg38",celltype_list)]

    hg38_resolution_list<-setNames(c("peaks_snn_res.0.05","peaks_snn_res.0.1","peaks_snn_res.0.05","peaks_snn_res.0.025",
        "peaks_snn_res.0.05","peaks_snn_res.0.01","peaks_snn_res.0.01","peaks_snn_res.0.025"),hg38_celltype_list)
    cell_order<-row.names(hg38_atac@meta.data)

    metadat_embedding<-lapply(hg38_celltype_list, function(celltype.x){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        embedding<-as.data.frame(atac_sub@reductions$umap@cell.embeddings)
        return(embedding)
        })

    metadat_subcluster<-lapply(hg38_celltype_list, function(celltype.x){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        seurat_subcluster<-data.frame(row.names=row.names(atac_sub@meta.data),seruat_subcluster=atac_sub@meta.data[,hg38_resolution_list[celltype.x]])
        return(seurat_subcluster)
        })

    embedding<-do.call("rbind",metadat_embedding)
    seurat_subcluster<-do.call("rbind",metadat_subcluster)
    hg38_atac<-AddMetaData(hg38_atac,embedding,col.name=c("subcluster_x","subcluster_y"))
    hg38_atac<-AddMetaData(hg38_atac,seurat_subcluster,col.name=c("seurat_subcluster"))
    hg38_atac$cluster_ID<-paste(hg38_atac$seurat_clusters,hg38_atac$seurat_subcluster,sep="_")
    as.data.frame(hg38_atac@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
    #na values are doublets or excluded indexes
    saveRDS(hg38_atac,"hg38_SeuratObject.Rds")


    #and mouse
    mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

    #adding all subclustering info back into main RDS object
    mm10_atac$seurat_subcluster<-"NA"
    mm10_atac$subcluster_x<-"NA"
    mm10_atac$subcluster_y<-"NA"

    #Assign clustering resolution based on clustering.pdf output
    celltype_list<-list.files(path="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/subcluster",pattern="clustering.pdf")
    mm10_celltype_list<-celltype_list[startsWith(prefix="mm10",celltype_list)]

    mm10_resolution_list<-setNames(c("peaks_snn_res.0.05","peaks_snn_res.0.05","peaks_snn_res.0.01","peaks_snn_res.0.025",
        "peaks_snn_res.0.01","peaks_snn_res.0.025","peaks_snn_res.0.01","peaks_snn_res.0.01","peaks_snn_res.0.01"),mm10_celltype_list)

    metadat_embedding<-lapply(mm10_celltype_list, function(celltype.x){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        embedding<-as.data.frame(atac_sub@reductions$umap@cell.embeddings)
        return(embedding)
        })

    metadat_subcluster<-lapply(mm10_celltype_list, function(celltype.x){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        seurat_subcluster<-data.frame(row.names=row.names(atac_sub@meta.data),seruat_subcluster=atac_sub@meta.data[,mm10_resolution_list[celltype.x]])
        return(seurat_subcluster)
        })

    embedding<-do.call("rbind",metadat_embedding)
    seurat_subcluster<-do.call("rbind",metadat_subcluster)
    mm10_atac<-AddMetaData(mm10_atac,embedding,col.name=c("subcluster_x","subcluster_y"))
    mm10_atac<-AddMetaData(mm10_atac,seurat_subcluster,col.name=c("seurat_subcluster"))
    mm10_atac$cluster_ID<-paste(mm10_atac$seurat_clusters,mm10_atac$seurat_subcluster,sep="_")
    as.data.frame(mm10_atac@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
    #na values are doublets or excluded indexes
    saveRDS(mm10_atac,"mm10_SeuratObject.Rds")
```

### Recoloring subclusters and plotting

Human

```R
    library(Signac)
    library(Seurat)
    library(SeuratWrappers)
    library(ggplot2)
    library(patchwork)
    library(cicero)
    library(cisTopic)
    library(GenomeInfoDb)
    set.seed(1234)
    library(Matrix)
    library(dplyr)
    library(ggrepel)
    library(RColorBrewer)
    library(palettetown)
    setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

    #Read in data and modify to monocle CDS file
    #read in RDS file.
    hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")

    #set main umap coordinates
    embedding<-as.data.frame(hg38_atac@reductions$umap@cell.embeddings)
    hg38_atac<-AddMetaData(hg38_atac,embedding,col.name=c("umap_x","umap_y"))

    #########Coloring and plotting human data#####################
    dat<-as.data.frame(hg38_atac@meta.data)
    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]

    #Seurat Clusters are Spectral
    clus_col<-setNames(sample(brewer.pal(n = 11, name = "Spectral"),length(unique(dat$seurat_clusters))),nm=sort(unique(dat$seurat_clusters)))
    #clus_col<-setNames("starmie" %>% ichooseyou(length(unique(dat$seurat_clusters))),nm=sort(unique(dat$seurat_clusters))) #maybe porygon?
    dat$cluster_col<-"NULL"
    dat$cluster_col<-clus_col[dat$seurat_clusters]


    sort(unique(dat$cluster_ID))
     #[1] "0_0" "0_1" "0_2" "0_3" "0_4" "0_5" "1_0" "1_1" "2_0" "2_1" "2_2" "3_0"
    #[13] "3_1" "4_0" "4_1" "4_2" "4_3" "5_0" "5_1" "6_0" "7_0" "7_1" "7_2"

    set_colors<-function(obj,i,pallet_list){
        x<-palette_list[i]
        palette<-names(palette_list)[i]
        subclus<-unique(obj@meta.data[obj@meta.data$seurat_clusters==x,]$cluster_ID)
        subclus<-subclus[!endsWith(subclus,"NA")]
        subclus_col<-setNames(palette %>% ichooseyou(length(subclus)),nm=unique(subclus)) #I love it.
        return(subclus_col)
    }

    palette_list<-c(charizard=0, pidgeot=1, weezing=2, rattata=3, noctowl=4, metapod=5, bulbasaur=6, blastoise=7) #these are from palettetown
    subclus_col<-unlist(lapply(1:length(palette_list),function(j) set_colors(obj=hg38_atac,i=j,pallet_list=pallet_list)))
    dat$subcluster_col<-"NULL"
    dat$subcluster_col<-subclus_col[dat$cluster_ID]
    dat<-dat[,c("subcluster_x","subcluster_y","cluster_col","subcluster_col")]
    hg38_atac<-AddMetaData(hg38_atac,dat,col.name=c("subcluster_x","subcluster_y","cluster_col","subcluster_col"))
    saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

    #Perform Plotting
    dat<-as.data.frame(hg38_atac@meta.data)
    dat$subcluster_x<-as.numeric(dat$subcluster_x)
    dat$subcluster_y<-as.numeric(dat$subcluster_y)

    label.df <- data.frame(seurat_clusters=unique(dat$seurat_clusters),label=unique(dat$seurat_clusters))
    label.df_2 <- dat %>% 
      group_by(seurat_clusters) %>% 
      summarize(umap_x = mean(umap_x), umap_y = mean(umap_y)) %>% left_join(label.df)

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=clus_col)+
    ggtitle("hg38")+ ggrepel::geom_label_repel(data = label.df_2, aes(label = label),fontface='bold') +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks = element_blank(),legend.position = "bottom")

    label.df <- data.frame(cluster_ID=unique(dat$cluster_ID),label=unique(dat$cluster_ID))
    label.df$seurat_clusters<-unlist(lapply(strsplit(label.df$cluster_ID,"_"),"[",1))
    label.df_3 <- dat %>% 
      group_by(cluster_ID) %>% 
      summarize(subcluster_x = mean(subcluster_x,na.rm=T), subcluster_y = mean(subcluster_y,na.rm=T)) %>% 
      merge(label.df,by="cluster_ID")

    plt_list<-ggplot(dat,aes(x=subcluster_x,y=subcluster_y,color=cluster_ID))+
    geom_point(alpha=0.05,size=0.5,shape=16)+ theme_bw()+ 
    ggrepel::geom_label_repel(data = label.df_3, aes(x=subcluster_x,y=subcluster_y,label = label), max.iter=10000,direction="both",size=2,force=5, fontface='bold') + 
    scale_color_manual(values=subclus_col) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
    facet_wrap(facets=vars(seurat_clusters),ncol=2) + coord_cartesian(clip = "off") 

    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 
    ggsave(plt,file="hg38_umap.subclus.pdf")
    #system("slack -F hg38_umap.subclus.pdf ryan_todo")

```
Mouse 
```R
    library(Signac)
    library(Seurat)
    library(SeuratWrappers)
    library(ggplot2)
    library(patchwork)
    library(cicero)
    library(cisTopic)
    library(GenomeInfoDb)
    set.seed(1234)
    library(Matrix)
    library(dplyr)
    library(ggrepel)
    library(RColorBrewer)
    library(palettetown)
    setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

    #Read in data and modify to monocle CDS file
    #read in RDS file.
    mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

    #set main umap coordinates
    embedding<-as.data.frame(mm10_atac@reductions$umap@cell.embeddings)
    mm10_atac<-AddMetaData(mm10_atac,embedding,col.name=c("umap_x","umap_y"))

    #########Coloring and plotting human data#####################
    dat<-as.data.frame(mm10_atac@meta.data)
    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]

    #Seurat Clusters are Spectral
    #clus_col<-setNames((brewer.pal(n = length(unique(dat$seurat_clusters)), name = "RdYlBu")),nm=sort(unique(dat$seurat_clusters)))
    clus_col<-setNames(sample(brewer.pal(n = 11, name = "Spectral"),length(unique(dat$seurat_clusters))),nm=sort(unique(dat$seurat_clusters)))
    dat$cluster_col<-"NULL"
    dat$cluster_col<-clus_col[dat$seurat_clusters]


    sort(unique(dat$cluster_ID))
     #[1] "0_0" "0_1" "0_2" "0_3" "0_4" "0_5" "1_0" "1_1" "2_0" "2_1" "2_2" "3_0"
    #[13] "3_1" "4_0" "4_1" "4_2" "4_3" "5_0" "5_1" "6_0" "7_0" "7_1" "7_2"

    set_colors<-function(obj,i,pallet_list){
        x<-palette_list[i]
        palette<-names(palette_list)[i]
        subclus<-unique(obj@meta.data[obj@meta.data$seurat_clusters==x,]$cluster_ID)
        subclus<-subclus[!endsWith(subclus,"NA")]
        subclus_col<-setNames(palette %>% ichooseyou(length(subclus)),nm=unique(subclus)) #I love it.
        return(subclus_col)
    }

    palette_list<-c(charizard=0, pidgeot=1, weezing=2, rattata=3, noctowl=4, metapod=5, bulbasaur=6, blastoise=7) #these are from palettetown
    subclus_col<-unlist(lapply(1:length(palette_list),function(j) set_colors(obj=mm10_atac,i=j,pallet_list=pallet_list)))
    dat$subcluster_col<-"NULL"
    dat$subcluster_col<-subclus_col[dat$cluster_ID]
    dat<-dat[,c("subcluster_x","subcluster_y","cluster_col","subcluster_col")]
    mm10_atac<-AddMetaData(mm10_atac,dat,col.name=c("subcluster_x","subcluster_y","cluster_col","subcluster_col"))
    saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

    #Perform Plotting
    dat<-as.data.frame(mm10_atac@meta.data)
    dat$subcluster_x<-as.numeric(dat$subcluster_x)
    dat$subcluster_y<-as.numeric(dat$subcluster_y)

    label.df <- data.frame(seurat_clusters=unique(dat$seurat_clusters),label=unique(dat$seurat_clusters))
    label.df_2 <- dat %>% 
      group_by(seurat_clusters) %>% 
      summarize(umap_x = mean(umap_x), umap_y = mean(umap_y)) %>% left_join(label.df)

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=clus_col)+
    ggtitle("mm10")+ ggrepel::geom_label_repel(data = label.df_2, aes(label = label),fontface='bold') +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks = element_blank(),legend.position = "bottom")

    label.df <- data.frame(cluster_ID=unique(dat$cluster_ID),label=unique(dat$cluster_ID))
    label.df$seurat_clusters<-unlist(lapply(strsplit(label.df$cluster_ID,"_"),"[",1))
    label.df_3 <- dat %>% 
      group_by(cluster_ID) %>% 
      summarize(subcluster_x = mean(subcluster_x,na.rm=T), subcluster_y = mean(subcluster_y,na.rm=T)) %>% 
      merge(label.df,by="cluster_ID")

    plt_list<-ggplot(dat,aes(x=subcluster_x,y=subcluster_y,color=cluster_ID))+
    geom_point(alpha=0.05,size=0.5,shape=16)+ theme_bw()+ 
    ggrepel::geom_label_repel(data = label.df_3, aes(x=subcluster_x,y=subcluster_y,label = label), max.iter=10000,direction="both",size=2,force=5, fontface='bold') + 
    scale_color_manual(values=subclus_col) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
    facet_wrap(facets=vars(seurat_clusters),ncol=2) + coord_cartesian(clip = "off") 

    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 
    ggsave(plt,file="mm10_umap.subclus.pdf")
    #system("slack -F mm10_umap.subclus.pdf ryan_todo")

```


## Cicero for Coaccessible Networks

Full Cicero Processing. Using CCANs to generate Gene Activity

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(monocle3)
library(cicero)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")


#Cicero processing function
cicero_processing<-function(object_input=hg38_atac,prefix="hg38"){

      #Generate CDS format from Seurat object
      #atac.cds <- as.cell_data_set(object_input,group_by="cluster_ID")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      #atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)
      #saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      atac.cicero<-readRDS(paste(prefix,"atac_cicero_cds.Rds",sep="_"))

      genome <- seqlengths(object_input) # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      #conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      #saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      conns<-readRDS(paste(prefix,"atac_cicero_conns.Rds",sep="_"))

      print("Generating CCANs")
      #ccans <- generate_ccans(conns) # generate ccans
      #saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      ccans<-readRDS(paste(prefix,"atac_cicero_ccans.Rds",sep="_"))

      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      return(object_input)
}

# generate unnormalized gene activity matrix
# gene annotation sample
annotation_generation<-function(ensdb_obj){
      annotations <- GetGRangesFromEnsDb(ensdb = ensdb_obj)
      pos <-as.data.frame(annotations,row.names=NULL)
      pos$chromosome<-paste0("chr",pos$seqnames)
      pos$gene<-pos$gene_id
      pos <- subset(pos, strand == "+")
      pos <- pos[order(pos$start),] 
      pos <- pos[!duplicated(pos$tx_id),] # remove all but the first exons per transcript
      pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS
      neg <-as.data.frame(annotations,row.names=NULL)
      neg$chromosome<-paste0("chr",neg$seqnames)
      neg$gene<-neg$gene_id
      neg <- subset(neg, strand == "-")
      neg <- neg[order(neg$start,decreasing=TRUE),] 
      neg <- neg[!duplicated(neg$tx_id),] # remove all but the first exons per transcript
      neg$end <- neg$end + 1 # make a 1 base pair marker of the TSS
      gene_annotation<- rbind(pos, neg)
      gene_annotation <- gene_annotation[,c("chromosome","start","end","gene_name")] # Make a subset of the TSS annotation columns containing just the coordinates and the gene name
      names(gene_annotation)[4] <- "gene" # Rename the gene symbol column to "gene"
      return(gene_annotation)
}


geneactivity_processing<-function(cds_input,conns_input,prefix,gene_annotation){
      atac.cds<- annotate_cds_by_site(cds_input, gene_annotation)
      unnorm_ga <- build_gene_activity_matrix(atac.cds, conns_input)
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
}


#hg38
hg38_annotation<-annotation_generation(ensdb_obj=EnsDb.Hsapiens.v86)
hg38_atac<-readRDS("hg38_SeuratObject.Rds")
conns<-as.data.frame(readRDS("hg38_atac_cicero_conns.Rds"))
#geneactivity_processing(cds_input=as.cell_data_set(hg38_atac,group_by="seurat_clusters"),conns_input=conns,prefix="hg38",gene_annotation=hg38_annotation)
cicero_gene_activities<-readRDS("hg38.unnorm_GA.Rds")  #Read in unnormalized GA
hg38_atac<-subset(hg38_atac,cells=which(colnames(hg38_atac) %in% colnames(cicero_gene_activities)))
hg38_atac[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
hg38_atac <- NormalizeData(object = hg38_atac,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(hg38_atac$nCount_GeneActivity))  # normalize
saveRDS(hg38_atac,"hg38_SeuratObject.PF.Rds") #this is limited to just cells passing filters (those with cluster IDs)

#mm10
mm10_annotation<-annotation_generation(ensdb_obj=EnsDb.Mmusculus.v79)
mm10_atac<-readRDS("mm10_SeuratObject.Rds")
conns<-as.data.frame(readRDS("mm10_atac_cicero_conns.Rds"))
#geneactivity_processing(cds_input=as.cell_data_set(mm10_atac,group_by="seurat_clusters"),conns_input=conns,prefix="mm10",gene_annotation=mm10_annotation)
cicero_gene_activities<-readRDS("mm10.unnorm_GA.Rds")  #Read in unnormalized GA
mm10_atac<-subset(mm10_atac,cells=which(colnames(mm10_atac) %in% colnames(cicero_gene_activities)))
cicero_gene_activities<-cicero_gene_activities[2:nrow(cicero_gene_activities),] #first feature is empy
mm10_atac[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
mm10_atac <- NormalizeData(object = mm10_atac,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(mm10_atac$nCount_GeneActivity))  # normalize
saveRDS(mm10_atac,"mm10_SeuratObject.PF.Rds")


```

## Public Data RNA Comparison
### Download data from Allen Brain-span
For human: https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
For mouse: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x

```bash
#Human download
cd /home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds
#Mouse download
cd /home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/87/14/8714d0a3-27d7-4a81-8c77-eebfd605a280/readme_mouse_10x.txt
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv 
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/trimmed_means.csv
wget http://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.RData
#Downloading Mouse whole brain data 
wget https://www.dropbox.com/s/kqsy9tvsklbu7c4/allen_brain.rds?dl=0
#Downloading Mouse Cerebrellum data
#downloaded using a curl command valid for 30 min. generated from https://singlecell.broadinstitute.org/single_cell/study/SCP795/a-transcriptomic-atlas-of-the-mouse-cerebellum#/
#file located here:
/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse/cerebellum/SCP795/other/cb_annotated_object.RDS

```

### Process Data into Seurat Object
Following https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

```R
library(Seurat)
library(ggplot2)
#Human
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex")
meta_data<-read.csv("metadata.csv",header=T)
row.names(meta_data)<-meta_data$sample_name
counts<-read.csv("matrix.csv",header=T,row.names=1)
brainspan <- CreateSeuratObject(counts = as.data.frame(t(counts)), project = "brainspain", min.cells = 3, min.features = 500, meta.data=meta_data)
saveRDS(brainspan, file = "allen_brainspan_humancortex.rds")
brainspan <- NormalizeData(brainspan, normalization.method = "LogNormalize", scale.factor = 10000)
brainspan <- FindVariableFeatures(brainspan, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(brainspan)
brainspan <- ScaleData(brainspan, features = all.genes)
brainspan <- RunPCA(brainspan, features = VariableFeatures(object = brainspan))
plt<-ElbowPlot(brainspan)
ggsave(plt,file="allen_brainspan_humancortex.elbowplot.pdf")
#system("slack -F allen_brainspan_humancortex.elbowplot.pdf ryan_todo")
brainspan <- FindNeighbors(brainspan, dims = 1:14)
brainspan <- FindClusters(brainspan, resolution = 0.5)
brainspan <- RunUMAP(brainspan, dims = 1:14)
plt<-DimPlot(brainspan, reduction = "umap",group.by=c("class_label","subclass_label"))
ggsave(plt,file="allen_brainspan_humancortex.dimplot.pdf",width=30)
#system("slack -F allen_brainspan_humancortex.dimplot.pdf ryan_todo")
saveRDS(brainspan, file = "allen_brainspan_humancortex.rds")

#Mouse
 setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse")
 meta_data<-read.csv("metadata.csv",header=T)
 row.names(meta_data)<-meta_data$sample_name
 counts<-read.csv("matrix.csv",header=T,row.names=1,nrows=100000) #this is 1million cells, i think 10% of that will be fine, hopefully they are evenly distributed
 brainspan <- CreateSeuratObject(counts = as.data.frame(t(counts)), project = "brainspain", min.cells = 3, min.features = 500, meta.data=meta_data)
 saveRDS(brainspan, file = "allen_brainspan_mouse.rds")
 brainspan <- NormalizeData(brainspan, normalization.method = "LogNormalize", scale.factor = 10000)
 brainspan <- FindVariableFeatures(brainspan, selection.method = "vst", nfeatures = 2000)
 all.genes <- rownames(brainspan)
 brainspan <- ScaleData(brainspan, features = all.genes)
 brainspan <- RunPCA(brainspan, features = VariableFeatures(object = brainspan))
 plt<-ElbowPlot(brainspan)
 ggsave(plt,file="allen_brainspan_mouse.elbowplot.pdf")
 #system("slack -F allen_brainspan_mouse.elbowplot.pdf ryan_todo")
 brainspan <- FindNeighbors(brainspan, dims = 1:15)
 brainspan <- FindClusters(brainspan, resolution = 0.5)
 brainspan <- RunUMAP(brainspan, dims = 1:15)
 plt<-DimPlot(brainspan, reduction = "umap",group.by=c("class_label","subclass_label"))
 ggsave(plt,file="allen_brainspan_mouse.dimplot.pdf",width=30)
 #system("slack -F allen_brainspan_mouse.dimplot.pdf ryan_todo")
 saveRDS(brainspan, file = "allen_brainspan_mouse.rds")




```

### Integration of ATAC and RNA for cell type identification

Follow this https://satijalab.org/seurat/articles/atacseq_integration_vignette.html#co-embedding-scrna-seq-and-scatac-seq-datasets-1

Retry mouse cluster id just straight up following https://satijalab.org/signac/articles/mouse_brain_vignette.html?

Human
```R
library(Seurat)
library(Signac)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")


prediction_transfer<-function(x,brainspan,feat_name,feat_col,obj){
    metadata_feat_name<-paste0("brainspan.predicted.",feat_name)
    metadata_col_name<-paste0("brainspan.predicted.",feat_col)
    Assay_Name<-paste0("predicted_",feat_name)
    named_vec<-setNames(x$predicted.id,nm=row.names(x)) #set features
    feat_assay<-CreateAssayObject(data = t(x[,2:ncol(x)])) #remove predicted id text
    col_df<-brainspan@meta.data[c(feat_name,feat_col)]
    col_df<-col_df[!duplicated(col_df),]
    cols<-setNames(col_df[,feat_col],nm=col_df[,feat_name])
    named_vec_color<-setNames(cols[named_vec],nm=names(named_vec)) #set named color vector
    obj <- AddMetaData(object=obj, metadata=named_vec,col.name=metadata_feat_name)
    obj <- AddMetaData(object=obj, metadata=named_vec_color,col.name=metadata_col_name)
    obj[[Assay_Name]]<-feat_assay
    return(obj)
}


predict_celltype<-function(object,brainspan,prefix){
    DefaultAssay(object)<-"GeneActivity"
    object<-ScaleData(object,features = row.names(object))
    transfer.anchors <- FindTransferAnchors(reference = brainspan,query = object,query.assay="GeneActivity",reduction = 'cca',features=VariableFeatures(object=brainspan))
    saveRDS(transfer.anchors,file=paste0(prefix,".transferanchors.rds"))
    transfer.anchors<-readRDS(file=paste0(prefix,".transferanchors.rds"))
    #predict labels for class and subclass
    predicted.labels.cluster <- TransferData(anchorset = transfer.anchors,refdata = brainspan$cluster_label,weight.reduction = object[["cistopic"]],dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.cluster,brainspan=brainspan,feat_name="cluster_label",feat_col="cluster_color",obj=object)
    predicted.labels.class <- TransferData(anchorset = transfer.anchors,refdata = brainspan$class_label,weight.reduction = "cca",dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.class,brainspan=brainspan,feat_name="class_label",feat_col="class_color",obj=object)
    predicted.labels.subclass <- TransferData(anchorset = transfer.anchors,refdata = brainspan$subclass_label,weight.reduction = "cca", dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.subclass,brainspan=brainspan,feat_name="subclass_label",feat_col="subclass_color",obj=object)
    #remove any metadata columns labelled "prediction"
    object@meta.data<-object@meta.data[!startsWith(colnames(object@meta.data),prefix="prediction.")]
    saveRDS(object,file=paste0(prefix,"_SeuratObject.PF.Rds"))
    return(object)

}

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
hg38_atac<-predict_celltype(object=hg38_atac,brainspan=brainspan.,prefix="hg38")
#update assay names for readability
hg38_atac<-RenameAssays(hg38_atac,
    predicted_class_label="allenbrainmap_class_prediction_values",
    predicted_subclass_label="allenbrainmap_subclass_preduction_values",
    predicted_cluster_label="allenbrainmap_cluserprediction_values")
saveRDS(hg38_atac,"hg38_SeuratObject.PF.Rds")

```

Mouse
```R
library(Seurat)
library(Signac)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")


prediction_transfer<-function(x,brainspan,feat_name,feat_col,obj){
    metadata_feat_name<-paste0("brainspan.predicted.",feat_name)
    metadata_col_name<-paste0("brainspan.predicted.",feat_col)
    Assay_Name<-paste0("predicted_",feat_name)
    named_vec<-setNames(x$predicted.id,nm=row.names(x)) #set features
    feat_assay<-CreateAssayObject(data = t(x[,2:ncol(x)])) #remove predicted id text
    col_df<-brainspan@meta.data[c(feat_name,feat_col)]
    col_df<-col_df[!duplicated(col_df),]
    cols<-setNames(col_df[,feat_col],nm=col_df[,feat_name])
    named_vec_color<-setNames(cols[named_vec],nm=names(named_vec)) #set named color vector
    obj <- AddMetaData(object=obj, metadata=named_vec,col.name=metadata_feat_name)
    obj <- AddMetaData(object=obj, metadata=named_vec_color,col.name=metadata_col_name)
    obj[[Assay_Name]]<-feat_assay
    return(obj)
}


predict_celltype<-function(object,brainspan,prefix){
    DefaultAssay(object)<-"GeneActivity"
    object<-ScaleData(object,features = row.names(object))
    transfer.anchors <- FindTransferAnchors(reference = brainspan,query = object,query.assay="GeneActivity",reduction = 'cca',features=VariableFeatures(object=brainspan))
    saveRDS(transfer.anchors,file=paste0(prefix,".transferanchors.rds"))
    transfer.anchors<-readRDS(file=paste0(prefix,".transferanchors.rds"))
    #predict labels for class and subclass
    predicted.labels.cluster <- TransferData(anchorset = transfer.anchors,refdata = brainspan$cluster_label,weight.reduction = object[["cistopic"]],dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.cluster,brainspan=brainspan,feat_name="cluster_label",feat_col="cluster_color",obj=object)
    predicted.labels.class <- TransferData(anchorset = transfer.anchors,refdata = brainspan$class_label,weight.reduction = "cca",dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.class,brainspan=brainspan,feat_name="class_label",feat_col="class_color",obj=object)
    predicted.labels.subclass <- TransferData(anchorset = transfer.anchors,refdata = brainspan$subclass_label,weight.reduction = "cca", dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer(x=predicted.labels.subclass,brainspan=brainspan,feat_name="subclass_label",feat_col="subclass_color",obj=object)
    #remove any metadata columns labelled "prediction"
    object@meta.data<-object@meta.data[!startsWith(colnames(object@meta.data),prefix="prediction.")]
    saveRDS(object,file=paste0(prefix,"_SeuratObject.PF.Rds"))
    return(object)

}

prediction_transfer_cerebellar<-function(x,ref,feat_name,obj){
    metadata_feat_name<-paste0("ref.predicted.",feat_name)
    Assay_Name<-paste0("predicted_",feat_name)
    named_vec<-setNames(x$predicted.id,nm=row.names(x)) #set features
    feat_assay<-CreateAssayObject(data = t(x[,2:ncol(x)])) #remove predicted id text
    col_df<-ref@meta.data[c(feat_name)]
    col_df<-col_df[!duplicated(col_df),]
    obj <- AddMetaData(object=obj, metadata=named_vec,col.name=metadata_feat_name)
    obj[[Assay_Name]]<-feat_assay
    return(obj)
}

#running cerebellar RNA data as well for whole brain mouse
predict_celltype_cerebellar<-function(object=mm10,ref=cerebellar,prefix="mm10"){
    DefaultAssay(object)<-"GeneActivity"
    object<-ScaleData(object,features = row.names(object))
    transfer.anchors <- FindTransferAnchors(reference = ref,query = object,query.assay="GeneActivity",reduction = 'cca',features=VariableFeatures(object=ref))
    saveRDS(transfer.anchors,file=paste0(prefix,".cerebellar.transferanchors.rds"))
    transfer.anchors<-readRDS(file=paste0(prefix,".cerebellar.transferanchors.rds"))
    #predict labels for cluster and subcluster
    predicted.labels.cluster <- TransferData(anchorset = transfer.anchors,refdata = ref$cluster,weight.reduction = object[["cistopic"]],dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer_cerebellar(x=predicted.labels.cluster,ref=ref,feat_name="cluster",obj=object)
    predicted.labels.subcluster<- TransferData(anchorset = transfer.anchors,refdata = ref$subcluster,weight.reduction = "cca",dims = 1:dim(object[["cistopic"]])[2])
    object<-prediction_transfer_cerebellar(x=predicted.labels.subcluster,ref=ref,feat_name="subcluster",obj=object)
    return(object)
}


#Mouse using smaller data set that is whole brain
    mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
    brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse/allen_brainspan_mouse.rds")
    brainspan. <- FindVariableFeatures(object = brainspan.,nfeatures = 5000)
    mm10_atac<-predict_celltype(object=mm10_atac,brainspan=brainspan.,prefix="mm10") #save RDS contained within function

#for mouse (whole brain) also going to integrate with cerebellar data set
    #mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
    cerebellar <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse/cerebellum/SCP795/other/cb_annotated_object.RDS")
    cerebellar<-UpdateSeuratObject(cerebellar) #update object
    cerebellar<-subset(cerebellar, downsample=2000) #downsample object so max number per cluster is 2000 cells
    cerebellar <- FindVariableFeatures(object = cerebellar,nfeatures = 5000)
    mm10_atac<-predict_celltype_cerebellar(object=mm10_atac,ref=cerebellar,prefix="mm10")
    #clean up prediction values to make the object more readable
    mm10_atac<-RenameAssays(mm10_atac,
        predicted_cluster="cerebellum_cluster_prediction_values",
        predicted_subcluster="cerebellum_subcluster_prediction_values",
        predicted_class_label="allenbrainmap_class_prediction_values",
        predicted_subclass_label="allenbrainmap_subclass_preduction_values",
        predicted_cluster_label="allenbrainmap_clusterprediction_values"
        )
    colnames(mm10_atac@meta.data)[which(colnames(mm10_atac@meta.data)=="ref.predicted.subcluster")]<-"cerebellum.predicted.subcluster"
    colnames(mm10_atac@meta.data)[which(colnames(mm10_atac@meta.data)=="ref.predicted.cluster")]<-"cerebellum.predicted.cluster"
saveRDS(mm10_atac,file="mm10_SeuratObject.PF.Rds")

```

### Confusion Matrices for Cell Type Identification

Human
```R
library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
library(ggdendro)
library(dendextend)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(circlize)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
hg38_atac<-subset(hg38_atac,cells=row.names(hg38_atac@meta.data[!endsWith(hg38_atac@meta.data$cluster_ID,suffix="NA"),]))

build_confusion_matrix<-function(obj,feat="brainspan.predicted.class_label"){
    #seurat_clusters
    confusion_matrix<-as.data.frame(table(obj$seurat_clusters,obj@meta.data[,which(colnames(obj@meta.data)==feat)]))
    confusion_matrix$Freq<-as.numeric(confusion_matrix$Freq)
    confusion_matrix<-reshape2::dcast(Var1~Var2,value.var="Freq",data=confusion_matrix,fill=0,fun.aggregate=sum)
    row.names(confusion_matrix)<-confusion_matrix[,1]
    confusion_matrix<-confusion_matrix[,2:ncol(confusion_matrix)]
    confusion_matrix<-t(apply(confusion_matrix, 1, function(x) x/sum(x)))
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    plt1<-Heatmap(confusion_matrix,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp=gpar(fontsize=7),
        column_names_rot=90,
        col=col_fun
    )
    pdf(paste0("hg38.",x,"confusion_mat.heatmap.pdf"),height=20,width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0("hg38.",x,"confusion_mat.heatmap.pdf")," ryan_todo"))

    #cluster_ID
    confusion_matrix<-as.data.frame(table(obj$cluster_ID,obj@meta.data[,which(colnames(obj@meta.data)==feat)]))
    confusion_matrix$Freq<-as.numeric(confusion_matrix$Freq)
    confusion_matrix<-reshape2::dcast(Var1~Var2,value.var="Freq",data=confusion_matrix,fill=0,fun.aggregate=sum)
    row.names(confusion_matrix)<-confusion_matrix[,1]
    confusion_matrix<-confusion_matrix[,2:ncol(confusion_matrix)]
    confusion_matrix<-t(apply(confusion_matrix, 1, function(x) x/sum(x)))
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    plt1<-Heatmap(confusion_matrix,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp=gpar(fontsize=7),
        column_names_rot=90,
        col=col_fun
    )
    pdf(paste0("hg38.",x,"confusion_mat.clusterID.heatmap.pdf"),height=20,width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0("hg38.",x,"confusion_mat.clusterID.heatmap.pdf")," ryan_todo"))
    return(confusion_matrix)
}

#hg38
for (x in c("brainspan.predicted.class_label","brainspan.predicted.subclass_label","brainspan.predicted.cluster_label")){
build_confusion_matrix(obj=hg38_atac,feat=x)
}

confusion_matrix<-build_confusion_matrix(obj=hg38_atac,feat="brainspan.predicted.cluster_label")
confusion_matrix<-t(confusion_matrix)

#cluster by all marker genes
sum_da_dend <- t(confusion_matrix) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="hg38.confusion.dend.rds") 

annot<-hg38_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(confusion_matrix),]
annot<-annot[match(colnames(confusion_matrix),annot$cluster_ID),]
annot_clus_col<-annot[!duplicated(annot$cluster_ID),]
col_fun = colorRamp2(c(0, 1), c("white", "red"))

top_ha<-columnAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
                        ),
                    show_legend = c(TRUE, TRUE, FALSE))


plt1<-Heatmap(confusion_matrix,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90,
    top_annotation=top_ha,
    cluster_columns=sum_da_dend,
    col=col_fun
)
pdf(paste0("hg38_mergeddataset_confusion_mat.heatmap.pdf"),height=20,width=20)
print(plt1)
dev.off()
system(paste0("slack -F ",paste0("hg38_mergeddataset_confusion_mat.heatmap.pdf")," ryan_todo"))


```

Mouse
```R
library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
library(ggdendro)
library(dendextend)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(circlize)
library(ggdendro)
library(dendextend)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
mm10_atac<-subset(mm10_atac,cells=row.names(mm10_atac@meta.data[!endsWith(mm10_atac@meta.data$cluster_ID,suffix="NA"),]))

build_confusion_matrix<-function(obj,feat="brainspan.predicted.class_label"){
    #seurat_clusters
    confusion_matrix<-as.data.frame(table(obj$seurat_clusters,obj@meta.data[,which(colnames(obj@meta.data)==feat)]))
    confusion_matrix$Freq<-as.numeric(confusion_matrix$Freq)
    confusion_matrix<-reshape2::dcast(Var1~Var2,value.var="Freq",data=confusion_matrix,fill=0,fun.aggregate=sum)
    row.names(confusion_matrix)<-confusion_matrix[,1]
    confusion_matrix<-confusion_matrix[,2:ncol(confusion_matrix)]
    confusion_matrix<-t(apply(confusion_matrix, 1, function(x) x/sum(x)))
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    plt1<-Heatmap(confusion_matrix,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp=gpar(fontsize=7),
        column_names_rot=90,
        col=col_fun
    )
    pdf(paste0("mm10.",x,"confusion_mat.heatmap.pdf"),height=20,width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0("mm10.",x,"confusion_mat.heatmap.pdf")," ryan_todo"))

    #cluster_ID
    confusion_matrix<-as.data.frame(table(obj$cluster_ID,obj@meta.data[,which(colnames(obj@meta.data)==feat)]))
    confusion_matrix$Freq<-as.numeric(confusion_matrix$Freq)
    confusion_matrix<-reshape2::dcast(Var1~Var2,value.var="Freq",data=confusion_matrix,fill=0,fun.aggregate=sum)
    row.names(confusion_matrix)<-confusion_matrix[,1]
    confusion_matrix<-confusion_matrix[,2:ncol(confusion_matrix)]
    confusion_matrix<-t(apply(confusion_matrix, 1, function(x) x/sum(x)))
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    plt1<-Heatmap(confusion_matrix,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp=gpar(fontsize=7),
        column_names_rot=90,
        col=col_fun
    )
    pdf(paste0("mm10.",x,"confusion_mat.clusterID.heatmap.pdf"),height=20,width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0("mm10.",x,"confusion_mat.clusterID.heatmap.pdf")," ryan_todo"))
}

#mm10
for (x in c("brainspan.predicted.class_label","brainspan.predicted.subclass_label","brainspan.predicted.cluster_label","cerebellum.predicted.cluster","cerebellum.predicted.subcluster")){
build_confusion_matrix(obj=mm10_atac,feat=x)
}

#Merge prediction matrices across data sets and select top cells
dat<-rbind(mm10_atac[["allenbrainmap_subclass_preduction_values"]]@data,mm10_atac[["cerebellum_cluster_prediction_values"]]@data)
dat<-data.frame(cellID=colnames(dat),celltype=row.names(dat)[apply(dat,2,which.max)])
dat$cluster<-mm10_atac$cluster_ID[names(mm10_atac$cluster_ID) %in% dat$cellID]
dat<-as.data.frame(table(dat$cluster,dat$celltype))
dat$Var2<-sub("prediction.score.","",dat$Var2)
confusion_matrix<-reshape2::dcast(Var1~Var2,value.var="Freq",data=dat,fill=0,fun.aggregate=sum)
row.names(confusion_matrix)<-confusion_matrix[,1]
confusion_matrix<-confusion_matrix[,2:ncol(confusion_matrix)]
confusion_matrix<-t(apply(confusion_matrix, 1, function(x) x/sum(x)))
confusion_matrix<-t(confusion_matrix)
col_fun = colorRamp2(c(0, 1), c("white", "red"))

#cluster by all marker genes
sum_da_dend <- t(confusion_matrix) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="mm10.confusion.dend.rds") 

annot<-mm10_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(confusion_matrix),]
annot<-annot[match(colnames(confusion_matrix),annot$cluster_ID),]
annot_clus_col<-annot[!duplicated(annot$cluster_ID),]

top_ha<-columnAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
                        ),
                    show_legend = c(TRUE, TRUE, FALSE))


plt1<-Heatmap(confusion_matrix,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90,
    top_annotation=top_ha,
    cluster_columns=sum_da_dend,
    col=col_fun
)
pdf(paste0("mm10_mergeddataset_confusion_mat.heatmap.pdf"),height=20,width=20)
print(plt1)
dev.off()
system(paste0("slack -F ",paste0("mm10_mergeddataset_confusion_mat.heatmap.pdf")," ryan_todo"))

```

### Add TF Motif Usage through ChromVAR

### ChromVar for Transcription Factor Motifs

Human
```R
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)

#lowerign cores to be used by chromvar to 10
library(BiocParallel)
register(MulticoreParam(3))

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
DefaultAssay(hg38_atac)<-"peaks"
  #Read in data and modify to monocle CDS file
  #read in RDS file.

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
x = JASPAR2020,
opts = list(species =9606, all_versions = FALSE))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
motif.matrix.hg38 <- CreateMotifMatrix(
features = granges(hg38_atac[["peaks"]]),
pwm = pfm,
genome = 'hg38',
use.counts = FALSE)

# Create a new Mofif object to store the results
motif.hg38 <- CreateMotifObject(
data = motif.matrix.hg38,
pwm = pfm)
# Add the Motif object to the assays and run ChromVar

#for human
hg38_atac <- SetAssayData(
object = hg38_atac,
assay = 'peaks',
slot = 'motifs',
new.data = motif.hg38)

hg38_atac <- RegionStats(object = hg38_atac, genome = BSgenome.Hsapiens.UCSC.hg38)
hg38_atac <- RunChromVAR( object = hg38_atac,genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(hg38_atac,file="hg38_SeuratObject.PF.Rds")


```

Mouse
```R
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)

#lowerign cores to be used by chromvar to 10
library(BiocParallel)
register(MulticoreParam(2))
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")

DefaultAssay(mm10_atac)<-"peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
x = JASPAR2020,
opts = list(species =9606, all_versions = FALSE))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)

motif.matrix.mm10 <- CreateMotifMatrix(
features = granges(mm10_atac[["peaks"]]),
pwm = pfm,
genome = 'mm10',
use.counts = FALSE)


motif.mm10 <- CreateMotifObject(
data = motif.matrix.mm10,
pwm = pfm)
# Add the Motif object to the assays and run ChromVar


#mouse
mm10_atac <- SetAssayData(
object = mm10_atac,
assay = 'peaks',
slot = 'motifs',
new.data = motif.mm10)

mm10_atac <- RegionStats(object = mm10_atac, genome = BSgenome.Mmusculus.UCSC.mm10)
mm10_atac <- RunChromVAR( object = mm10_atac,genome = BSgenome.Mmusculus.UCSC.mm10)
saveRDS(mm10_atac,file="mm10_SeuratObject.PF.Rds")
```

### Differential Gene Activity through Subclusters

```R
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(ggdendro)
library(dendextend)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(circlize)
library(TFBSTools)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
hg38_atac$cluster_ID<-paste(hg38_atac$seurat_clusters,hg38_atac$seurat_subcluster,sep="_")
hg38_atac$celltype<-"NA"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==0,]$celltype<-"ExN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==1,]$celltype<-"Oligo"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==2,]$celltype<-"Astro"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==3,]$celltype<-"iN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==4,]$celltype<-"iN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==5,]$celltype<-"Micro.PVM"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==6,]$celltype<-"OPC"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==7,]$celltype<-"NonN"
hg38_atac$celltype_col<-"NA"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="iN",]$celltype_col<-"#cc3366"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="ExN",]$celltype_col<-"#3399cc"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Oligo",]$celltype_col<-"#66cc99"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Astro",]$celltype_col<-"#99cc99"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Micro.PVM",]$celltype_col<-"#ff6633"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="OPC",]$celltype_col<-"#ffcc99"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="NonN",]$celltype_col<-"#808080"

saveRDS(hg38_atac,"hg38_SeuratObject.PF.Rds")


hg38_atac<-subset(hg38_atac,cells=which(!endsWith(hg38_atac@meta.data$cluster_ID,"NA")))

#set clusters to test
cluster_to_test<-unique(hg38_atac$seurat_clusters)
#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group,assay.="GeneActivity",latent.vars.="nCount_GeneActivity"){
    da_ga_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = latent.vars.,
        only.pos=T,
        assay=assay.
        )
    da_ga_tmp$da_region<-row.names(da_ga_tmp)
    da_ga_tmp$enriched_group<-c(i)
    da_ga_tmp$compared_group<-c("all_other_cells")
    return(da_ga_tmp)
  }

#Gene Activity
da_ga<-list() #set up an empty list for looping through
n.cores=8 #Perform parallel application of DA test
da_ga<-lapply(
    cluster_to_test,
    FUN=da_one_v_rest,
    obj=hg38_atac,
    group="seurat_clusters",
    assay.="GeneActivity")
da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
write.table(da_ga_df,file="hg38.onevrest.da_ga.txt",sep="\t",col.names=T,row.names=T,quote=F)

#Chromvar TFs
da_ga<-list() #set up an empty list for looping through
n.cores=10 #Perform parallel application of DA test
da_ga<-mclapply(
    cluster_to_test,
    FUN=da_one_v_rest,
    obj=hg38_atac,
    group="seurat_clusters",
    assay.="chromvar", latent.vars.="nCount_peaks",
    mc.cores=n.cores)
da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
#To convert JASPAR ID TO TF NAME
da_ga_df$tf_name <- unlist(lapply(unlist(lapply(da_ga_df$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
write.table(da_ga_df,file="hg38.onevrest.da_chromvar.txt",sep="\t",col.names=T,row.names=T,quote=F)

#Plot out top ga for each cluster
da_ga<-read.csv(file="hg38.onevrest.da_ga.txt",head=T,sep="\t",row.names=NULL)
da_ga$gene_name<-da_ga$da_region
da_ga<-da_ga[complete.cases(da_ga),]

da_ga$label<-""
for (x in unique(da_ga$enriched_group)){
selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:3))$da_region
da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_ga<-as.data.frame(t(as.data.frame(hg38_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,hg38_atac$seurat_clusters) #group by rows to seurat clusters
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame

sum_ga<-t(scale(sum_ga))
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

#cluster by all marker genes
sum_da_dend <- t(sum_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="hg38.geneactivity.dend.rds") 

########## RUNNING #########

sum_ga<-sum_ga[row.names(sum_ga) %in% unique(da_ga$label),]

annot<-hg38_atac@meta.data[,c("celltype","seurat_clusters","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$seurat_clusters),]
annot<-annot[annot$seurat_clusters %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$seurat_clusters),]
sum_ga_plot<-t(sum_ga)

annot_clus_col<-annot[!duplicated(annot$seurat_clusters),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$seurat_clusters),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$seurat_clusters) #due to nonunique colors present
                        ))

bottom_ha<-columnAnnotation(foo = anno_mark(at = 1:ncol(sum_ga_plot), labels = colnames(sum_ga_plot)))

colfun=colorRamp2(quantile(unlist(sum_ga_plot), probs=c(0.5,0.90,0.95)),magma(3))
plt1<-Heatmap(sum_ga_plot,
    cluster_rows=sum_da_dend,
    left_annotation=side_ha,
    col=colfun,
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

pdf("hg38.geneactivity.heatmap.pdf",height=20,width=20)
plt1
dev.off()
#system("slack -F hg38.geneactivity.heatmap.pdf ryan_todo")

```


## Plotting heatmap using brainspan given markers

```R
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(ggdendro)
library(dendextend)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(circlize)
library(TFBSTools)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
sum_da_dend<-readRDS(file="hg38.geneactivity.dend.rds")
dat_ga<-as.data.frame(t(as.data.frame(hg38_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,hg38_atac$seurat_clusters) #group by rows to seurat clusters 
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame
sum_ga<-t(sum_ga)

markers_limited<-   as.data.frame(rbind(
    c("Mark","L1-6","GABAergic","GAD1"),
    c("Mark","L1-6","GABAergic","LHX6"),
    c("Mark","L1-6","GABAergic","ADARB2"),
    c("Mark","L1-6","GABAergic","LAMP5"),
    c("Mark","L1-6","GABAergic","VIP"),
    c("Mark","L1-6","GABAergic","SVIL"),
    c("Mark","L1-6","GABAergic","NXPH2"),
    c("Mark","L1-6","GABAergic","SST"),
    c("Mark","L1-6","Glutamatergic","SLC17A7"),
    c("Mark","L1-6","Glutamatergic","HTR1F"),
    c("Mark","L1-6","Glutamatergic","FEZF2"),
    c("Mark","L1-6","Glutamatergic","RORB"),
    c("Mark","L1-6","Glutamatergic","PCP4"),
    c("Mark","L1-6","Glutamatergic","POU3F2"),
    c("Mark","L1-6","Glutamatergic","CUX2"),
    c("Mark","L1-6","Glutamatergic","CALB1"),
    c("Mark","L1-6","Glutamatergic","RASGRF2"),
    c("Mark","L1-6","Astrocyte","AQP4"),
    c("Mark","L1-6","Astrocyte","ALDH1L1"),
    c("Mark","L1-6","Astrocyte","GFAP"),
    c("Mark","L1-6","Astrocyte","MFGE8"), #CAT
    c("Mark","L1-6","OPC","PDGFRA"),
    c("Mark","L1-6","OPC","CSPG4"),
    c("Mark","L1-6","Oligodendrocyte","OPALIN"),
    c("Mark","L1-6","Oligodendrocyte","PROX1"),
    c("Mark","L1-6","Oligodendrocyte","OLIG1"),
    c("Mark","L1-6","Oligodendrocyte","MOBP"),
    c("Mark","L1-6","Oligodendrocyte","ENPP4"), #CAT
    c("Mark","L1-6","Oligodendrocyte","LDLRAP1"), #CAT
    c("Mark","L1-6","Microglia","FYB"),
    c("Mark","L1-6","Microglia","C1QA"),
    c("Mark","L1-6","Microglia","C1QC"),
    c("Mark","L1-6","Endothelia","FLT1"),
    c("Mark","L1-6","Endothelia","KDR"),
    c("Mark","L1-6","Astrocyte","SLC1A2")
    ))
colnames(markers_limited)<-c("celltype","layer","marker_1","marker_2") #i decided to ignore layer markers so made them all L1-6

sum_ga_sub<-sum_ga[match(markers_limited$marker_2,row.names(sum_ga),nomatch=0),]
markers_limited<-markers_limited[match(row.names(sum_ga_sub),markers_limited$marker_2,nomatch=0),]

sum_ga_sub<-t(scale(t(sum_ga_sub),center=F,scale=T))

annot<-hg38_atac@meta.data[,c("celltype","cluster_col","seurat_clusters","celltype_col")] 
annot<-annot[!duplicated(annot$seurat_clusters),]
annot<-annot[annot$seurat_clusters %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$seurat_clusters),]
annot_clus_col<-annot[!duplicated(annot$seurat_clusters),]

top_ha<-columnAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters)))
                        ),
                    show_legend = c(TRUE, TRUE))


colfun=colorRamp2(quantile(unlist(sum_ga_sub), probs=c(0.5,0.90,0.95)),magma(3))

plt1<-Heatmap(sum_ga_sub,
    clustering_distance_columns="euclidean",
    clustering_distance_rows="euclidean",
    row_order=1:nrow(sum_ga_sub),
    column_split=annot$celltype,
    row_split=paste(markers_limited$marker_1),
    cluster_row_slices=F,
    top_annotation=top_ha,
    left_annotation=rowAnnotation(df=data.frame(celltype_mark=markers_limited$marker_1)),
    #left_annotation=rowAnnotation(df=data.frame(celltype_mark=markers$celltype,layers_mark=markers$layer)),
    col=colfun,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp=gpar(fontsize=3)
)

pdf("hg38.geneactivity.markers.heatmap.pdf",height=5)
plt1
dev.off()
#system("slack -F hg38.geneactivity.markers.heatmap.pdf ryan_todo")



########################Plot out top TF for each cluster###################
da_tf<-read.csv(file="hg38.onevrest.da_chromvar.txt",head=T,sep="\t",row.names=NULL)
da_tf$gene_name<-da_tf$da_region
da_tf<-da_tf[complete.cases(da_tf),]
da_tf<-da_tf[!endsWith(da_tf$enriched_group,"NA"),]

da_tf$label<-""
for (x in unique(da_tf$enriched_group)){
selc_genes<-as.data.frame(da_tf %>% filter(enriched_group==x) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:3))$tf_name
da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$label<- da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$tf_name
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_tf<-as.data.frame(t(as.data.frame(hg38_atac[["chromvar"]]@data)))
sum_tf<-split(dat_tf,hg38_atac$seurat_clusters) #group by rows to seurat clusters
sum_tf<-lapply(sum_tf,function(x) apply(x,2,mean)) #take average across group
sum_tf<-do.call("rbind",sum_tf) #condense to smaller data frame

sum_tf<-t(scale(sum_tf))
sum_tf<-sum_tf[,!endsWith(colnames(sum_tf),"NA")] #remove NA (doublet cells)

#clustered by all marker genes ga
sum_da_dend<-readRDS(file="hg38.geneactivity.dend.rds") 


sum_tf<-sum_tf[row.names(sum_tf) %in% unique(da_tf[da_tf$label!="",]$da_region),]
row.names(sum_tf)<-da_tf[match(row.names(sum_tf),da_tf$da_region,nomatch=0),]$tf_name
annot<-hg38_atac@meta.data[,c("celltype","seurat_clusters","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$seurat_clusters),]
annot<-annot[annot$seurat_clusters %in% colnames(sum_tf),]
annot<-annot[match(colnames(sum_tf),annot$seurat_clusters),]
sum_tf_plot<-t(sum_tf)

annot_clus_col<-annot[!duplicated(annot$seurat_clusters),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters)))
                        ))

bottom_ha<-columnAnnotation(foo = anno_mark(at = 1:ncol(sum_tf_plot), labels = colnames(sum_tf_plot)))

colfun=colorRamp2(quantile(unlist(sum_tf_plot), probs=c(0.5,0.90,0.95)),cividis(3))
plt1<-Heatmap(sum_tf_plot,
    cluster_rows=sum_da_dend,
    left_annotation=side_ha,
    col=colfun,
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

plt1<-draw(plt1)


pdf("hg38.tf.heatmap.pdf",height=20,width=20)
draw(plt1)
dev.off()
#system("slack -F hg38.tf.heatmap.pdf ryan_todo")

#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)

motif_order<-names(hg38_atac@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(hg38_atac@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = hg38_atac,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())

ggsave(plt,file="hg38.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
#system("slack -F hg38.tf.heatmap.motif.pdf ryan_todo")


```

```R
library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
library(ggdendro)
library(dendextend)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(circlize)
library(viridis)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
mm10_atac@meta.data$celltype<-"NA"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="0",]$celltype<-"Gran"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="1",]$celltype<-"iN"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="2",]$celltype<-"ExN"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="3",]$celltype<-"Astro"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="4",]$celltype<-"Oligo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="5",]$celltype<-"Oligo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="6",]$celltype<-"Endo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="7",]$celltype<-"Endo"

mm10_atac$celltype_col<-"NA"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="iN",]$celltype_col<-"#cc3366"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="ExN",]$celltype_col<-"#3399cc"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Oligo",]$celltype_col<-"#66cc99"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Astro",]$celltype_col<-"#99cc99"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Endo",]$celltype_col<-"#A52A2A"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Gran",]$celltype_col<-"#4B0082"

saveRDS(mm10_atac,"mm10_SeuratObject.PF.Rds")


#Clusters to test
cluster_to_test<-unique(mm10_atac$seurat_clusters)
#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group,assay.="GeneActivity",latent.vars.="nCount_GeneActivity"){
    da_ga_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = latent.vars.,
        only.pos=T,
        assay=assay.
        )
    da_ga_tmp$da_region<-row.names(da_ga_tmp)
    da_ga_tmp$enriched_group<-c(i)
    da_ga_tmp$compared_group<-c("all_other_cells")
    return(da_ga_tmp)
  }

da_ga<-list() #set up an empty list for looping through

n.cores=10 #Perform parallel application of DA test
da_ga<-mclapply(
    cluster_to_test,
    FUN=da_one_v_rest,
    obj=mm10_atac,
    group="seurat_clusters",
    assay.="GeneActivity",
    mc.cores=n.cores)

da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
write.table(da_ga_df,file="mm10.onevrest.da_ga.txt",sep="\t",col.names=T,row.names=T,quote=F)


da_ga<-list() #set up an empty list for looping through

n.cores=5 #Perform parallel application of DA test
da_ga<-mclapply(
    cluster_to_test,
    FUN=da_one_v_rest,
    obj=mm10_atac,
    group="seurat_cluster",
    assay.="chromvar", latent.vars.="nCount_peaks",
    mc.cores=n.cores)

da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
#To convert JASPAR ID TO TF NAME
da_ga_df<-da_ga_df[startsWith(row.names(da_ga_df),"MA"),]
da_ga_df$tf_name <- unlist(lapply(unlist(lapply(da_ga_df$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
write.table(da_ga_df,file="mm10.onevrest.da_chromvar.txt",sep="\t",col.names=T,row.names=T,quote=F)


#Plot out top ga for each cluster
da_ga<-read.csv(file="mm10.onevrest.da_ga.txt",head=T,sep="\t",row.names=NULL)
da_ga$gene_name<-da_ga$da_region
da_ga<-da_ga[complete.cases(da_ga),]
da_ga<-da_ga[!endsWith(da_ga$enriched_group,"NA"),]

da_ga$label<-""
for (x in unique(da_ga$enriched_group)){
selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% dplyr::arrange(p_val_adj) %>% dplyr::slice(1:3))$da_region
da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_ga<-as.data.frame(t(as.data.frame(mm10_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,mm10_atac$seurat_clusters) #group by rows to seurat clusters
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame

sum_ga<-t(scale(sum_ga))
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

#cluster by all marker genes
sum_da_dend <- t(sum_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="mm10.geneactivity.dend.rds") 


sum_ga<-sum_ga[row.names(sum_ga) %in% unique(da_ga$label),]

annot<-mm10_atac@meta.data[,c("celltype","cluster_col","seurat_clusters","celltype_col")]
annot<-annot[!is.na(annot$cluster_col),]
annot<-annot[!duplicated(annot$seurat_clusters),]
annot<-annot[annot$seurat_clusters %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$seurat_clusters),]
sum_ga_plot<-t(sum_ga)

annot_clus_col<-annot[!duplicated(annot$seurat_clusters),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters)))
                        ))

bottom_ha<-columnAnnotation(foo = anno_mark(at = 1:ncol(sum_ga_plot), labels = colnames(sum_ga_plot)))

colfun=colorRamp2(quantile(unlist(sum_ga_plot), probs=c(0.5,0.90,0.95)),magma(3))
plt1<-Heatmap(sum_ga_plot,
    cluster_rows=sum_da_dend,
    left_annotation=side_ha,
    col=colfun,
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

pdf("mm10.geneactivity.heatmap.pdf",height=4,width=20)
plt1
dev.off()
#system("slack -F mm10.geneactivity.heatmap.pdf ryan_todo")
############################################################
#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)############
#genes from brain map and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7918299/
sum_da_dend<-readRDS(file="mm10.geneactivity.dend.rds")
dat_ga<-as.data.frame(t(as.data.frame(mm10_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,mm10_atac$seurat_clusters) #group by rows to seurat clusters 
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame
sum_ga<-t(sum_ga)
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
      sep="", collapse=" ")
}

markers_limited<-as.data.frame(rbind(
    c("Mark","L1-6","GABAergic","GAD1"),
    c("Mark","L1-6","GABAergic","LHX6"),
    c("Mark","L1-6","GABAergic","ADARB2"),
    c("Mark","L1-6","GABAergic","LAMP5"),
    c("Mark","L1-6","GABAergic","VIP"),
    c("Mark","L1-6","GABAergic","SVIL"),
    c("Mark","L1-6","GABAergic","NXPH2"),
    c("Mark","L1-6","GABAergic","SST"),
    c("Mark","L1-6","Glutamatergic","SLC17A7"),
    c("Mark","L1-6","Glutamatergic","HTR1F"),
    c("Mark","L1-6","Glutamatergic","FEZF2"),
    c("Mark","L1-6","Glutamatergic","RORB"),
    c("Mark","L1-6","Glutamatergic","PCP4"),
    c("Mark","L1-6","Glutamatergic","POU3F2"),
    c("Mark","L1-6","Glutamatergic","CUX2"),
    c("Mark","L1-6","Glutamatergic","CALB1"),
    c("Mark","L1-6","Glutamatergic","RASGRF2"),
    c("Mark","Cerebellar","Granule","GABRA6"),
    c("Mark","Cerebellar","Granule","ZIC1"),
    c("Mark","Cerebellar","Granule","NEUROD1"),
    c("Mark","Cerebellar","Granule","ETV1"),
    c("Mark","Cerebellar","Granule","NFIA"),
    c("Mark","L1-6","Astrocyte","AQP4"),
    c("Mark","L1-6","Astrocyte","ALDH1L1"),
    c("Mark","L1-6","Astrocyte","GFAP"),
    c("Mark","L1-6","Astrocyte","SLC1A2"),
    c("Mark","L1-6","OPC","PDGFRA"),
    c("Mark","L1-6","OPC","CSPG4"),
    c("Mark","L1-6","Oligodendrocyte","OPALIN"),
    c("Mark","L1-6","Oligodendrocyte","PROX1"),
    c("Mark","L1-6","Oligodendrocyte","OLIG1"),
    c("Mark","L1-6","Oligodendrocyte","MOBP"),
    c("Mark","L1-6","Microglia","FYB"),
    c("Mark","L1-6","Microglia","C1QA"),
    c("Mark","L1-6","Microglia","C1QC"),
    c("Mark","L1-6","Endothelia","FLT1"),
    c("Mark","L1-6","Endothelia","KDR")

    ))

colnames(markers_limited)<-c("celltype","layer","marker_1","marker_2")

markers_limited$marker_2<-unlist(lapply(markers_limited$marker_2,simpleCap)) # correct case
sum_ga_sub<-sum_ga[match(markers_limited$marker_2,row.names(sum_ga),nomatch=0),]
markers_limited<-markers_limited[match(row.names(sum_ga_sub),markers_limited$marker_2,nomatch=0),]

sum_ga_sub<-t(scale(t(sum_ga_sub),center=F,scale=T))


annot<-mm10_atac@meta.data[,c("celltype","cluster_col","seurat_clusters","celltype_col")]
annot<-annot[!(annot$cluster_col=="NA"),]
annot<-annot[!duplicated(annot$seurat_clusters),]
annot<-annot[annot$seurat_clusters %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$seurat_clusters,nomatch=0),]
annot_clus_col<-annot[!duplicated(annot$seurat_clusters),]
annot_clus_col<-annot_clus_col[complete.cases(annot_clus_col),]
sum_ga_sub<-sum_ga_sub[,colnames(sum_ga_sub) %in% annot$seurat_clusters]

celltype=setNames(unique(annot$celltype_col),unique(annot$celltype))
cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters)))
celltype<-celltype[!is.na(celltype)]
cluster<-cluster[!is.na(cluster)]
cluster<-cluster[!is.na(names(cluster))]
top_ha<-columnAnnotation(df= data.frame(cluster=annot$seurat_clusters, celltype=annot$celltype), 
                col=list(celltype=celltype,cluster=cluster),
                    show_legend = c(TRUE, TRUE))

colfun=colorRamp2(quantile(unlist(sum_ga_sub), probs=c(0.5,0.90,0.95)),magma(3))

plt1<-Heatmap(sum_ga_sub,
    clustering_distance_columns="euclidean",
    clustering_distance_rows="euclidean",
    row_order=1:nrow(sum_ga_sub),
    column_split=annot$celltype,
    row_split=markers_limited$marker_1,
    cluster_row_slices=F,
    top_annotation=top_ha,
    left_annotation=rowAnnotation(df=data.frame(celltype_mark=markers_limited$marker_1)),
    #left_annotation=rowAnnotation(df=data.frame(celltype_mark=markers$celltype,layers_mark=markers$layer)),
    col=colfun,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp=gpar(fontsize=3)
)


pdf("mm10.geneactivity.markers.heatmap.pdf",height=5)
plt1
dev.off()
#system("slack -F mm10.geneactivity.markers.heatmap.pdf ryan_todo")



########################Plot out top TF for each cluster###################
da_tf<-read.csv(file="mm10.onevrest.da_chromvar.txt",head=T,sep="\t",row.names=NULL)
da_tf$gene_name<-da_tf$da_region
da_tf<-da_tf[complete.cases(da_tf),]
da_tf<-da_tf[!endsWith(da_tf$enriched_group,"NA"),]


da_tf$label<-""
for (x in unique(da_tf$enriched_group)){
selc_genes<-as.data.frame(da_tf %>% filter(enriched_group==x) %>% dplyr::arrange(p_val_adj) %>% dplyr::slice(1:3))$tf_name
da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$label<- da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$tf_name
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_tf<-as.data.frame(t(as.data.frame(mm10_atac[["chromvar"]]@data)))
sum_tf<-split(dat_tf,mm10_atac$cluster_ID) #group by rows to seurat clusters
sum_tf<-lapply(sum_tf,function(x) apply(x,2,mean)) #take average across group
sum_tf<-do.call("rbind",sum_tf) #condense to smaller data frame

sum_tf<-t(scale(sum_tf))
sum_tf<-sum_tf[,!endsWith(colnames(sum_tf),"NA")] #remove NA (doublet cells)

#clustered by all marker genes ga
sum_da_dend<-readRDS(file="mm10.geneactivity.dend.rds") 


sum_tf<-sum_tf[row.names(sum_tf) %in% unique(da_tf[da_tf$label!="",]$da_region),]
row.names(sum_tf)<-da_tf[match(row.names(sum_tf),da_tf$da_region,nomatch=0),]$tf_name

annot<-mm10_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(sum_tf),]
annot<-annot[match(colnames(sum_tf),annot$cluster_ID),]
sum_tf_plot<-t(sum_tf)

annot_clus_col<-annot[!duplicated(annot$cluster_ID),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
                        ))

#bottom_ha<-columnAnnotation(foo = anno_mark(at = 1:ncol(sum_tf_plot), labels = colnames(sum_tf_plot)))

colfun=colorRamp2(quantile(unlist(sum_tf_plot), probs=c(0.5,0.90,0.95)),cividis(3))
plt1<-Heatmap(sum_tf_plot,
    cluster_rows=sum_da_dend,
    left_annotation=side_ha,
    col=colfun,
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

plt1<-draw(plt1)


pdf("mm10.tf.heatmap.pdf",height=20,width=20)
draw(plt1)
dev.off()
#system("slack -F mm10.tf.heatmap.pdf ryan_todo")

#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)

motif_order<-names(mm10_atac@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(mm10_atac@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = mm10_atac,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())

ggsave(plt,file="mm10.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
#system("slack -F mm10.tf.heatmap.motif.pdf ryan_todo")

```

### Marker Plotting for cell type refinement


```R
library(Seurat)
library(Signac)
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")

marker_list<-NULL
marker_list[["Oligodendrocyte"]]<-c('Mobp','Cldn1','Prox1','Olig1')
marker_list[["Polydendrocyte"]]<-c('Cspg4','Pdgfra')
marker_list[["Astrocyte"]]<-c('Gfap','Glul','Slc1a2','Agt')
marker_list[["Ex_Neuron"]]<-c('Satb2','Cux2','Rorb','Pcp4','Foxp2','Col19a1','Slc17a7')
marker_list[["Microglia"]]<-c('C1qa','C1qc','Cxcr1')
marker_list[["Inhib_Neuron"]]<-c('Cnr1','Bcl11b','Gad1','Dlx1','Dlx2')
saveRDS(marker_list,"grosscelltype_markerlist.rds")

#Function to plot cell type markers
#For now just using the gross celltype marker list

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(parallel)

region_check<-function(i,object=hg38_atac){
    return(nrow(as.data.frame(LookupGeneCoords(object=object,gene=toupper(gene_list[i])))))
    }

marker_plot<-function(j,k=celltype_name,l=hg38_atac,m="hg38_",n="seurat_clusters"){
    plt<-CoveragePlot(
        object = l,
        region = j,
        group.by=n,
        extend.upstream = 1000,
        extend.downstream = 1000,
        ncol = 1
    )
    pdf(paste0("./marker_sets/",m,k,"_",j,"_genebody_accessibility.pdf")) 
    print(plt)
    dev.off()
}

dir.create("marker_sets")

hg38_atac<-readRDS(file="hg38_SeuratObject.PF.Rds")
marker_list<-readRDS("grosscelltype_markerlist.rds")


for (x in 1:length(marker_list)){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    gene_list<-toupper(gene_list)
    celltype_name<-names(marker_list)[x]
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,,n="peaks_snn_res.0.01",mc.cores=5)
}

for (i in list.files(path="./marker_sets",pattern="markerset_hg38")){
    system(paste0("slack -F ./marker_sets/",i," ryan_todo"))
}

#Now setting up for mm10 
region_check<-function(i,object=mm10_atac){
    return(nrow(as.data.frame(LookupGeneCoords(object=object,gene=gene_list[i]))))
    }

marker_plot<-function(j,k=celltype_name,l=mm10_atac,m="mm10_",n="seurat_clusters"){
    plt<-CoveragePlot(
        object = l,
        region = j,
        group.by=n,
        extend.upstream = 1000,
        extend.downstream = 1000,
        ncol = 1,
        tile=T,
        tile.size=500,
        tile.cells=downsamp
    )
    pdf(paste0("./marker_sets/",m,k,"_",j,"_genebody_accessibility.pdf"))
    print(plt)
    dev.off()
}

mm10_atac<-readRDS(file="mm10_SeuratObject.PF.Rds")
marker_list<-readRDS("grosscelltype_markerlist.rds")
summary(mm10_atac@meta.data$seurat_clusters) #setting tile cells to 99 to downsample
downsamp=30

for (x in 1:length(marker_list)){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    celltype_name<-names(marker_list)[x]
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,mc.cores=20) 
}

#For concatenating cell types into a single scrollable pdf
celltype=`ls mm10*genebody_accessibility.pdf | awk '{split($1,a,"_");print a[2]}' - | uniq`
for i in $celltype ; do convert `echo mm10_${i}_*genebody_accessibility.pdf` markerset_mm10_${i}.pdf; done

celltype=`ls hg38*genebody_accessibility.pdf | awk '{split($1,a,"_");print a[2]}' - | uniq`
for i in $celltype ; do convert `echo hg38_${i}_*genebody_accessibility.pdf` markerset_hg38_${i}.pdf; done

```

### Updating Marker List with Higher Resolution

Using Allen Brain Span dendrogram to determine cell types in mouse data.
http://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.RData

### Plotting of subcluster marker sets


```R
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(parallel)


region_check<-function(i,object=hg38_atac){
    return(nrow(as.data.frame(LookupGeneCoords(object=object,gene=toupper(gene_list[i])))))
    }

marker_plot<-function(j,k=celltype_name,l=hg38_atac,m="hg38_",n="seurat_clusters",outdir){
    plt<-CoveragePlot(
        object = l,
        region = j,
        group.by=n,
        extend.upstream = 1000,
        extend.downstream = 1000,
        ncol = 1
    )
    pdf(paste0(outdir,m,k,"_",j,"_genebody_accessibility.pdf"))
    print(plt)
    dev.off()
}

#Subcell type labelling of hg38 subcluster 4 (excitatory neurons)
dir.create("./subcluster/hg38_4_markers")
hg38_atac<-readRDS(file="./subcluster/_hg38_4_SeuratObject.Rds")
marker_list<-readRDS("hg38_markerlist.rds")
outdirec="./subcluster/hg38_4_markers/"
#Subset marker list to excitatory neuron subclusters
for (x in c("L5ET","L6CT","L56NP","L4IT","L6b","L56ITCAR3")){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    gene_list<-toupper(gene_list)
    celltype_name<-x
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,outdir=outdirec,n="peaks_snn_res.0.2",mc.cores=5)
}


#Subcell type labelling of hg38 subcluster 6 (inhibitory neurons)
dir.create("./subcluster/hg38_6_markers")
hg38_atac<-readRDS(file="./subcluster/_hg38_6_SeuratObject.Rds")
marker_list<-readRDS("hg38_markerlist.rds")
outdirec="./subcluster/hg38_6_markers/"
#Subset marker list to inhibitory neuron subclusters
for (x in c("VLMC","PVALB","VIP","PAX6","LAMP5","SST")){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    gene_list<-toupper(gene_list)
    celltype_name<-names(marker_list)[x]
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,outdir=outdirec,n="peaks_snn_res.0.2",mc.cores=20)
}


celltype=`ls hg38*genebody_accessibility.pdf | awk '{split($1,a,"_");print a[2]}' - | uniq`
for i in $celltype ; do convert `echo hg38_${i}_*genebody_accessibility.pdf` markerset_hg38_${i}.pdf; done



```

## Identify Marker Genes with specificity to subclusters
Based on seurat tutorial https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1

```R
library(Signac)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(seriation)
library(viridis)
library(circlize)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(ggdendro)
library(dendextend)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")


#Grab top overlapping TFs
topTFs <- function(markers_list,celltype, padj.cutoff = 1e-2,ga=NA,motifs=NA) {
    ctmarkers_ga<- dplyr::filter(
      ga, GeneActivity.group == celltype) %>% 
      arrange(-GeneActivity.auc)
    ctmarkers_ga$gene<-toupper(ctmarkers_ga$gene)

    if(is.data.frame(motifs)){
        ctmarkers_motif <- dplyr::filter(
      motifs, chromvar.group == celltype) %>% 
      arrange(-chromvar.auc)}

    if(is.data.frame(motifs) && is.data.frame(ga)){    
      top_tfs <- inner_join(
        x = ctmarkers_ga [,c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )} else {
        top_tfs<-ctmarkers_ga [,c(2, 11, 6, 7)]
      }

  auc_colnames<-grep(".auc$",colnames(top_tfs))
  top_tfs$avg_auc <-  rowMeans(top_tfs[auc_colnames])
  top_tfs <- arrange(top_tfs, -avg_auc)
  top_tfs$celltype<-celltype
  return(top_tfs)
}

#Identify top markers
Identify_Marker_TFs<-function(x,group_by.="predicted.id",assay.="RNA",prefix){
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by., assay = 'data', seurat_assay = assay.)
    colnames(markers) <- paste(assay., colnames(markers),sep=".")
    if (assay. == "chromvar") {
      motif.names <- markers[,paste0(assay.,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="peaks")
    } else {
    markers$gene <- markers[,paste0(assay.,".feature")]
    }
    return(markers) 
}

#Average markers across groups
average_features<-function(x=hg38_atac,features=da_tf_markers$motif.feature,assay="chromvar",group_by.="predicted.id"){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    dat_motif<-x[[assay]]@data[features,]
    dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
    sum_motif<-split(dat_motif,x@meta.data[,group_by.]) #group by rows to seurat clusters
    sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take average across group
    sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame

    sum_motif<-t(scale(sum_motif))
    sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

#Make a heatmap of aligned multiple modalities
plot_top_TFs<-function(x=stromal,tf_markers=da_tf_markers,prefix="stromal",group_by.="cluster_ID",CHROMVAR=TRUE,GA=TRUE){
    if(CHROMVAR){
    tf_motif<-average_features(x=x,features=tf_markers$chromvar.feature,assay="chromvar",group_by.=group_by.)
    tf_motif<-tf_motif[row.names(tf_motif) %in% tf_markers$chromvar.feature,]
    row.names(tf_motif)<-tf_markers[tf_markers$chromvar.feature %in% row.names(tf_motif),]$gene
    }
    if(GA){
    row.names(x[["GeneActivity"]]@data)<-toupper(row.names(x[["GeneActivity"]]@data))#change geneactivity names to upper
    tf_ga<-average_features(x=x,features=tf_markers$gene,assay="GeneActivity",group_by.=group_by.)
    tf_ga<-tf_ga[row.names(tf_ga) %in% tf_markers$gene,]
    }
    if(GA && CHROMVAR){
    markers_list<-Reduce(intersect, list(row.names(tf_motif),row.names(tf_ga)))
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
    }

    #dend_col <- t(tf_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  

    #set up heatmap seriation and order by RNA
    #o = seriate(max(tf_ga) - tf_ga, method = "BEA_TSP")
    o = seriate(tf_ga, method = "Heatmap")
    saveRDS(o,file=paste0(prefix,".geneactivity.dend.rds")) 
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o,1),]$GeneActivity.auc)

    if(CHROMVAR){
    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.80,0.95)),cividis(3))
    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="peaks",motifs = motif_list[get_order(o,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

    }
    if(GA){
    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.80,0.95)),magma(3))
    }

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_ga)), labels =row.names(tf_ga),labels_gp=gpar(fontsize=6)))

    annot<-x@meta.data[,c("celltype","cluster_col","seurat_clusters","celltype_col")]
    annot<-annot[!duplicated(annot$seurat_clusters),]
    annot<-annot[annot$seurat_clusters%in% colnames(tf_ga),]
    annot<-annot[match(colnames(tf_ga),annot$seurat_clusters),]

    col_ha<-columnAnnotation(df= data.frame(cluster=annot$seurat_clusters,celltype=annot$celltype),
                    col=list(
                        cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                        celltype=setNames(unique(annot$celltype_col),unique(as.character(annot$celltype)))
                            ))

  if(GA==TRUE ){
      ga_auc<-Heatmap(side_ha_ga,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          top_annotation=col_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)
  }
  if(CHROMVAR==TRUE ){
      motif_auc<-Heatmap(side_ha_motif,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          top_annotation=col_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
  }
  if(GA==TRUE &&(CHROMVAR==FALSE)){
          ga_plot<-Heatmap(tf_ga,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          top_annotation=col_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
        right_annotation=gene_ha)

  }

  if(all(CHROMVAR,GA)){
      plt1<-draw(ga_auc+ga_plot+motif_auc+motif_plot)
  } else if(CHROMVAR){
      plt1<-draw(motif_auc+motif_plot)
  } else {
      plt1<-draw(ga_auc+ga_plot)
  }


    pdf(paste0(prefix,".tf.heatmap.pdf"),height=20)
    print(plt1)
    dev.off()

    system(paste0("slack -F ",prefix,".tf.heatmap.pdf ryan_todo"))
    if(CHROMVAR){
    system(paste0("slack -F ",prefix,".tf.heatmap.motif.pdf ryan_todo"))
    }
}

#Final wrapper function
run_top_TFs<-function(obj=hg38_atac,prefix="hg38_atac",i="cluster_ID",marker_number=3,CHROMVAR.=TRUE){
    if(CHROMVAR.==TRUE){
        markers<-lapply(c("GeneActivity","chromvar"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay))
        names(markers)<-c("GeneActivity","chromvar")
        clusters<-unique(obj@meta.data[,i])
        markers_out<-do.call("rbind",lapply(clusters, function(j) head(topTFs(markers_list=markers,celltype=j,ga=markers$GeneActivity,motifs=markers$chromvar),n=marker_number))) #grab top 5 TF markers per celltype
    } else {
        markers<-lapply(c("GeneActivity"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay))
        names(markers)<-c("GeneActivity")
        clusters<-unique(obj@meta.data[,i])
        markers_out<-do.call("rbind",lapply(clusters, function(j) head(topTFs(markers_list=markers,celltype=j,ga=markers$GeneActivity),n=marker_number))) #grab top 5 TF markers per celltype
    }
    dim(markers_out)
    markers_out<-markers_out[!duplicated(markers_out$gene),]
    dim(markers_out)
    saveRDS(markers_out,file=paste0(prefix,"_celltype_TF_markers.RDS"))
    da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
    plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i,CHROMVAR=CHROMVAR.,GA=TRUE)
}


#hg38 TF markers
hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
hg38_atac<-subset(hg38_atac,subcluster_x!="NA")
run_top_TFs(obj=hg38_atac,prefix="hg38_TF_revised",i="seurat_clusters",marker_number=5)

#Marker Genes per celltype
#for(j in unique(hg38_atac$celltype)){
#    hg38_sub<-subset(hg38_atac,celltype==j)
#    if(length(unique(hg38_sub$cluster_ID))>1){
#    run_top_TFs(obj=hg38_sub,prefix=paste0("hg38_markergenes_",j),i="cluster_ID",marker_number=25,CHROMVAR.=FALSE)
#    }
#}

#mm10 markers
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
mm10_atac<-subset(mm10_atac,subcluster_x!="NA")
run_top_TFs(obj=mm10_atac,prefix="mm10_TF_revised",i="seurat_clusters",marker_number=3)

#Marker Genes per celltype
#for(j in unique(mm10_atac$celltype)){
#    mm10_sub<-subset(mm10_atac,celltype==j)
#    if(length(unique(mm10_sub$cluster_ID))>1){
#    run_top_TFs(obj=mm10_sub,prefix=paste0("mm10_markergenes_",j),i="cluster_ID",marker_number=10,CHROMVAR.=FALSE)
#    }
#}



```


## Plotting with Marker Genes

```R
    library(Signac)
    library(Seurat)
    library(SeuratWrappers)
    library(ggplot2)
    library(patchwork)
    library(cicero)
    library(cisTopic)
    library(GenomeInfoDb)
    set.seed(1234)
    library(Matrix)
    library(dplyr)
    library(ggrepel)
    library(RColorBrewer)
    setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

plot_markers<-function(obj=mm10_atac,gene_name="Gad1"){
    feat_data<-setNames(obj[["GeneActivity"]]@data[gene_name,],nm=colnames(obj[["GeneActivity"]]@data))
    seurat_x<-setNames(colnames(obj[["GeneActivity"]]@data),obj[["GeneActivity"]]@data[gene_name,])
    obj<-AddMetaData(obj,feat_data,col.name="FEAT")
    obj<-AddMetaData(obj,obj@reductions$umap@cell.embeddings,col.name=c("umap_x","umap_y"))
    dat<-as.data.frame(obj@meta.data)

    dat$subcluster_x<-as.numeric(dat$subcluster_x)
    dat$subcluster_y<-as.numeric(dat$subcluster_y)

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=FEAT))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+ggtitle(gene_name)+
    ggtitle(gene_name) +scale_color_gradient(low="white",high="red",na.value="red",limits=c(0,quantile(dat$FEAT,0.99)))+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks = element_blank(),legend.position = "bottom")

    plt_list<-ggplot(dat,aes(x=subcluster_x,y=subcluster_y,color=FEAT))+
    geom_point(alpha=0.05,size=0.5,shape=16)+ theme_bw()+ scale_color_gradient(low="white",high="red",na.value="red",limits=c(0,quantile(dat$FEAT,0.9)))+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
    facet_wrap(facets=vars(seurat_clusters),ncol=2) + coord_cartesian(clip = "off") 
    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 
    ggsave(plt,file="feat.png")
    #system("slack -F feat.png ryan_todo")   
}

mm10_atac<-readRDS(file="mm10_SeuratObject.PF.Rds")
gene_list<-c("Gad1","Gad2","Slc32a1","Slc17a7","Lamp5","Ndnf","Sncg","Vip","Sst","Chodl","Pvalb","Cux2","Rorb","Fezf2","Sulf1","Fam84b","Sla2","Foxp2","Nxph4","Aqp4","Mbp","Cldn5","Ctss","C1qa") #from https://celltypes.brain-map.org/rnaseq/mouse_ctx-hpf_10x?selectedVisualization=Heatmap&colorByFeature=Cell+Type&colorByFeatureValue=Gad1
gene_list<-c("Gabra6") #evan macosko kozareva mouse cerebellar cortex https://www.nature.com/articles/s41586-021-03220-z
gene_list<-gene_list[unlist(lapply(gene_list,function(x) x %in% row.names(mm10_atac[["GeneActivity"]]@data)))]
lapply(gene_list,function(x) plot_markers(gene_name=x))

hg38_atac<-readRDS(file="hg38_SeuratObject.PF.Rds")
gene_list<-c('RORB',
'LINC00507',
'GLRA3',
'PTPN3',
'CCDC68',
'LINC01202',
'OTOGL',
'TNNT2',
'LNX2',
'THEMIS',
'TNFAIP6',
'FEZF2',
'KLK7',
'OPALIN',
'FTH1P3',
'FGFR3',
'AQP1',
'LAMP5',
'PAX6',
'PVALB',
'SST',
'GGTLC3',
'COL15A1',
'VIP',
'EXPH5',
'SMOC1',
'IGDCC3',
'BMP2',
'CRABP1',
'AARD',
'TYROBP',
'CD74',
'PDGFRA',
'COL20A1',
'NES',
'FREM2',
'PDGFRA')

#https://celltypes.brain-map.org/rnaseq/human_m1_10x?selectedVisualization=Heatmap&colorByFeature=Gene+Expression&colorByFeatureValue=GAD1
gene_list<-gene_list[unlist(lapply(gene_list,function(x) x %in% row.names(hg38_atac[["GeneActivity"]]@data)))]
lapply(gene_list,function(x) plot_markers(obj=hg38_atac,gene_name=x))
```

## Output final metadata tables
```R
library(Signac)
library(Seurat)
library(patchwork)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")

write.table(as.data.frame(mm10_atac@meta.data),file="mm10_PF.metadata.tsv",sep="\t",quote=F,col.names=T,row.names=T)
#system("slack -F mm10_PF.metadata.tsv ryan_todo")

write.table(as.data.frame(hg38_atac@meta.data),file="hg38_PF.metadata.tsv",sep="\t",quote=F,col.names=T,row.names=T)
#system("slack -F hg38_PF.metadata.tsv ryan_todo")
```

## Output Tab separated 3D clustering for Blender Plot
```R
library(Signac)
library(Seurat)
library(patchwork)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")

hg38_atac <- RunUMAP(object = hg38_atac, reduction = 'cistopic', reduction.name="umap3d",dims = 1:ncol(hg38_atac@reductions$cistopic),n.components=3)
mm10_atac <- RunUMAP(object = mm10_atac, reduction = 'cistopic', reduction.name="umap3d",dims = 1:ncol(mm10_atac@reductions$cistopic),n.components=3)

hg38_out<-cbind(hg38_atac$celltype,
    colnames(hg38_atac),
    as.data.frame(hg38_atac@reductions$umap3d@cell.embeddings),
    hg38_atac$celltype_col)
write.table(hg38_out,"hg38_3d.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
#system("slack -F hg38_3d.umap.tsv ryan_todo")


mm10_out<-cbind(mm10_atac$celltype,
    colnames(mm10_atac),
    as.data.frame(mm10_atac@reductions$umap3d@cell.embeddings),
    mm10_atac$celltype_col)
write.table(mm10_out,"mm10_3d.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
#system("slack -F mm10_3d.umap.tsv ryan_todo")

###Repeat with subcluster coloring
hg38_out<-cbind(hg38_atac$cluster_ID,
    colnames(hg38_atac),
    as.data.frame(hg38_atac@reductions$umap3d@cell.embeddings),
    hg38_atac$subcluster_col)
write.table(hg38_out,"hg38_3d.subclus.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
#system("slack -F hg38_3d.subclus.umap.tsv ryan_todo")


###Repeat with subcluster coloring
mm10_out<-cbind(mm10_atac$cluster_ID,
    colnames(mm10_atac),
    as.data.frame(mm10_atac@reductions$umap3d@cell.embeddings),
    mm10_atac$subcluster_col)
write.table(mm10_out,"mm10_3d.subclus.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
#system("slack -F mm10_3d.subclus.umap.tsv ryan_todo")
```


## Save final seurat files for UCSC cellbrowser and set up.

To view data in an interactive way, I converted the final seurat object for UCSCs cell browser.

I then hosted locally to test.
Install cellbrowser via this (guide.)[https://cellbrowser.readthedocs.io/en/master/installation.html]

```R
#Downloaded seurat objects off clusters
#Prepare seurat data as cellbrowser file

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
library(Signac)
library(Seurat)
library(Matrix)
library(R.utils)
library(dplyr)


#modified function from https://github.com/satijalab/seurat-wrappers/blob/master/R/cellbrowser.R just for export
ExportToCellbrowser <- function(
object,
reduction="umap",
assay.name="GeneActivity",
dir,
dataset.name,
marker.file,
marker.n=5,
cluster.field="seurat_clusters",
gene_list,
skip.expr.mat=FALSE) {

  reducNames = reduction
  # see https://satijalab.org/seurat/essential_commands.html
    idents <- Idents(object = object)
    cellOrder <- row.names(object@meta.data)
    counts <-  object@assays[[assay.name]]@data
    counts<-counts[,cellOrder]
    genes <- rownames(x = counts)
    dr<-object@reductions[[reduction]]
    meta.fields<-object@meta.data
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)}

  # Export expression matrix
    # we have to write the matrix to an mtx file
    matrixPath <- file.path(dir, "matrix.mtx")
    genesPath <- file.path(dir, "features.tsv")
    barcodesPath <- file.path(dir, "barcodes.tsv")
if(!skip.expr.mat){
    message("Writing expression matrix to ", matrixPath)
    message("This may take a couple of minutes...")
    writeMM(counts, matrixPath)
    # easier to load if the genes file has at least two columns. Even though seurat objects
    write.table(as.data.frame(cbind(rownames(counts), rownames(counts))), file=genesPath, sep="\t", row.names=F, col.names=F, quote=F)
    write(colnames(counts), file = barcodesPath)

    message("Gzipping expression matrix")
    gzip(matrixPath,overwrite=T,remove=T)
    gzip(genesPath,overwrite=T,remove=T)
    gzip(barcodesPath,overwrite=T,remove=T)
}
    matrixOutPath <- "matrix.mtx.gz"

    #Export embeddings
    df <- dr@cell.embeddings
    if (ncol(x = df) > 2) {
      warning('Embedding ', embedding, ' has more than 2 coordinates, taking only the first 2')
      df <- df[, 1:2]
    }
    colnames(x = df) <- c("x", "y")
    df <- data.frame(cellId = rownames(x = df), df, check.names = FALSE)
    fname <- file.path(dir,"umap.coords.tsv")
    message("Writing embeddings to ", fname)
    write.table(df[cellOrder, ], sep="\t", file=fname, quote = FALSE, row.names = FALSE)
  embeddings.conf <- sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',reduction)
  coords.string <- paste0("[",paste(embeddings.conf, collapse = ",\n"),"]")

  # Export metadata
  enum.fields=colnames(object@meta.data)[colnames(object@meta.data) %in% c("cellid","pcr_idx","tn5_idx",cluster.field,"seurat_clusters","cluster_ID","predicted.id","celltype")]
  fname <- file.path(dir, "meta.tsv")
  message("Writing meta data to ", fname)
  out.meta<-as.matrix(object@meta.data[cellOrder, enum.fields])
  write.table(out.meta, sep = "\t", file = fname, quote = FALSE, row.names = FALSE)
  enum.string <- paste0("[",paste(paste0('"', enum.fields, '"'), collapse = ", "), "]")

  # Export markers
    file <- file.path("markers.tsv")
    fname <- file.path(dir, file)
      da_ga<-read.table(file=marker.file)
        da_ga$label<-""
        for (x in unique(da_ga$enriched_group)){
        selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% dplyr::arrange(p_val_adj) %>% dplyr::slice(1:marker.n))$da_region
        da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
        }
    da_ga<-da_ga[da_ga$label!="",]
      message("Writing top ", marker.n, " cluster markers to ", fname)
    da_ga$cluster<-unlist(lapply(strsplit(da_ga$enriched_group,"_"),"[",1))
      da_ga<-da_ga[c("cluster","label","p_val_adj","enriched_group")]
      colnames(da_ga)<-c("cluster","symbol","Adjusted p Value","Subcluster Enrichment")
      write.table(x = da_ga, file = fname, quote = FALSE, sep = "\t", col.names = T,row.names=F)
    markers.string <- sprintf('markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',file)

#Export Gene list
  file <- file.path("quickGenes.csv")
  fname <- file.path(dir, file)
  write.table(gene_list,sep=",",file = fname, quote = FALSE, row.names = FALSE)

  config <- '
# This is a bare-bones cellbrowser config file auto-generated from R.
# Look at https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf
# for a full file that shows all possible options
name="%s"
shortLabel="%1$s"
exprMatrix="%s"
tags = ["sciatac"]
meta="meta.tsv"
# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"
geneIdType="auto"
# file with gene,description (one per line) with highlighted genes, called "Dataset Genes" in the user interface
quickGenesFile="quickGenes.csv"
clusterField="%s"
labelField="%s"
enumFields=%s
%s
coords=%s
geneLabel="Gene Activity"
unit="GA"'

  config <- sprintf(
    config,
    dataset.name,
    matrixOutPath,
    cluster.field,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  confPath = file.path(dir, "cellbrowser.conf")
  message("Writing cellbrowser config to ", confPath)
  cat(config, file = confPath)
  message("Prepared cellbrowser directory ", dir)
}


hg38_atac<-readRDS(file="hg38_SeuratObject.PF.Rds")

marker_list<-readRDS("grosscelltype_markerlist.rds")
hg38_marker_list<-toupper(unname(unlist(marker_list)))

ExportToCellbrowser(
object=hg38_atac,
reduction="umap",
assay.name="GeneActivity",
dir="hg38_cb",
dataset.name="hg38",
marker.file="hg38.onevrest.da_ga.txt",
marker.n=5,
cluster.field="seurat_clusters",
gene_list=hg38_marker_list,
skip.expr.mat=FALSE)

mm10_atac<-readRDS(file="mm10_SeuratObject.PF.Rds")
mm10_marker_list<-unname(unlist(marker_list))

ExportToCellbrowser(
object=mm10_atac,
reduction="umap",
assay.name="GeneActivity",
dir="mm10_cb",
dataset.name="mm10",
marker.file="mm10.onevrest.da_ga.txt",
marker.n=5,
cluster.field="seurat_clusters",
gene_list=mm10_marker_list,
skip.expr.mat=TRUE)
```

Locally building and hosting the cell browser.
Using windows wsl2 ubuntu environment

```bash
#install cellbrowser via conda
conda install -c bioconda ucsc-cell-browser

cd /mnt/c/Users/mulqueen/Documents/sciDROP/hg38_cb
cbBuild -o hg38_cellbrowser -p 8888

cd /mnt/c/Users/mulqueen/Documents/sciDROP/mm10_cb
cbBuild -o mm10_cellbrowser -p 8888
#Go to http://localhost:8888/ to use the interactive viewer

```

## Single Cell ATAC Comparison Across Adult Mouse Brains

### Set up SRA-toolkit
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

```bash
cd /home/groups/CEDAR/mulqueen/src
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.0-mac64/bin
vdb-config --prefetch-to-cwd
```

### Downloading Raw scATAC data Comparisons using SRA-toolkit

sciATAC Mouse Brain
https://www.ncbi.nlm.nih.gov/sra?term=SRX9850743
Downloading the files which contains a standard genomic Read 1 and Read 2. Read 3 (technical read) has the read assigned cellID barcode.
Going to add read 3 into the fastq read name and pipe these through out standard scitools processing.

Note R3 is properly assigned, corrected cellIDs, so no fuzzy matching, or white-listed barcode concatenation is necessary.

```bash
sra="SRR13437232"
ref_dir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"
mkdir ${ref_dir}/sciATAC
prefetch $sra --type fastq -O ${ref_dir}/sciATAC/${sra}
fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 20 -m 50G -O ${ref_dir}/sciATAC/${sra} ${sra}
gzip ${ref_dir}/${sra}/*fastq &

```

Generating a script for bam setup and output. Note that this does hold all the barcode sequences in memory, which isn't ideal. 
barcode_sciatac.py
```python
import gzip
from Bio import SeqIO
import sys

def change_id(x):
    global i
    x.id=barc[i]+":"+str(i)
    i+=1
    return(x)

cwd="/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciATAC/SRR13437232"
fq1=sys.argv[1] #read argument
fq3=sys.argv[2] #index argument
barc=[str(record.seq) for record in SeqIO.parse(gzip.open(fq3,"rt"), "fastq")] #make list of barcodes
i=0
SeqIO.write((change_id(x=record) for record in SeqIO.parse(gzip.open(fq1,"rt"), 'fastq')), sys.stdout, "fastq")
```

Running the python script. Then compressing files, then performing alignment on the bam files.

```bash
python ./barcode_sciatac.py SRR13437232_1.fastq.gz SRR13437232_3.fastq.gz > SRR13437232_1.barc.fastq &
python ./barcode_sciatac.py SRR13437232_2.fastq.gz SRR13437232_3.fastq.gz > SRR13437232_2.barc.fastq &
```

sciMAP Mouse Brain
https://www.ncbi.nlm.nih.gov/sra?term=SRX9850744
```bash
sra="SRR13437233"
ref_dir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"
mkdir ${ref_dir}/sciMAP
prefetch $sra --type fastq -O ${ref_dir}/sciMAP/${sra}
fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 20 -m 50G -O ${ref_dir}/sciMAP/${sra} ${sra}
gzip ${ref_dir}/${sra}/*fastq &
````

Generating a script for bam setup and output, just as above.
barcode_scimap.py
```python
import gzip
from Bio import SeqIO
import sys

def change_id(x):
    global i
    x.id=barc[i]+"."+str(i)
    i+=1
    return(x)

fq1=sys.argv[1] #read argument
fq3=sys.argv[2] #index argument
barc=[str(record.seq) for record in SeqIO.parse(gzip.open(fq3,"rt"), "fastq")] #make list of barcodes
i=0
SeqIO.write((change_id(x=record) for record in SeqIO.parse(gzip.open(fq1,"rt"), 'fastq')), sys.stdout, "fastq")
```

Now running this python script for both fastq 1 and fastq 2.

```bash
python ./barcode_scimap.py SRR13437233_1.fastq.gz SRR13437233_3.fastq.gz > SRR13437233_1.barc.fastq &
python ./barcode_scimap.py SRR13437233_2.fastq.gz SRR13437233_3.fastq.gz > SRR13437233_2.barc.fastq &
```

snATAC 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2668124
```bash
sra="SRR6768122"
ref_dir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"
mkdir ${ref_dir}/snATAC
prefetch $sra --type fastq -X 50G -O ${ref_dir}/snATAC/${sra}
fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 20 -m 50G -O ${ref_dir}/snATAC/${sra}/${sra} ${ref_dir}/snATAC/${sra}/${sra}.sra
gzip ${ref_dir}/${sra}/*fastq &
#download list of cells passing filter
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2668nnn/GSM2668124/suppl/GSM2668124_p56.nchrM.merge.sel_cell.qc.txt.gz
gzip -d GSM2668124_p56.nchrM.merge.sel_cell.qc.txt.gz
#make filter cell list
awk '{print $1'} GSM2668124_p56.nchrM.merge.sel_cell.qc.txt > GSM2668124_filtcells.list.txt
```

Adjusting fastq readname format to jive with scitools.

Note the script changes relative to sci-platform. Indexes are already in the fastq file, just need to rearrange format. "PLease note that the barcode information was integrated into the read name."
barcode_snatac.py
```python
import gzip
from Bio import SeqIO
import sys

def change_id(x):
    barc=x.description.split(" ")[1].split(":")[0]
    number_out=x.description.split(" ")[0].split(".")[2]
    x.id=barc+"."+number_out
    return(x)

cwd="/home/groups/CEDAR/mulqueen/mouse_brain_ref/snATAC/SRR6768122/SRR6768122"
fq1=sys.argv[1] #read argument
#fq1="SRR6768122.sra_1.fastq.gz"
SeqIO.write((change_id(x=record) for record in SeqIO.parse(gzip.open(fq1,"rt"), 'fastq')), sys.stdout, "fastq")
```

Now running this python script for both fastq 1 and fastq 2.

```bash
python ./barcode_snatac.py SRR6768122.sra_1.fastq.gz > SRR6768122.sra_1.barc.fastq &
python ./barcode_snatac.py SRR6768122.sra_2.fastq.gz > SRR6768122.sra_2.barc.fastq &
```


ddscATAC
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123581
Use this to parse reads? https://github.com/buenrostrolab/dscATAC_analysis_code/blob/master/tag_validation/parse_oligos.py

```bash
mkdir ddscATAC
ref_dir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"
mkdir ${ref_dir}/ddscATAC
for sra in SRR8310661 SRR8310662 SRR8310663 SRR8310664 SRR8310665 SRR8310666 SRR8310667 SRR8310668; do
prefetch $sra --type fastq -X 50G -O ${ref_dir}/ddscATAC/${sra} &
done &

for sra in SRR8310661 SRR8310662 SRR8310663 SRR8310664 SRR8310665 SRR8310666 SRR8310667 SRR8310668; do
fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 20 -m 50G -O ${ref_dir}/ddscATAC/${sra}/${sra} ${ref_dir}/ddscATAC/${sra}/${sra}.sra &
done &

gzip -r ${ref_dir}/ddscATAC/*/*/*fastq
```

Using a custom script to parse through the oligoes. Following the barcode dscATAC_dsciATAC_bed_structures.xlsx logic file provided. Fastqs are already split by PCR barcode, so going to read into memory for speed.
barcode_dscatac.py
```python
import gzip
from Bio import SeqIO
import sys
from fuzzysearch import find_near_matches
from Bio.SeqIO.QualityIO import FastqGeneralIterator

idx_list=[]

me_seq="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
idx_list=['AACCACA','AACGGTG','AACGTAA','AACTCTT','AAGCAGC','AAGGTTC','AAGTGCG','AAGTTAT','AATGATT','AATGGCC','AATTCCA','AATTGGT','ACACGCG','ACAGCTT','ACAGTAC','ACCATGC','ACCGGCT','ACCTACC','ACGAGAA','ACGTATT','ACGTTGG','ACTACGA','ACTCAAT','AGACCAT','AGACTTC','AGAGACC','AGATAGG','AGCAACG','AGCCGCC','AGCGAAT','AGCTGAG','AGGACGT','AGGATAC','AGGCATG','AGTAAGC','AGTTTCT','ATAGGCA','ATAGTTG','ATATAAC','ATGAATA','ATGAGCT','ATGGTGT','ATGTTCC','ATTATTC','ATTCACG','ATTCGTT','ATTGCCT','CACATGA','CACGCCA','CACGGAC','CACTTCT','CAGAATT','CAGAGAG','CAGGCGG','CATACGC','CATCAGT','CATCTTA','CATGTAT','CCAAGCT','CCACTTG','CCAGTCA','CCATAAT','CCGAACC','CCGCGAT','CCGGTTT','CCTATGT','CCTCCTT','CCTTAGG','CGACACT','CGAGGTC','CGAGTGG','CGATGCA','CGATTAC','CGCAATC','CGCCTAA','CGCGCTT','CGCGGCG','CGGAGGA','CGGATCT','CGGCCAG','CGGCTGC','CGGTACG','CGTACAA','CGTAGCC','CGTGATA','CGTTTGA','CTAACTC','CTAAGAA','CTAGAGC','CTATTCG','CTCATTT','CTCTTGC','CTGCGCC','CTGGCAT','CTTACCG','CTTCATC','CTTGCGA','CTTGTCC','GAACCGT','GAATATG','GAATCAA','GACAATA','GAGAGGT','GAGCGTG','GAGCTAA','GAGGACA','GAGTTGC','GATAGAC','GATCACC','GCACAGC','GCAGTGT','GCCTCGT','GCCTTTG','GCGACTC','GCGCACG','GCGTAGA','GCTAATT','GCTCCAA','GCTTTAT','GGAAGTT','GGACGAC','GGAGCCT','GGCAGGC','GGCGGAA','GGCGTCC','GGTAACA','GGTCGTA','GGTGTTT','GGTTAGT','GGTTCAC','GTAATAC','GTCCTTC','GTCGGTT','GTGCATT','GTGGCGC','GTGGTAG','GTGTCCA','GTGTGTC','GTTAGGA','GTTGATG','TAACGCC','TAAGAGG','TAAGGTA','TACCGAA','TACGCAT','TACTTTC','TAGTACC','TAGTGTT','TATACTT','TATGTGC','TCAAGAC','TCAGCAA','TCATACA','TCCAGTT','TCCGCTC','TCCTGGC','TCCTTAA','TCGACAG','TCGCGCA','TCGGATG','TCGGCGT','TCGTTCT','TCTGAAC','TCTTGTA','TGAATCC','TGACCGC','TGAGATT','TGAGGAG','TGATTGT','TGCGAGC','TGCGTTG','TGCTACT','TGGAAGG','TGGACCA','TGGCAAC','TGGCCTT','TGGTGAA','TGTAGTG','TGTCGCT','TGTTTAG','TTAAGCG','TTACAGA','TTATCAT','TTCCTCT','TTCGTAC','TTCTGCA','TTGAGGC','TTGGACT','TTGGTTA','TTGTAAG','TTTCCTA','TTTGGTC'] 

constant_1="TATGCATGAC"
constant_2="AGTCACTGAG" #i think this is actually reversed in the sequence
linker1 = "CCTAGTCGCGTAGAC" #this is from the parse_oligoes.py, not sure if it is in this sequence though

#Function to find a single mismatch against list of given indexes
def one_mismatch_find(bc):
    matched=False
    matched_idx=[find_near_matches(bc, l, max_deletions=0, max_insertions=0, max_substitutions=1) for l in idx_list]
    idx_match_list=[idx for idx, x in enumerate(matched_idx) if len(x)>0]
    if len(idx_match_list)>0:
        matched_idx=idx_list[idx_match_list[0]]
        bc=matched_idx
        matched=True
    return([bc,matched])

#Function to see if read fits expected format, and change the read id to include the barcode
def change_id(title1,seq1,qual1):
    bc1_matched=False
    bc2_matched=False
    bc3_matched=False
    tn5_adapt_cut=find_near_matches(me_seq, seq1, max_deletions=0, max_insertions=0, max_substitutions=1)
    if len(tn5_adapt_cut)>0:
        if (len(str(seq1))-tn5_adapt_cut[0].end >= 14) and (tn5_adapt_cut[0].start >= (7*3)+len(constant_1)+len(constant_2)):
            out_seq=seq1[tn5_adapt_cut[0].start+len(me_seq):] #get output sequence and qual scores
            idx=str(seq1[:tn5_adapt_cut[0].start]) #pull out index as string
            bc3=idx[len(idx)-7:len(idx)+1] #get 7bp bc3
            idx=idx[:len(idx)-len(bc3)-len(constant_2)] #trim
            bc2=idx[len(idx)-7:len(idx)+1] #get 7bp bc2
            bc1=idx[:7] #pull bc1 from the front since there is variable phasing bases between
            if bc1 in idx_list:
                bc1_matched=True
            else:
                bc=one_mismatch_find(bc1)
                bc1_matched=bc[1]
                bc1=bc[0]
            if bc2 in idx_list:
                bc2_matched=True
            else:
                bc=one_mismatch_find(bc2)
                bc2_matched=bc[1]
                bc2=bc[0]
            if bc3 in idx_list:
                bc3_matched=True
            else:
                bc=one_mismatch_find(bc3)
                bc3_matched=bc[1]
                bc3=bc[0]
            if bc1_matched and bc2_matched and bc3_matched:
                number_out=title1.split(" ")[0].split(".")[2]
                barc=bc1+bc2+bc3
                seq1=seq1[tn5_adapt_cut[0].end:]
                qual1=qual1[tn5_adapt_cut[0].end:]
                title1=barc+"."+number_out
            else:
                pass
        else:
            pass
    else:
        pass
    return((title1, seq1, qual1))


fq1=sys.argv[1] #read argument fq1="SRR8310661.sra_1.fastq.gz"
fq2=sys.argv[2] #read argument fq2="SRR8310661.sra_2.fastq.gz"

#open fastq files, correct barcode read names based on fuzzy matching (1 substitution allowed per barcode) then out fastq 1 and 2 which pass filter with new read name
with gzip.open(fq1, "rt") as handle1:
    with gzip.open(fq2, "rt") as handle2:
        with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
            with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
                for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):
                    out=change_id(title1, seq1, qual1)
                    if out[0].startswith("SRR") is False:
                        outfile_fq1.write("@%s\n%s\n+\n%s\n" % (out[0], out[1], out[2]))
                        outfile_fq2.write("@%s\n%s\n+\n%s\n" % (out[0], seq2, qual2))

```


Running python script as batch script
ddscatc_barcordes_sbatch.slurm
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-7
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=20gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --
array_in=("SRR8310661" "SRR8310662" "SRR8310663" "SRR8310664" "SRR8310665" "SRR8310666" "SRR8310667" "SRR8310668") 
i=${array_in[$SLURM_ARRAY_TASK_ID]}

srun python /home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC/barcode_ddscatac.py \
/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC/${i}/${i}/${i}.sra_1.fastq.gz \
/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC/${i}/${i}/${i}.sra_2.fastq.gz

```

```bash
ref_dir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"
fastq_in=`ls ${ref_dir}/ddscATAC/*/*/*fastq`; 
for i in $fastq_in; do gzip -f $i & done &
```


10x ATAC v1 Chemistry
https://www.10xgenomics.com/resources/datasets/fresh-cortex-from-adult-mouse-brain-p-50-1-standard-1-2-0
For 10x libraries I will use standard cellranger analysis pipeline.

```bash
mkdir 10x_atac_v1
cd 10x_atac_v1

# Input Files
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fastqs.tar

# Output Files
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_summary.csv

tar -xvf atac_v1_adult_brain_fresh_5k_fastqs.tar

cd /home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1/atac_v1_adult_brain_fresh_5k_fastqs

#10x v1 atac chemistry
tenxatacv1_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1"
tenx_bam="atac_v1_adult_brain_fresh_5k_possorted_bam.bam"
((samtools view -H ${tenxatacv1_outdir}/${tenx_bam}) & \
(samtools view ${tenxatacv1_outdir}/${tenx_bam} | \
        awk '{OFS = "\t"} {for(i=1;i<=NF;i++)  if ($i ~ /^CB:Z:/) {$1=$i"."NR; gsub("CB:Z:","",$1); print}}')) | samtools view -b - | samtools sort -@5 -n - > tenxv1_mus.bam &


```


10x ATAC v2 Chemistry
https://www.10xgenomics.com/resources/datasets/8k-adult-mouse-cortex-cells-atac-v2-chromium-x-2-standard
```bash
mkdir 10x_atac_v2
cd 10x_atac_v2

# Input Files
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs.tar
# Output Files
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_analysis.tar.gz
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_possorted_bam.bam
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_possorted_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_summary.csv


#Convert 10x bam files back to fastq then realign to same genome
##First put cell barcodes into read name just as with the others

#10x v2 atac chemistry
tenxatacv2_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v2"
tenx_bam="8k_mouse_cortex_ATACv2_nextgem_Chromium_X_possorted_bam.bam"
((samtools view -H ${tenxatacv2_outdir}/${tenx_bam}) & \
(samtools view ${tenxatacv2_outdir}/${tenx_bam} | \
        awk '{OFS = "\t"} {for(i=1;i<=NF;i++)  if ($i ~ /^CB:Z:/) {$1=$i"."NR; gsub("CB:Z:","",$1); print}}')) | samtools view -b - | samtools sort -@5 -n - > tenxv2_mus.bam &

```

s3ATAC
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5289637

For this I used the bam file on SRA since it is in the same reference genome. Downloaded via AWS bucket.

```bash
#metadata
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5289nnn/GSM5289637/suppl/GSM5289637_s3atac.mm10.metadata.csv.gz
```

Once all fastq files have their cell-barcode identifier within their read name, use scitools for fastq alignments and processing.

```bash
scitools="/home/groups/oroaklab/src/scitools/scitools-dev/scitools"
sciDROP_70k_dir="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_70k"
sciatac_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciATAC/SRR13437232"
scimap_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciMAP/SRR13437233"
snatac_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/snATAC/SRR6768122/SRR6768122"
ddscATAC_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC"
tenxatacv1_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1"
tenxatacv2_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v2"
s3atac_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/s3atac"
ref_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"

#Alignment
##sciatac alignment
$scitools fastq-align -m 5G -n -t 20 -r 20 mm10 \
$sciatac_outdir/sciATAC_mus \
$sciatac_outdir/SRR13437232_1.barc.fastq.gz \
$sciatac_outdir/SRR13437232_2.barc.fastq.gz &

##scimap alignment
$scitools fastq-align -m 5G -n -t 20 -r 20 mm10 \
$scimap_outdir/sciMAP_mus \
$scimap_outdir/SRR13437233_1.barc.fastq.gz \
$scimap_outdir/SRR13437233_2.barc.fastq.gz &

##snATAC alignment
$scitools fastq-align -m 5G -n -t 20 -r 20 mm10 \
$snatac_outdir/snATAC_mus \
$snatac_outdir/SRR6768122.sra_1.barc.fastq.gz \
$snatac_outdir/SRR6768122.sra_2.barc.fastq.gz &

##ddscATAC alignment
for i in "SRR8310661" "SRR8310662" "SRR8310663" "SRR8310664" "SRR8310665" "SRR8310666" "SRR8310667" "SRR8310668"; do
$scitools fastq-align -m 5G -n -t 10 -r 10 mm10 \
$ddscATAC_outdir/ddscATAC_mus_${i} \
$ddscATAC_outdir/${i}/${i}/${i}.sra_1.barc.fastq.gz \
$ddscATAC_outdir/${i}/${i}/${i}.sra_2.barc.fastq.gz ; done &

##Some bam files how read names with the wrong delimter, so now I'm just going to corret (replace a "." with a ":" so it is scitools compatible)
#sciATAC
((samtools view -H $sciatac_outdir/sciATAC_mus.nsrt.bam) && (samtools view $sciatac_outdir/sciATAC_mus.nsrt.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $sciatac_outdir/sciATAC_mus.RG.nsrt.bam &
#sciMAP
((samtools view -H $scimap_outdir/sciMAP_mus.nsrt.bam) && (samtools view $scimap_outdir/sciMAP_mus.nsrt.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $scimap_outdir/sciMAP_mus.RG.nsrt.bam &
#snATAC
((samtools view -H $snatac_outdir/snATAC_mus.nsrt.bam) && (samtools view $snatac_outdir/snATAC_mus.nsrt.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $snatac_outdir/snATAC_mus.RG.nsrt.bam &
#10xv1
((samtools view -H $tenxatacv1_outdir/tenxv1_mus.bam) && (samtools view $tenxatacv1_outdir/tenxv1_mus.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $tenxatacv1_outdir/tenxv1_mus.RG.bam &
#10xv2
((samtools view -H $tenxatacv2_outdir/tenxv2_mus.bam) && (samtools view $tenxatacv2_outdir/tenxv2_mus.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $tenxatacv2_outdir/tenxv2_mus.RG.bam &
#ddscATAC_outdir
for i in "SRR8310661" "SRR8310662" "SRR8310663" "SRR8310664" "SRR8310665" "SRR8310666" "SRR8310667" "SRR8310668"; do
    ((samtools view -H $ddscATAC_outdir/ddscATAC_mus_${i}.nsrt.bam ) && (samtools view $ddscATAC_outdir/ddscATAC_mus_${i}.nsrt.bam |  awk 'OFS="\t" {gsub(/[.]/,":",$1);print}')) | samtools view -b > $ddscATAC_outdir/ddscATAC_mus_${i}.nsrt.RG.bam & done &

```
Project bams for each method for scaled library complexity

```bash
$scitools bam-project -X ${sciatac_outdir}/sciATAC_mus.RG.nsrt.bam &
$scitools bam-project -X ${scimap_outdir}/sciMAP_mus.RG.nsrt.bam &
#prefilter snatac before continuing
$scitools bam-filter -L ${snatac_outdir}/GSM2668124_filtcells.list.txt -O ${snatac_outdir}/snATAC_mus.cellidfilt ${snatac_outdir}/snATAC_mus.RG.nsrt.bam &
$scitools bam-project -X ${snatac_outdir}/snATAC_mus.cellidfilt.filt.bam &
$scitools bam-project -X ${tenxatacv1_outdir}/tenxv1_mus.RG.bam &
$scitools bam-project -X ${tenxatacv2_outdir}/tenxv2_mus.RG.bam &
for i in "SRR8310661" "SRR8310662" "SRR8310663" "SRR8310664" "SRR8310665" "SRR8310666" "SRR8310667" "SRR8310668"; do
$scitools bam-project -X $ddscATAC_outdir/ddscATAC_mus_${i}.nsrt.RG.bam; done & 
$scitools bam-project -X ${s3atac_outdir}/s3atac.mm10.bam  &


#Project read count for sciDROP mm10 across tech comparisons (towards end of code). 
sciDROP_70k_dir="/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_70k"
awk 'OFS="\t" {a=substr($2,1,10);print $2,a}' mm10.complexity.txt > arbitrary_barcsplit.annot #will be used to split the big bam into ~30 smaller bams, so more memory managable
$scitools bam-split -A arbitrary_barcsplit.annot mm10.bam & wait 
for i in mm10.C*.bam; do
    $scitools bam-project -X $i ; done &

```

Remove PCR duplicate reads
```bash
#Remove duplicates based on cellID, chromosome and start sites per read
##using the name sorted (-n) barcode based removal of duplicates
#Note: by default this also filters out cells with <1000 unique reads. This is the same way sciDROP data was processed.

#back to removal of duplicate reads
$scitools bam-rmdup -n -t 4 ${snatac_outdir}/snATAC_mus.cellidfilt.filt.bam & #done
$scitools bam-rmdup -t 4 ${tenxatacv1_outdir}/tenxv1_mus.RG.bam & #done
$scitools bam-rmdup -t 4 ${tenxatacv2_outdir}/tenxv2_mus.RG.bam & #done
$scitools bam-rmdup -n -t 4 ${s3atac_outdir}/s3atac.mm10.bam & #done
$scitools bam-rmdup -n -t 4 ${sciatac_outdir}/sciATAC_mus.RG.nsrt.bam & #done
$scitools bam-rmdup -n -t 4 ${scimap_outdir}/sciMAP_mus.RG.nsrt.bam & #done
for i in "SRR8310661" "SRR8310662" "SRR8310663" "SRR8310664" "SRR8310665" "SRR8310666" "SRR8310667" "SRR8310668"; do
$scitools bam-rmdup -n -t 4 $ddscATAC_outdir/ddscATAC_mus_${i}.nsrt.RG.bam ; done & #done

```

Filter all bam files

```bash
$scitools bam-filter -N 1000 -C ${sciatac_outdir}/sciATAC_mus.RG.complexity.txt ${sciatac_outdir}/sciATAC_mus.RG.bbrd.q10.bam &
$scitools bam-filter -N 1000 -C ${scimap_outdir}/sciMAP_mus.RG.complexity.txt ${scimap_outdir}/sciMAP_mus.RG.bbrd.q10.bam &
#snatac is prefiltered
for i in 1 2 3 4 5 6 7 8; do
$scitools bam-filter -N 1000 -C ${ddscATAC_outdir}/ddscATAC_mus_SRR831066${i}.RG.complexity.txt ${ddscATAC_outdir}/ddscATAC_mus_SRR831066${i}.RG.bbrd.q10.bam & done & 
$scitools bam-filter -N 1000 -C ${tenxatacv1_outdir}/tenxv1_mus.RG.complexity.txt ${tenxatacv1_outdir}/tenxv1_mus.RG.bbrd.q10.bam &
$scitools bam-filter -N 1000 -C ${tenxatacv2_outdir}/tenxv2_mus.RG.complexity.txt ${tenxatacv2_outdir}/tenxv2_mus.RG.bbrd.q10.bam &
$scitools bam-filter -N 1000 -C ${s3atac_outdir}/s3atac.mm10.complexity.txt ${s3atac_outdir}/s3atac.mm10.bbrd.q10.bam &
#scidrop is prefiltered
```
Merge all bam files and generate a tabix file for fragments.

```bash
$scitools bam-merge -m 10G ${ref_outdir}/all_methods_merged_mm10.bbrd.q10.bam \
${sciatac_outdir}/sciATAC_mus.RG.bbrd.q10.filt.bam \
${scimap_outdir}/sciMAP_mus.RG.bbrd.q10.filt.bam \
${snatac_outdir}/snATAC_mus.cellidfilt.filt.bbrd.q10.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310661.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310662.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310663.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310664.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310665.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310666.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310667.RG.bbrd.q10.filt.bam \
${ddscATAC_outdir}/ddscATAC_mus_SRR8310668.RG.bbrd.q10.filt.bam \
${tenxatacv1_outdir}/tenxv1_mus.RG.bbrd.q10.filt.bam \
${tenxatacv2_outdir}/tenxv2_mus.RG.bbrd.q10.filt.bam \
${s3atac_outdir}/s3atac.mm10.bbrd.q10.filt.bam \
${sciDROP_70k_dir}/mm10.bbrd.q10.bam 

#Make fragment files
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"

#mouse processing
input_bam=${ref_outdir}/all_methods_merged_mm10.bbrd.q10.bam
$scitools bam-tssenrich -X all_methods_merged_mm10.bbrd.q10.bam mm10 &


output_name=${input_bam::-4}
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 20G -T . --parallel=10 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz ;
zcat all_methods_merged_mm10.bbrd.q10.fragments.tsv.gz | awk 'OFS="\t" { if($1 !~ /M|Y|L|K|G|Un|Random|Alt|random/) {print $0}}' | $bgzip > $output_name.fragments.tsv.gz ;
$tabix -p bed $output_name.fragments.tsv.gz ; 

```

Call peaks by read pileups
```bash
$scitools callpeaks -f mm10 ${ref_outdir}/all_methods_merged_mm10.bbrd.q10.bam &

```
Make a counts matrix and a list of fragment sizes per cell ID.

```bash
input_bed="/home/groups/CEDAR/mulqueen/mouse_brain_ref/all_methods_merged_mm10.bbrd.q10.500.bed"
input_bam="/home/groups/CEDAR/mulqueen/mouse_brain_ref/all_methods_merged_mm10.bbrd.q10.bam"
scitools="/home/groups/oroaklab/src/scitools/scitools-dev/scitools"
ref_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref"

#Make counts matrix on merged peaks per technology
$scitools atac-counts \
$input_bam \
$input_bed &

#make an annotation file of bam source per cell id and scitools wrapper for samtools isize
samtools view all_methods_merged_mm10.bbrd.q10.bam | awk 'OFS="\t" {split($1,a,":");print a[1],a[3]}' | sort --parallel=10 -T . -S 2G | uniq > bam_id.annot ;
$scitools isize -A ${ref_outdir}/all_methods.tech.annot ${ref_outdir}/all_methods_merged_mm10.bbrd.q10.bam &

```

## Get author-defined metadata per cell type

```bash
scimap_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciMAP/SRR13437233"
tenxatacv1_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1"
tenxatacv2_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v2"
ddscATAC_outdir="/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC"

#sciMAP and sciATAC data
cd $scimap_outdir
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164849/suppl/GSE164849_sciMAP.metadata.csv.gz

#snATAC

#10x v1
cd $tenxatacv1_outdir
wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_analysis.tar.gz
tar -xvf atac_v1_adult_brain_fresh_5k_analysis.tar.gz
#clusters in /home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1/analysis/clustering/graphclust/clusters.csv

#10x v2
cd $tenxatacv1_outdir
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/8k_mouse_cortex_ATACv2_nextgem_Chromium_X/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_analysis.tar.gz
tar -xvf 8k_mouse_cortex_ATACv2_nextgem_Chromium_X_analysis.tar.gz
#clusters in /home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v2/analysis/clustering/graphclust/clusters.csv

#ddATAC
cd $ddscATAC_outdir
wget https://github.com/buenrostrolab/dscATAC_analysis_code/blob/master/mousebrain/data/mousebrain-master_dataframe.rds
#clusters in dat<-readRDS("/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC/mousebrain-master_dataframe.rds")
```

### Set up Seurat Object

```R
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")

# extract gene annotations from EnsDb
mm10_annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(mm10_annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(mm10_annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead
genome(mm10_annotation) <- "mm10"
seqlevelsStyle(mm10_annotation) <- 'UCSC'


#function to read in sparse matrix format from atac-count
read_in_sparse<-function(x){ #x is character file prefix followed by .bbrd.q10.500.counts.sparseMatrix.values.gz
    IN<-as.matrix(read.table(paste0(x,".counts.sparseMatrix.values.gz")))
    IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
    COLS<-read.table(paste0(x,".counts.sparseMatrix.cols.gz"))
    colnames(IN)<-COLS$V1
    ROWS<-read.table(paste0(x,".counts.sparseMatrix.rows.gz"))
    row.names(IN)<-ROWS$V1
    writeMM(IN,file=paste0(x,".counts.mtx")) #this is to generate counts matrices in scrublet friendly format
    return(IN)
}

mm10_counts<-read_in_sparse("all_methods_merged_mm10.bbrd.q10.500") # make mm10 counts matrix from sparse matrix


#write out as MM format
#Read in fragment path for coverage plots
mm10_fragment.path="/home/groups/CEDAR/mulqueen/mouse_brain_ref/all_methods_merged_mm10.bbrd.q10.fragments.tsv.gz"

#Generate ChromatinAssay Objects
mm10_chromatinassay <- CreateChromatinAssay(
  counts = mm10_counts,
  #genome="mm10",
  min.cells = 1,
  annotation=mm10_annotation,
  sep=c("_","_")#,
  #fragments=mm10_fragment.path
)


#Create Seurat Objects
mm10_atac <- CreateSeuratObject(
  counts = mm10_chromatinassay,
  assay = "peaks"
)

#Meta.data to be updated after clustering
bamid_annot<-read.table("bam_id.annot",header=F)
bamid_readable<-c("BAMID=1"="sciATAC",
    "BAMID=2"="sciMAP",
    "BAMID=3"="snATAC",
    "BAMID=4"="ddscATAC",
    "BAMID=5"="ddscATAC",
    "BAMID=6"="ddscATAC",
    "BAMID=7"="ddscATAC",
    "BAMID=8"="ddscATAC",
    "BAMID=9"="ddscATAC",
    "BAMID=10"="ddscATAC",
    "BAMID=11"="ddscATAC",
    "BAMID=12"="tenxv1",
    "BAMID=13"="tenxv2",
    "BAMID=14"="s3ATAC",
    "BAMID=15"="sciDROP") #this is based on the order of bam-merge scitools call above
bamid<-setNames(nm=bamid_annot$V1,bamid_annot$V2)
bamid_reads<-setNames(nm=bamid_annot$V1,bamid_readable[bamid_annot$V2])
mm10_atac<-AddMetaData(mm10_atac,bamid,col.name="bam_id")
mm10_atac<-AddMetaData(mm10_atac,bamid_reads,col.name="tech")
tech_out<-cbind(colnames(mm10_atac),mm10_atac@meta.data$tech)
write.table(tech_out,file="all_methods.tech.annot",sep="\t",quote=F,col.names=F,row.names=F)

#saving unprocessed SeuratObjects
saveRDS(mm10_atac,file="allmethods_merged_SeuratObject.Rds")

mm10_atac <- RunTFIDF(mm10_atac)
mm10_atac <- FindTopFeatures(mm10_atac, min.cutoff = 20)
mm10_atac <- RunSVD(mm10_atac)
mm10_atac <- RunUMAP(mm10_atac, dims = 2:50, reduction = 'lsi')

plt<-DimPlot(mm10_atac, group.by = 'tech', pt.size = 0.1)
ggsave(plt,file="all_methods.tech.pdf")
system(paste0("slack -F ","all_methods.tech.pdf"," ryan_todo"))

library(harmony,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/lib_backup_210125")
#Correcting bias with harmony
pdf("mm10.harmony.convergence.pdf")
harm_mat<-HarmonyMatrix(mm10_atac@reductions$lsi@cell.embeddings, mm10_atac@meta.data$tech,do_pca=FALSE,nclust=14,plot_convergence=T)
dev.off()
#system("slack -F mm10.harmony.convergence.pdf ryan_todo")
mm10_atac@reductions$harmony<-CreateDimReducObject(embeddings=as.matrix(harm_mat),assay="peaks",key="topic_")
mm10_atac<-RunUMAP(mm10_atac, reduction = "harmony",dims=1:ncol(mm10_atac@reductions$harmony))
mm10_atac <- FindNeighbors(object = mm10_atac,reduction = 'harmony')
mm10_atac <- FindClusters(object = mm10_atac,verbose = TRUE,resolution=0.05)

#Add scimap/sciatac metadata
scimap_meta<-read.csv("/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciMAP/SRR13437233/GSE164849_sciMAP.metadata.csv")
scimap_celltype<-setNames(nm=scimap_meta$Barcode,scimap_meta$Celltype)
mm10_atac<-AddMetaData(mm10_atac,scimap_celltype,col.name="celltype")
plt<-DimPlot(mm10_atac,group.by = 'celltype', pt.size = 0.1)
ggsave(plt,file="all_methods.scimapcelltype.harmony.pdf")
system(paste0("slack -F ","all_methods.scimapcelltype.harmony.pdf"," ryan_todo"))

#Add FRIP based on the merged data peakset
frip<-read.table("all_methods_merged_mm10.bbrd.q10.500.fracOnTarget.values",col.names=c("cellID","FRIP"))
frip_in<-setNames(nm=frip$cellID,frip$FRIP)
mm10_atac<-AddMetaData(mm10_atac,frip_in,col.name="FRIP")
saveRDS(mm10_atac,file="allmethods_merged_SeuratObject.Rds")


#Add FRIP based on the merged data peakset
tss<-read.table("all_methods_merged_mm10.bbrd.q10.TSSenrich.value",col.names=c("cellID","TSSEnrichment"))
tss_in<-setNames(nm=tss$cellID,tss$TSSEnrichment)
mm10_atac<-AddMetaData(mm10_atac,tss_in,col.name="TSSenrichment")
saveRDS(mm10_atac,file="allmethods_merged_SeuratObject.Rds")
```
Generate library complexity per technology from cell projections

```R
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(palettetown)
library(plyr)
library(dplyr)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")
mm10_atac<-readRDS("allmethods_merged_SeuratObject.Rds")

tech_col<-c("snATAC"="grey","sciATAC"="lightcyan4","ddscATAC"="honeydew4","sciMAP"="lightslategrey","sciDROP"="red","s3ATAC"="azure3","tenxv1"="cadetblue","tenxv2"="cadetblue2")
tech_order<-c("snATAC","sciATAC","tenxv1","tenxv2","ddscATAC","sciMAP","s3ATAC","sciDROP")

#Projected complexity models per cell formatted as #cellid int slope
#Vmax=1/V2 ;Km=V3/V2
#test_uniq = (($Vmax*$test_count)/($Km+$test_count))

tenxv1<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v1/tenxv1_mus.RG.read_projections/model.txt"),tech="tenxv1")
tenxv2<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/10x_atac_v2/tenxv2_mus.RG.read_projections/model.txt"),tech="tenxv2")
sciatac<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciATAC/SRR13437232/sciATAC_mus.RG.nsrt.read_projections/model.txt"),tech="sciATAC")
scimap<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/sciMAP/SRR13437233/sciMAP_mus.RG.nsrt.read_projections/model.txt"),tech="sciMAP")
snatac<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/snATAC/SRR6768122/SRR6768122/snATAC_mus.cellidfilt.filt.read_projections/model.txt"),tech="snATAC")
s3atac<-cbind(read.table("/home/groups/CEDAR/mulqueen/mouse_brain_ref/s3atac/s3atac.mm10.read_projections/model.txt"),tech="s3ATAC")
scidrop<-cbind(do.call("rbind", 
    lapply(c('CAGAGAGGAA', 'CAGAGAGGAC', 'CAGAGAGGAG', 'CAGAGAGGAT', 'CAGAGAGGCA', 'CAGAGAGGCC', 'CAGAGAGGCG', 'CAGAGAGGCT', 'CAGAGAGGGA', 'CAGAGAGGGC', 'CAGAGAGGGG', 'CAGAGAGGGT', 'CAGAGAGGTA', 'CAGAGAGGTC', 'CAGAGAGGTG', 'CAGAGAGGTT', 'CTCTCTACAA', 'CTCTCTACAC', 'CTCTCTACAG', 'CTCTCTACAT', 'CTCTCTACCA', 'CTCTCTACCC', 'CTCTCTACCG', 'CTCTCTACCT', 'CTCTCTACGA', 'CTCTCTACGC', 'CTCTCTACGG', 'CTCTCTACGT', 'CTCTCTACTA', 'CTCTCTACTC', 'CTCTCTACTG', 'CTCTCTACTT'), function(x) read.table(paste0("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/sciDROP_70k/mm10.",x,".read_projections/model.txt")))),tech="sciDROP")
ddscatac<-cbind(do.call("rbind",
    lapply(1:8,function(x) read.table(paste0("/home/groups/CEDAR/mulqueen/mouse_brain_ref/ddscATAC/ddscATAC_mus_SRR831066",x,".nsrt.RG.read_projections/model.txt")))),tech="ddscATAC")
compl<-rbind(tenxv1,tenxv2,sciatac,scimap,snatac,s3atac,scidrop,ddscatac)
compl<-compl[!duplicated(compl$V1),]
compl<-compl[compl$V1 %in% row.names(dat),]
f<-function(x,V2,V3) ((1/V2)*x)/((V3/V2)+x) #determine unique read count given a x read effort, using modelled projection

f_read_in<-function(i){
return(cbind(cellid=compl$V1,tech=compl$tech,read_effort=i,read_uniq=unlist(lapply(1:nrow(compl), function(j) f(x=i,V2=compl$V2[j],V3=compl$V3[j])))))
}

out<-lapply(c(5000,10000,20000,50000,100000), f_read_in)
out<-as.data.frame(do.call("rbind",out))
out$read_uniq<-as.numeric(out$read_uniq)
out$tech <- factor(out$tech, levels=tech_order)
out$read_effort <- factor(out$read_effort, levels=c(5000,10000,20000,50000,100000))

plt<-ggplot(out,aes(x=read_effort,color=tech,fill=tech,y=read_uniq))+geom_boxplot(alpha=0.8,fill="white",outlier.shape=NA)+scale_color_manual(values=tech_col)+theme_minimal()+scale_y_continuous(breaks=c(seq(0,120000,10000),c(5000,10000,20000,50000,100000)),limits=c(0,120000))
ggsave(plt,file="complexity_tech.pdf")
#system("slack -F complexity_tech.pdf ryan_todo")


#Adding projected read counts per cell to metadata
out_df<-as.data.frame(do.call("cbind",lapply(unique(out$read_effort), function(x) as.data.frame(cbind(out[out$read_effort==x,]$cellid,out[out$read_effort==x,]$read_uniq)))))
out_df<-out_df[,c(1,which(seq(1,ncol(out_df)) %% 2 == 0))]
colnames(out_df)<-c("cellid",paste0("effort_",unique(out$read_effort)))
out_df<-out_df[out_df$cellid %in% row.names(dat),]
row.names(out_df)<-out_df$cellid
mm10_atac<-AddMetaData(mm10_atac,out_df)
saveRDS(mm10_atac,file="allmethods_merged_SeuratObject.Rds")

```

QC Plot comparisons

```R
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(palettetown)
library(plyr)
library(dplyr)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")
mm10_atac<-readRDS("allmethods_merged_SeuratObject.Rds")
#table(mm10_atac$tech)
#ddscATAC   s3ATAC  sciATAC  sciDROP   sciMAP   snATAC   tenxv1   tenxv2 
#    6611      920     4638    38606     8206     3034     4663     9167 

tech_col<-c("snATAC"="grey","sciATAC"="lightcyan4","ddscATAC"="honeydew4","sciMAP"="lightslategrey","sciDROP"="red","s3ATAC"="azure3","tenxv1"="cadetblue","tenxv2"="cadetblue2")
tech_order<-c("snATAC","sciATAC","tenxv1","tenxv2","ddscATAC","sciMAP","s3ATAC","sciDROP")

#Plot integration
plt<-DimPlot(mm10_atac,group.by = 'tech', pt.size = 0.1)+scale_fill_manual(values=c(NA,NA,NA,NA,NA))+scale_color_manual(values=tech_col)
ggsave(plt,file="all_methods.tech.harmony.pdf")
system(paste0("slack -F ","all_methods.tech.harmony.pdf"," ryan_todo"))

#FRIP boxplot
Idents(mm10_atac)<-"tech"
dat<-as.data.frame(mm10_atac@meta.data)
dat$tech <- factor(dat$tech, levels=tech_order)
plt<-ggplot(dat,aes(x=tech,y=FRIP,color=tech))+geom_boxplot(alpha=0.8,fill="white",outlier.shape=NA)+scale_color_manual(values=tech_col)+ylim(c(0,1))+theme_minimal() #+geom_jitter(alpha=0.1,size=0.1)
ggsave(plt,file="frip_tech.pdf")
#system("slack -F frip_tech.pdf ryan_todo")

#TSS boxplot
Idents(mm10_atac)<-"tech"
dat<-as.data.frame(mm10_atac@meta.data)
dat$tech <- factor(dat$tech, levels=tech_order)
plt<-ggplot(dat,aes(x=tech,y=TSSenrichment,color=tech))+geom_boxplot(alpha=0.8,fill="white",outlier.shape=NA)+ylim(c(0,25))+scale_color_manual(values=tech_col)+theme_minimal() #+geom_jitter(alpha=0.1,size=0.1)
ggsave(plt,file="tss_tech.pdf")
#system("slack -F tss_tech.pdf ryan_todo")

#Stats
library(dplyr)
detach(package:plyr)
dat<-data.frame(mm10_atac@meta.data)
out<-as.data.frame(dat %>% group_by(tech) %>% summarize(mean_tss=mean(TSSenrichment),median_tss=median(TSSenrichment),sd_tss=sd(TSSenrichment),mean_frip=mean(FRIP),median_frip=median(FRIP),sd_frip=sd(FRIP),count=n()))
row.names(out)<-out$tech
out<-out[tech_order,]
out$median_tss/out[out$tech=="sciDROP",]$median_tss
out$median_frip/out[out$tech=="sciDROP",]$median_frip

# tech     mean_tss median_tss sd_tss mean_frip median_frip sd_frip count
#  <chr>       <dbl>      <dbl>  <dbl>     <dbl>       <dbl>   <dbl> <int>
#1 ddscATAC     5.33       5.04   1.70     0.577       0.583  0.0846  6611
#2 s3ATAC       4.17       4.00   1.42     0.411       0.41   0.0744   920
#3 sciATAC      8.95       8.15   4.42     0.591       0.595  0.152   4638
#4 sciDROP      6.03       5.92   1.86     0.549       0.575  0.0902 38606
#5 sciMAP       3.04       2.70   1.42     0.293       0.285  0.0783  8206
#6 snATAC       6.33       6.02   2.10     0.493       0.492  0.0987  3034
#7 tenxv1       9.06       7.96   3.72     0.715       0.721  0.0868  4663
#8 tenxv2       6.12       5.73   3.39     0.488       0.521  0.174   9167


#frip
test_out<-lapply(unique(dat$tech), function(i) setNames(wilcox.test(dat[dat$tech=="sciDROP",]$FRIP,dat[dat$tech==i,]$FRIP)$p.value,nm=i))
names(test_out)<-unique(dat$tech)

#tss
kruskal.test(TSSenrichment ~ tech, data = dat)
pairwise.wilcox.test(dat$TSSenrichment, dat$tech,
                 p.adjust.method = "bonferroni")

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction

#data:  dat$TSSenrichment and dat$tech

#        ddscATAC s3ATAC  sciATAC sciDROP sciMAP  snATAC  tenxv1
#s3ATAC  < 2e-16  -       -       -       -       -       -
#sciATAC < 2e-16  < 2e-16 -       -       -       -       -
#sciDROP < 2e-16  < 2e-16 < 2e-16 -       -       -       -
#sciMAP  < 2e-16  < 2e-16 < 2e-16 < 2e-16 -       -       -
#snATAC  < 2e-16  < 2e-16 < 2e-16 4.1e-06 < 2e-16 -       -
#tenxv1  < 2e-16  < 2e-16 0.018   < 2e-16 < 2e-16 < 2e-16 -
#tenxv2  < 2e-16  < 2e-16 < 2e-16 4.1e-14 < 2e-16 5.4e-14 < 2e-16

#P value adjustment method: bonferroni

pairwise.wilcox.test(dat$FRIP, dat$tech,
                 p.adjust.method = "bonferroni")

#effort_5000
dat$effort_5000<-as.numeric(dat$effort_5000)
pairwise.wilcox.test(dat$effort_5000, dat$tech,
                 p.adjust.method = "bonferroni",alternative="greater")
```

## Reviewer Responses

### Correlation between pseudobulked samples
https://theislab.github.io/scib-reproducibility/

Using gene body count for correlation
```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(Matrix)
library(rliger)
library(SeuratWrappers)
library(parallel)
library(corrplot)
library(patchwork)
library(GenomicRanges)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")
mm10_atac<-readRDS("allmethods_merged_SeuratObject.Rds")
mm10_fragment.path="/home/groups/CEDAR/mulqueen/mouse_brain_ref/all_methods_merged_mm10.bbrd.q10.fragments.tsv.gz" #adding this is now for GA and plotting purposes.
Fragments(mm10_atac)<-CreateFragmentObject(path=mm10_fragment.path)

#from https://github.com/stuart-lab/signac/blob/HEAD/R/utilities.R
CollapseToLongestTranscript <- function(ranges) {
  range.df <- data.table::as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#from https://github.com/stuart-lab/signac/blob/HEAD/R/utilities.R
Extend <- function(
  x,
  upstream = 0,
  downstream = 0,
  from.midpoint = FALSE
    ) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}
#doing feature setting following signac gene activity calculation
#filter to protein coding
#subset to which protein coding gene is longest that has the same name
#extend 2kb upstream for promoter inclusion
feat=mm10_atac@assays$peaks@annotation[mm10_atac@assays$peaks@annotation$gene_biotype=="protein_coding",]
feat<-mclapply(unique(feat$gene_name),function(x) CollapseToLongestTranscript(feat[feat$gene_name==x,]),mc.cores=10) #collapse to longest transcripts
feat<-unlist(as(feat, "GRangesList"))
feat<-setNames(feat,feat$gene_name)#set row names as gene names
feat<-feat[feat@ranges@width<500000,]#filter extra long transcripts
transcripts <- Extend(x = feat,upstream = 2000,downstream = 0)# extend to include promoters

feat_split<-split(transcripts, rep_len(1:300, length(transcripts)))
#parallelize gene count to speed up feature matrix generation

split_gene_count<-function(x){
    FeatureMatrix(fragments = Fragments(mm10_atac),cells=Cells(mm10_atac),
                              features= feat_split[[x]],
                              verbose = TRUE,
                              process_n=20000)
}

mm10_atac_counts<-mclapply(1:length(feat_split),split_gene_count,mc.cores=10)
x<-do.call("rbind",mm10_atac_counts)
mm10_atac_counts<-x
saveRDS(mm10_atac_counts,file="allmethods_merged_SeuratObject_counts.Rds")
mm10_atac_counts<-readRDS(file="allmethods_merged_SeuratObject_counts.Rds")
mm10_atac[['GeneCount']] <- CreateAssayObject(counts = mm10_atac_counts)

mm10_atac <- NormalizeData(
  object = mm10_atac,
  assay = 'GeneCount',
  normalization.method = 'LogNormalize',
  scale.factor = median(mm10_atac$nCount_GeneCount)
)

DefaultAssay(mm10_atac)<-"GeneCount"
mm10_atac <- FindVariableFeatures(mm10_atac)
mm10_atac <- ScaleData(mm10_atac, split.by = "tech", do.center = FALSE)

mm10_atac_genes<-AverageExpression(mm10_atac,assay="GeneCount",group.by="tech")
cor_out<-cor(mm10_atac_genes$GeneCount,method="spearman")

#testing significance in pairwise manner

for(x in 1:ncol(mm10_atac_genes$GeneCount)){
    for(y in 1:ncol(mm10_atac_genes$GeneCount)){
        if(x!=y){
        print(cor.test(x=mm10_atac_genes$GeneCount[,x],y=mm10_atac_genes$GeneCount[,y],method="spearman"))
        }
    }
}
#all have p-value < 2.2e-16

pdf("technology_correlations.pdf")
#corrplot(cor_out,type='upper',col = rev(COL2('RdBu', 100)))
corrplot(cor_out,type='lower',method='number',col = rev(COL2('RdBu', 100)))
dev.off()
#system("slack -F technology_correlations.pdf ryan_todo")
```

Plotting on neuronal signal marker
```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(Matrix)
library(rliger)
library(SeuratWrappers)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")
mm10_atac<-readRDS("allmethods_merged_SeuratObject.Rds")
mm10_fragment.path="/home/groups/CEDAR/mulqueen/mouse_brain_ref/all_methods_merged_mm10.bbrd.q10.fragments.tsv.gz" #adding this is now for GA and plotting purposes.
Fragments(mm10_atac)<-CreateFragmentObject(path=mm10_fragment.path)

#Slc17a7
plt<-CoveragePlot(mm10_atac,group.by="tech",region="Slc17a7",extend.upstream=50000,extend.downstream=50000,tile=TRUE)
ggsave(plt,file="all_methods_Slc17a7_covplot.pdf")
#system("slack -F all_methods_Slc17a7_covplot.pdf ryan_todo")
```

## Modifying Batch Correction across mouse brain sets

Changing batch (method) correction across methods for the mouse brain from Harmony to LIGER
based on https://www.nature.com/articles/s41592-021-01336-8/figures/4
Retrying LIGER replicating scIB's large peak analysis instead:
https://github.com/theislab/scib-reproducibility/blob/main/notebooks/data_preprocessing/mouse_brain_atac/preprocessing_large_dataset.ipynb

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(Matrix)
library(rliger)
library(SeuratWrappers)
setwd("/home/groups/CEDAR/mulqueen/mouse_brain_ref")
mm10_atac<-readRDS("allmethods_merged_SeuratObject.Rds")
mm10_atac_counts<-readRDS(file="allmethods_merged_SeuratObject_counts.Rds")
mm10_atac_txsciatac<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard/mm10_SeuratObject.PF.Rds")

#table(mm10_atac$tech)
#ddscATAC   s3ATAC  sciATAC  sciDROP   sciMAP   snATAC   tenxv1   tenxv2 
#    6611      920     4638    38606     8206     3034     4663     9167 

mm10_atac<-AddMetaData(mm10_atac,mm10_atac_txsciatac$celltype,col.name="txsciatac_celltype")

#following scIB filtering
mm10_atac<-BinarizeCounts(mm10_atac) # binarize peak accessibility
FindTopFeatures(mm10_atac,assay="peaks",min.cutoff=50)# minimum number of cells sharing a feature min_cells = 50 #all features pass so not doign proper subsetting function call
mm10_atac <- subset(x = mm10_atac,subset = nFeature_peaks  > 500) # set a minimum number of cells to keep min_features = 500
mm10_atac <- FindVariableFeatures(mm10_atac,nfeatures=50000) # 150000 most variable windows used in scIB
mm10_atac<-SetAssayData(mm10_atac,assay="peaks",slot="scale.data",new.data=as.matrix(mm10_atac@assays$peaks@data[mm10_atac@assays$peaks@var.features,]))
table()
#mm10_atac <- ScaleData(mm10_atac, split.by = "tech", do.center = FALSE) #scale data, maybe in the future just write @counts to @data
mm10_atac <- RunOptimizeALS(mm10_atac, k = 20, lambda = 5, split.by = "tech")
saveRDS(mm10_atac,"allmethods_merged_SeuratObject.LIGER.genecount.Rds")


mm10_atac.liger <- RunQuantileNorm(mm10_atac, split.by = "tech")
saveRDS(mm10_atac.liger,"allmethods_merged_SeuratObject.LIGER.genecount.Rds")

mm10_atac.liger <- RunUMAP(object = mm10_atac.liger, reduction = 'iNMF', dims = 1:ncol(mm10_atac.liger@reductions$iNMF@cell.embeddings))

plt<-DimPlot(mm10_atac.liger,group.by=c("tech","txsciatac_celltype"))
ggsave(plt,file="allmethods_merged_liger_umap.pdf",width=10)
#system("slack -F allmethods_merged_liger_umap.pdf ryan_todo")

write.table(mm10_atac@meta.data,file="mm10_brain_all_methods.metadata.tsv",sep="\t",col.names=T,row.names=T)
#system("slack -F mm10_brain_all_methods.metadata.tsv ryan_todo")
```


Final R session information with all packages loaded.

```R
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

```R
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/groups/CEDAR/mulqueen/src/miniconda3/lib/libopenblasp-r0.3.17.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   grid      stats4    parallel  stats     graphics  grDevices
 [8] utils     datasets  methods   base     

other attached packages:
 [1] corrplot_0.92                      rliger_1.0.0                      
 [3] cowplot_1.1.1                      plyr_1.8.8                        
 [5] R.utils_2.12.2                     R.oo_1.25.0                       
 [7] R.methodsS3_1.8.2                  motifmatchr_1.12.0                
 [9] chromVAR_1.12.0                    seriation_1.4.1                   
[11] forcats_0.5.2                      stringr_1.5.0                     
[13] purrr_1.0.1                        readr_2.1.3                       
[15] tidyr_1.2.1                        tibble_3.1.8                      
[17] tidyverse_1.3.2                    BSgenome.Mmusculus.UCSC.mm10_1.4.0
[19] circlize_0.4.15                    reshape2_1.4.4                    
[21] viridis_0.6.2                      viridisLite_0.4.1                 
[23] dendextend_1.16.0                  ggdendro_0.1.23                   
[25] palettetown_0.1.1.90000            RColorBrewer_1.1-3                
[27] clustree_0.5.0                     ggraph_2.1.0                      
[29] ggrepel_0.9.2                      dplyr_1.0.9                       
[31] BSgenome.Hsapiens.UCSC.hg38_1.4.3  BSgenome_1.58.0                   
[33] rtracklayer_1.50.0                 TFBSTools_1.28.0                  
[35] JASPAR2020_0.99.10                 cicero_1.8.1                      
[37] Gviz_1.34.1                        monocle_2.18.0                    
[39] DDRTree_0.1.5                      irlba_2.3.5.1                     
[41] VGAM_1.1-7                         SeuratWrappers_0.3.0              
[43] harmony_1.0                        Rcpp_1.0.9                        
[45] cisTopic_0.3.0                     ComplexHeatmap_2.6.2              
[47] patchwork_1.1.2                    Matrix_1.5-3                      
[49] EnsDb.Mmusculus.v79_2.99.0         EnsDb.Hsapiens.v86_2.99.0         
[51] ensembldb_2.14.1                   AnnotationFilter_1.14.0           
[53] GenomicFeatures_1.42.3             AnnotationDbi_1.52.0              
[55] Biobase_2.50.0                     GenomicRanges_1.42.0              
[57] GenomeInfoDb_1.26.7                SeuratObject_4.1.3                
[59] Seurat_4.3.0                       Signac_1.5.0                      
[61] Biostrings_2.58.0                  XVector_0.30.0                    
[63] IRanges_2.24.1                     S4Vectors_0.28.1                  
[65] BiocGenerics_0.36.1                ggplot2_3.4.0                  

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                  Hmisc_4.7-2                
  [3] ica_1.0-3                   RcppRoll_0.3.0             
  [5] Rsamtools_2.6.0             foreach_1.5.2              
  [7] lmtest_0.9-40               crayon_1.5.2               
  [9] MASS_7.3-58.1               nlme_3.1-161               
 [11] backports_1.4.1             reprex_2.0.2               
 [13] qlcMatrix_0.9.7             rlang_1.0.6                
 [15] readxl_1.4.1                ROCR_1.0-11                
 [17] ca_0.71.1                   limma_3.46.0               
 [19] BiocParallel_1.24.1         rjson_0.2.21               
 [21] CNEr_1.26.0                 bit64_4.0.5                
 [23] glue_1.6.2                  pheatmap_1.0.12            
 [25] poweRlaw_0.70.6             text2vec_0.6.3             
 [27] sctransform_0.3.5           spatstat.sparse_3.0-0      
 [29] spatstat.geom_3.0-3         haven_2.5.1                
 [31] tidyselect_1.2.0            SummarizedExperiment_1.20.0
 [33] fitdistrplus_1.1-8          XML_3.99-0.13              
 [35] zoo_1.8-11                  GenomicAlignments_1.26.0   
 [37] xtable_1.8-4                magrittr_2.0.3             
 [39] cli_3.4.1                   zlibbioc_1.36.0            
 [41] rstudioapi_0.14             miniUI_0.1.1.1             
 [43] sp_1.6-0                    rpart_4.1.19               
 [45] fastmatch_1.1-3             shiny_1.7.4                
 [47] xfun_0.36                   askpass_1.1                
 [49] clue_0.3-63                 cluster_2.1.4              
 [51] caTools_1.18.2              TSP_1.2-1                  
 [53] tidygraph_1.2.2             doSNOW_1.0.20              
 [55] KEGGREST_1.30.1             biovizBase_1.38.0          
 [57] listenv_0.9.0               TFMPvalue_0.0.9            
 [59] lda_1.4.2                   png_0.1-8                  
 [61] future_1.30.0               withr_2.5.0                
 [63] lsa_0.73.3                  bitops_1.0-7               
 [65] slam_0.1-50                 ggforce_0.4.1              
 [67] cellranger_1.1.0            GSEABase_1.52.1            
 [69] sparsesvd_0.2-2             pracma_2.4.2               
 [71] pillar_1.8.1                GlobalOptions_0.1.2        
 [73] cachem_1.0.6                fs_1.5.2                   
 [75] hdf5r_1.3.8                 GetoptLong_1.0.5           
 [77] vctrs_0.5.1                 ellipsis_0.3.2             
 [79] generics_0.1.3              tools_4.0.3                
 [81] foreign_0.8-84              feather_0.3.5              
 [83] mlapi_0.1.1                 munsell_0.5.0              
 [85] tweenr_2.0.2                DelayedArray_0.16.3        
 [87] fastmap_1.1.0               compiler_4.0.3             
 [89] HSMMSingleCell_1.10.0       abind_1.4-5                
 [91] httpuv_1.6.8                plotly_4.10.1              
 [93] GenomeInfoDbData_1.2.4      gridExtra_2.3              
 [95] riverplot_0.10              lattice_0.20-45            
 [97] deldir_1.0-6                snow_0.4-4                 
 [99] utf8_1.2.2                  later_1.3.0                
[101] BiocFileCache_1.14.0        jsonlite_1.8.4             
[103] scales_1.2.1                docopt_0.7.1               
[105] graph_1.68.0                pbapply_1.7-0              
[107] lazyeval_0.2.2              promises_1.2.0.1           
[109] doParallel_1.0.17           latticeExtra_0.6-30        
[111] goftest_1.2-3               spatstat.utils_3.0-1       
[113] reticulate_1.27             checkmate_2.1.0            
[115] Rtsne_0.16                  dichromat_2.0-0.1          
[117] uwot_0.1.14                 igraph_1.3.5               
[119] survival_3.5-0              htmltools_0.5.4            
[121] memoise_2.0.1               VariantAnnotation_1.36.0   
[123] lgr_0.4.4                   graphlayouts_0.8.4         
[125] arrow_5.0.0.2               digest_0.6.31              
[127] assertthat_0.2.1            RhpcBLASctl_0.21-247.1     
[129] mime_0.12                   rappdirs_0.3.3             
[131] densityClust_0.3.2          registry_0.5-1             
[133] RSQLite_2.2.20              future.apply_1.10.0        
[135] remotes_2.4.2               data.table_1.14.6          
[137] blob_1.2.3                  fastICA_1.2-3              
[139] Formula_1.2-4               googledrive_2.0.0          
[141] Cairo_1.5-12.2              ProtGenerics_1.22.0        
[143] RCurl_1.98-1.9              broom_1.0.2                
[145] hms_1.1.2                   modelr_0.1.10              
[147] colorspace_2.0-3            base64enc_0.1-3            
[149] BiocManager_1.30.19         shape_1.4.6                
[151] nnet_7.3-18                 mclust_6.0.0               
[153] RANN_2.6.1                  ggseqlogo_0.1              
[155] fansi_1.0.3                 tzdb_0.3.0                 
[157] parallelly_1.34.0           SnowballC_0.7.0            
[159] R6_2.5.1                    ggridges_0.5.4             
[161] lifecycle_1.0.3             googlesheets4_1.0.1        
[163] curl_5.0.0                  leiden_0.4.3               
[165] RcppAnnoy_0.0.20            iterators_1.0.14           
[167] spatstat.explore_3.0-5      htmlwidgets_1.6.1          
[169] polyclip_1.10-4             biomaRt_2.46.3             
[171] timechange_0.2.0            seqLogo_1.56.0             
[173] rvest_1.0.3                 globals_0.16.2             
[175] rsparse_0.5.1               openssl_2.0.5              
[177] htmlTable_2.4.1             spatstat.random_3.0-1      
[179] progressr_0.13.0            lubridate_1.9.0            
[181] codetools_0.2-18            matrixStats_0.63.0         
[183] GO.db_3.12.1                FNN_1.1.3.1                
[185] gtools_3.9.4                prettyunits_1.1.1          
[187] dbplyr_2.2.1                gtable_0.3.1               
[189] DBI_1.1.3                   tensor_1.5                 
[191] httr_1.4.4                  KernSmooth_2.23-20         
[193] RcisTarget_1.11.10          stringi_1.7.12             
[195] progress_1.2.2              farver_2.1.1               
[197] annotate_1.68.0             DT_0.27                    
[199] xml2_1.3.3                  combinat_0.0-8             
[201] AUCell_1.13.3               interp_1.1-3               
[203] float_0.3-0                 scattermore_0.8            
[205] bit_4.0.5                   jpeg_0.1-10                
[207] MatrixGenerics_1.2.1        spatstat.data_3.0-0        
[209] gargle_1.2.1                pkgconfig_2.0.3            
[211] DirichletMultinomial_1.32.0 knitr_1.41                 
```