---
title: sciDROP
layout: analysis
permalink: /scidrop/
category: alternative
---


## Processing for sciDROP

![sciDROP Overview](/assets/images/sciDROP.png){:width="80%"}


This notebook details the processing of the "20K" and "70K" loaded mouse brain and human cortex samples. It begins with scitools wrapper functions for intial alignment to a concatenated mouse and human genome, following with splitting of reads and realignment to separate human and mouse genomes. It then follows the established scitools formation of a counts matrix and Signac processing.

```bash
#libraries were generated as two separate lanes of a NovaSeq S4 flowcell.
#bcl2fastq was run prior to the run transfer
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

        table(dat$i7_idx_namez`)
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
            system("slack -F 70k.gem_count.pdf ryan_todo")

        dat$condition<-"Mix"
        dat[dat$tn5_row %in% c("A","B"),]$condition<-"Human"
        dat[dat$tn5_row %in% c("C","D"),]$condition<-"Mouse"

        ggplot(dat,aes(x=perc_uniq,y=log10(uniq_reads_q10),color=species_call))+geom_point()+theme_bw()+ylim(c(0,7))+xlim(c(0,100))
        ggsave("complexity.pdf")
        system("slack -F complexity.pdf ryan_todo")

        library(patchwork)
        plt_barnyard<-ggplot(dat[dat$condition=="Mix",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Barnyard")
        plt_human<-ggplot(dat[dat$condition=="Human",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Human")
        plt_mouse<-ggplot(dat[dat$condition=="Mouse",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Mouse")
        plt<-plt_barnyard+plt_human+plt_mouse
        ggsave(plt,file="barnyard.pdf",width=20)
        system("slack -F barnyard.pdf ryan_todo")

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
        system("slack -F 20k.gem_count.pdf ryan_todo")

        dat$condition<-"Mix"
        dat[dat$tn5_row %in% c("A","B"),]$condition<-"Human"
        dat[dat$tn5_row %in% c("C","D"),]$condition<-"Mouse"

        ggplot(dat,aes(x=perc_uniq,y=log10(uniq_reads_q10),color=species_call))+geom_point()+theme_bw()+ylim(c(0,7))+xlim(c(0,100))
        ggsave("complexity.pdf")
        system("slack -F complexity.pdf ryan_todo")

        library(patchwork)
        plt_barnyard<-ggplot(dat[dat$condition=="Mix",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Barnyard")
        plt_human<-ggplot(dat[dat$condition=="Human",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Human")
        plt_mouse<-ggplot(dat[dat$condition=="Mouse",],aes(x=hg38_count,y=mm10_count,color=species_call,alpha=0.1))+geom_point()+theme_bw()+ylim(c(0,200000))+xlim(c(0,200000))+ggtitle("Mouse")
        plt<-plt_barnyard+plt_human+plt_mouse
        ggsave(plt,file="barnyard.pdf",width=20)
        system("slack -F barnyard.pdf ryan_todo")

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
system("slack -F mm10.hg38.complexity.2d.pdf ryan_todo")

#hard coded these numbers just because i had them for a meeting
#cell_count_75k<-data.frame(count=c(10040,18985,170,12663-(141+19),19530-(484+50),141+19,484+50),names=c("by_h","by_m","by_mix","hum","mus","scrub_h","scrub_m"),loading=c("75k"))
#cell_count_20k<-data.frame(count=c(3001,6572,24,3703-(27),5841-(449),27,449),names=c("by_h","by_m","by_mix","hum","mus","scrub_h","scrub_m"),loading=c("20k"))
#cell_count<-rbind(cell_count_75k,cell_count_20k)
#plt<-ggplot(cell_count,aes(x=loading,y=count,fill=factor(names,levels=rev(c("mus","hum","by_h","by_m","by_mix","scrub_h","scrub_m")))))+geom_bar(position="stack",stat="identity")
#ggsave(plt,file="mm10.hg38.cellcount.pdf")
#system("slack -F mm10.hg38.cellcount.pdf ryan_todo")

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

system("slack -F hg38_atac_model_selection.pdf ryan_todo")
system("slack -F mm10_atac_model_selection.pdf ryan_todo")

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
system("slack -F hg38.umap.i7idx.pdf ryan_todo")

plt<-DimPlot(mm10_atac,group.by=c('seurat_clusters','predicted_doublets',"pcr_idx"))
ggsave(plt,file="mm10.umap.i7idx.pdf",width=10)
system("slack -F mm10.umap.i7idx.pdf ryan_todo")

plt<-FeaturePlot(hg38_atac,feature=c('doublet_scores'))
ggsave(plt,file="hg38.umap.scrub.pdf")
system("slack -F hg38.umap.scrub.pdf ryan_todo")

plt<-FeaturePlot(mm10_atac,feature=c('doublet_scores'))
ggsave(plt,file="mm10.umap.scrub.pdf")
system("slack -F mm10.umap.scrub.pdf ryan_todo")


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
system("slack -F hg38.harmony.convergence.pdf ryan_todo")
hg38_atac@reductions$harmony<-CreateDimReducObject(embeddings=as.matrix(harm_mat),assay="peaks",key="topic_")
hg38_atac<-RunUMAP(hg38_atac, reduction = "harmony",dims=1:ncol(hg38_atac@reductions$harmony))
hg38_atac <- FindNeighbors(object = hg38_atac,reduction = 'harmony')
hg38_atac <- FindClusters(object = hg38_atac,verbose = TRUE,resolution=0.05)

plt<-DimPlot(hg38_atac,group.by=c("pcr_idx","seurat_clusters"))
ggsave(plt,file="hg38.umap.i7idx.harm.pdf",width=15)
system("slack -F hg38.umap.i7idx.harm.pdf ryan_todo")

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
system("slack -F mm10.umap.i7idx.harm.pdf ryan_todo")

plt1<-DimPlot(hg38_atac,group.by="seurat_clusters")
plt2<-DimPlot(mm10_atac,group.by="seurat_clusters")
plt<-plt1+plt2
ggsave(plt,file="hg38_mm10_seurat.clusters.pdf")
system("slack -F hg38_mm10_seurat.clusters.pdf ryan_todo")

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
set.seed(1234)
library(dplyr)
library(ggrepel)


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
    object_input<-RunUMAP(object_input,reduction="cistopic",dims=1:n_topics)    #finally recluster
    #Clustering with multiple resolutions to account for different celltype complexities
    object_input <- FindNeighbors(object = object_input, reduction = 'cistopic', dims = 1:n_topics)
    object_input <- FindClusters(object = object_input,resolution=0.1)
    object_input <- FindClusters(object = object_input,verbose = TRUE,resolution=0.2)
    object_input <- FindClusters(object = object_input,verbose = TRUE,resolution=0.5)
    object_input <- FindClusters(object = object_input,verbose = TRUE,resolution=0.9)
    
    saveRDS(object_input,paste0(outname,".75k.subset.SeuratObject.Rds"))
    plt<-DimPlot(object_input,group.by=c('peaks_snn_res.0.1','peaks_snn_res.0.2','peaks_snn_res.0.5','peaks_snn_res.0.9'))
    ggsave(plt,file=paste(outname,"clustering.pdf",sep="."))
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

    hg38_resolution_list<-setNames(c("peaks_snn_res.0.2","peaks_snn_res.0.2","peaks_snn_res.0.1","peaks_snn_res.0.1",
        "peaks_snn_res.0.5","peaks_snn_res.0.2","peaks_snn_res.0.1","peaks_snn_res.0.2"),hg38_celltype_list)
    cell_order<-row.names(hg38_atac@meta.data)

    for(celltype.x in hg38_celltype_list){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        embedding_order<-row.names(atac_sub@reductions$umap@cell.embeddings)
        row_order<-match(cell_order,embedding_order,nomatch=0)
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$seurat_subcluster<-as.character(unlist(atac_sub@meta.data[which(colnames(atac_sub@meta.data)==hg38_resolution_list[celltype.x])]))
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$subcluster_x<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,1])
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$subcluster_y<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,2])
    }

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

    mm10_resolution_list<-setNames(c("peaks_snn_res.0.2","peaks_snn_res.0.2","peaks_snn_res.0.2","peaks_snn_res.0.1",
        "peaks_snn_res.0.2","peaks_snn_res.0.2","peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1"),mm10_celltype_list)
    cell_order<-row.names(mm10_atac@meta.data)

    for(celltype.x in mm10_celltype_list){
        outname<-strsplit(celltype.x,split="[.]")[[1]][1]
        atac_sub<-readRDS(paste0("./subcluster/",outname,".75k.subset.SeuratObject.Rds"))
        embedding_order<-row.names(atac_sub@reductions$umap@cell.embeddings)
        row_order<-match(cell_order,embedding_order,nomatch=0)
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$seurat_subcluster<-as.character(unlist(atac_sub@meta.data[which(colnames(atac_sub@meta.data)==mm10_resolution_list[celltype.x])]))
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_x<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,1])
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_y<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,2])
    }

    as.data.frame(mm10_atac@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
    #na values are doublets or excluded indexes
    saveRDS(mm10_atac,"mm10_SeuratObject.Rds")

```

### Recoloring subclusters

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

    #Read in data and modify to monocle CDS file
    #read in RDS file.
    hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
    hg38_atac$cluster_ID<-paste(hg38_atac$seurat_clusters,hg38_atac$seurat_subcluster,sep="_")

    #########Coloring and plotting human data#####################
    dat<-as.data.frame(hg38_atac@meta.data)
    cell_order<-row.names(dat)
    embedding_order<-row.names(hg38_atac@reductions$umap@cell.embeddings)
    row_order<-match(cell_order,embedding_order)
    dat$umap_x<-as.numeric(hg38_atac@reductions$umap@cell.embeddings[row_order,1])
    dat$umap_y<-as.numeric(hg38_atac@reductions$umap@cell.embeddings[row_order,2])
    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]
    dat$seurat_clusters<-as.character(dat$seurat_clusters)
   
    dat$cluster_col<-"NULL"
    dat[dat$seurat_clusters=="0",]$cluster_col<-"#1f78b4"
    dat[dat$seurat_clusters=="1",]$cluster_col<-"#b2df8a"
    dat[dat$seurat_clusters=="2",]$cluster_col<-"#FDE725"
    dat[dat$seurat_clusters=="3",]$cluster_col<-"#E31A1C"
    dat[dat$seurat_clusters=="4",]$cluster_col<-"#6a3d9a"
    dat[dat$seurat_clusters=="5",]$cluster_col<-"#B15928"
    dat[dat$seurat_clusters=="6",]$cluster_col<-"#ff7f00"
    dat[dat$seurat_clusters=="7",]$cluster_col<-"#228b22"

    seurat_clus_col<-setNames(unique(dat$cluster_col),unique(dat$seurat_clusters))
    sort(unique(dat$cluster_ID))
    # [1] "0_0"  "0_1"  "0_2"  "0_3"  "0_4"  "0_5"  "0_6"  "0_7"  "0_8"  "0_NA"
    #[11] "1_0"  "1_1"  "1_2"  "1_3"  "1_NA" "2_0"  "2_1"  "2_2"  "3_0"  "3_1"
    #[21] "3_2"  "3_3"  "3_4"  "3_5"  "3_6"  "3_7"  "3_NA" "4_0"  "4_1"  "4_2"
    #[31] "4_3"  "4_NA" "5_0"  "5_1"  "5_2"  "5_3"  "5_4"  "5_NA" "6_0"  "6_NA"
    dat$subcluster_col<-"NULL"
    dat[dat$cluster_ID=="0_0",]$subcluster_col<-"#820933"
    dat[dat$cluster_ID=="0_1",]$subcluster_col<-"#D84797"
    dat[dat$cluster_ID=="0_2",]$subcluster_col<-"#D2FDFF"
    dat[dat$cluster_ID=="0_3",]$subcluster_col<-"#26FFE6"
    dat[dat$cluster_ID=="0_4",]$subcluster_col<-"#001427"
    dat[dat$cluster_ID=="0_5",]$subcluster_col<-"#F4D58D"
    dat[dat$cluster_ID=="0_6",]$subcluster_col<-"#C68866"
    dat[dat$cluster_ID=="0_7",]$subcluster_col<-"#CCC9E7"
    dat[dat$cluster_ID=="0_8",]$subcluster_col<-"#138A87"

    dat[dat$cluster_ID=="1_0",]$subcluster_col<-"#A9F5A4"
    dat[dat$cluster_ID=="1_1",]$subcluster_col<-"#0DA802"
    dat[dat$cluster_ID=="1_2",]$subcluster_col<-"#1C3B1A"

    dat[dat$cluster_ID=="2_0",]$subcluster_col<-"#B38A1B"
    dat[dat$cluster_ID=="2_1",]$subcluster_col<-"#6B0E44"
    dat[dat$cluster_ID=="2_2",]$subcluster_col<-"#E0BB55"

    dat[dat$cluster_ID=="3_0",]$subcluster_col<-"#DBAD6A"
    dat[dat$cluster_ID=="3_1",]$subcluster_col<-"#007EA7"
    dat[dat$cluster_ID=="3_2",]$subcluster_col<-"#C191A1"
    dat[dat$cluster_ID=="3_3",]$subcluster_col<-"#1818b2"

    dat[dat$cluster_ID=="4_0",]$subcluster_col<-"#FAA916"
    dat[dat$cluster_ID=="4_1",]$subcluster_col<-"#5E503F"
    dat[dat$cluster_ID=="4_2",]$subcluster_col<-"#EDFFAB"
    dat[dat$cluster_ID=="4_3",]$subcluster_col<-"#FFE1EA"
    dat[dat$cluster_ID=="4_4",]$subcluster_col<-"#E952DE"
    dat[dat$cluster_ID=="4_5",]$subcluster_col<-"#848FA5"
    dat[dat$cluster_ID=="4_6",]$subcluster_col<-"#832161"
    dat[dat$cluster_ID=="4_7",]$subcluster_col<-"#B0A3D4"
    dat[dat$cluster_ID=="4_8",]$subcluster_col<-"#7fff00"

    dat[dat$cluster_ID=="5_0",]$subcluster_col<-"#230903"
    dat[dat$cluster_ID=="5_1",]$subcluster_col<-"#FFBFB7"
    dat[dat$cluster_ID=="5_2",]$subcluster_col<-"#FFD447"

    dat[dat$cluster_ID=="6_0",]$subcluster_col<-"#80552B"

    dat[dat$cluster_ID=="7_0",]$subcluster_col<-"#355d39"
    dat[dat$cluster_ID=="7_1",]$subcluster_col<-"#95b391"
    dat[dat$cluster_ID=="7_2",]$subcluster_col<-"#f7f2f2"


    seurat_subclus_col<-setNames(unique(dat$subcluster_col),unique(dat$cluster_ID))

    dat$subcluster_x<-as.numeric(dat$subcluster_x)
    dat$subcluster_y<-as.numeric(dat$subcluster_y)

    label.df <- data.frame(seurat_clusters=unique(dat$seurat_clusters),label=unique(dat$seurat_clusters))
    label.df_2 <- dat %>% 
      group_by(seurat_clusters) %>% 
      summarize(umap_x = mean(umap_x), umap_y = mean(umap_y)) %>% left_join(label.df)

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=seurat_clus_col)+
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
    scale_color_manual(values=subcluster_col) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
    facet_wrap(facets=vars(seurat_clusters),ncol=2) + coord_cartesian(clip = "off") 

    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 
    ggsave(plt,file="hg38_umap.subclus.pdf")
    system("slack -F hg38_umap.subclus.pdf ryan_todo")

    hg38_atac@meta.data$cluster_col<-"NA"
    hg38_atac@meta.data[match(dat$cellid,row.names(hg38_atac@meta.data),nomatch=0),]$cluster_col<-dat$cluster_col
    hg38_atac@meta.data$subcluster_col<-"NA"
    hg38_atac@meta.data[match(dat$cellid,row.names(hg38_atac@meta.data),nomatch=0),]$subcluster_col<-dat$subcluster_col
    saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

    #########Coloring and plotting mouse data#####################

    mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
    mm10_atac$cluster_ID<-paste(mm10_atac$seurat_clusters,mm10_atac$seurat_subcluster,sep="_")

    dat<-as.data.frame(mm10_atac@meta.data)
    cell_order<-row.names(dat)
    embedding_order<-row.names(mm10_atac@reductions$umap@cell.embeddings)
    row_order<-match(cell_order,embedding_order)
    dat$umap_x<-as.numeric(mm10_atac@reductions$umap@cell.embeddings[row_order,1])
    dat$umap_y<-as.numeric(mm10_atac@reductions$umap@cell.embeddings[row_order,2])
    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]
    dat$seurat_clusters<-as.character(dat$seurat_clusters)

    dat$cluster_col<-"NA"
    dat[dat$seurat_clusters=="0",]$cluster_col<-"#A7ACD9" 
    dat[dat$seurat_clusters=="1",]$cluster_col<-"#2CA58D" 
    dat[dat$seurat_clusters=="2",]$cluster_col<-"#BA3B46" 
    dat[dat$seurat_clusters=="3",]$cluster_col<-"#76B041" 
    dat[dat$seurat_clusters=="4",]$cluster_col<-"#639cd9" 
    dat[dat$seurat_clusters=="5",]$cluster_col<-"#FCDC4D" 
    dat[dat$seurat_clusters=="6",]$cluster_col<-"#EB5A0C" 
    dat[dat$seurat_clusters=="7",]$cluster_col<-"#A7754D" 

    seurat_clus_col<-setNames(unique(dat$cluster_col),unique(dat$seurat_clusters))
    sort(unique(dat$cluster_ID))
#[1] "0_0"  "0_1"  "0_2"  "0_3"  "1_0"  "1_1"  "1_10" "1_11" "1_12" "1_13"
#[11] "1_14" "1_2"  "1_3"  "1_4"  "1_5"  "1_6"  "1_7"  "1_8"  "1_9"  "2_0"
#[21] "2_1"  "2_2"  "2_3"  "2_4"  "2_5"  "2_6"  "2_7"  "3_0"  "3_1"  "3_2"
#[31] "3_3"  "3_4"  "3_5"  "4_0"  "4_1"  "4_2"  "4_3"  "5_0"  "5_1"  "5_2"
#[41] "5_3"  "5_4"  "5_5"  "6_0"  "6_1"  "7_0"  "7_1"  "7_2"  "7_3"  "7_4"


    dat$subcluster_col<-"NA"
    dat[dat$cluster_ID=="0_0",]$subcluster_col<-"#5EC6F2"
    dat[dat$cluster_ID=="0_1",]$subcluster_col<-"#F2BC5E"
    dat[dat$cluster_ID=="0_2",]$subcluster_col<-"#5EF28B"
    dat[dat$cluster_ID=="0_3",]$subcluster_col<-"#217A3C"

    dat[dat$cluster_ID=="1_0",]$subcluster_col<-"#D46F63"
    dat[dat$cluster_ID=="1_1",]$subcluster_col<-"#68A882"
    dat[dat$cluster_ID=="1_2",]$subcluster_col<-"#496E51"
    dat[dat$cluster_ID=="1_3",]$subcluster_col<-"#CF5A52"
    dat[dat$cluster_ID=="1_4",]$subcluster_col<-"#F5CCE8"
    dat[dat$cluster_ID=="1_5",]$subcluster_col<-"#EC9DED"
    dat[dat$cluster_ID=="1_6",]$subcluster_col<-"#12355B"
    dat[dat$cluster_ID=="1_7",]$subcluster_col<-"#420039"
    dat[dat$cluster_ID=="1_8",]$subcluster_col<-"#FF6B61"
    dat[dat$cluster_ID=="1_9",]$subcluster_col<-"#6A35F0"
    dat[dat$cluster_ID=="1_10",]$subcluster_col<-"#F05135"
    dat[dat$cluster_ID=="1_11",]$subcluster_col<-"#D4F035"
    dat[dat$cluster_ID=="1_12",]$subcluster_col<-"#35F06A"
    dat[dat$cluster_ID=="1_13",]$subcluster_col<-"#D5C4FF"
    dat[dat$cluster_ID=="1_13",]$subcluster_col<-"#077827"
    dat[dat$cluster_ID=="1_14",]$subcluster_col<-"#F502B8"

    dat[dat$cluster_ID=="2_0",]$subcluster_col<-"#A9A587"
    dat[dat$cluster_ID=="2_1",]$subcluster_col<-"#6E171E"
    dat[dat$cluster_ID=="2_2",]$subcluster_col<-"#7871AA"
    dat[dat$cluster_ID=="2_3",]$subcluster_col<-"#533B4D"
    dat[dat$cluster_ID=="2_4",]$subcluster_col<-"#F63E02"
    dat[dat$cluster_ID=="2_5",]$subcluster_col<-"#92BFB1"
    dat[dat$cluster_ID=="2_6",]$subcluster_col<-"#F5E102"
    dat[dat$cluster_ID=="2_7",]$subcluster_col<-"#7A005C"

    dat[dat$cluster_ID=="3_0",]$subcluster_col<-"#47624F"
    dat[dat$cluster_ID=="3_1",]$subcluster_col<-"#A67C72"
    dat[dat$cluster_ID=="3_2",]$subcluster_col<-"#76D91A"
    dat[dat$cluster_ID=="3_3",]$subcluster_col<-"#55251D"
    dat[dat$cluster_ID=="3_4",]$subcluster_col<-"#A17C37"
    dat[dat$cluster_ID=="3_5",]$subcluster_col<-"#74AB3A"

    dat[dat$cluster_ID=="4_0",]$subcluster_col<-"#0975E8"
    dat[dat$cluster_ID=="4_1",]$subcluster_col<-"#535D69"
    dat[dat$cluster_ID=="4_2",]$subcluster_col<-"#80DED9"
    dat[dat$cluster_ID=="4_3",]$subcluster_col<-"#C6B9CD"

    dat[dat$cluster_ID=="5_0",]$subcluster_col<-"#857A4D"
    dat[dat$cluster_ID=="5_1",]$subcluster_col<-"#E8DFB5"
    dat[dat$cluster_ID=="5_2",]$subcluster_col<-"#B39615"
    dat[dat$cluster_ID=="5_3",]$subcluster_col<-"#E55812"
    dat[dat$cluster_ID=="5_4",]$subcluster_col<-"#93032E"
    dat[dat$cluster_ID=="5_5",]$subcluster_col<-"#FF5E5B"

    dat[dat$cluster_ID=="6_0",]$subcluster_col<-"#F57631"
    dat[dat$cluster_ID=="6_1",]$subcluster_col<-"#CF8259"

    dat[dat$cluster_ID=="7_0",]$subcluster_col<-"#9C8938"
    dat[dat$cluster_ID=="7_1",]$subcluster_col<-"#D4CB22"
    dat[dat$cluster_ID=="7_2",]$subcluster_col<-"#9E6D18"
    dat[dat$cluster_ID=="7_3",]$subcluster_col<-"#5C6F9E"
    dat[dat$cluster_ID=="7_4",]$subcluster_col<-"#B22CF5"

  

    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]

    dat$subcluster_x<-as.numeric(dat$subcluster_x)
    dat$subcluster_y<-as.numeric(dat$subcluster_y)
    subcluster_col<-setNames(unique(dat$subcluster_col),unique(dat$cluster_ID))

    label.df <- data.frame(seurat_clusters=unique(dat$seurat_clusters),label=unique(dat$seurat_clusters))
    label.df_2 <- dat %>% 
      group_by(seurat_clusters) %>% 
      summarize(umap_x = mean(umap_x), umap_y = mean(umap_y)) %>% left_join(label.df)

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=seurat_clus_col)+
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
    ggrepel::geom_label_repel(data = label.df_3, aes(x=subcluster_x,y=subcluster_y,label = label), max.iter=10000,direction="both",size=2,force=5, fontface='bold')+
    scale_color_manual(values=subcluster_col) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none",strip.background = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
    facet_wrap(facets=vars(seurat_clusters),ncol=2) + coord_cartesian(clip = "off") 

    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 
    ggsave(plt,file="mm10_umap.subclus.pdf")
    system("slack -F mm10_umap.subclus.pdf ryan_todo")

    mm10_atac@meta.data$cluster_col<-"NA"
    mm10_atac@meta.data[match(dat$cellid,row.names(mm10_atac@meta.data),nomatch=0),]$cluster_col<-dat$cluster_col
    mm10_atac@meta.data$subcluster_col<-"NA"
    mm10_atac@meta.data[match(dat$cellid,row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_col<-dat$subcluster_col
    saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

```
## Cicero for Coaccessible Networks

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
      atac.cds <- as.cell_data_set(object_input,group_by="cluster_ID")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      atac.cicero<-readRDS(paste(prefix,"atac_cicero_cds.Rds",sep="_"))

      genome <- seqlengths(object_input) # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      
      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      
      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      return(object_input)
  }

  hg38_atac<-readRDS("hg38_SeuratObject.Rds")
  hg38_atac<-cicero_processing(object_input=hg38_atac,prefix="hg38")
  saveRDS(hg38_atac,"hg38_SeuratObject.GA.Rds")
  hg38_atac<-readRDS("hg38_SeuratObject.GA.Rds")

  mm10_atac<-readRDS("mm10_SeuratObject.Rds")
  mm10_atac<-cicero_processing(object_input=mm10_atac,prefix="mm10")
  saveRDS(mm10_atac,"mm10_SeuratObject.GA.Rds")
  mm10_atac<-readRDS("mm10_SeuratObject.GA.Rds")
  
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

    hg38_annotation<-annotation_generation(ensdb_obj=EnsDb.Hsapiens.v86)
    mm10_annotation<-annotation_generation(ensdb_obj=EnsDb.Mmusculus.v79)

  geneactivity_processing<-function(cds_input,conns_input,prefix,gene_annotation){
      atac.cds<- annotate_cds_by_site(cds_input, gene_annotation)
      unnorm_ga <- build_gene_activity_matrix(atac.cds, conns_input)
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
  }

  #hg38
  conns<-as.data.frame(readRDS("hg38_atac_cicero_conns.Rds"))
  geneactivity_processing(cds_input=as.cell_data_set(hg38_atac,group_by="seurat_clusters"),conns_input=conns,prefix="hg38",gene_annotation=hg38_annotation)
  cicero_gene_activities<-readRDS("hg38.unnorm_GA.Rds")  #Read in unnormalized GA
  hg38_atac<-subset(hg38_atac,cells=which(colnames(hg38_atac) %in% colnames(cicero_gene_activities)))
  hg38_atac[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
  hg38_atac <- NormalizeData(object = hg38_atac,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(hg38_atac$nCount_GeneActivity))  # normalize
  saveRDS(hg38_atac,"hg38_SeuratObject.PF.Rds") #this is limited to just cells passing filters (those with cluster IDs)

  #mm10
  conns<-as.data.frame(readRDS("mm10_atac_cicero_conns.Rds"))
  geneactivity_processing(cds_input=as.cell_data_set(mm10_atac,group_by="seurat_clusters"),conns_input=conns,prefix="mm10",gene_annotation=mm10_annotation)
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
#Mouse download
cd /home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/87/14/8714d0a3-27d7-4a81-8c77-eebfd605a280/readme_mouse_10x.txt
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv 
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/trimmed_means.csv
#Downloading Mouse whole brain data 
wget https://www.dropbox.com/s/kqsy9tvsklbu7c4/allen_brain.rds?dl=0

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
system("slack -F allen_brainspan_humancortex.elbowplot.pdf ryan_todo")
brainspan <- FindNeighbors(brainspan, dims = 1:14)
brainspan <- FindClusters(brainspan, resolution = 0.5)
brainspan <- RunUMAP(brainspan, dims = 1:14)
plt<-DimPlot(brainspan, reduction = "umap",group.by=c("class_label","subclass_label"))
ggsave(plt,file="allen_brainspan_humancortex.dimplot.pdf",width=30)
system("slack -F allen_brainspan_humancortex.dimplot.pdf ryan_todo")
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
 system("slack -F allen_brainspan_mouse.elbowplot.pdf ryan_todo")
 brainspan <- FindNeighbors(brainspan, dims = 1:15)
 brainspan <- FindClusters(brainspan, resolution = 0.5)
 brainspan <- RunUMAP(brainspan, dims = 1:15)
 plt<-DimPlot(brainspan, reduction = "umap",group.by=c("class_label","subclass_label"))
 ggsave(plt,file="allen_brainspan_mouse.dimplot.pdf",width=30)
 system("slack -F allen_brainspan_mouse.dimplot.pdf ryan_todo")
 saveRDS(brainspan, file = "allen_brainspan_mouse.rds")




```

### Integration of ATAC and RNA for cell type identification

Retry mouse cluster id just straight up following https://satijalab.org/signac/articles/mouse_brain_vignette.html?

```R
library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

predict_celltype<-function(object,brainspan,prefix){
    transfer.anchors <- FindTransferAnchors(reference = brainspan,query = object,query.assay="GeneActivity",reduction = 'cca')
    saveRDS(transfer.anchors,file=paste0(prefix,".transferanchors.rds"))
    #predict labels for class and subclass
    predicted.labels.class <- TransferData(anchorset = transfer.anchors,refdata = brainspan$class_label,weight.reduction = "cca",dims = 1:30)
    predicted.labels.subclass <- TransferData(anchorset = transfer.anchors,refdata = brainspan$subclass_label,weight.reduction = "cca", dims = 1:30)
    object <- AddMetaData(object = object, metadata = predicted.labels.class)
    object <- AddMetaData(object = object, metadata = predicted.labels.subclass)
    saveRDS(object,file=paste0(prefix,"_SeuratObject.PF.Rds"))

    feat<-colnames(object@meta.data)[which(grepl("prediction.score",colnames(object@meta.data)))]
    feat<-feat[feat !="prediction.score.max"]
    plt<-FeaturePlot(object,features=feat,order=T)#plot feature plots
    ggsave(plt,file=paste0(prefix,"_predicted.umap.png"),width=20,height=30)
    system(paste0("slack -F ",prefix,"_predicted.umap.png ryan_todo"))

    #plot merged cluster heatmap
    #Generate Heatmap of TF ChromVar score and TF modules
    pred<-object@meta.data[,feat]
    colnames(pred)<-substring(colnames(pred),18) #remove "prediction.score" for readability
    #Combine over subclusters
    pred<-split(pred,paste(object$seurat_clusters,object$seurat_subcluster,sep="_")) #group by rows to seurat clusters
    pred<-lapply(pred,function(x) apply(x,2,mean)) #take average across group
    pred<-do.call("rbind",pred) #condense to smaller data frame
    pred<-pred[!endsWith(row.names(pred),"NA"),] #remove doublets not assigned a subcluster

    plt1<-Heatmap(pred,column_split=c(rep("class",3),rep("subclass",20)),clustering_distance_rows="maximum",clustering_distance_columns="maximum")
    pdf(paste0(prefix,".complexheatmap.pdf"),height=20)
    plt1
    dev.off()
    system(paste0("slack -F ",prefix,".complexheatmap.pdf ryan_todo"))

}

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
predict_celltype(object=hg38_atac,brainspan=brainspan.,prefix="hg38")

#Mouse using smaller data set that is whole brain
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")
brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse/allen_brainspan_mouse.rds")
brainspan. <- FindVariableFeatures(object = brainspan.,nfeatures = 5000)
predict_celltype(object=mm10_atac,brainspan=brainspan.,prefix="mm10")


```

### Add TF Motif Usage through ChromVAR

### ChromVar for Transcription Factor Motifs


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
  register(MulticoreParam(10))

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.PF.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.PF.Rds")

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

  motif.matrix.mm10 <- CreateMotifMatrix(
    features = granges(mm10_atac[["peaks"]]),
    pwm = pfm,
    genome = 'mm10',
    use.counts = FALSE)

  # Create a new Mofif object to store the results
  motif.hg38 <- CreateMotifObject(
    data = motif.matrix.hg38,
    pwm = pfm)

  motif.mm10 <- CreateMotifObject(
    data = motif.matrix.mm10,
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


#Clusters to test
hg38_atac<-subset(hg38_atac,cells=which(!endsWith(hg38_atac@meta.data$cluster_ID,"NA")))
cluster_to_test<-unique(hg38_atac$cluster_ID)
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
    obj=hg38_atac,
    group="cluster_ID",
    assay.="GeneActivity",
    mc.cores=n.cores)

da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
write.table(da_ga_df,file="hg38.onevrest.da_ga.txt",sep="\t",col.names=T,row.names=T,quote=F)

da_ga<-list() #set up an empty list for looping through

n.cores=10 #Perform parallel application of DA test
da_ga<-mclapply(
    cluster_to_test,
    FUN=da_one_v_rest,
    obj=hg38_atac,
    group="cluster_ID",
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
da_ga<-da_ga[!endsWith(da_ga$enriched_group,"NA"),]

da_ga$label<-""
for (x in unique(da_ga$enriched_group)){
selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:3))$da_region
da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_ga<-as.data.frame(t(as.data.frame(hg38_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,hg38_atac$cluster_ID) #group by rows to seurat clusters
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame

sum_ga<-t(scale(sum_ga))
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

#cluster by all marker genes
sum_da_dend <- t(sum_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="hg38.geneactivity.dend.rds") 


sum_ga<-sum_ga[row.names(sum_ga) %in% unique(da_ga$label),]

annot<-hg38_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$cluster_ID),]
sum_ga_plot<-t(sum_ga)

annot_clus_col<-annot[!duplicated(annot$cluster_ID),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
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
system("slack -F hg38.geneactivity.heatmap.pdf ryan_todo")

saveRDS(hg38_atac,"hg38_SeuratObject.PF.Rds")

###########Plotting heatmap using brainspan given markers#############################

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
sum_da_dend<-readRDS(file="hg38.geneactivity.dend.rds")
dat_ga<-as.data.frame(t(as.data.frame(hg38_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,hg38_atac$cluster_ID) #group by rows to seurat clusters 
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame
sum_ga<-t(sum_ga)
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

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
    c("Mark","L1-6","Endothelia","KDR"),
    c("Mark","L1-6","Astrocyte","SLC1A2")
    ))
colnames(markers_limited)<-c("celltype","layer","marker_1","marker_2") #i decided to ignore layer markers so made them all L1-6

sum_ga_sub<-sum_ga[match(markers_limited$marker_2,row.names(sum_ga),nomatch=0),]
markers_limited<-markers_limited[match(row.names(sum_ga_sub),markers_limited$marker_2,nomatch=0),]

sum_ga_sub<-t(scale(t(sum_ga_sub),center=F,scale=T))

annot<-hg38_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$cluster_ID),]
annot_clus_col<-annot[!duplicated(annot$cluster_ID),]

top_ha<-columnAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
                        ),
                    show_legend = c(TRUE, TRUE, FALSE))


colfun=colorRamp2(quantile(unlist(sum_ga_plot), probs=c(0.5,0.90,0.95)),magma(3))

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
system("slack -F hg38.geneactivity.markers.heatmap.pdf ryan_todo")



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
sum_tf<-split(dat_tf,hg38_atac$cluster_ID) #group by rows to seurat clusters
sum_tf<-lapply(sum_tf,function(x) apply(x,2,mean)) #take average across group
sum_tf<-do.call("rbind",sum_tf) #condense to smaller data frame

sum_tf<-t(scale(sum_tf))
sum_tf<-sum_tf[,!endsWith(colnames(sum_tf),"NA")] #remove NA (doublet cells)

#clustered by all marker genes ga
sum_da_dend<-readRDS(file="hg38.geneactivity.dend.rds") 


sum_tf<-sum_tf[row.names(sum_tf) %in% unique(da_tf[da_tf$label!="",]$da_region),]
row.names(sum_tf)<-da_tf[match(row.names(sum_tf),da_tf$da_region,nomatch=0),]$tf_name
annot<-hg38_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
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
system("slack -F hg38.tf.heatmap.pdf ryan_todo")

#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)

motif_order<-names(hg38_atac@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(hg38_atac@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = hg38_atac,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())

ggsave(plt,file="hg38.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
system("slack -F hg38.tf.heatmap.motif.pdf ryan_todo")


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
#mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Unknown",]$celltype_col<-"#ff6633"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Endo",]$celltype_col<-"#A52A2A"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Gran",]$celltype_col<-"#4B0082"

saveRDS(mm10_atac,"mm10_SeuratObject.PF.Rds")

#Clusters to test
cluster_to_test<-unique(mm10_atac$cluster_ID)
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
    group="cluster_ID",
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
    group="cluster_ID",
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
sum_ga<-split(dat_ga,mm10_atac$cluster_ID) #group by rows to seurat clusters
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame

sum_ga<-t(scale(sum_ga))
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

#cluster by all marker genes
sum_da_dend <- t(sum_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="mm10.geneactivity.dend.rds") 


sum_ga<-sum_ga[row.names(sum_ga) %in% unique(da_ga$label),]

annot<-mm10_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$cluster_ID),]
sum_ga_plot<-t(sum_ga)

annot_clus_col<-annot[!duplicated(annot$cluster_ID),]

side_ha<-rowAnnotation(df= data.frame(celltype=annot$celltype, cluster=annot$seurat_clusters, subcluster=annot$cluster_ID),
                col=list(
                    celltype=setNames(unique(annot$celltype_col),unique(annot$celltype)),
                    cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters))),
                    subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID) #due to nonunique colors present
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
system("slack -F mm10.geneactivity.heatmap.pdf ryan_todo")
############################################################
#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)############
sum_da_dend<-readRDS(file="mm10.geneactivity.dend.rds")
dat_ga<-as.data.frame(t(as.data.frame(mm10_atac[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,mm10_atac$cluster_ID) #group by rows to seurat clusters 
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame
sum_ga<-t(sum_ga)
sum_ga<-sum_ga[,!endsWith(colnames(sum_ga),"NA")] #remove NA (doublet cells)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
      sep="", collapse=" ")
}

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
    c("Mark","L1-6","Endothelia","KDR"),
    c("Mark","L1-6","PyrCA1","Crym"),
    c("Mark","L1-6","PyrCA1","Neurod6"),
    c("Mark","L1-6","PyrCA1","Gria1"),
    c("Mark","L1-6","PyrCA1","Cpne6"),
    c("Mark","L1-6","PyrCA1","Tspan13"),
    c("Mark","L1-6","PyrCA1","Grm5")
    ))

colnames(markers_limited)<-c("celltype","layer","marker_1","marker_2")

markers_limited$marker_2<-unlist(lapply(markers_limited$marker_2,simpleCap)) # correct case
sum_ga_sub<-sum_ga[match(markers_limited$marker_2,row.names(sum_ga),nomatch=0),]
markers_limited<-markers_limited[match(row.names(sum_ga_sub),markers_limited$marker_2,nomatch=0),]

sum_ga_sub<-t(scale(t(sum_ga_sub),center=F,scale=T))


annot<-mm10_atac@meta.data[,c("celltype","cluster_ID","subcluster_col","cluster_col","seurat_clusters","seurat_subcluster","celltype_col")]
annot<-annot[!(annot$subcluster_col=="NA"),]
annot<-annot[!duplicated(annot$cluster_ID),]
annot<-annot[annot$cluster_ID %in% colnames(sum_ga),]
annot<-annot[match(colnames(sum_ga),annot$cluster_ID,nomatch=0),]
annot_clus_col<-annot[!duplicated(annot$cluster_ID),]
annot_clus_col<-annot_clus_col[complete.cases(annot_clus_col),]
sum_ga_sub<-sum_ga_sub[,colnames(sum_ga_sub) %in% annot$cluster_ID]

celltype=setNames(unique(annot$celltype_col),unique(annot$celltype))
cluster=setNames(unique(annot$cluster_col),unique(as.character(annot$seurat_clusters)))
subcluster=setNames(annot_clus_col$subcluster_col,annot_clus_col$cluster_ID)
celltype<-celltype[!is.na(celltype)]
cluster<-cluster[!is.na(cluster)]
subcluster<-subcluster[!is.na(subcluster)]
cluster<-cluster[!is.na(names(cluster))]
top_ha<-columnAnnotation(df= data.frame(cluster=annot$seurat_clusters, celltype=annot$celltype, subcluster=annot$cluster_ID), 
                col=list(celltype=celltype,cluster=cluster,subcluster=subcluster),
                    show_legend = c(TRUE, TRUE,FALSE))  

762.6284 pt 113.3457 pt 
176.6231 pt 239.3815 pt
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
system("slack -F mm10.geneactivity.markers.heatmap.pdf ryan_todo")



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
system("slack -F mm10.tf.heatmap.pdf ryan_todo")

#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)

motif_order<-names(mm10_atac@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(mm10_atac@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = mm10_atac,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())

ggsave(plt,file="mm10.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
system("slack -F mm10.tf.heatmap.motif.pdf ryan_todo")

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
system("slack -F hg38_3d.umap.tsv ryan_todo")


mm10_out<-cbind(mm10_atac$celltype,
    colnames(mm10_atac),
    as.data.frame(mm10_atac@reductions$umap3d@cell.embeddings),
    mm10_atac$celltype_col)
write.table(mm10_out,"mm10_3d.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
system("slack -F mm10_3d.umap.tsv ryan_todo")

###Repeat with subcluster coloring
hg38_out<-cbind(hg38_atac$cluster_ID,
    colnames(hg38_atac),
    as.data.frame(hg38_atac@reductions$umap3d@cell.embeddings),
    hg38_atac$subcluster_col)
write.table(hg38_out,"hg38_3d.subclus.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
system("slack -F hg38_3d.subclus.umap.tsv ryan_todo")


###Repeat with subcluster coloring
mm10_out<-cbind(mm10_atac$cluster_ID,
    colnames(mm10_atac),
    as.data.frame(mm10_atac@reductions$umap3d@cell.embeddings),
    mm10_atac$subcluster_col)
write.table(mm10_out,"mm10_3d.subclus.umap.tsv",sep="\t",col.names=F,row.names=F,quote=F)
system("slack -F mm10_3d.subclus.umap.tsv ryan_todo")
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


<!---

## Cortex Layering 

### Generating Monocle trajectory of excitatory neurons


```python
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#Read in data and modify to monocle CDS file
#read in RDS file.
hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")

#Subset hg38 to just excitatory neurons in the cortex and apply a trajectory, trying to get layer information here
hg38_ex_neurons<-subset(hg38_atac,idents="2") #cluster 2 is excitatory neurons


monocle_processing<-function(prefix, seurat_input){
    atac.cds <- as.cell_data_set(seurat_input)
    atac.cds <- cluster_cells(cds = atac.cds, reduction_method = "UMAP") 
    #Read in cds from cicero processing earlier and continue processing
    atac.cds<- learn_graph(atac.cds, 
                           use_partition = F, 
                           learn_graph_control=list(
                               minimal_branch_len=10,
                               orthogonal_proj_tip=F,
                               prune_graph=T))
    #plot to see nodes for anchoring
    plt1<-plot_cells(
                    cds = atac.cds,
                    show_trajectory_graph = TRUE,
                    label_leaves=T,
                    label_branch_points=F,
                    label_roots=T)
    #Also make a plot of just node names for easier identification
    root_nodes<-as.data.frame(t(atac.cds@principal_graph_aux$UMAP$dp_mst))
    root_nodes$label<-row.names(root_nodes)
    plt2<-ggplot(
        root_nodes,
        aes(x=UMAP_1,y=UMAP_2))+
        geom_text(aes(label=label),size=3)+
        theme_bw()
    plt<-(plt1+plt2)
    ggsave(plt,file=paste(prefix,"trajectory.pdf",sep="_"),width=20)
    system(paste0("slack -F ",paste(prefix,"trajectory.pdf",sep="_")," ryan_todo"))
    return(atac.cds)
}


hg38_ex_neurons.cicero<-monocle_processing(seurat_input=hg38_ex_neurons,prefix="hg38_excNeurons")

#Setting starting node as Y_21
hg38_ex_neurons.cicero<-order_cells(hg38_ex_neurons.cicero,reduction_method="UMAP",root_pr_nodes="Y_21")


#variance over pseudotime
#pr_graph_test <- principalGraphTest(hg38_ex_neurons.cicero, k=6, cores=1)
#dplyr::add_rownames(pr_graph_test) %>%
#    dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval)) %>% head(3)

#Now replotting with pseudotime
pdf("hg38_excNeurons.trajectory.pseudotime.pdf")
plot_cells(
  cds = hg38_ex_neurons.cicero,
  show_trajectory_graph = TRUE,
color_cells_by = "pseudotime"
)
dev.off()
system("slack -F hg38_excNeurons.trajectory.pseudotime.pdf ryan_todo")

saveRDS(hg38_ex_neurons.cicero,"hg38_excNeurons.monocle.cds.Rds")

#Append pseudotime to meta data of seurat object
hg38_excNeurons_pseudotime<-as.data.frame(hg38_ex_neurons.cicero@principal_graph_aux@listData$UMAP$pseudotime)
colnames(hg38_excNeurons_pseudotime)<-c("pseudotime")
hg38_excNeurons_pseudotime$cellID<-row.names(hg38_excNeurons_pseudotime)

hg38_ex_neurons$pseudotime<-hg38_excNeurons_pseudotime[match(hg38_excNeurons_pseudotime$cellID,hg38_ex_neurons$cellID),]$pseudotime

saveRDS(hg38_ex_neurons,"hg38_excNeurons.SeuratObject.Rds")

```