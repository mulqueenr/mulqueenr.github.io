---
title: sciDROP
layout: analysis
permalink: /scidrop/
category: alternative
---


## Processing for sciDROP

![sciDROP Overview](/assets/images/sciDROP.png){:width="80%"}


This notebook details the processing of the "20K" and "70K" loaded mouse brain and human cortex samples. It begins with scitools wrapper functions for intial alignment to a concatenated mouse and human genome, following with splitting of reads and realignment to separate human and mouse genomes. It then follows the established scitools formation of a counts matrix and Signac processing.

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 

## Combine 10% Sampling of 20K Experiment with current run data

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 

## Calculate collision rate from barnyard experiment

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 

## Split out species from barnyard experiments

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 


### Tabix fragment file generation

Tabix file format is a tab separated multicolumn data structure.

| Column Number | Name | Description |
|:--------|:-------:|:--------|
|1 |chrom |  Reference genome chromosome of fragment |
|2 |chromStart | Adjusted start position of fragment on chromosome. |
|3 |chromEnd   | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
|4 |barcode | The 10x (or sci) cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment. |
|5 |duplicateCount |The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeq instrument are excluded from this count. |

{% capture summary %} Code {% endcapture %} {% capture details %}  

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

{% endcapture %} {% include details.html %} 


# sciDROP Full Processing

### Generating Seurat Objects

Using R v4.0.0 and Signac v1.0

{% capture summary %} Code {% endcapture %} {% capture details %}  


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

{% endcapture %} {% include details.html %} 


### Perform Scrublet on Data to Ensure Single-cells

Code from tutorial here.[https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb]

{% capture summary %} Code {% endcapture %} {% capture details %}  
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
{% endcapture %} {% include details.html %} 


### Add library complexity data to RDS files.

{% capture summary %} Code {% endcapture %} {% capture details %}  

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

{% endcapture %} {% include details.html %} 

## Running Initial Clustering of Cells
Using CisTopic for Dimensionality reduction and UMAP for projection.

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
{% endcapture %} {% include details.html %} 

## Correcting for systematic bias with harmony

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
```
{% endcapture %} {% include details.html %} 

### Subclustering 

Going to exclude cells in subclustering that are identified by scrublet as potential doublets.

{% capture summary %} Code {% endcapture %} {% capture details %}  


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
    atac_sub<-subset(atac_sub,subset=pcr_idx %in% c("CAGAGAGG","CTCTCTAC"))
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
    sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10,22:28),nCores=8,addModels=FALSE)
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


setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

#Read in data and modify to monocle CDS file
#read in RDS file.
hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
dir.create("subcluster")

#UMAP Projection and clustering on selected cistopic model
clustering_loop<-function(topicmodel_list.=topicmodel_list,object_input=hg38_atac,celltype.x,topiccount_list.=topic_count_list){
    #set up subset object again
    atac_sub<-subset(object_input,subset=(seurat_clusters==celltype.x & predicted_doublets=="False")) 
    #select_topic
    models_input<-readRDS(paste("./subcluster/",species,celltype.x,".CisTopicObject.Rds",sep="_"))
    cisTopicObject<-cisTopic::selectModel(models_input,select=topiccount_list.[celltype.x],keepModels=F)
    
    #perform UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
    dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
    row.names(dims)<-colnames(topic_df)
    colnames(dims)<-c("x","y")
    dims$cellID<-row.names(dims)
    dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")

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
    umap_dims<-as.data.frame(as.matrix(dims[2:3]))
    colnames(umap_dims)<-c("UMAP_1","UMAP_2")
    row.names(umap_dims)<-dims$cellID
    cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
    atac_sub@reductions$cistopic<-cistopic_obj
    atac_sub@reductions$umap<-cistopic_umap
    #finally recluster
    n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic"))
    #Clustering with multiple resolutions to account for different celltype complexities
    atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics)
    atac_sub <- FindClusters(object = atac_sub,resolution=0.1)
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.2)
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.5)
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.9)
    
    saveRDS(atac_sub,paste("./subcluster/",species,celltype.x,"SeuratObject.Rds",sep="_"))
    plt<-DimPlot(atac_sub,group.by=c('peaks_snn_res.0.1','peaks_snn_res.0.2','peaks_snn_res.0.5','peaks_snn_res.0.9'))
    ggsave(plt,file=paste("./subcluster/",species,celltype.x,"clustering.pdf",sep="_"))
}

####################################
### Processing ###
#hg38
    #Set up variables
    species="hg38"
    celltype_list<-unique(hg38_atac$seurat_clusters)

    #Running cistopic subclustering on all identified cell types
    #lapply(celltype_list, function(i){atac_sub<-subset(hg38_atac,subset=seurat_clusters==i); cistopic_generation(x=atac_sub,celltype.x=i,species="hg38")})

    topicmodel_list<-paste("./subcluster/",species,celltype_list,".CisTopicObject.Rds",sep="_")

    #determine model count to use for each cell type
    for (i in paste("./subcluster/",species,celltype_list,"_model_selection.pdf",sep="_")){system(paste0("slack -F ",i," ryan_todo"))}

    #selecting topics based on derivative, making a named vector, note some of these seem to suggest a topic count over 28, but were artificially capped
    topic_count_list<-setNames(c(27,28,25,26,28,28,26),celltype_list)
    names(topic_count_list)<-celltype_list

    #Running clustering loop
    for (i in celltype_list){clustering_loop(object_input=hg38_atac,celltype.x=i)}

    #selecting resolution by plots
    for (i in celltype_list){system(paste0("slack -F ",paste("./subcluster/",species,i,"clustering.pdf",sep="_")," ryan_todo"))}
    
    #adding all subclustering info back into main RDS object
    hg38_atac$seurat_subcluster<-"NA"
    hg38_atac$subcluster_x<-"NA"
    hg38_atac$subcluster_y<-"NA"

    #Assign clustering resolution based on clustering.pdf output
    resolution_list<-setNames(c("peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.2","peaks_snn_res.0.2",
        "peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.2"),celltype_list)
    cell_order<-row.names(hg38_atac@meta.data)

    for(celltype.x in celltype_list){
        atac_sub<-readRDS(paste("./subcluster/",species,celltype.x,"SeuratObject.Rds",sep="_"))
        embedding_order<-row.names(atac_sub@reductions$umap@cell.embeddings)
        row_order<-match(cell_order,embedding_order,nomatch=0)
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$seurat_subcluster<-as.character(unlist(atac_sub@meta.data[which(colnames(atac_sub@meta.data)==resolution_list[celltype.x])]))
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$subcluster_x<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,1])
        hg38_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(hg38_atac@meta.data),nomatch=0),]$subcluster_y<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,2])
    }

    as.data.frame(hg38_atac@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
    #na values are doublets  
    saveRDS(hg38_atac,"hg38_SeuratObject.Rds")

    #mm10
    #Set up variables
    species="mm10"
    celltype_list<-unique(mm10_atac$seurat_clusters)

    #Running cistopic subclustering on all identified cell types
    #####DOING CUSTERS 1:4 to START, CLEARING MEMORY THEN DOING AGAIN
    #lapply(celltype_list[4:length(celltype_list)], function(i){atac_sub<-subset(mm10_atac,subset=seurat_clusters==i); cistopic_generation(x=atac_sub,celltype.x=i,species="mm10")})
    
    #determine model count to use for each cell type
    for (i in paste("./subcluster/",species,celltype_list,"_model_selection.pdf",sep="_")){system(paste0("slack -F ",i," ryan_todo"))}
    topicmodel_list<-paste("./subcluster/",species,celltype_list,".CisTopicObject.Rds",sep="_")

    #selecting topics based on derivative, making a named vector 
    topic_count_list<-setNames(c(24,24,24,27,25,27,24,27,28),celltype_list)
    names(topic_count_list)<-celltype_list

    #Running clustering loop
    for (i in celltype_list){clustering_loop(object_input=mm10_atac,celltype.x=i)}

    #selecting resolution by plots 
    for (i in celltype_list){system(paste0("slack -F ",paste("./subcluster/",species,i,"clustering.pdf",sep="_")," ryan_todo"))}

    #adding all subclustering info back into main RDS object
    mm10_atac$seurat_subcluster<-"NA"
    mm10_atac$subcluster_x<-"NA"
    mm10_atac$subcluster_y<-"NA"

    #Assign clustering resolution based on clustering.pdf output
    resolution_list<-setNames(c("peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1",
        "peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1","peaks_snn_res.0.1"),celltype_list)
    cell_order<-row.names(mm10_atac@meta.data)

    for(celltype.x in celltype_list){
        atac_sub<-readRDS(paste("./subcluster/",species,celltype.x,"SeuratObject.Rds",sep="_"))
        embedding_order<-row.names(atac_sub@reductions$umap@cell.embeddings)
        row_order<-match(cell_order,embedding_order,nomatch=0)
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$seurat_subcluster<-as.character(unlist(atac_sub@meta.data[which(colnames(atac_sub@meta.data)==resolution_list[celltype.x])]))
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_x<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,1])
        mm10_atac@meta.data[match(row.names(atac_sub@meta.data),row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_y<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,2])
    }

    as.data.frame(mm10_atac@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
    #na values are doublets  
    saveRDS(mm10_atac,"mm10_SeuratObject.Rds")

```

{% endcapture %} {% include details.html %} 


### Recoloring subclusters

{% capture summary %} Code {% endcapture %} {% capture details %}  


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
    dat<-dat[dat$predicted_doublets=="False",]
    dat$seurat_clusters<-as.character(dat$seurat_clusters)
   
    dat$cluster_col<-"NULL"
    dat[dat$seurat_clusters=="0",]$cluster_col<-"#1f78b4"
    dat[dat$seurat_clusters=="1",]$cluster_col<-"#b2df8a"
    dat[dat$seurat_clusters=="2",]$cluster_col<-"#FDE725"
    dat[dat$seurat_clusters=="3",]$cluster_col<-"#E31A1C"
    dat[dat$seurat_clusters=="4",]$cluster_col<-"#6a3d9a"
    dat[dat$seurat_clusters=="5",]$cluster_col<-"#B15928"
    dat[dat$seurat_clusters=="6",]$cluster_col<-"#ff7f00"

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
    dat[dat$cluster_ID=="1_2",]$subcluster_col<-"#ED1169"
    dat[dat$cluster_ID=="1_3",]$subcluster_col<-"#1C3B1A"

    dat[dat$cluster_ID=="2_0",]$subcluster_col<-"#B38A1B"
    dat[dat$cluster_ID=="2_1",]$subcluster_col<-"#6B0E44"
    dat[dat$cluster_ID=="2_2",]$subcluster_col<-"#E0BB55"

    dat[dat$cluster_ID=="3_0",]$subcluster_col<-"#FAA916"
    dat[dat$cluster_ID=="3_1",]$subcluster_col<-"#5E503F"
    dat[dat$cluster_ID=="3_2",]$subcluster_col<-"#EDFFAB"
    dat[dat$cluster_ID=="3_3",]$subcluster_col<-"#FFE1EA"
    dat[dat$cluster_ID=="3_4",]$subcluster_col<-"#E952DE"
    dat[dat$cluster_ID=="3_5",]$subcluster_col<-"#848FA5"
    dat[dat$cluster_ID=="3_6",]$subcluster_col<-"#832161"
    dat[dat$cluster_ID=="3_7",]$subcluster_col<-"#B0A3D4"

    dat[dat$cluster_ID=="4_0",]$subcluster_col<-"#DBAD6A"
    dat[dat$cluster_ID=="4_1",]$subcluster_col<-"#007EA7"
    dat[dat$cluster_ID=="4_2",]$subcluster_col<-"#C191A1"
    dat[dat$cluster_ID=="4_3",]$subcluster_col<-"#02010A"

    dat[dat$cluster_ID=="5_0",]$subcluster_col<-"#230903"
    dat[dat$cluster_ID=="5_1",]$subcluster_col<-"#FFBFB7"
    dat[dat$cluster_ID=="5_2",]$subcluster_col<-"#FFD447"
    dat[dat$cluster_ID=="5_3",]$subcluster_col<-"#493657"
    dat[dat$cluster_ID=="5_4",]$subcluster_col<-"#B1B695"

    dat[dat$cluster_ID=="6_0",]$subcluster_col<-"#80552B"

    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]

    seurat_subclus_col<-setNames(unique(dat$subcluster_col),unique(dat$cluster_ID))

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=seurat_clus_col)+
    ggtitle("hg38")+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "bottom")

    plt_list<-ggplot(dat,aes(x=as.numeric(subcluster_x),y=as.numeric(subcluster_y),color=cluster_ID))+
    geom_point(alpha=0.05,size=0.5,shape=16)+scale_color_manual(values=seurat_subclus_col)+
    theme_bw()+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none")+
    facet_wrap(facets=vars(seurat_clusters),ncol=2)
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
    dat<-dat[dat$predicted_doublets=="False",]
    dat$seurat_clusters<-as.character(dat$seurat_clusters)

    dat$cluster_col<-"NA"
    dat[dat$seurat_clusters=="0",]$cluster_col<-"#A7ACD9" #granule
    dat[dat$seurat_clusters=="1",]$cluster_col<-"#2CA58D" #granule
    dat[dat$seurat_clusters=="2",]$cluster_col<-"#BA3B46" #inhib
    dat[dat$seurat_clusters=="3",]$cluster_col<-"#76B041" #opc and oligo
    dat[dat$seurat_clusters=="4",]$cluster_col<-"#639cd9" #exN
    dat[dat$seurat_clusters=="5",]$cluster_col<-"#FCDC4D" #astro
    dat[dat$seurat_clusters=="6",]$cluster_col<-"#EB5A0C" #inhib
    dat[dat$seurat_clusters=="7",]$cluster_col<-"#A7754D" #micro
    dat[dat$seurat_clusters=="8",]$cluster_col<-"#342056" #exN

    seurat_clus_col<-setNames(unique(dat$cluster_col),unique(dat$seurat_clusters))
    sort(unique(dat$cluster_ID))
# [1] "0_0" "1_0" "1_1" "1_2" "1_3" "1_4" "1_5" "1_6" "1_7" "1_8" "2_0" "2_1"
#[13] "2_2" "2_3" "2_4" "2_5" "3_0" "3_1" "3_2" "3_3" "4_0" "4_1" "4_2" "4_3"
#[25] "4_4" "4_5" "5_0" "5_1" "5_2" "5_3" "5_4" "5_5" "5_6" "6_0" "6_1" "7_0"
#[37] "7_1" "7_2" "7_3" "7_4" "8_0"


    dat$subcluster_col<-"NA"
    dat[dat$cluster_ID=="0_0",]$subcluster_col<-"#5EC6F2"
    
    dat[dat$cluster_ID=="1_0",]$subcluster_col<-"#D46F63"
    dat[dat$cluster_ID=="1_1",]$subcluster_col<-"#68A882"
    dat[dat$cluster_ID=="1_2",]$subcluster_col<-"#496E51"
    dat[dat$cluster_ID=="1_3",]$subcluster_col<-"#CF5A52"
    dat[dat$cluster_ID=="1_4",]$subcluster_col<-"#F5CCE8"
    dat[dat$cluster_ID=="1_5",]$subcluster_col<-"#EC9DED"
    dat[dat$cluster_ID=="1_6",]$subcluster_col<-"#12355B"
    dat[dat$cluster_ID=="1_7",]$subcluster_col<-"#420039"
    dat[dat$cluster_ID=="1_8",]$subcluster_col<-"#FF6B61"

    dat[dat$cluster_ID=="2_0",]$subcluster_col<-"#A9A587"
    dat[dat$cluster_ID=="2_1",]$subcluster_col<-"#6E171E"
    dat[dat$cluster_ID=="2_2",]$subcluster_col<-"#7871AA"
    dat[dat$cluster_ID=="2_3",]$subcluster_col<-"#533B4D"
    dat[dat$cluster_ID=="2_4",]$subcluster_col<-"#F63E02"
    dat[dat$cluster_ID=="2_5",]$subcluster_col<-"#92BFB1"

    dat[dat$cluster_ID=="3_0",]$subcluster_col<-"#47624F"
    dat[dat$cluster_ID=="3_1",]$subcluster_col<-"#8BA672"
    dat[dat$cluster_ID=="3_2",]$subcluster_col<-"#76D91A"
    dat[dat$cluster_ID=="3_3",]$subcluster_col<-"#55251D"

    dat[dat$cluster_ID=="4_0",]$subcluster_col<-"#0975E8"
    dat[dat$cluster_ID=="4_1",]$subcluster_col<-"#CAD5E0"
    dat[dat$cluster_ID=="4_2",]$subcluster_col<-"#535D69"
    dat[dat$cluster_ID=="4_3",]$subcluster_col<-"#AEECEF"
    dat[dat$cluster_ID=="4_4",]$subcluster_col<-"#80DED9"
    dat[dat$cluster_ID=="4_5",]$subcluster_col<-"#C6B9CD"

    dat[dat$cluster_ID=="5_0",]$subcluster_col<-"#857A4D"
    dat[dat$cluster_ID=="5_1",]$subcluster_col<-"#E8DFB5"
    dat[dat$cluster_ID=="5_2",]$subcluster_col<-"#B39615"
    dat[dat$cluster_ID=="5_3",]$subcluster_col<-"#E55812"
    dat[dat$cluster_ID=="5_4",]$subcluster_col<-"#93032E"
    dat[dat$cluster_ID=="5_5",]$subcluster_col<-"#FF5E5B"
    dat[dat$cluster_ID=="5_6",]$subcluster_col<-"#5FBB97"

    dat[dat$cluster_ID=="6_0",]$subcluster_col<-"#F57631"
    dat[dat$cluster_ID=="6_1",]$subcluster_col<-"#CF8259"

    dat[dat$cluster_ID=="7_0",]$subcluster_col<-"#9C8938"
    dat[dat$cluster_ID=="7_1",]$subcluster_col<-"#D4CB22"
    dat[dat$cluster_ID=="7_2",]$subcluster_col<-"#9E6D18"
    dat[dat$cluster_ID=="7_3",]$subcluster_col<-"#5C6F9E"
    dat[dat$cluster_ID=="7_4",]$subcluster_col<-"#B22CF5"

    dat[dat$cluster_ID=="8_0",]$subcluster_col<-"#6A35F0"
  

    dat<-dat[!endsWith(dat$cluster_ID,"NA"),]

    seurat_subclus_col<-setNames(unique(dat$subcluster_col),unique(dat$cluster_ID))

    plt1<-ggplot(dat,aes(x=umap_x,y=umap_y,color=seurat_clusters))+
    geom_point(alpha=0.1,size=0.5,shape=16)+
    theme_bw()+scale_color_manual(values=seurat_clus_col)+
    ggtitle("mm10")+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "bottom")


    plt_list<-ggplot(dat,aes(x=as.numeric(subcluster_x),y=as.numeric(subcluster_y),color=cluster_ID))+
    geom_point(alpha=0.05,size=0.5,shape=16)+scale_color_manual(values=seurat_subclus_col)+
    theme_bw()+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none")+
    facet_wrap(facets=vars(seurat_clusters),ncol=2)
    plt<-plt1+plt_list+plot_layout(width=c(8,3)) 

    ggsave(plt,file="mm10_umap.subclus.pdf")
    system("slack -F mm10_umap.subclus.pdf ryan_todo")

    mm10_atac@meta.data$cluster_col<-"NA"
    mm10_atac@meta.data[match(dat$cellid,row.names(mm10_atac@meta.data),nomatch=0),]$cluster_col<-dat$cluster_col
    mm10_atac@meta.data$subcluster_col<-"NA"
    mm10_atac@meta.data[match(dat$cellid,row.names(mm10_atac@meta.data),nomatch=0),]$subcluster_col<-dat$subcluster_col
    saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

```
{% endcapture %} {% include details.html %} 


## Cicero for Coaccessible Networks

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
      atac.cds <- as.cell_data_set(object_input,group_by="seurat_clusters")

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
  hg38_atac[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
  hg38_atac <- NormalizeData(object = hg38_atac,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(hg38_atac$nCount_peaks))  # normalize
  saveRDS(hg38_atac,"hg38_SeuratObject.Rds")


  #mm10
  conns<-as.data.frame(readRDS("mm10_atac_cicero_conns.Rds"))
  geneactivity_processing(cds_input=as.cell_data_set(mm10_atac,group_by="seurat_clusters"),conns_input=conns,prefix="mm10",gene_annotation=mm10_annotation)
  cicero_gene_activities<-readRDS("mm10.unnorm_GA.Rds")  #Read in unnormalized GA
  cicero_gene_activities<-cicero_gene_activities[2:nrow(cicero_gene_activities),] #first feature is empy
  mm10_atac[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
  mm10_atac <- NormalizeData(object = mm10_atac,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(mm10_atac$nCount_peaks))  # normalize
  saveRDS(mm10_atac,"mm10_SeuratObject.Rds")


```
{% endcapture %} {% include details.html %} 

## Public Data RNA Comparison
### Download data from Allen Brain-span
For human: https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
For mouse: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x

{% capture summary %} Code {% endcapture %} {% capture details %}  
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
{% endcapture %} {% include details.html %} 


### Process Data into Seurat Object
Following https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
{% capture summary %} Code {% endcapture %} {% capture details %}  
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
{% endcapture %} {% include details.html %} 

### Integration of ATAC and RNA for cell type identification

Retry mouse cluster id just straight up following https://satijalab.org/signac/articles/mouse_brain_vignette.html?
{% capture summary %} Code {% endcapture %} {% capture details %}  

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
    saveRDS(object,file=paste0(prefix,"_SeuratObject.Rds"))

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

hg38_atac<-readRDS("hg38_SeuratObject.Rds")
brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
predict_celltype(object=hg38_atac,brainspan=brainspan.,prefix="hg38")

#Mouse using smaller data set that is whole brain
mm10_atac<-readRDS("mm10_SeuratObject.Rds")
brainspan. <- readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_mouse/allen_brainspan_mouse.rds")
brainspan. <- FindVariableFeatures(object = brainspan.,nfeatures = 5000)
predict_celltype(object=mm10_atac,brainspan=brainspan.,prefix="mm10")


```

{% endcapture %} {% include details.html %} 

### Add TF Motif Usage through ChromVAR

### ChromVar for Transcription Factor Motifs

{% capture summary %} Code {% endcapture %} {% capture details %}  

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
  register(MulticoreParam(5))

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.Rds")

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
  saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

  #mouse
   mm10_atac <- SetAssayData(
    object = mm10_atac,
    assay = 'peaks',
    slot = 'motifs',
    new.data = motif.mm10)

  mm10_atac <- RegionStats(object = mm10_atac, genome = BSgenome.Mmusculus.UCSC.mm10)
  mm10_atac <- RunChromVAR( object = mm10_atac,genome = BSgenome.Mmusculus.UCSC.mm10)
  saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

```
{% endcapture %} {% include details.html %} 

### Differential Gene Activity through Subclusters

{% capture summary %} For hg38 {% endcapture %} {% capture details %}  

```R
  library(JASPAR2020)
  library(BSgenome.Hsapiens.UCSC.hg38)
  #library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
  #library(TFBSTools)

library(ggdendro)
library(dendextend)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(circlize)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

hg38_atac<-readRDS("hg38_SeuratObject.Rds")
hg38_atac$cluster_ID<-paste(hg38_atac$seurat_clusters,hg38_atac$seurat_subcluster,sep="_")
hg38_atac$celltype<-"NA"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==0,]$celltype<-"ExN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==1,]$celltype<-"Oligo"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==2,]$celltype<-"Astro"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==3,]$celltype<-"Micro.PVM"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==4,]$celltype<-"iN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==5,]$celltype<-"iN"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters==6,]$celltype<-"OPC"
hg38_atac$celltype_col<-"NA"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="iN",]$celltype_col<-"#cc3366"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="ExN",]$celltype_col<-"#3399cc"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Oligo",]$celltype_col<-"#66cc99"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Astro",]$celltype_col<-"#99cc99"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="Micro.PVM",]$celltype_col<-"#ff6633"
hg38_atac@meta.data[hg38_atac@meta.data$celltype=="OPC",]$celltype_col<-"#ffcc99"

saveRDS(hg38_atac,"hg38_SeuratObject.Rds")


#Clusters to test
cluster_to_test<-names(table(hg38_atac$cluster_ID)[table(hg38_atac$cluster_ID)>50 & !endsWith(names(table(hg38_atac$cluster_ID)),"NA")])
#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group,assay.="GeneActivity"){
    da_ga_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
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
    assay.="chromvar",
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

saveRDS(hg38_atac,"hg38_SeuratObject.Rds")

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
{% endcapture %} {% include details.html %} 

{% capture summary %} For mm10 {% endcapture %} {% capture details %}  
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

mm10_atac<-readRDS("mm10_SeuratObject.Rds")
mm10_atac@meta.data$celltype<-"NA"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="0",]$celltype<-"Gran"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="1",]$celltype<-"iN"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="2",]$celltype<-"ExN"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="3",]$celltype<-"Oligo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="4",]$celltype<-"iN"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="5",]$celltype<-"Astro"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="6",]$celltype<-"Endo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="7",]$celltype<-"Endo"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters=="8",]$celltype<-"Unknown"

mm10_atac$celltype_col<-"NA"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="iN",]$celltype_col<-"#cc3366"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="ExN",]$celltype_col<-"#3399cc"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Oligo",]$celltype_col<-"#66cc99"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Astro",]$celltype_col<-"#99cc99"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Unknown",]$celltype_col<-"#ff6633"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Endo",]$celltype_col<-"#A52A2A"
mm10_atac@meta.data[mm10_atac@meta.data$celltype=="Gran",]$celltype_col<-"#4B0082"

saveRDS(mm10_atac,"mm10_SeuratObject.Rds")

#Clusters to test
cluster_to_test<-names(table(mm10_atac$cluster_ID)[table(mm10_atac$cluster_ID)>50 & !endsWith(names(table(mm10_atac$cluster_ID)),"NA")])
#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group,assay.="GeneActivity"){
    da_ga_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T,
        assay=assay.
        )
    da_ga_tmp$da_region<-row.names(da_ga_tmp)
    da_ga_tmp$enriched_group<-c(i)
    da_ga_tmp$compared_group<-c("all_other_cells")
    return(da_ga_tmp)
  }

da_ga<-list() #set up an empty list for looping through

n.cores=5 #Perform parallel application of DA test
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
    assay.="chromvar",
    mc.cores=n.cores)

da_ga_df<-do.call("rbind",da_ga) #Merge the final data frame from the list for 1vrest DA
#To convert JASPAR ID TO TF NAME
da_ga_df$tf_name <- unlist(lapply(unlist(lapply(da_ga_df$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
write.table(da_ga_df,file="mm10.onevrest.da_chromvar.txt",sep="\t",col.names=T,row.names=T,quote=F)



#Plot out top ga for each cluster
da_ga<-read.csv(file="mm10.onevrest.da_ga.txt",head=T,sep="\t",row.names=NULL)
da_ga$gene_name<-da_ga$da_region
da_ga<-da_ga[complete.cases(da_ga),]
da_ga<-da_ga[!endsWith(da_ga$enriched_group,"NA"),]

da_ga$label<-""
for (x in unique(da_ga$enriched_group)){
selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:3))$da_region
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

###Subset by columns and markers
#sum_ga_sub<-sum_ga[match(markers$marker_2,row.names(sum_ga),nomatch=0),]
#sum_ga_sub<-sum_ga_sub[!duplicated(sum_ga_sub),]
#markers<-markers[match(row.names(sum_ga_sub),markers$marker_2,nomatch=0),]

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
top_ha<-columnAnnotation(df= data.frame(cluster=annot$seurat_clusters, celltype=annot$celltype), #subcluster=annot$cluster_ID,
                col=list(celltype=celltype,cluster=cluster,subcluster=subcluster),
                    show_legend = c(TRUE, TRUE,FALSE))  


colfun=colorRamp2(quantile(unlist(sum_ga_plot), probs=c(0.5,0.90,0.95)),magma(3))

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


pdf("mm10.tf.heatmap.pdf",height=20,width=20)
plt1
dev.off()
system("slack -F mm10.tf.heatmap.pdf ryan_todo")


#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)


motif_order<-names(hg38_atac@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(hg38_atac@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = hg38_atac,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())
ggsave(plt,file="mm10.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
system("slack -F mm10.tf.heatmap.motif.pdf ryan_todo")


feat<-c("Gad1","Lhx6","Adarb2","Slc17a7","Fezf2","Bcl11b","Cux2","Reln","Atoh1","Calb2","Ppp1r17","Gabra6","Prox1","Dsp")
plt<-FeaturePlot(mm10_atac,features=feat,order=T)
ggsave(plt,file="test.png",width=10*as.integer(length(feat)/3),height=10*as.integer(length(feat)/3),limitsize=F)
system("slack -F test.png ryan_todo")
```
{% endcapture %} {% include details.html %} 

### Marker Plotting for cell type refinement

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(Seurat)
library(Signac)
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
hg38_atac<-readRDS("hg38_SeuratObject.Rds")

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

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
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

mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
marker_list<-readRDS("./marker_sets/grosscelltype_markerlist.rds")
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
{% endcapture %} {% include details.html %} 

### Plotting of subcluster marker sets

{% capture summary %} Code {% endcapture %} {% capture details %}  

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

{% endcapture %} {% include details.html %} 


### Looks like there is a bias in the data going to integrate with harmony
{% capture summary %} For hg38 {% endcapture %} {% capture details %}  

```R
library(Seurat)
library(Signac)
library(ggplot2)

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")
hg38_atac<-readRDS("hg38_SeuratObject.Rds")
mm10_atac<-readRDS("mm10_SeuratObject.Rds")

plt<-DimPlot(hg38_atac,group.by="pcr_idx")
ggsave(plt,file="hg38.pcridx.png")

plt<-DimPlot(mm10_atac,group.by="pcr_idx")
ggsave(plt,file="mm10.pcridx.png")

hg38_1<-subset(hg38_atac,pcr_idx=="CAGAGAGG")
hg38_2<-subset(hg38_atac,pcr_idx=="CGTACTAG")
hg38_3<-subset(hg38_atac,pcr_idx=="CTCTCTAC")

anchors <- FindIntegrationAnchors(
  object.list = c(hg38_1,hg38_2,hg38_3),
  anchor.features = rownames(hg38_1),
  assay = c("peaks","peaks","peaks"),
  k.filter = NA
)


# integrate data and create a new merged object
integrated <- IntegrateData(
  anchorset = anchors,
  weight.reduction = hg38_1[['cistopic']],
  dims = 2:30,
  preserve.order = TRUE
)

###TO BE UPDATED

# we now have a "corrected" TF-IDF matrix, and can run LSI again on this corrected matrix
integrated <- RunSVD(
  object = integrated,
  n = 30,
  reduction.name = 'integratedLSI'
)

integrated <- RunUMAP(
  object = integrated,
  dims = 2:30,
  reduction = 'integratedLSI'
)

###

```

{% endcapture %} {% include details.html %} 


# Haven't run following code. Commented it out for now.

<!--
## Comparison of sciDROP adult mouse brain reads with available data sets
 TO BE DONE

```R
setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/201107_sciDROP_Barnyard")

dat_10x<-read.csv("/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/s3ATAC_AdultMusBrain_Comparison/10x_atac_v1_adult_brain_fresh_5k_singlecell.csv")
dat_10x<-dat_10x[c("barcode","passed_filters")]
barcused_10x<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/s3ATAC_AdultMusBrain_Comparison/10x_atac_v1_adult_brain_fresh_5k_barcodesused.txt")
dat_10x<-dat_10x[dat_10x$barcode %in% barcused_10x$V1,]
colnames(dat_10x)<-c("cellID","uniq_reads")
dat_10x$assay<-"10x_scATAC"

dat_biorad<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/s3ATAC_AdultMusBrain_Comparison/mousebrain-master_dataframe.rds")
dat_biorad_projected<-dat_biorad[c("DropBarcode","librarySize")]
colnames(dat_biorad_projected)<-c("cellID","uniq_reads")
dat_biorad_projected$assay<-"dscATAC_projected"

dat_biorad<-dat_biorad[c("DropBarcode","uniqueNuclearFrags")]
colnames(dat_biorad)<-c("cellID","uniq_reads")
dat_biorad$assay<-"dscATAC"

dat_ren<-read.csv("/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/s3ATAC_AdultMusBrain_Comparison/GSM2668124_p56.nchrM.merge.sel_cell.qc.csv",header=F)
colnames(dat_ren)<-c("cellID","uniq_reads","perc_uniq","frip")
dat_ren<-dat_ren[c("cellID","uniq_reads")]
dat_ren$assay<-"snATAC"

dat_sciatac<-read.table("/home/groups/oroaklab/adey_lab/projects/spatial/wholebrain/wholebrain.complexity.txt",header=F)
colnames(dat_sciatac)<-c("rowname_carryover","cellID","tot_reads","uniq_reads","perc_uniq")
dat_sciatac_condition<-read.table("/home/groups/oroaklab/adey_lab/projects/spatial/wholebrain/annots_wholebrain/wholebrain_FreshvFrozen.annot")
colnames(dat_sciatac_condition)<-c("cellID","condition")
dat_sciatac<-merge(dat_sciatac,dat_sciatac_condition,by="cellID")
dat_sciatac<-dat_sciatac[dat_sciatac$condition=="Frozen",] #filter by frozen condition
dat_sciatac_dims<-read.table("/home/groups/oroaklab/adey_lab/projects/spatial/wholebrain/wholebrain.bbrd.q10.filt.500.counts.cistopic.UMAP.dims") 
colnames(dat_sciatac_dims)<-c("cellID","x","y") #filter by cells used in analysis (dims file are those which underwent clustering)
dat_sciatac<-dat_sciatac[c("cellID","uniq_reads")]
dat_sciatac<-dat_sciatac[dat_sciatac$cellID %in% dat_sciatac_dims$cellID,]
dat_sciatac$assay<-"sciatac"


dat_s3<-read.table("mm10.complexity.txt",header=F)
colnames(dat_s3)<-c("rowname_carryover","cellID","tot_reads","uniq_reads","perc_uniq")
dat_s3_pcrplate<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data/barnyard_pcrplate.annot",header=F,col.names=c("cellID","pcr_plate"))
dat_s3<-merge(dat_s3,dat_s3_pcrplate,by="cellID")
#summarize percent unique per plate
library(dplyr)

dat_s3<-dat_s3[dat_s3$uniq_reads>=10000,]

data.frame(dat_s3 %>% group_by(pcr_plate) %>% summarise(mean_compl=mean(perc_uniq),
                                                        median_compl=median(perc_uniq),
                                                        sd=sd(perc_uniq),
                                                        mean_uniq=mean(uniq_reads),
                                                        median_uniq=median(uniq_reads),
                                                        sd_uniq=sd(uniq_reads),
                                                        cell_count=n()))

#  pcr_plate mean_compl median_compl        sd mean_uniq median_uniq   sd_uniq cell_count
#1         B   71.48840       71.820  3.668856  46299.46     33502.0  39659.36 381
#2         C   36.31242       27.895 16.913973 264133.51    181248.5 269722.56 298
#3         D   66.16967       67.690  7.614098  43358.67     26490.0  48451.74 273 
       
dat_s3_b<-dat_s3[dat_s3$pcr_plate=="B",c("cellID","uniq_reads")]
dat_s3_b$assay<-"s3ATAC_Plate_B"
dat_s3_c<-dat_s3[dat_s3$pcr_plate=="C",c("cellID","uniq_reads")]
dat_s3_c$assay<-"s3ATAC_Plate_C"

dat_10x<-dat_10x[dat_10x$uniq_reads>=1000,]

dat<-rbind(dat_10x,dat_ren,dat_biorad,dat_sciatac,dat_s3_c)
dat$assay  = factor(dat$assay, levels=c("snATAC", "10x_scATAC","dscATAC","sciatac","s3ATAC_Plate_C"))

library(ggplot2)
ggplot(dat,aes(x=as.factor(assay),y=log10(uniq_reads),color=as.factor(assay)))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust=1))
ggsave("adultmusbrain_atacprotocol_comparisons.svg")
ggsave("adultmusbrain_atacprotocol_comparisons.pdf")
ggsave("adultmusbrain_atacprotocol_comparisons.png")
i="s3ATAC_Plate_C"
for (j in unique(dat$assay)){
	ttemp<-t.test(dat[dat$assay==i,]$uniq_reads,dat[dat$assay==j,]$uniq_reads,alternative="greater")
	print(paste(i,j,ttemp$p.value,ttemp$estimate[1],ttemp$estimate[2],ttemp$alternative,ttemp$method,sep=","))
	}
#[1] "s3ATAC_Plate_C,10x_scATAC,3.13338819899197e-40,264133.506711409,22715.9052708283,greater,Welch Two Sample t-test"
#[1] "s3ATAC_Plate_C,snATAC,5.7106611661251e-42,264133.506711409,15459.4017798286,greater,Welch Two Sample t-test"
#[1] "s3ATAC_Plate_C,dscATAC,1.5654951763642e-37,264133.506711409,34045.7355797055,greater,Welch Two Sample t-test"
#[1] "s3ATAC_Plate_C,sciatac,9.87044162982429e-43,264133.506711409,12263.8299081767,greater,Welch Two Sample t-test"
#[1] "s3ATAC_Plate_C,s3ATAC_Plate_C,0.5,264133.506711409,264133.506711409,greater,Welch Two Sample t-test"
```


## Plotting and updating metadata



```python
#renaming annot for simplified annotation file making
#rename processing_ processing. *annot

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)

annot.files<-list.files(pattern="annot$")
annot.files<-annot.files[startsWith(annot.files,prefix="barnyard")]
annot_append<-read.table(annot.files[1],col.names=c("cellID",strsplit(annot.files[1],"[.]")[[1]][1]))
for (i in 2:length(annot.files)){
  tmp<-read.table(annot.files[i],col.names=c("cellID",strsplit(annot.files[i],"[.]")[[1]][1]))
  annot_append<-merge(annot_append,tmp,by="cellID")
}

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
hg38_atac@meta.data$cellID<-row.names(hg38_atac@meta.data)
hg38_atac@meta.data$i5_idx_seq<-substr(hg38_atac@meta.data$cellID,9,16)
hg38_atac@meta.data$i7_idx_seq<-substr(hg38_atac@meta.data$cellID,0,8)
annot<-as.data.frame(hg38_atac@meta.data)
annot<-merge(annot,annot_append,by="cellID",all.x=T)
compl_hg38<-read.table("hg38.complexity.txt")
colnames(compl_hg38)<-c("row_carryover","cellID","total_reads","uniq_reads","perc_uniq")
annot<-merge(annot,compl_hg38,by="cellID")
hg38_fracontarget<-read.table("hg38.bbrd.q10.500.fracOnTarget.values")
colnames(hg38_fracontarget)<-c("cellID","frac_reads_in_peaks")
annot<-merge(annot,hg38_fracontarget,by="cellID")
row.names(annot)<-annot$cellID
hg38_atac@meta.data<-annot


mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
mm10_atac@meta.data$cellID<-row.names(mm10_atac@meta.data)
mm10_atac@meta.data$i5_idx_seq<-substr(mm10_atac@meta.data$cellID,9,16)
mm10_atac@meta.data$i7_idx_seq<-substr(mm10_atac@meta.data$cellID,0,8)
annot<-as.data.frame(mm10_atac@meta.data)
annot<-merge(annot,annot_append,by="cellID",all.x=T)
compl_mm10<-read.table("mm10.complexity.txt")
colnames(compl_mm10)<-c("row_carryover","cellID","total_reads","uniq_reads","perc_uniq")
annot<-merge(annot,compl_mm10,by="cellID")
mm10_fracontarget<-read.table("mm10.bbrd.q10.500.fracOnTarget.values")
colnames(mm10_fracontarget)<-c("cellID","frac_reads_in_peaks")
annot<-merge(annot,mm10_fracontarget,by="cellID")
row.names(annot)<-annot$cellID
mm10_atac@meta.data<-annot


#Filtering cells to those with fraction of reads in peaks >0.2
mm10_atac <- subset(mm10_atac, subset = frac_reads_in_peaks >= 0.2)
hg38_atac <- subset(hg38_atac, subset = frac_reads_in_peaks >= 0.2)

saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")
saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

```

## Performing cisTopic and UMAP



```python


## Plotting and writing out cell table summaries


```python
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(parallel)

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

dat_hg38<-hg38_atac@meta.data
dat_mm10<-mm10_atac@meta.data

library(ggplot2)

plt<-ggplot(dat_hg38,aes(x=barnyard_pcrplate,y=log10(uniq_reads)))+
geom_jitter()+
geom_boxplot(outlier.shape=NA)+
theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="hg38_complexity_boxplot.pdf")
ggsave(plt,file="hg38_complexity_boxplot.png")

plt<-ggplot(dat_mm10,aes(x=barnyard_pcrplate,y=log10(uniq_reads)))+
geom_jitter()+
geom_boxplot(outlier.shape=NA)+
theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="mm10_complexity_boxplot.pdf")
ggsave(plt,file="mm10_complexity_boxplot.png")


plt<-DimPlot(hg38_atac,group.by=c('seurat_clusters',
                                  'barnyard_platevolume',
                                  'barnyard_polymerase',
                                  'barnyard_storagecondition',
                                  'barnyard_pcrplate'
                                 ))

pdf("hg38.umap.pdf",width=20)
plt
dev.off()
ggsave(plt,file="hg38.umap.png",width=20)




plt<-DimPlot(mm10_atac,group.by=c('seurat_clusters',
                                  'barnyard_platevolume',
                                  'barnyard_polymerase',
                                  'barnyard_storagecondition',
                                  'barnyard_pcrplate'))



pdf("mm10.umap.pdf",width=20)
plt
dev.off()
ggsave(plt,file="mm10.umap.png",width=20)

dat_hg38$UMAP_1<-hg38_atac@reductions$umap@cell.embeddings[,1]
dat_hg38$UMAP_2<-hg38_atac@reductions$umap@cell.embeddings[,2]

dat_mm10$UMAP_1<-mm10_atac@reductions$umap@cell.embeddings[,1]
dat_mm10$UMAP_2<-mm10_atac@reductions$umap@cell.embeddings[,2]

write.table(dat_hg38,file="hg38_cell_summary.txt",col.names=T,row.names=T,sep="\t",quote=F)
write.table(dat_mm10,file="mm10_cell_summary.txt",col.names=T,row.names=T,sep="\t",quote=F)
```

## Adding gene activity matrix to ATAC processing
This section is not fully developed.

### Cicero for co-accessible sites


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
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

#Processing hg38 cells
# convert to CellDataSet format and make the cicero object
hg38_atac.cds <- as.cell_data_set(x = hg38_atac,group_by="seurat_clusters")
hg38_atac.cicero <- make_cicero_cds( hg38_atac.cds, reduced_coordinates = reducedDims(hg38_atac.cds)$UMAP)
genome <- seqlengths(hg38_atac) # get the chromosome sizes from the Seurat object
genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
conns <- run_cicero(hg38_atac.cicero, genomic_coords = genome.df, sample_num = 10) # run cicero
saveRDS(conns,"hg38_cicero_conns.Rds")
ccans <- generate_ccans(conns) # generate ccans
saveRDS(ccans,"hg38_cicero_ccans.Rds")
links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
Links(hg38_atac) <- links
saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

#Processing mm10 cells
mm10_atac.cds <- as.cell_data_set(x = mm10_atac,group_by="seurat_clusters")
mm10_atac.cicero <- make_cicero_cds( mm10_atac.cds, reduced_coordinates = reducedDims(mm10_atac.cds)$UMAP)
genome <- seqlengths(mm10_atac) # get the chromosome sizes from the Seurat object
genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
conns <- run_cicero(mm10_atac.cicero, genomic_coords = genome.df, sample_num = 10) # run cicero
saveRDS(conns,"mm10_cicero_conns.Rds")
ccans <- generate_ccans(conns) # generate ccans
saveRDS(ccans,"mm10_cicero_ccans.Rds")
links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
Links(mm10_atac) <- links
saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

```

## Plotting Coverage Plots



```python
#Making a set of marker lists for human cortex and mouse whole brain

#Human cortex marker list is from https://cells.ucsc.edu/?ds=allen-celltypes%2Fhuman-cortex&meta=clusterlabel
#And stored in https://docs.google.com/spreadsheets/d/15YwK2lWh018IyxgQsRxnrINxBFlG3XtU5kE8gArYxJ4/edit?usp=sharing
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data/marker_sets")

marker_list<-NULL
marker_list[["VLMC"]]<-c('ITIH5', 'ABCA9', 'DCN', 'ATP1A2', 'SAT1', 'NEAT1', 'SLC7A11', 'ABCA8', 'APOD', 'NKD1', 'PTN', 'FBLN1', 'LAMC3', 'PDGFRB', 'CXCL12', 'AKR1C2', 'SVIL', 'ABCA6', 'COLEC12', 'CYP1B1', 'UACA', 'AHNAK', 'TPCN1', 'SLC1A3', 'ADAM33', 'LAMA2', 'COL18A1', 'CYTL1', 'SLC6A1', 'EDN3', 'ARHGAP29', 'COL4A5', 'DAB2', 'BGN', 'LEPR', 'PTGDS', 'OGN', 'WDR86', 'PRELP', 'ISLR', 'PDGFRA', 'MRC2', 'COL15A1', 'SPRY4', 'SELENBP1', 'RHOJ', 'KIAA1755', 'SASH1', 'LHFP', 'RPS12P5', 'KANK2')
marker_list[["Astrocyte"]]<-c('ADGRV1', 'SLC1A3', 'ATP1A2', 'C1orf61', 'FGFR3', 'SLC1A2', 'PTGDS', 'GLUL', 'NDRG2', 'PRODH', 'NTM', 'COL5A3', 'CST3', 'MACF1', 'PTPRZ1', 'ATP1B2', 'SLC25A18', 'ZBTB20', 'MT3', 'MT2A', 'ATP13A4', 'PON2', 'AQP4', 'DTNA', 'GRAMD3', 'NCAN', 'ACSS1', 'MSI2', 'GJA1', 'SLCO1C1', 'SLC7A11', 'PAPLN', 'F3', 'RYR3', 'PSD2', 'ETNPPL', 'MYO10', 'RAPGEF3', 'TTYH1', 'SOX2', 'EMX2', 'HIF3A', 'DST', 'APOE', 'PDGFRB', 'SFXN5', 'AASS', 'FNBP1', 'SLC4A4', 'GRIN2C', 'PLXNB1')
marker_list[["L5ET"]]<-c('COL24A1', 'ADRA1A', 'COL5A2', 'CRYM', 'GRIK2', 'DGKH', 'SULF2', 'ATP6V1C2', 'VAT1L', 'NEFM', 'SPHKAP', 'NRP1', 'SLC26A4', 'COL21A1', 'TOX', 'PDE1C', 'NTNG1', 'IQCA1', 'BCL11B', 'ANXA4', 'ARAP2', 'ADAMTSL3', 'FAM126A', 'LRP2', 'LOC101927745', 'SLC5A8', 'DOCK4', 'NEFH', 'PRUNE2', 'GRB14', 'GRM8', 'MYO16', 'FAM84B', 'ANKRD18A', 'TRDMT1', 'PCSK6', 'ANO4', 'SSTR2', 'ESRRG', 'MYLIP', 'VASH2', 'KLHL32', 'RYR3', 'ZNF189', 'SLC26A4-AS1', 'EPM2A', 'CNR1', 'ST3GAL1', 'LOC105374239', 'ANKRD34B', 'SEPT4')
marker_list[["L6CT"]]<-c('SEMA3E', 'TLE4', 'EGFEM1P', 'HS3ST4', 'SYT6', 'RYR3', 'ENPP7P4', 'DLC1', 'SEMA3A', 'MEIS2', 'COL24A1', 'DIP2A', 'PCDH17', 'RXFP1', 'LRP8', 'ANKRD18A', 'FOXP2', 'NRP1', 'SEMA5A', 'LRP1B', 'LOC105374199', 'ABCC9', 'SLC24A2', 'VWA2', 'SYNJ2', 'FAM95C', 'SLC8A1-AS1', 'ITGA11', 'SERPINE2', 'KLHL5', 'NFIA', 'CDH6', 'SPARCL1', 'CPE', 'CTGF', 'TENM1', 'ANKRD26P3', 'SEZ6L', 'GRIK3', 'MCC', 'SYNPO2', 'ZFHX3', 'SULF1', 'FILIP1', 'DGKG', 'LUZP2', 'ANKRD20A5P', 'ANKRD20A7P', 'DMD', 'PCLO', 'TMEFF2')
marker_list[["L4IT"]]<-c('NTNG1', 'HTR2A', 'VAV3', 'PCP4', 'TRPC3', 'CDR1', 'PCSK1', 'BHLHE22', 'PLEKHH2', 'ITM2A', 'BTBD3', 'SORL1', 'CUX2', 'VSTM2A', 'SLC17A6', 'POU6F2', 'MEF2A', 'SEMA6D', 'RORA', 'LOC105377864', 'LRRK2', 'NEFM', 'GPR26', 'KCNH8', 'FSTL1', 'MET', 'ACVR1C', 'LOC105375817', 'KCNS1', 'ND4', 'EGR1', 'RET', 'CD74', 'TIAM1', 'PTTG1', 'PHACTR2', 'TADA1', 'SLC38A11', 'SYT2', 'TMEM215', 'LOC101927653', 'RHBDL3', 'PRKCA', 'GRAMD3', 'ILDR2', 'BMP8A', 'ITPR2', 'ND2', 'TRAF5', 'RALB', 'AIFM3')
marker_list[["L56NP"]]<-c('HTR2C', 'CRYM', 'TSHZ2', 'PCP4', 'FGFR1', 'ROBO3', 'NPSR1-AS1', 'ENPP7P8', 'SEMA3E', 'TLE4', 'IFNG-AS1', 'ENPP7P7', 'RPS3AP34', 'COL11A1', 'ALCAM', 'PCDH17', 'ENPP7P4', 'NXPH2', 'CDH6', 'MYLIP', 'DIP2A', 'MEIS2', 'CALCRL', 'LRMP', 'KCNT2', 'FABP7', 'PHLDB2', 'SORCS2', 'ZNF385D', 'CD36', 'ITGA8', 'SEL1L3', 'DOCK4', 'ANKRD20A5P', 'SPOCK1', 'TLL1', 'SLC24A2', 'KLHL5', 'MYHAS', 'BCL11B', 'KCNIP1', 'NR4A3', 'LUZP2', 'COL9A1', 'CCDC80', 'TOX', 'TRHDE', 'VWC2L', 'FAM89A', 'ABCC9', 'LINC01279')
marker_list[["L6b"]]<-c('TLE4', 'CTGF', 'RXFP1', 'SEMA3E', 'ZFHX3', 'LOC100996635', 'ENPP7P7', 'EGFEM1P', 'NFIA', 'SEMA3A', 'DLC1', 'KLHL5', 'LOC401134', 'HS3ST4', 'DGKH', 'CDH9', 'NR4A2', 'PCSK5', 'ENPP7P8', 'RYR3', 'SEMA3D', 'ADD3', 'PCDH17', 'OLFML2B', 'HPCAL1', 'MDFIC', 'ITPR2', 'COL19A1', 'TMEFF2', 'SEMA3C', 'SEZ6L', 'PQLC3', 'ENPP7P4', 'CDH11', 'DOCK4', 'MGST1', 'TENM1', 'DPP4', 'TNIK', 'GLCE', 'SLC1A2', 'FILIP1', 'VWA5A', 'ARSJ', 'ROS1', 'LRP8', 'DGKG', 'CPLX3', 'SLC24A2', '43891', 'LRP1B')
marker_list[["PVALB"]]<-c('ERBB4', 'ZNF385D', 'TAC1', 'GAD1', 'SLIT2', 'KCNC2', 'GAD2', 'PVALB', 'KLHL5', 'KCNAB3', 'PTPRM', 'C8orf4', 'SAT1', 'SLC6A1', 'RGS5', 'LHX6', 'BTBD11', 'ANK1', 'RAB3IP', 'ZNF536', 'SOX6', 'SPARCL1', 'TENM1', 'LANCL1', 'LRP8', 'MYO16', 'KCNAB1', 'PAM', 'ZNF804A', 'NEAT1', 'NXPH1', 'OSBPL3', 'SLC9A9', 'SULF1', 'WIF1', 'CRHBP', 'RND3', 'LGI2', 'SEPT4', 'SERPINI1', 'ANO4', 'LOC105369345', 'KLF12', 'NPPC', 'RUNX2', 'DOCK11', 'SPOCK3', 'PLCE1', 'GRIA4', 'KCNS3', 'ARX')
marker_list[["VIP"]]<-c('CXCL14', 'SYNPR', 'RGS12', 'ERBB4', 'SCG2', 'VIP', 'GAD1', 'CALB2', 'ADARB2', 'THSD7A', 'CNR1', 'PROX1', 'ARL4C', 'DLX6-AS1', 'SLC6A1', 'DNER', 'NR2F2', 'FXYD6', 'ATP1B2', 'TIMP2', 'COL21A1', 'DOCK10', 'TAC3', 'CHRNA2', 'IGF1', 'PTHLH', 'SLC24A3', 'AP1S2', 'RGS16', 'DLX1', 'SHISA8', 'ZNF536', 'ADAM33', 'ADRA1A', 'ASIC4', 'SLC10A4', 'ROBO1', 'CRH', 'PCDH8', 'CNTNAP4', 'ENTPD3', 'HLA-A', 'SERINC1', 'GRIK2', 'ROBO2', 'PLD5', 'SEZ6', 'LOC105373643', 'HTRA1', 'KCNT2', 'PTPRE')
marker_list[["IT"]]<-c('ENC1', 'PTPRK', 'VSNL1', 'LINC00507', 'FAM13A', 'LOC101927745', 'CCK', 'TNNT2', 'ARAP2', 'LMF1', 'EPB41L2', 'ART3', 'HTR2A', 'ANKRD18B', 'LOC105374973', 'LOC101929974', 'LOC105373009', 'ARL4A', 'OTOGL', 'CRYM', 'RFX3', 'ADAMTS3', 'SLC38A11', 'NWD2', 'PTPN3', 'COL5A2', 'COL24A1', 'ID2', 'DMD', 'MYO5C', 'RFTN1', 'LOC102724736', 'LRRIQ1', 'ATP2B4', 'IQCA1', 'HSD3BP4''LOC105370610', 'RXFP1', 'LOC102724834', 'SYN3', 'ADCYAP1', 'DPYD', 'DGKA', 'LOC101929293', 'SMAD3', 'MEF2A', 'CNGB1', 'MYH9', 'TRABD2A', 'CA10', 'SYNJ2')
marker_list[["L56ITCAR3"]]<-c('RGS12', 'POSTN', 'NR4A2', 'SYNPR', 'SMYD1', 'ITGB8', 'FRAS1', 'ANXA1', 'OLFML2B', 'ITGA11', 'SPOCK1', 'ARAP2', 'GNB4', 'ABCC3', 'CCK', 'HSD3BP4', 'ENPP7P8', 'NWD2', 'THEMIS', 'MCTP2', 'FAP', 'SEMA6D', 'PDLIM5', 'EGFEM1P', 'GFRA1', 'MAMDC2', 'CD52', 'RXFP1', 'IL7R', 'GPR63', 'GNG7', 'CRABP1', 'ANXA5', 'TPMT', 'DGKH', 'B2M', 'SYNJ2', 'DCSTAMP', 'LOC105371828', 'LOC100421670', 'GALNT14', 'PRLR', 'HLA-L', 'MYOM2', 'SFMBT2', 'MREG', 'ATP10A', 'CD63', 'GSG1L', 'GAS2L3', 'SOCS2')
marker_list[["PAX6"]]<-c('CXCL14', 'RELN', 'CNR1', 'WIF1', 'AP1S2', 'GRIK2', 'CRH', 'CCK', 'ADARB2', 'DDR2', 'SP8', 'RGS12', 'HMBOX1', 'SLC6A1', 'ENTPD3', 'FXYD6', 'GAD1', 'SEZ6L', 'SPOCK3', 'ZNF385D', 'NR2F2', 'NXPH1', 'DOCK10', 'PNOC', 'TIMP2', 'SORCS3', 'MYO16', 'SCG2', 'ROBO2', 'PLS3', 'ZNF536', 'PCDH8', 'NECAB2', 'DNER', 'FSTL5', 'PAX6', 'CDH4', 'NPAS3', 'SERINC1', 'CPLX3', 'RBMS3', 'DLX6-AS1', 'SPOCK1', 'ARL4C', 'RAB3C', 'COL16A1', 'PTPRM', 'WWP1', 'LOC102724124', 'IGF1', 'SYNPR')
marker_list[["LAMP5"]]<-c('SLC6A1', 'GAD2', 'DNER', 'KIT', 'LAMP5', 'SV2C', 'GAD1', 'FSTL5', 'FXYD6', 'PTPRT', 'AP1S2', 'DOCK10', 'ADARB2', 'GRIK1', 'GRIA4', 'ZNF536', 'CNTNAP4', 'HAPLN1', 'NXPH2', 'ERBB4', 'CXCL14', 'RAB3C', 'GRIN3A', 'MYO16', 'ARL4C', 'PTCHD4', 'NXPH1', 'ATP1B2', 'FREM1', 'RAB3IP', 'SPOCK1', 'EYA4', 'KCNC2', 'NECAB2', 'CACNA2D1', 'DLX6-AS1', 'FAT1', 'PMEPA1', 'RELN', 'GRIK2', 'ADRA1A', 'PCP4L1', 'MTSS1', 'CA2', 'SPHKAP', 'BCL11B', 'TPD52L1', 'IGF1', 'SOX2-OT', 'PTPRM', 'ATP11C')
marker_list[["Oligodendrocyte"]]<-c('PLP1', 'MBP', 'PTGDS', 'CERCAM', 'TF', 'MOBP', 'MOG', 'CARNS1', 'LOC101929249', 'SCD', 'SLC44A1', 'CLDND1', 'CLDN11', 'PXK', 'RNASE1', 'CNTN2', 'PLEKHH1', 'ABCA2', 'PHLDB1', 'CRYAB', 'SPP1', 'ENPP2', 'QKI', 'TMEM144', 'APOD', 'MYRF', 'CNDP1', 'APLP1', 'TMEM63A', 'PPP1R14A', 'OPALIN', 'GPRC5B', 'MAG', 'UGT8', 'ST18', 'CAPN3', 'LAMP2', 'QDPR', 'EDIL3', 'HAPLN2', 'ERMN', 'TTYH2', 'NDRG1', 'PMP22', 'LINC00844', 'MAL', 'SOX2-OT', 'ABCA8', 'DBNDD2', 'ANLN', 'S100B')
marker_list[["SST"]]<-c('SST', 'GRIK1', 'GAD1', 'SYNPR', 'SPOCK3', 'LOC105372768', 'PAWR', 'LRP8', 'NMU', 'RAB3B', 'GRIN3A', 'LHX6', 'TRB', 'ROBO2', 'CDH13', 'KCNA3', 'KLHL5', 'NTM', 'GAD2', 'COL25A1', 'MAFB', 'CORT', 'STXBP6', 'SLC6A1', 'PNOC', 'ARX', 'TIMP2', 'GRIK2', 'FXYD6', 'GRIK1-AS2', 'SOX6', 'MAF', 'TRHDE', 'NXPH1', 'PLCH1', 'CDH9', 'GRIK3', 'PAM', 'ELFN1', 'FLT3', 'DLX1', 'PTPRM', 'NPAS3', 'CRHBP', 'RAB3IP', 'RNASET2', 'TMEM47', 'SPARCL1', 'ANK1', 'ADCYAP1R1', 'RBP4')
marker_list[["OPC"]]<-c('PTPRZ1', 'VCAN', 'OLIG1', 'PDGFRA', 'BCAN', 'NTM', 'APOD', 'COL9A1', 'C1orf61', 'PCDH15', 'EPN2', 'HIP1R', 'COL20A1', 'LOC105379054', 'ZBTB20', 'PLEKHH2', 'SCD5', 'CST3', 'OLIG2', 'PTN', 'SOX6', 'HTRA1', 'GPNMB', 'SNX22', 'SMOC1', 'SEMA5A', 'SULF2', 'LUZP2', 'NAV1', 'CCDC50', 'CSPG5', 'SPRY4', 'ADGRG1', 'LHFPL3', 'COL9A2', 'NKAIN4', 'BCAS1', 'KAT2B', 'SDC3', 'FERMT1', 'QKI', 'CSPG4', 'LPPR1', 'DSCAM', 'KANK1', 'COL11A1', 'FGFR1', 'MYT1', 'MEGF11', 'LRP4', 'PHLDA1')
marker_list[["Pericyte"]]<-c('PDGFRB', 'IGFBP7', 'ATP1A2', 'NOTCH3', 'ITIH5', 'SLC6A12', 'HIGD1B', 'RGS5', 'DCN', 'SLC6A1', 'PTN', 'IFITM3', 'CALD1', 'GNG11', 'BGN', 'FN1', 'IFITM1', 'SLC19A1', 'NDUFA4L2', 'PLXDC1', 'B2M', 'DLC1', 'TXNIP', 'COLEC12', 'PRELP', 'ADIRF', 'SPARC', 'SLC12A7', 'SLC6A13', 'CD9', 'EPS8', 'ISYNA1', 'LGALS1', 'IFITM2', 'EPAS1', 'MYO1B', 'STOM', 'HSPB1', 'EMP2', 'HES1', 'CYBA', 'PTGDS', 'MCAM', 'MYL12A', 'PDE7B', 'CD63', 'NR2F2', 'ARHGAP29', 'FRZB', 'TNS2', 'GJC1')
marker_list[["Endothelial"]]<-c('HLA-B', 'ABCG2', 'CLDN5', 'FLT1', 'HLA-E', 'B2M', 'A2M', 'ADGRF5', 'XAF1', 'HLA-C', 'ABCB1', 'SLC2A1', 'PODXL', 'EMP2', 'IFITM3', 'CLEC3B', 'EPAS1', 'IFI27', 'ID3', 'PTPRB', 'ENG', 'IFI44', 'MFSD2A', 'NEAT1', 'GPR85', 'ID1', 'VIM', 'SRGN', 'COBLL1', 'IFI44L', 'ITM2A', 'GPCPD1', 'DUSP1', 'HLA-H', 'LIMS2', 'MT2A', 'TM4SF1', 'GIMAP7', 'IGFBP7', 'VWF', 'SLC7A5', 'STOM', 'ITGA6', 'LMO2', 'PDGFB', 'SLC38A5', 'ESAM', 'KLF2', 'SLC7A1', 'HLA-A', 'HSPB1')
marker_list[["Microglia"]]<-c('LPAR6', 'CD74', 'ADAM28', 'P2RY12', 'CSF1R', 'CX3CR1', 'C3', 'RNASET2', 'SAT1', 'A2M', 'MAF', 'ITGAX', 'BHLHE41', 'BLNK', 'SRGAP2', 'HLA-DRB1', 'APBB1IP', 'ZFP36L2', 'LAPTM5', 'P2RY13', 'IFNGR1', 'PTPRC', 'CSF3R', 'ZFP36L1', 'CYBA', 'USP53', 'SFMBT2', 'C1QC', 'CSF2RA', 'LINC01268', 'SAMSN1', 'HLA-DRA', 'RCSD1', 'ADAP2', 'FYB', 'LTC4S', 'TYROBP', 'CD53', 'ARHGAP24', 'IFI16', 'C1QB', 'SLCO2B1', 'SLC1A3', 'MS4A7', 'GPR34', 'AIF1', 'TAL1', 'ITPR2', 'DOCK8', 'HLA-DPA1', 'MEF2A')
saveRDS(marker_list,"hg38_markerlist.rds")

#Human cortex marker list is from https://cells.ucsc.edu/?ds=aging-brain
#And stored in https://docs.google.com/spreadsheets/d/15YwK2lWh018IyxgQsRxnrINxBFlG3XtU5kE8gArYxJ4/edit?usp=sharing
marker_list<-NULL
marker_list[["SMC"]]<-c('Acta2', 'Tagln', 'Myl9', 'Tpm2', 'Mylk', 'Myh11', 'Cald1', 'Tpm1', 'Crip1', 'Mustn1', 'Flna', 'Pln', 'Sncg', 'Gm13889', 'Csrp1', 'Lmod1', 'Ppp1r14a', 'Des', 'Rgs4', 'Palld', 'Fxyd1', 'Nexn', 'Rasl11a', 'Aspn', 'Filip1l', 'Map1b', 'Pde3a', 'Igfbp7', 'Gm13861', 'Ccnd2', 'Snhg18', 'Lgals1', 'Gja4', 'Gucy1a1', 'Sparcl1', 'Errfi1', 'Gper1', 'Lbh', 'Hspb2', 'Crispld2', 'Mfge8', 'Cryab', 'Bgn', 'Cnn1', 'Higd1b', 'Rrad', 'Rasl12', 'Wtip', 'Notch3', 'Nrip2', 'Perp')
marker_list[["mNeur"]]<-c('Meg3', 'Snhg11', 'Rtn1', 'Bex2', 'Ly6h', 'Celf4', 'Syt1', 'Snap25', 'Pcsk1n', 'Map1b', 'Stmn1', 'Stmn3', 'Ahi1', '6330403K07Rik', 'Reln', 'Stmn2', 'Pcp4', 'Gng3', 'Uchl1', 'Zcchc18', 'Bcl11a', 'Peg3', 'Nap1l5', 'Pcsk2', 'Resp18', 'Atp1b1', 'Rufy3', 'Aplp1', 'Scg2', 'Gria2', 'Scg5', 'Nhlh2', 'Sgip1', 'Gm1673', 'Gpm6a', 'Nrxn3', 'Syn2', 'Lhx1os', 'Mapt', 'Cacna2d2', 'Lhx1', 'Tubb3', 'Elavl3', 'Cck', 'Fxyd6', 'Map7d2', 'Negr1', 'Cit', 'Thy1', 'Fxyd7', 'Cbarp')
marker_list[["Hb_EC"]]<-c('Hba-a2', 'Hba-a1', 'Hbb-bt', 'Hbb-bs', 'Alas2', 'Cldn5', 'Ly6c1', 'Ly6a', 'Itm2a', 'Flt1', 'Cxcl12', 'Spock2', 'Egfl7', 'Pltp', 'Slco1a4', 'Snca', 'Slc9a3r2', 'Bsg', 'Hspb1', 'Pglyrp1', 'Ctla2a', 'Apold1', 'Ifitm3', 'Bpgm', 'Id1', 'Car4', 'Klf2', 'Igfbp7', 'Crip1', 'Abhd2', '8430408G22Rik', 'Rgs5', 'Foxq1', 'Cyr61', 'Cd24a', 'Slc25a37', 'Tpm1', 'Rec114', 'Hspa1a', 'Ptn', 'Ftl1', 'Pam', 'Ube2c', 'Cald1', 'Vtn', 'Higd1b', 'Jund', 'Sparcl1', 'Tprgl', 'Jun', 'Nfkbia')
marker_list[["imNeur"]]<-c('Sox11', 'Stmn1', 'Sox4', 'Ccnd2', 'Tubb3', 'Meis2', 'Stmn2', 'Marcksl1', 'Cd24a', 'Tubb2b', 'Rtn1', 'Pfn2', 'Dcx', 'Dlx6os1', 'Map1b', 'Nrxn3', 'Hmgn2', 'Ncam1', 'Nnat', 'Pbx1', 'Celf4', 'Meg3', 'Elavl3', 'Dlx1', 'Igfbpl1', 'Stmn4', 'Gm1673', 'Bcl11a', 'Gm17750', 'Stmn3', 'Cdk5r1', 'H3f3b', '6330403K07Rik', 'Serf1', 'Basp1', 'Gad2', 'Hist3h2ba', 'Ptprs', 'Dlx2', 'Mpped2', 'Gng3', 'Islr2', 'Ccdc88a', 'Fxyd6', 'Auts2', 'Arx', 'Dlx5', 'Tiam2', 'Sp9', 'Foxg1', 'Etv1')
marker_list[["NRP"]]<-c('Sox11', 'Hmgn2', 'Stmn1', 'Ccnd2', 'Pclaf', 'Marcksl1', 'Sox4', 'Cd24a', 'Marcks', 'Mdk', 'Cenpf', 'Meis2', 'Pantr1', 'Gm1673', 'Pfn2', 'Ube2c', 'Dlx1', 'Rtn1', 'Hist1h2ap', 'Tubb3', 'Elavl3', 'Ncam1', 'Dlx2', 'Bex2', 'Dcx', 'Nnat', 'Bcl11a', 'Gm17750', 'Nrxn3', 'Foxg1', 'Arx', 'Ptprs', 'Pou3f2', 'Igfbpl1', 'Tubb2b', 'Fxyd6', 'Ppp1r14b', 'Serf1', 'Pax6', 'Pbx1', 'Elavl2', 'Map1b', 'Ccnb1', 'Mfap2', 'Dlx6os1', 'Trim2', 'Stmn3', 'Hes5', 'Sox2', 'Fabp5', 'Mpped2')
marker_list[["OLG"]]<-c('Mgp', 'Igfbp6', 'Cldn11', 'Nov', 'Igfbp5', 'Rspo3', 'Nbl1', 'Lmo4', 'Slc47a1', 'Rbp1', '1500015O10Rik', 'Efemp1', 'Ptgds', 'Emp3', 'Prg4', 'Atp1b1', 'Col1a2', 'Emb', 'Apod', 'Ogn', 'Cryab', 'Bgn', 'Pcolce', 'Perp', 'Foxp2', 'Cald1', 'Gpc3', 'S100b', 'Cd74', 'Prrx2', 'Alcam', 'Ank2', 'Aspa', 'Asgr1', 'Nupr1', 'Thbs1', 'Serping1', 'Scara3', 'Timp2', 'Gjb6', 'Gjb2', 'Foxc2', 'Spp1', 'Foxd1', 'Gstm5', 'Penk', 'Zic1', 'Cped1', 'Dapl1', 'Fstl1', 'Col1a1')
marker_list[["EC"]]<-c('Cldn5', 'Flt1', 'Ly6c1', 'Ly6a', 'Itm2a', 'Cxcl12', 'Egfl7', 'Spock2', 'Pltp', 'Slco1a4', 'Pglyrp1', 'Id1', 'Slc9a3r2', 'Ifitm3', 'Bsg', 'Apold1', 'Hspb1', 'Ctla2a', 'Klf2', 'Car4', 'Foxq1', 'Crip1', 'Igfbp7', '8430408G22Rik', 'Hspa1a', 'Edn1', 'Abhd2', 'Cyr61', 'Jun', 'Ptn', 'Nfkbia', 'Jund', 'Junb', 'Tpm1', 'Gkn3', 'Stmn2', 'Fos', 'Rgs5', 'Trim47', 'Rbpms', 'Tprgl', 'Cdkn1c', 'Pam', 'Sparcl1', 'Usp53', 'Slc38a3', 'Stmn1', 'Ctnnbip1', 'Gm26532', 'Sept4', 'B430010I23Rik')
marker_list[["AC"]]<-c('Aldoc', 'Slc1a3', 'Slc1a2', 'Plpp3', 'Ntm', 'Atp1a2', 'Mt2', 'Mt3', 'Gpr37l1', 'Ndrg2', 'Gm3764', 'Gja1', 'Clu', 'Atp1b2', 'Ntsr2', 'Mt1', 'Bcan', 'Cldn10', 'Prdx6', 'Htra1', 'Prnp', 'Dclk1', 'Gstm5', 'Cxcl14', 'Ttyh1', 'Apoe', 'Pla2g7', 'Acsl3', 'Ntrk2', 'Dbi', 'Slc6a11', 'Gpm6b', 'Mmd2', 'Ptprz1', 'Gria2', 'Mfge8', 'Sparcl1', 'Slc4a4', 'F3', 'Cpe', 'Tspan7', 'Gstm1', 'Gpm6a', 'Rorb', 'Ckb', 'Dtna', 'Slc6a1', 'Cspg5', 'Scd2', 'Slc7a10', 'Ptn')
marker_list[["TNC"]]<-c('6330403K07Rik', 'S100a6', 'Prdx6', 'Ntrk2', 'Fbxo2', 'Dbi', 'Tmem47', 'Nnat', 'Mt3', 'Cpe', 'Gpm6b', 'Mt2', 'Mlc1', 'Igfbp5', 'Fxyd1', 'Mt1', 'Scd2', 'Id4', 'Sox9', 'Aldoc', 'Slc1a3', 'Gfap', 'Apoe', 'Thrsp', 'Gstm1', 'Mdk', 'Pcsk1n', 'Meg3', 'Zcchc18', 'Myoc', 'Rax', 'S100a1', 'Ttyh1', 'Mgst1', 'Six3', 'Dtna', 'Riiad1', 'Ndrg2', 'Chchd10', 'Pbx1', 'Map1b', 'Mia', 'Pygb', 'Kctd14', 'Gprc5b', 'Slc1a2', 'Scg5', 'Ank2', 'Atp1b2', 'Apoc1', 'Nrxn2')
marker_list[["EPC"]]<-c('Ccdc153', 'Tmem212', 'Rarres2', 'Tppp3', 'Riiad1', 'Gm19935', 'Nnat', 'Chchd10', 'Mia', 'Cd24a', 'Mt3', 'Dbi', 'Clu', 'S100b', 'Calml4', 'Mt2', 'Hspa2', '1500015O10Rik', 'Gm14964', 'Sox9', 'Tmem47', 'Mlc1', 'Prdx6', 'Ntrk2', 'Fam213a', 'Ramp1', 'Mt1', 'Cpe', 'Scd2', 'Mgst1', 'Gm1673', 'Lmo4', 'Prelp', 'Ppil6', 'Gstm1', 'Id4', 'Gm973', 'Ttll3', 'Plekhb1', 'Dalrd3', 'Cfap100', 'Zcchc18', 'S100a6', 'Trp53bp2', 'Aqp4', 'Sox2', 'Mro', 'Ckb', 'Wdr60', 'Tspan15', 'Spint2')
marker_list[["MG"]]<-c('Ctss', 'C1qc', 'Hexb', 'C1qa', 'C1qb', 'Tyrobp', 'Selplg', 'Trem2', 'Tmem119', 'Cx3cr1', 'Fcer1g', 'Lgmn', 'Csf1r', 'Laptm5', 'Rgs10', 'Ctsd', 'Fcrls', 'P2ry12', 'Fcgr3', 'Olfml3', 'Ctsz', 'Gpr34', 'Lpcat2', 'Pld4', 'Grn', 'Unc93b1', 'Aif1', 'Ly86', 'Vsir', 'Siglech', 'Ccl4', 'Ccl3', 'Trf', 'Irf8', 'Cd53', 'Bin1', 'Hexa', 'Atf3', 'Junb', 'Spi1', 'Marcks', 'Cd83', 'Hk2', 'Ltc4s', 'Cd37', 'Coro1a', 'Ptgs1', 'Tgfbr1', 'Slc15a3', 'Lair1', 'Hpgds')
marker_list[["MNC"]]<-c('Cd52', 'Lyz2', 'Plac8', 'Fcer1g', 'Msrb1', 'Tyrobp', 'Coro1a', 'Alox5ap', 'Cybb', 'Lgals3', 'Lsp1', 'Lst1', 'S100a6', 'Ifitm6', 'Ccl6', 'Ifitm3', 'Rac2', 'Cytip', 'Pim1', 'Lcp1', 'Hp', 'Stk17b', 'Emp3', 'Ifi27l2a', 'Spi1', 'Itgal', 'Ms4a6c', 'S100a4', 'Wfdc17', 'Il1b', 'S100a8', 'S100a9', 'Samsn1', 'Gmfg', 'Slpi', 'Ftl1', 'Cebpb', 'Napsa', 'Ly6c2', 'Ear2', 'Retnlg', 'Wfdc21', 'Cd53', 'Mcemp1', 'G0s2', 'Plbd1', 'Cd74', 'AB124611', 'Rps11', 'Sirpb1c', 'Gpr141')
marker_list[["PC"]]<-c('Vtn', 'Higd1b', 'Ndufa4l2', 'Rgs5', 'Cald1', 'Myl9', 'Cox4i2', 'Sept4', 'Pdgfrb', 'Kcnj8', 'Ifitm1', 'Rgs4', 'Ptn', 'Atp1a2', 'P2ry14', 'Rarres2', 'Atp13a5', 'Abcc9', 'Mgp', 'Sod3', 'Slc19a1', 'Slc6a20a', 'Nbl1', 'Art3', 'Zic1', 'Ppp1r14a', 'Lhfp', 'Gper1', 'Gucy1b1', 'Notch3', 'Bgn', 'Pitpnc1', 'Car4', 'Tbx3os1', 'Gja4', 'Gm13861', 'Perp', 'Phlda1', 'Il34', 'Uchl1', 'Slc38a11', 'Carmn', 'Gpx8', 'Mfge8', 'Pth1r', 'Ecm2', 'Gucy1a1', 'Igfbp7', 'Pde4d', 'Ddit4l', 'G0s2')
marker_list[["MAC"]]<-c('Pf4', 'Lyz2', 'Ms4a7', 'Dab2', 'Mrc1', 'Fcer1g', 'C1qc', 'Ctsc', 'Wfdc17', 'C1qa', 'C1qb', 'Tyrobp', 'F13a1', 'Maf', 'Stab1', 'Trf', 'Apoe', 'Ccl7', 'Ms4a6c', 'Ccl24', 'Ccl2', 'Lst1', 'Hexa', 'Clec4n', 'Cd14', 'Cbr2', 'Ccl12', 'Cd209f', 'Cybb', 'Fcrls', 'Grn', 'Ptpn18', 'Ifi207', 'Cxcl2', 'Cd68', 'Ms4a6b', 'Fcgr2b', 'Marcksl1', 'Lgmn', 'Ftl1', 'Csf1r', 'Ctss', 'Lgals1', 'Folr2', 'Igfbp4', 'Ltc4s', 'C5ar1', 'Fcgr3', 'Tslp', 'Tgfbi', 'Cd163')
saveRDS(marker_list,"mm10_markerlist.rds")

#Gross cell type list (for both organisms)
marker_list<-NULL
marker_list[["Oligodendrocyte"]]<-c('Mobp','Cldn1','Prox1','Olig1')
marker_list[["Polydendrocyte"]]<-c('Cspg4','Pdgfra')
marker_list[["Astrocyte"]]<-c('Gfap','GluI','Slc1a2','Agt')
marker_list[["Ex_Neuron"]]<-c('Satb2','Cux2','Rorb','Pcp4','Foxp2','Col19a1','Slc17a7')
marker_list[["Microglia"]]<-c('C1qa','C1qc','Cxcr1')
marker_list[["Inhib_Neuron"]]<-c('Cnr1','Bcl11b','Gad1','Dlx1','Dlx2')
saveRDS(marker_list,"grosscelltype_markerlist.rds")


```


```python
#Function to plot cell type markers
#For now just using the gross celltype marker list
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

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
        ncol = 1,
        tile=T,
        tile.size=500,
        tile.cells=downsamp
    )
    pdf(paste0("./marker_sets/",m,k,"_",j,"_genebody_accessibility.pdf"))
    print(plt)
    dev.off()
}


hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
marker_list<-readRDS("./marker_sets/grosscelltype_markerlist.rds")
summary(hg38_atac@meta.data$seurat_clusters) #setting tile cells to 90 to downsample
downsamp=90

for (x in 1:length(marker_list)){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    gene_list<-toupper(gene_list)
    celltype_name<-names(marker_list)[x]
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,mc.cores=20)
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

mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")
marker_list<-readRDS("./marker_sets/grosscelltype_markerlist.rds")
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

### Differential Accessibillity on Clusters



```python

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(parallel)
library(ggplot2)
library(ggrepel)
library(dplyr)

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")



#Perform One vs. rest DA enrichment

write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

#set up an empty list for looping through
hg38_da_peaks<-list()
mm10_da_peaks<-list()

#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group){
    da_peaks_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T
        )
    da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
    closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
    da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
    da_peaks_tmp$enriched_group<-c(i)
    da_peaks_tmp$compared_group<-c("all_other_cells")
    return(da_peaks_tmp)
  }


da_one_v_one<-function(i,obj,group,j_list){
    i<-as.character(i)
    da_tmp_2<-list()
    for (j in j_list){
        if ( i != j){
        da_peaks_tmp <- FindMarkers(
            object = obj,
            ident.1 = i,
            ident.2 = j,
            group.by = group,
            test.use = 'LR',
            latent.vars = 'nCount_peaks',
            only.pos=T
            )
        da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
        closest_genes <- ClosestFeature(hg38_atac,da_peaks_tmp$da_region)
        da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
        da_peaks_tmp$enriched_group<-c(i)
        da_peaks_tmp$compared_group<-c(j)
        da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
        }
    }
    return(da_tmp_2)
  }

#Perform parallel application of DA test
library(parallel)
n.cores=length(unique(hg38_atac@meta.data$seurat_clusters))
hg38_da_peaks<-mclapply(
    unique(hg38_atac@meta.data$seurat_clusters),
    FUN=da_one_v_rest,
    obj=hg38_atac,
    group="seurat_clusters",
    mc.cores=n.cores)

n.cores=length(unique(mm10_atac@meta.data$seurat_clusters))
mm10_da_peaks<-mclapply(
    unique(mm10_atac@meta.data$seurat_clusters),
    FUN=da_one_v_rest,
    obj=mm10_atac,
    group="seurat_clusters",
    mc.cores=n.cores)

#Merge the final data frame from the list for 1vrest DA
hg38_da_peaks<-do.call("rbind",hg38_da_peaks)
mm10_da_peaks<-do.call("rbind",mm10_da_peaks)

write("Outputting One v Rest DA Table.", stderr())
write.table(hg38_da_peaks,file="hg38.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)
write.table(mm10_da_peaks,file="mm10.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)


dat<-read.table("hg38.onevrest.da_peaks.txt",header=T,sep="\t")
dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% slice(1:2) #grabbing top 2 most significant peaks to label
plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(dat=dat_select,aes(label=gene_name,size=-distance),force=3)+theme_bw()
pdf("hg38_da_peaks.pdf")
plt
dev.off()

dat<-read.table("mm10.onevrest.da_peaks.txt",header=T,sep="\t")
dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% slice(1:2) #grabbing top 2 most significant peaks to label
plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(dat=dat_select,aes(label=gene_name,size=-distance),force=3)+theme_bw()
pdf("mm10_da_peaks.pdf")
plt
dev.off()

#Empty list to rerun for 1v1 comparisons
hg38_da_peaks<-list()
mm10_da_peaks<-list()
    
n.cores=length(unique(hg38_atac@meta.data$seurat_clusters))
hg38_da_peaks<-mclapply(
    unique(hg38_atac@meta.data$seurat_clusters),
    FUN=da_one_v_one,
    obj=hg38_atac,
    group="seurat_clusters",
    j_list=do.call("as.character",list(unique(hg38_atac@meta.data$seurat_clusters))),
    mc.cores=n.cores)

n.cores=length(unique(mm10_atac@meta.data$seurat_clusters))
mm10_da_peaks<-mclapply(
    unique(mm10_atac@meta.data$seurat_clusters),
    FUN=da_one_v_one,
    obj=mm10_atac,
    group="seurat_clusters",
    j_list=do.call("as.character",list(unique(mm10_atac@meta.data$seurat_clusters))),
    mc.cores=n.cores)

#Merge the final data frame from the list for 1v1 DA
hg38_da_peaks<-do.call("rbind",do.call("rbind",hg38_da_peaks))
mm10_da_peaks<-do.call("rbind",do.call("rbind",mm10_da_peaks))

write("Outputting One v One DA Table.", stderr())
write.table(hg38_da_peaks,file="hg38.onevone.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)
write.table(mm10_da_peaks,file="mm10.onevone.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)


#write("Generating a sanity check PDF for enrichment. Using Top 5 sites per annotation.", stderr())
#grab top N peaks per cell grouping
#N<-as.numeric(as.character(5))
#hg38_da_peaks_topN<-do.call("rbind",lapply(unique(hg38_da_peaks$enriched_group),function(x) head(hg38_da_peaks[hg38_da_peaks$enriched_group==x,][order(hg38_da_peaks[hg38_da_peaks$enriched_group==x,]$p_val,na.last=T),],N)))
#mm10_da_peaks_topN<-do.call("rbind",lapply(unique(mm10_da_peaks$enriched_group),function(x) head(mm10_da_peaks[mm10_da_peaks$enriched_group==x,][order(mm10_da_peaks[mm10_da_peaks$enriched_group==x,]$p_val,na.last=T),],N)))

#hg38_da_peaks_topN<-row.names(hg38_da_peaks_topN)
#mm10_da_peaks_topN<-row.names(mm10_da_peaks_topN)

#pdf("hg38.DApeaks.violinplot.pdf",height=100,width=100)
#VlnPlot(object = hg38_atac,features = hg38_da_peaks_topN,ncol = N,group.by="peaks_snn_res.0.2")
#dev.off()

#pdf("mm10.DApeaks.violinplot.pdf",height=100,width=100)
#VlnPlot(object = mm10_atac,features = mm10_da_peaks_topN,ncol = N,group.by="peaks_snn_res.0.2")
#dev.off()

#pdf("hg38.DApeaks.featureplot.pdf",height=100,width=100)
#FeaturePlot(object = hg38_atac,features = hg38_da_peaks_topN,ncol = N,pt.size = 1,max.cutoff=1)
#dev.off()

#pdf("mm10.DApeaks.featureplot.pdf",height=100,width=100)
#FeaturePlot(object = mm10_atac,features = mm10_da_peaks_topN,ncol = N,pt.size = 1,max.cutoff=1)
#dev.off()

#write("Performing GREAT on all enriched sites per annotation group", stderr())

#skipped this because rGREAT couldn't be installed in 3.6.1 right now
#library(rGREAT)

#format data as bed file
#write("Preparing Background Set as all called peaks.", stderr())
#hg38_bg_bed<-do.call("rbind",strsplit(unlist(hg38_atac@assays$peaks@counts@Dimnames[1]),"[-]"))
#hg38_bg_bed<-as.data.frame(hg38_bg_bed)
#colnames(hg38_bg_bed)<-c("chr","start","end")
#hg38_bg_bed$start<-as.numeric(as.character(hg38_bg_bed$start))
#hg38_bg_bed$end<-as.numeric(as.character(hg38_bg_bed$end))

#mm10_bg_bed<-do.call("rbind",strsplit(unlist(mm10_atac@assays$peaks@counts@Dimnames[1]),"[-]"))
#mm10_bg_bed<-as.data.frame(mm10_bg_bed)
#colnames(mm10_bg_bed)<-c("chr","start","end")
#mm10_bg_bed$start<-as.numeric(as.character(mm10_bg_bed$start))
#mm10_bg_bed$end<-as.numeric(as.character(mm10_bg_bed$end))

#write("Beginning loop through all annotation groups.", stderr())
#for (i in unique(hg38_da_peaks$enriched_group)){
#hg38_bed<-do.call("rbind",strsplit(hg38_da_peaks[hg38_da_peaks$enriched_group==i,]$da_region,"-"))
#hg38_bed<-as.data.frame(hg38_bed)
#colnames(hg38_bed)<-c("chr","start","end")
#hg38_bed$start<-as.numeric(as.character(hg38_bed$start))
#hg38_bed$end<-as.numeric(as.character(hg38_bed$end))
#write(paste("Using",nrow(hg38_bed), "DA peaks from",i), stderr())
#job = submitGreatJob(hg38_bed,hg38_bg_bed,species="hg38",request_interval=30)
#tb = getEnrichmentTables(job, ontology = c("GO Molecular Function", "GO Biological Process","GO Cellular Component"))
#tb = getEnrichmentTables(job, category = c("GO","Phenotype","Genes"))
#pdf(paste0("hg38_DApeaks_",i,".GeneAssociation.pdf"))
#plotRegionGeneAssociationGraphs(job)
#dev.off()
#for (j in 1:length(names(tb))){
#  write(paste("Outputting DA GREAT Analysis for", i, as.character(names(tb))[j]), stderr())
#  tabl_name<-gsub(" ","",as.character(names(tb))[j])
#  write.table(as.data.frame(tb[[j]]),file=paste("hg38_DApeaks_",i,tabl_name,"txt",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
#  }
#}

#tb<-NULL
#for (i in unique(mm10_da_peaks$enriched_group)){
#mm10_bed<-do.call("rbind",strsplit(mm10_da_peaks[mm10_da_peaks$enriched_group==i,]$da_region,"-"))
#mm10_bed<-as.data.frame(mm10_bed)
#colnames(mm10_bed)<-c("chr","start","end")
#mm10_bed$start<-as.numeric(as.character(mm10_bed$start))
#mm10_bed$end<-as.numeric(as.character(mm10_bed$end))
#write(paste("Using",nrow(mm10_bed), "DA peaks from",i), stderr())
#job = submitGreatJob(mm10_bed,mm10_bg_bed,species="mm10",request_interval=30)
#tb = getEnrichmentTables(job, ontology = c("GO Molecular Function", "GO Biological Process","GO Cellular Component"))
#tb = getEnrichmentTables(job, category = c("GO","Phenotype","Genes"))
#pdf(paste0("mm10_DApeaks_",i,".GeneAssociation.pdf"))
#plotRegionGeneAssociationGraphs(job)
#dev.off()
#for (j in 1:length(names(tb))){
#  write(paste("Outputting DA GREAT Analysis for", i, as.character(names(tb))[j]), stderr())
#  tabl_name<-gsub(" ","",as.character(names(tb))[j])
#  write.table(as.data.frame(tb[[j]]),file=paste("mm10_DApeaks_",i,tabl_name,"txt",sep="."),sep="\t",col.names=T,row.names=T,quote=F)
#  }
#}
```

## Labeling Cell Types


```python

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(parallel)
library(ggplot2)
library(ggrepel)
library(dplyr)

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

#Cell types deteremined by TF motif and marker genebody accessibility
hg38_atac@meta.data$celltype<-"NA"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("2"),]$celltype<-"excitatory_neuron"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("3","4"),]$celltype<-"inhibitory_neuron"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("0","1"),]$celltype<-"oligodendrocytes"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("7"),]$celltype<-"polydendrocytes"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("5"),]$celltype<-"microglia"
hg38_atac@meta.data[hg38_atac@meta.data$seurat_clusters %in% c("6"),]$celltype<-"astrocyte"
saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

mm10_atac@meta.data$celltype<-"NA"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters %in% c("0","2"),]$celltype<-"excitatory_neuron"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters %in% c("1"),]$celltype<-"inhibitory_neuron"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters %in% c("4"),]$celltype<-"oligodendrocytes_polydendrocytes"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters %in% c("5"),]$celltype<-"microglia"
mm10_atac@meta.data[mm10_atac@meta.data$seurat_clusters %in% c("3"),]$celltype<-"astrocyte"
saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

#Also adding TSS enrichment values to meta data
hg38_tss<-read.table("hg38.bbrd.q10.tss_reads.value",header=F)
colnames(hg38_tss)<-c("cellID","TSS_enrichment")
mm10_tss<-read.table("mm10.bbrd.q10.tss_reads.value",header=F)
colnames(mm10_tss)<-c("cellID","TSS_enrichment")

hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

hg38_annot<-hg38_atac@meta.data
mm10_annot<-mm10_atac@meta.data

hg38_annot<-merge(hg38_annot,hg38_tss,by="cellID")
mm10_annot<-merge(mm10_annot,mm10_tss,by="cellID")

row.names(hg38_annot)<-hg38_annot$cellID
row.names(mm10_annot)<-mm10_annot$cellID

hg38_atac@meta.data<-hg38_annot
mm10_atac@meta.data<-mm10_annot

saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")
saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")

```

# Subclustering and Trajectory Analysis

For Human cells, going to generate a trajectory through excitatory neurons to link to cortical layer markers

For mouse cells, going to subcluster inhibitory neurons and classify

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

### Running Cicero and gene activity on Exc. Neurons 


```python

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(EnsDb.Hsapiens.v86)


setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#Read in data and modify to monocle CDS file
#read in RDS file.
hg38_ex_neurons.cds<-readRDS("hg38_excNeurons.monocle.cds.Rds")
excneuron<-readRDS("hg38_excNeurons.SeuratObject.Rds")

#Processing hg38 cells
# convert to CellDataSet format and make the cicero object
hg38_ex_neurons.cicero <- make_cicero_cds( hg38_ex_neurons.cds, reduced_coordinates = reducedDims(hg38_ex_neurons.cds)$UMAP)
genome <- seqlengths(excneuron) # get the chromosome sizes from the Seurat object
genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
conns <- run_cicero(hg38_ex_neurons.cicero, genomic_coords = genome.df, sample_num = 10) # run cicero
saveRDS(conns,"hg38_excneuron_cicero_conns.Rds")
ccans <- generate_ccans(conns) # generate ccans
saveRDS(ccans,"hg38_excneuron_cicero_ccans.Rds")
links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
Links(excneuron) <- links
saveRDS(excneuron,file="hg38_exNeuron.SeuratObject.RDS")

# generate unnormalized gene activity matrix
#gene annotation sample

hg38_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
pos <-as.data.frame(hg38_annotations,row.names=NULL)
pos$chromosome<-paste0("chr",pos$seqnames)
pos$gene<-pos$gene_id
pos <- subset(pos, strand == "+")
pos <- pos[order(pos$start),] 
pos <- pos[!duplicated(pos$tx_id),] # remove all but the first exons per transcript
pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS

neg <-as.data.frame(hg38_annotations,row.names=NULL)
neg$chromosome<-paste0("chr",neg$seqnames)
neg$gene<-neg$gene_id
neg <- subset(neg, strand == "-")
neg <- neg[order(neg$start,decreasing=TRUE),] 
neg <- neg[!duplicated(neg$tx_id),] # remove all but the first exons per transcript
neg$end <- neg$end + 1 # make a 1 base pair marker of the TSS
gene_annotation<- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation <- gene_annotation[,c("chromosome","start","end","gene_name")]


# Rename the gene symbol column to "gene"
names(gene_annotation)[4] <- "gene"

hg38_ex_neurons.cds<- annotate_cds_by_site(hg38_ex_neurons.cds, gene_annotation)
unnorm_exNeuron_ga <- build_gene_activity_matrix(hg38_ex_neurons.cds, conns)
saveRDS(unnorm_exNeuron_ga,"hg38_exNeuron.unnorm_GA.Rds")

unnorm_exNeuron_ga<-readRDS("hg38_exNeuron.unnorm_GA.Rds")
excneuron<-readRDS("hg38_excNeurons.SeuratObject.Rds")

# add the gene activity matrix to the Seurat object as a new assay and normalize it
excneuron[['GA']]<-CreateAssayObject(counts = unnorm_exNeuron_ga)
excneuron <- NormalizeData(
  object = excneuron,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(excneuron$nFeature_peaks)
)

saveRDS(excneuron,"hg38_excNeurons.SeuratObject.Rds")
```

### Add ChromVar Motif Analysis To Exc. Neurons


```python
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#read in RDS file.
excneuron<-readRDS("hg38_exNeuron.SeuratObject.RDS")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(excneuron),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)


# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
excneuron <- SetAssayData(
  object = excneuron,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)

excneuron <- RegionStats(object = excneuron, 
    genome = BSgenome.Hsapiens.UCSC.hg38)

excneuron <- RunChromVAR(
  object = excneuron,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(excneuron,"hg38_exNeuron.SeuratObject.RDS")
excneuron<-readRDS("hg38_exNeuron.SeuratObject.RDS")
tfList <- getMatrixByID(JASPAR2020, ID=row.names(excneuron@assays$chromvar@data)) 
tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
dat<-excneuron@assays$chromvar@data
row.names(dat)<-tfList

dat<-data.frame(t(dat))
dat$cellID<-row.names(dat)
pseudospace<-excneuron@meta.data[c("cellID","pseudotime")]
dat<-merge(dat,pseudospace,by="cellID")

library(ComplexHeatmap)

cell_row_order=order(dat$pseudotime)
dat<-dat[!(colnames(dat) %in% c("cellID","pseudotime"))]
dat<-dat[head(order(apply(dat,2,var),decreasing=T),50)]

pdf("excneuron_tf_layering.pdf")
Heatmap(as.matrix(dat),
        row_order=cell_row_order,
        show_row_names=F,
        column_km=4,
        column_names_gp = gpar(fontsize = 2)
)
dev.off()
system("slack -F excneuron_tf_layering.pdf ryan_todo")
```

### Plot Gene Activity and TF Motifs through Pseudo"space" 


```python

##############New Session #####################
#Plotting layer markers
#Sourced from sciMAP ATAC paper and Tasic et al.
marker_list<-NULL
#marker_list[["L2_3"]]<-c("CALB1","RASGRF2","CUX2","PTGS2")
#marker_list[["L2_3_4"]]<-c("CUX1","POU3F2","RORB")
#marker_list[["L5"]]<-c("FEZF2","PCP4","RBP4","UCMA")
#marker_list[["L5_6"]]<-c("BCL11B","GRIK3")
#marker_list[["L6"]]<-c("HTR1F","NFIA","SLA","SYT6","TLE4","NTSR1")
#marker_list[["L4"]]<-c("SCNN1A-TG3","CTXN3","ARF5")
#marker_list[["L2"]]<-c("NGB")
#marker_list[["L5a"]]<-c("HSD11B1","FAM5C","SYT17","BATF3")
#marker_list[["L5b"]]<-c("SAMD3","CDH1","CHRNA6")
#marker_list[["L6a"]]<-c("MGP","SLA","CAR12","SYT17")
#Just the CAT curated list
marker_list<-c("CALB1","RASGRF2","CUX1","POU3F2","RORB","FEZF2","PCP4","BCL11B","GRIK3","HTR1F","NFIA","SLA","SYT6","TLE4")

library(ComplexHeatmap)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#Read in data and modify to monocle CDS file
#read in RDS file.

excneuron<-readRDS("hg38_exNeuron.SeuratObject.RDS") #using this for pseudospace meta data

marker_list<-unlist(marker_list)
cicero_gene_activities<-as.data.frame(cicero_gene_activities[row.names(cicero_gene_activities) %in% marker_list,])
dat<-data.frame(t(cicero_gene_activities))
dat$cellID<-row.names(dat)
dat<-merge(dat,excneuron@meta.data[c("cellID","pseudospace")],by="cellID")

library(reshape2)
dat.m<-melt(dat,id.vars=c("cellID","pseudospace"))
ggplot(dat.m,aes(x=pseudospace,y=value,color=as.factor(variable)))+geom_point()+geom_smooth()+theme_bw()
ggsave("test_excneuron_ALL.pdf")

cell_row_order=order(dat$pseudospace)
dat<-dat[colnames(dat) %in% marker_list]

library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "black"))
dat_bin<-data.frame(do.call("cbind",
        lapply(1:ncol(dat), 
               function(x) dat[,x]<-ifelse(dat[,x]>colMeans(dat)[x],1,0)
              )
       ))
colnames(dat_bin)<-colnames(dat)
row.names(dat_bin)<-row.names(dat)
dat<-dat_bin
#trying binary matrix again
dat<-dat>0
pdf("test_excneuron_layering.pdf")
Heatmap(as.matrix(dat),
        row_order=cell_row_order,
        show_row_names=F,
        col=col_fun,
        clustering_distance_columns="binary"
)
dev.off()



###Use detect genes for thresholding?



```

## Interneuron Clustering  


```python
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
library(EnsDb.Mmusculus.v79)
library(Matrix)


setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#Read in data and modify to monocle CDS file
#read in RDS file.
mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

#Subclustering inhibitory neurons with cisTopic and UMAP
mm10_inhib_neurons<-subset(mm10_atac,celltype=="inhibitory_neuron") #cluster 1 is inhibitory neurons

mm10_cistopic_counts_frmt<-mm10_inhib_neurons$peaks@counts
row.names(mm10_cistopic_counts_frmt)<-sub("-", ":", row.names(mm10_cistopic_counts_frmt))
mm10_atac_cistopic<-cisTopic::createcisTopicObject(mm10_cistopic_counts_frmt)
mm10_atac_cistopic_models<-cisTopic::runWarpLDAModels(mm10_atac_cistopic,topic=c(5,10,20:30,40,50,55,60:70),nCores=27,addModels=FALSE)
saveRDS(mm10_atac_cistopic_models,file="mm10_inhibitoryNeurons.CisTopicObject.Rds")
mm10_atac_cistopic_models<-readRDS("mm10_inhibitoryNeurons.CisTopicObject.Rds")

pdf("mm10_inhibitoryNeurons_model_selection.pdf")
par(mfrow=c(3,3))
mm10_atac_cistopic_models <- selectModel(mm10_atac_cistopic_models, type='maximum')
mm10_atac_cistopic_models <- selectModel(mm10_atac_cistopic_models, type='perplexity')
mm10_atac_cistopic_models <- selectModel(mm10_atac_cistopic_models, type='derivative')
dev.off()

#Loop through cistopic models
cistopic_loop<-function(topic_number,object_input=mm10_inhib_neurons,models_input=mm10_atac_cistopic_models){
    #select_topic
    cisTopicObject<-cisTopic::selectModel(models_input,select=topic_number,keepModels=F)
    
    #perform UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
    dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
    row.names(dims)<-colnames(topic_df)
    colnames(dims)<-c("x","y")
    dims$cellID<-row.names(dims)
    dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")

    #add cell embeddings
    cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
    colnames(cell_embeddings)<-cisTopicObject@cell.names
    n_topics<-nrow(cell_embeddings)
    row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
    cell_embeddings<-as.data.frame(t(cell_embeddings))
    
    #add feature loadings
    feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
    row.names(feature_loadings)<-paste0("topic_",1:n_topics)
    feature_loadings<-as.data.frame(t(feature_loadings))
    
    #combine with seurat object
    cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
    umap_dims<-as.data.frame(as.matrix(dims[2:3]))
    colnames(umap_dims)<-c("UMAP_1","UMAP_2")
    row.names(umap_dims)<-dims$cellID
    cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
    object_input@reductions$cistopic<-cistopic_obj
    object_input@reductions$umap<-cistopic_umap
    #finally recluster
    n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))
    plt<-DimPlot(object_input,size=0.1)+ggtitle(as.character(topic_number))
    return(plt)
}

library(parallel)
library(patchwork)
plt_list<-mclapply(mm10_atac_cistopic_models@calc.params$runWarpLDAModels$topic,FUN=cistopic_loop,mc.cores=5)
ggsave(wrap_plot(plt_list),file="mm10_inhibitoryneurons.modelchoice.pdf",height=20)


#set topics based on derivative
#selected topics subject to change
mm10_selected_topic=28
mm10_cisTopicObject<-cisTopic::selectModel(mm10_atac_cistopic_models,select=mm10_selected_topic,keepModels=T)
saveRDS(mm10_cisTopicObject,file="mm10_inhibitoryNeurons.CisTopicObject.Rds")

#perform UMAP on topics
mm10_topic_df<-as.data.frame(mm10_cisTopicObject@selected.model$document_expects)
row.names(mm10_topic_df)<-paste0("Topic_",row.names(mm10_topic_df))
mm10_dims<-as.data.frame(uwot::umap(t(mm10_topic_df),n_components=2))
row.names(mm10_dims)<-colnames(mm10_topic_df)
colnames(mm10_dims)<-c("x","y")
mm10_dims$cellID<-row.names(mm10_dims)
mm10_dims<-merge(mm10_dims,mm10_inhib_neurons@meta.data,by.x="cellID",by.y="row.names")

#add cell embeddings
mm10_cell_embeddings<-as.data.frame(mm10_cisTopicObject@selected.model$document_expects)
colnames(mm10_cell_embeddings)<-mm10_cisTopicObject@cell.names
mm10_n_topics<-nrow(mm10_cell_embeddings)
row.names(mm10_cell_embeddings)<-paste0("topic_",1:mm10_n_topics)
mm10_cell_embeddings<-as.data.frame(t(mm10_cell_embeddings))

#add feature loadings
mm10_feature_loadings<-as.data.frame(mm10_cisTopicObject@selected.model$topics)
row.names(mm10_feature_loadings)<-paste0("topic_",1:mm10_n_topics)
mm10_feature_loadings<-as.data.frame(t(mm10_feature_loadings))

#combine with seurat object
mm10_cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(mm10_cell_embeddings),loadings=as.matrix(mm10_feature_loadings),assay="peaks",key="topic_")
mm10_umap_dims<-as.data.frame(as.matrix(mm10_dims[2:3]))
colnames(mm10_umap_dims)<-c("UMAP_1","UMAP_2")
row.names(mm10_umap_dims)<-mm10_dims$cellID
mm10_cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(mm10_umap_dims),assay="peaks",key="UMAP_")
mm10_inhib_neurons@reductions$cistopic<-mm10_cistopic_obj
mm10_inhib_neurons@reductions$umap<-mm10_cistopic_umap

#finally recluster
mm10_n_topics<-ncol(Embeddings(mm10_inhib_neurons,reduction="cistopic"))

mm10_inhib_neurons <- FindNeighbors(
  object = mm10_inhib_neurons,
  reduction = 'cistopic',
  dims = 1:mm10_n_topics
)
mm10_inhib_neurons <- FindClusters(
  object = mm10_inhib_neurons,
  verbose = TRUE,
  resolution=0.1
)
mm10_inhib_neurons <- FindClusters(
  object = mm10_inhib_neurons,
  verbose = TRUE,
  resolution=0.2
)
mm10_inhib_neurons <- FindClusters(
  object = mm10_inhib_neurons,
  verbose = TRUE,
  resolution=0.5
)
mm10_inhib_neurons <- FindClusters(
  object = mm10_inhib_neurons,
  verbose = TRUE,
  resolution=0.9
)

###save Seurat files

saveRDS(mm10_inhib_neurons,file="mm10_inhibitoryNeurons.SeuratObject.Rds")

plt<-DimPlot(mm10_inhib_neurons,group.by=c('peaks_snn_res.0.1','peaks_snn_res.0.2','peaks_snn_res.0.5','peaks_snn_res.0.9',
                                  'barnyard_platevolume',
                                  'barnyard_storagecondition',
                                  'barnyard_pcrplate'))


ggsave(plt,file="mm10.inhibitoryneuron.umap.pdf",width=20)

pdf("mm10.inhibitoryneuron.umap.pdf",width=20)
plt
dev.off()
ggsave(plt,file="mm10.inhibitoryneuron.umap.png",width=20)

#Plotting topics, to be added to the notebook
library(ggplot2)
library(patchwork)

plt<-FeaturePlot(inhib,
pt.size=1,
features=colnames(inhib@reductions$cistopic@cell.embeddings))
ggsave(plt,file="inhibitoryneurons_topics.pdf",width=30,height=30)

########################
#New session#
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
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

mm10_inhib_neurons<-readRDS(file="mm10_inhibitoryNeurons.SeuratObject.Rds")

#set up an empty list for looping through
mm10_da_peaks<-list()

#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group){
    da_peaks_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T
        )
    da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
    closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
    da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
    da_peaks_tmp$enriched_group<-c(i)
    da_peaks_tmp$compared_group<-c("all_other_inhibneurons")
    return(da_peaks_tmp)
  }

#Perform parallel application of DA test
library(parallel)

n.cores=length(unique(mm10_inhib_neurons@meta.data$final_clusters))
n.cores=3
mm10_da_peaks<-mclapply(
    unique(mm10_inhib_neurons@meta.data$final_clusters),
    FUN=da_one_v_rest,
    obj=mm10_inhib_neurons,
    group="final_clusters",
    mc.cores=n.cores)

mm10_da_peaks<-do.call("rbind",mm10_da_peaks)

write("Outputting One v Rest DA Table.", stderr())
write.table(mm10_da_peaks,file="mm10.inhibNeurons.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

as.data.frame(mm10_da_peaks %>% group_by(enriched_group) %>% head(n=5,wt=p_val_adj))

#Checking if genes are properly listed
mm10_genes <- data.frame(GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),row.names=NULL)
unique(mm10_genes[mm10_genes$gene_name %in% c("Pvalb","Sst","Reln","Calb1","Calb2","Cck","Npy","Nos1","Vip"),]$gene_name)

mm10_da_peaks[mm10_da_peaks$gene_name %in% c("Pvalb","Sst","Reln","Calb1","Calb2","Cck","Npy","Nos1","Vip"),]


#Coverage plots of different marker genes
region_check<-function(i,object=mm10_inhib_neurons){
    return(nrow(as.data.frame(LookupGeneCoords(object=object,gene=gene_list[i]))))
    }

marker_plot<-function(j,k=celltype_name,l=mm10_inhib_neurons,m="mm10_inhibneurons_",n="final_clusters"){
    plt<-CoveragePlot(
        object = l,
        region = j,
        group.by=n,
        extend.upstream = 1000,
        extend.downstream = 1000,
        ncol = 1,
        tile=T,
        tile.size=500
    )
    pdf(paste0("./marker_sets/",m,k,"_",j,"_genebody_accessibility.pdf"))
    print(plt)
    dev.off()
}

marker_list<-c("Pvalb","Sst","Reln","Calb1","Calb2","Cck","Npy","Nos1","Vip")


for (x in 1:length(marker_list)){
    gene_list<-marker_list[[x]]
    gene_list<-gene_list[as.numeric(unlist(lapply(1:length(gene_list),FUN=region_check)))==1]
    celltype_name<-"inhibitory_neuron_classicmarkers"
    mclapply(gene_list,FUN=marker_plot,k=celltype_name,mc.cores=20) 
}

#########################New Session #####################################
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")

#Read in data and modify to monocle CDS file
#read in RDS file.
inhibneuron<-readRDS("mm10_inhibitoryNeurons.SeuratObject.Rds")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = T)
) ######species should be 10090??, this is listed in signac tutorial

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(inhibneuron),
  pwm = pfm,
  genome = 'mm10',
  use.counts = FALSE
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
inhibneuron <- SetAssayData(
  object = inhibneuron,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)

inhibneuron <- RegionStats(object = inhibneuron, 
    genome = BSgenome.Mmusculus.UCSC.mm10
)

inhibneuron <- RunChromVAR(
  object = inhibneuron,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(inhibneuron,"mm10_inhibitoryNeurons.SeuratObject.Rds")


#replacing naming scheme with TF names and plotting
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)
library(ggplot2)
library(parallel)
library(dplyr)
library(ggrepel)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data")
inhibneuron<-readRDS("mm10_inhibitoryNeurons.SeuratObject.Rds")

tfList <- getMatrixByID(JASPAR2020, ID=row.names(inhibneuron@assays$chromvar@data)) 
tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
row.names(inhibneuron@assays$chromvar@data)<-tfList

DefaultAssay(inhibneuron) <- 'chromvar'

plt <- FeaturePlot(
  object = inhibneuron,
  features = row.names(inhibneuron@assays$chromvar@data),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

ggsave(plt,file="inhibitoryNeuron_TF.pdf",width=25,height=500,limitsize=F)


DefaultAssay(inhibneuron) <- 'chromvar'

#set up an empty list for looping through
mm10_da_peaks<-list()

#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group){
    da_peaks_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T
        )
    da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
    da_peaks_tmp$enriched_group<-c(i)
    da_peaks_tmp$compared_group<-c("all_other_cells")
    return(da_peaks_tmp)
  }


da_one_v_one<-function(i,obj,group,j_list){
    i<-as.character(i)
    da_tmp_2<-list()
    for (j in j_list){
        if ( i != j){
        da_peaks_tmp <- FindMarkers(
            object = obj,
            ident.1 = i,
            ident.2 = j,
            group.by = group,
            test.use = 'LR',
            latent.vars = 'nCount_peaks',
            only.pos=T
            )
        da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
        da_peaks_tmp$enriched_group<-c(i)
        da_peaks_tmp$compared_group<-c(j)
        da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
        }
    }
    return(da_tmp_2)
  }

#Perform parallel application of DA test


n.cores=length(unique(inhibneuron@meta.data$final_clusters))
mm10_tf<-mclapply(
    unique(inhibneuron@meta.data$final_clusters),
    FUN=da_one_v_rest,
    obj=inhibneuron,
    group="final_clusters",
    mc.cores=n.cores)

#Merge the final data frame from the list for 1vrest DA
mm10_tf<-do.call("rbind",mm10_tf)

write("Outputting One v Rest DA Table.", stderr())
write.table(mm10_tf,file="mm10.inhibneuron.onevrest.tffactors.txt",sep="\t",col.names=T,row.names=T,quote=F)

dat<-read.table("mm10.inhibneuron.onevrest.tffactors.txt",header=T,sep="\t")

dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% slice(1:2) #grabbing top 2 most significant peaks to label
plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+
geom_point(aes(alpha=0.1))+
geom_label_repel(dat=dat_select,aes(label=da_region),force=3)+
theme_bw()

ggsave(plt,file="mm10_inhibtneuron.tffactors.pdf")


#Empty list to rerun for 1v1 comparisons
mm10_tf<-list()
 
n.cores=length(unique(inhibneuron@meta.data$final_clusters))
mm10_tf<-mclapply(
    unique(inhibneuron@meta.data$final_clusters),
    FUN=da_one_v_one,
    obj=inhibneuron,
    group="final_clusters",
    j_list=do.call("as.character",list(unique(inhibneuron@meta.data$final_clusters))),
    mc.cores=n.cores)

#Merge the final data frame from the list for 1v1 DA
mm10_tf<-do.call("rbind",do.call("rbind",mm10_tf))

write("Outputting One v One DA Table.", stderr())
write.table(mm10_tf,file="mm10.inhibneurons.onevone.tffactors.txt",sep="\t",col.names=T,row.names=T,quote=F)


```

https://satijalab.org/seurat/v3.0/integration.html

Integrating the two data sets. First making seurat objects.


```python
    library(Seurat)
    library(cicero)
    setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200218_s3ATAC_brainbarnyard_reprocessing")
    hg38_ga<-readRDS("hg38_unnorm_GA.Rds")
    mm10_ga<-readRDS("mm10_unnorm_GA.Rds")
    mm10_atac<-readRDS("mm10_SeuratObject.Rds")
    hg38_atac<-readRDS("hg38_SeuratObject.Rds")

    #limit seurat object cells to match those of the filtered GA cells
 hg38_atac<-subset(hg38_atac,cells=which(hg38_atac$cellID %in% colnames(hg38_ga)))
 mm10_atac<-subset(mm10_atac,cells=which(mm10_atac$cellID %in% colnames(mm10_ga)))
 row.names(mm10_ga)<-toupper(row.names(mm10_ga))
 
 hg38_atac[["ACTIVITY"]] <- CreateAssayObject(counts = hg38_ga)
 mm10_atac[["ACTIVITY"]] <- CreateAssayObject(counts = mm10_ga)
DefaultAssay(hg38_atac)<-"ACTIVITY"
DefaultAssay(mm10_atac)<-"ACTIVITY"

atac.list<-merge(hg38_atac,mm10_atac)
atac.list <- SplitObject(atac.list, split.by = "called_species")
atac.list <- lapply(X = atac.list, FUN = function(x) {
    x <- NormalizeData(x, assay="ACTIVITY",verbose = FALSE)
    x <- FindVariableFeatures(x, assay="ACTIVITY",verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = atac.list,nfeatures=5000)
atac.list <- lapply(X = atac.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = TRUE)
})
anchors <- FindIntegrationAnchors(object.list = atac.list, reference = c(1, 2), reduction = "cca",
    dims = 1:10)
atac.integrated <- IntegrateData(anchorset = anchors, dims = 1:10)
atac.integrated <- ScaleData(atac.integrated, verbose = FALSE)
atac.integrated <- RunPCA(atac.integrated, verbose = FALSE)
atac.integrated <- RunUMAP(atac.integrated, dims = 1:10)
pdf("integrated_species.umap.pdf")
DimPlot(atac.integrated, group.by = "called_species")
dev.off()



    # remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0,
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
unnorm_ga2 <- unnorm_ga
cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2),
                                                    num_genes)
    
```

Further code, not yet implemented


```python

    
    
    
    
	#using umap coordinates stored in the Seurat object
	umap_coords<-Embeddings(hg38_atac,"umap")

	#running cicero aggregation function
	cicero_cds<-cicero::make_cicero_cds(cds,k=50,reduced_coordinates=umap_coords)

	#generate list of chromosome lengths for bounding
	library(BSgenome.Hsapiens.UCSC.hg38)
	chr_len_df<-as.data.frame(seqlengths(BSgenome.Hsapiens.UCSC.hg38))
	chr_list<-row.names(chr_len_df)[!grepl("_",row.names(chr_len_df))]
	chr_len_df<-chr_len_df[row.names(chr_len_df)%in%chr_list,]
	chr_len_df<-data.frame(cbind(chr_list,chr_len_df))
	colnames(chr_len_df)<-c("V1","V2")
	chr_len_df$V2<-as.numeric(as.character(chr_len_df$V2))

	conns<-run_cicero(cicero_cds,chr_len_df)
	saveRDS(conns,"hg38_conns.Rds")

	#plot_connections(conns, "chr2", 9773451, 9848598,
	#                 gene_model = gene_anno,
	#                 coaccess_cutoff = .25,
	#                 connection_width = .5,
	#                 collapseTranscripts = "longest" )

	CCAN_assigns <- generate_ccans(conns)
	saveRDS(CCAN_assigns,"hg38_CCANS.Rds")

	# Download the GTF associated with this data (hg38) from ensembl and load it
	# using rtracklayer

	# download and unzip
	temp <- tempfile()
	download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/homo_sapiens/Homo_sapiens.GRCh37.65.gtf.gz", temp)
	gene_anno <- rtracklayer::readGFF(temp)
	unlink(temp)

	# rename some columns to match requirements
	gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
	gene_anno$gene <- gene_anno$gene_id
	gene_anno$transcript <- gene_anno$transcript_id
	gene_anno$symbol <- gene_anno$gene_name

	#### Add a column for the pData table indicating the gene if a peak is a promoter ####
	# Create a gene annotation set that only marks the transcription start sites of
	# the genes. We use this as a proxy for promoters.
	# To do this we need the first exon of each transcript
	pos <- subset(gene_anno, strand == "+")
	pos <- pos[order(pos$start),]

	# remove all but the first exons per transcript
	pos <- pos[!duplicated(pos$transcript),]

	# make a 1 base pair marker of the TSS
	pos$end <- pos$start + 1
	neg <- subset(gene_anno, strand == "-")
	neg <- neg[order(neg$start, decreasing = TRUE),]

	# remove all but the first exons per transcript
	neg <- neg[!duplicated(neg$transcript),]
	neg$start <- neg$end - 1

	gene_annotation_sub <- rbind(pos, neg)

	# Make a subset of the TSS annotation columns containing just the coordinates
	# and the gene name
	gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

	# Rename the gene symbol column to "gene"
	names(gene_annotation_sub)[4] <- "gene"

	cds <- annotate_cds_by_site(cds, gene_annotation_sub)
	saveRDS(cds,"hg38_cds.Rds")
	#### Generate gene activity scores ####
	# generate unnormalized gene activity matrix
	unnorm_ga <- build_gene_activity_matrix(cds, conns)
	saveRDS(unnorm_ga,"hg38_unnorm_GA.Rds")
```

### Cicero on mouse cells
```{r echo=TRUE, eval=FALSE}
library(Signac)
	library(cicero)
	library(monocle)
	library(monocle3)
	set.seed(123)
	setwd("C:/Users/grima_000/Desktop/s3_ATAC")

	#Read in data and modify to monocle CDS file
	#Signac seurat objects not yet supported
	mm10_atac<-readRDS("mm10_SeuratObject.Rds")
	data<-mm10_atac@assays$peaks@counts
	fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
	cds <- new_cell_data_set(as.matrix(data), cell_metadata = mm10_atac@meta.data, gene_metadata =fData)

	#using umap coordinates stored in the Seurat object
	umap_coords<-Embeddings(mm10_atac,"umap")

	#running cicero aggregation function
	cicero_cds<-cicero::make_cicero_cds(cds,k=50,reduced_coordinates=umap_coords)

	#generate list of chromosome lengths for bounding
	library(BSgenome.Mmusculus.UCSC.mm10)
	chr_len_df<-as.data.frame(seqlengths(BSgenome.Mmusculus.UCSC.mm10))
	chr_list<-row.names(chr_len_df)[!grepl("_",row.names(chr_len_df))]
	chr_len_df<-chr_len_df[row.names(chr_len_df)%in%chr_list,]
	chr_len_df<-data.frame(cbind(chr_list,chr_len_df))
	colnames(chr_len_df)<-c("V1","V2")
	chr_len_df$V2<-as.numeric(as.character(chr_len_df$V2))

	conns<-run_cicero(cicero_cds,chr_len_df)
	saveRDS(conns,"mm10_conns.Rds")

	#plot_connections(conns, "chr2", 9773451, 9848598,
	#                 gene_model = gene_anno,
	#                 coaccess_cutoff = .25,
	#                 connection_width = .5,
	#                 collapseTranscripts = "longest" )

	CCAN_assigns <- generate_ccans(conns)
	saveRDS(CCAN_assigns,"mm10_CCANS.Rds")

	#done to here

	# Download the GTF associated with this data (mm10) from ensembl and load it
	# using rtracklayer

	# download and unzip
	temp <- tempfile()
	download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz", temp)
	gene_anno <- rtracklayer::readGFF(temp)
	unlink(temp)

	# rename some columns to match requirements
	gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
	gene_anno$gene <- gene_anno$gene_id
	gene_anno$transcript <- gene_anno$transcript_id
	gene_anno$symbol <- gene_anno$gene_name

	#### Add a column for the pData table indicating the gene if a peak is a promoter ####
	# Create a gene annotation set that only marks the transcription start sites of
	# the genes. We use this as a proxy for promoters.
	# To do this we need the first exon of each transcript
	pos <- subset(gene_anno, strand == "+")
	pos <- pos[order(pos$start),]

	# remove all but the first exons per transcript
	pos <- pos[!duplicated(pos$transcript),]

	# make a 1 base pair marker of the TSS
	pos$end <- pos$start + 1
	neg <- subset(gene_anno, strand == "-")
	neg <- neg[order(neg$start, decreasing = TRUE),]

	# remove all but the first exons per transcript
	neg <- neg[!duplicated(neg$transcript),]
	neg$start <- neg$end - 1

	gene_annotation_sub <- rbind(pos, neg)

	# Make a subset of the TSS annotation columns containing just the coordinates
	# and the gene name
	gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

	# Rename the gene symbol column to "gene"
	names(gene_annotation_sub)[4] <- "gene"

	cds <- annotate_cds_by_site(cds, gene_annotation_sub)
	saveRDS(cds,"mm10_cds.Rds")
	#### Generate gene activity scores ####
	# generate unnormalized gene activity matrix
	unnorm_ga <- build_gene_activity_matrix(cds, conns)
	saveRDS(unnorm_ga,"mm10_unnorm_GA.Rds")

```


# Loading in cicero output for cluster analysis
```{r echo=TRUE, eval=FALSE}

	export R_LIBS_USER='/home/groups/oroaklab/src/R/R-3.5.1/library_arsn'

	setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/191205_s3ATAC_brainbarnyard/200114_Reprocessing")

	library(Signac)
	library(Seurat)
	library(GenomeInfoDb)
	library(ggplot2)
	set.seed(1234)
	library(EnsDb.Hsapiens.v86)
	library(EnsDb.Mmusculus.v79)

	hg38_atac<-readRDS(file="hg38_SeuratObject.Rds")
	mm10_atac<-readRDS(file="mm10_SeuratObject.Rds")

	cicero_processing_dir<-"/home/groups/oroaklab/adey_lab/projects/sciWGS/191205_s3ATAC_brainbarnyard/200114_Reprocessing/s3_ATAC_Cicero_processing"
	hg38_unnorm_GA<-readRDS(file=paste(cicero_processing_dir,"hg38_unnorm_GA.Rds",sep="/"))
	mm10_unnorm_GA<-readRDS(file=paste(cicero_processing_dir,"mm10_unnorm_GA.Rds",sep="/"))


	hg38_atac[["unnorm_GA"]]<-CreateAssayObject(counts=as.data.frame(hg38_unnorm_GA))
	mm10_atac[["unnorm_GA"]]<-CreateAssayObject(counts=as.data.frame(mm10_unnorm_GA))

	hg38_atac <- NormalizeData(object = hg38_atac,assay = 'unnorm_GA',normalization.method = 'LogNormalize',scale.factor = median(hg38_atac$nCount_unnorm_GA))
	mm10_atac <- NormalizeData(object = mm10_atac,assay = 'unnorm_GA',normalization.method = 'LogNormalize',scale.factor = median(hg38_atac$nCount_unnorm_GA))

	DefaultAssay(hg38_atac) <- 'unnorm_GA'
	DefaultAssay(mm10_atac) <- 'unnorm_GA'

	#look at cluster marker gene activity scores
	hg38_atac.markers <- FindAllMarkers(hg38_atac,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	write.table(hg38_atac.markers,"hg38_atac.GA.markers.txt",col.names=T,quote=F,sep="\t")
	mm10_atac.markers <- FindAllMarkers(mm10_atac,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	write.table(mm10_atac.markers,"mm10_atac.GA.markers.txt",col.names=T,quote=F,sep="\t")

	pdf("hg38_GA_markerfeats.pdf",width=100,height=100)
	FeaturePlot(
	  object = hg38_atac,
	  features = row.names(hg38_atac.markers),
	  pt.size = 2,
	  max.cutoff = 'q95',
	)
	dev.off()


	pdf("mm10_GA_markerfeats.pdf",width=100,height=100)
	FeaturePlot(
	  object = mm10_atac,
	  features = row.names(mm10_atac.markers),
	  pt.size = 2,
	  max.cutoff = 'q95',
	)
	dev.off()

	saveRDS(hg38_atac,file="hg38_SeuratObject.Rds")
	saveRDS(mm10_atac,file="mm10_SeuratObject.Rds")

	#################################################################
	#Coembedding of ATAC GA and RNA
	#################################################################
	#https://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html
	export R_LIBS_USER='/home/groups/oroaklab/src/R/R-3.5.1/library_arsn'
	library(Signac)
	library(Seurat)
	library(GenomeInfoDb)
	library(EnsDb.Hsapiens.v86)
	library(ggplot2)
	set.seed(1234)
	orgo_atac<-readRDS(file="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/180703_sciATAC_Organoids/NatureLetterData/annot/180707_sciATAC_Organoids.filt.nodigi.bbrd.q10.filt.500.counts.cistopic.k2000.pg.geneactivitySeurat.Rds")
	orgo_rna<-readRDS(file="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids/orgo_rna.Rds")
	# Load the pre-processed scRNA-seq data

	transfer.anchors <- FindTransferAnchors(
	  reference = orgo_rna,
	  query = orgo_atac,
	  reduction = 'cca',
	)

	predicted.labels <- TransferData(
	  anchorset = transfer.anchors,
	  refdata = orgo_rna$RNA_snn_res.0.5,
	  weight.reduction = orgo_atac[['lsi']],
	)

	orgo_atac <- AddMetaData(object = orgo_atac, metadata = predicted.labels)
	DefaultAssay(orgo_atac) <- 'peaks'

	plot1 <- DimPlot(
	  object = orgo_rna,
	  group.by = 'RNA_snn_res.0.5',
	  label = TRUE,
	  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
	
	plot2 <- DimPlot(
	  object = orgo_atac,
	  group.by = 'predicted.id',
	  reduction="umap",
	  label = TRUE,
	  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

	pdf("integrated.data.pdf")
	CombinePlots(list(plot1,plot2), ncol = 2)
	dev.off()
```

```

--->
