---
title: gccACT-seq Wafergen Test
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /gccACT_wafergen/
category: mdanderson
---

## Initial analysis of Nextseq 2000 Run
Ran P1 100 cycle kit with 100% allocated to a gccACT wafergen test (split 50|50 for SKBR3 and MDA-MB-231)

Stored here:
```bash
/volumes/seq/flowcells/MDA/nextseq2000/2023/20230629_Mariam_HiC_MDA231_SKor3/

#make symbolic link to working directory
cd /volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest
ln -s /volumes/seq/flowcells/MDA/nextseq2000/2023/20230629_Mariam_HiC_MDA231_SKor3/

#located here:
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/20230629_Mariam_HiC_MDA231_SKor3
```

```bash
#run transfered to /volumes/seq/tmp
cd /volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/20230629_Mariam_HiC_MDA231_SKor3
bcl2fastq -R /volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/20230629_Mariam_HiC_MDA231_SKor3/230628_VH00219_441_AACN3HKM5 \
-o /volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest \
-r 4 \
-p 10 \
-w 4 \
--ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls \
--create-fastq-for-index-reads
```

Split out gccACT reads via python script

```python
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd

#set up plate
n7="/volumes/lab/users/wet_lab/protocols/WaferDT/Barcodes/barcode_txt/wafer_v2_N7_72.txt" #copied to home directory as backup
s5="/volumes/lab/users/wet_lab/protocols/WaferDT/Barcodes/barcode_txt/wafer_v2_S5_72.txt" #copied to home directory as backup
n7=pd.read_table(n7)
s5=pd.read_table(s5)
plate_i7=n7["RC_N7"]
plate_i5=s5["RC_S5_Hiseq4000_nextseq"]

plate_i7=[i.strip() for i in plate_i7]
plate_i5=[i.strip() for i in plate_i5]


fq1=sys.argv[1] 
fq2=sys.argv[2] 
idx3=sys.argv[3] 
idx4=sys.argv[4] 
#fq1="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R1_001.fastq.gz"
#fq2="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R2_001.fastq.gz"
#idx3="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_I1_001.fastq.gz"
#idx4="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_I2_001.fastq.gz"

i=0
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1:
	with gzip.open(fq2, "rt") as handle2:
		with gzip.open(idx3, "rt") as handle3:
			with gzip.open(idx4, "rt") as handle4:
				with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
					with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
						for (title1, seq1, qual1), (title2, seq2, qual2), (title3,seq3,qual3), (title3,seq4,qual4) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2),FastqGeneralIterator(handle3),FastqGeneralIterator(handle4)):
							if seq3[:8] in plate_i7 and seq4[:8] in plate_i5:
								i+=1									
								readname=seq3[:8]+seq4[:8]
								outfile_fq1.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq1, qual1))
								outfile_fq2.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq2, qual2))

```
Running fastq splitter
This can be made a lot faster using something like a hash table, or limiting search space more.

```bash
python ~/src/plate2_fastqsplitter.py \
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R1_001.fastq.gz \
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R2_001.fastq.gz \
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_I1_001.fastq.gz \
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_I2_001.fastq.gz

#then gzip
for i in *fastq; do gzip $i & done &
```
Got 116,212,259 reads total from assignment
Moving fastq.gz files to a project directory

Processing will be in: 
```bash
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest
```

# Processing Samples via Copykit to start
## Alignment
```bash
#set up variables and directory
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
fq1="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R1_001.barc.fastq.gz"
fq2="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/Undetermined_S0_L001_R2_001.barc.fastq.gz"

#Map reads with BWA Mem
bwa mem -t 50 $ref $fq1 $fq2 | samtools view -b - > $dir/230701_wafergen_gccact.bam
```


## Split out single-cells
```bash
#set up cell directory
mkdir $dir/cells
```
Split by readname bam field into cells subdir and add metadata column

```bash
samtools view $dir/230701_wafergen_gccact.bam | awk -v dir=$dir 'OFS="\t" {split($1,a,":"); print $0,"XM:Z:"a[1] > "./cells/"a[1]".230701_wafergen_gccact.sam"}' &
 #split out bam to cell level sam
```

## Remove sam files that have less than 100000 reads

```bash
wc -l *sam | sort -k1,1n - | head -n -1 | awk '{if($1<100000) print $2}' | xargs rm
#586 cells returned
```

## Add header to each sam and convert to bam and sort
Using parallel to save time
```bash

add_header() {
	ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
	dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
	samtools view -bT $ref $1 | samtools sort -T $dir -o ${1::-4}.bam -
}
export -f add_header
sam_in=`ls *sam`
parallel --jobs 50 add_header ::: $sam_in
```

## Mark duplicate reads
```bash
#name sort, fix mates, sort by position, mark dup
sort_and_markdup() {
  samtools sort -T . -n -o - $1 | samtools fixmate -m - -| samtools sort -T . -o - - | samtools markdup -s - ${1::-4}.rmdup.bam 2> ${1::-4}.rmdup.stats.txt
}
export -f sort_and_markdup

bam_in=`ls *bam`
parallel --jobs 50 sort_and_markdup ::: $bam_in
```

## Run Fastqc on everything and clean up

```bash
bam_in=`ls *rmdup.bam`
parallel --jobs 30 fastqc ::: $bam_in

mkdir $dir/cells/fastqc
mv *fastqc* $dir/cells/fastqc
mv *stats.txt $dir/cells/fastqc
rm -rf *gccact.bam #only keep duplicate marked bams

#run multiqc to aggregate
multiqc . 

```

## Run CopyKit for WGS portion
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/cells")

tumor2 <- runVarbin("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
tumor2 <- findAneuploidCells(tumor2)

# Mark low-quality cells for filtering
tumor2 <- findOutliers(tumor2)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(tumor2, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

tumor<-tumor2

# kNN smooth profiles
tumor <- knnSmooth(tumor)


# Create a umap embedding 
tumor <- runUmap(tumor)
k_clones<-findSuggestedK(tumor) #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
tumor  <- findClusters(tumor,k_subclones=16)#output from k_clones

pdf("subclone.umap.pdf")
plotUmap(tumor, label = 'subclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
tumor <- calcConsensus(tumor)
tumor <- runConsensusPhylo(tumor)

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.pdf")
plotHeatmap(tumor, label = 'subclones',order='hclust')
dev.off()

saveRDS(tumor,file="/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
tumor<-readRDS("/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
clone_out<-data.frame(bam=paste0(row.names(tumor@colData),".bam"),clone=tumor@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("clone_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}

#based on CNV profiles, clone c3 is skbr3 and all others are mda-mb-231.
#going to split and rerun subclone analysis and heatmaping
# Remove cells marked as low-quality and/or aneuploid from the copykit object
skbr3<- tumor[,SummarizedExperiment::colData(tumor)$subclones == "c3"]
mdamb231 <- tumor[,SummarizedExperiment::colData(tumor)$subclones != "c3"]


# Create a umap embedding 
skbr3 <- runUmap(skbr3) ; mdamb231<- runUmap(mdamb231)
k_clones_skbr3<-findSuggestedK(skbr3); k_clones_mdamb231<-findSuggestedK(mdamb231); #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
skbr3 <- findClusters(skbr3,k_subclones=2)#output from k_clones
mdamb231 <- findClusters(mdamb231,k_subclones=4)#output from k_clones

pdf("skbr3_subclone.umap.pdf")
plotUmap(skbr3, label = 'subclones')
dev.off()

pdf("mdamb231_subclone.umap.pdf")
plotUmap(mdamb231, label = 'subclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
skbr3 <- calcConsensus(skbr3) ; mdamb231<- calcConsensus(mdamb231)
skbr3 <- runConsensusPhylo(skbr3); mdamb231 <- runConsensusPhylo(mdamb231)

# Plot a copy number heatmap with clustering annotation
pdf("skbr3_subclone.heatmap.pdf")
plotHeatmap(skbr3, label = 'subclones',order='hclust')
dev.off()

pdf("mdamb231_subclone.heatmap.pdf")
plotHeatmap(mdamb231, label = 'subclones',order='hclust')
dev.off()

```

### Count of WGS and GCC Reads
```bash
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
#count reads based on paired alignments
readtype_count() {
	#print out reads on same chromosome that are >= 1000 bp apart (distal)
	distal=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))>=1000) print $0}' | wc -l)
	#print out reads on different chromosomes
	trans=$(samtools view -F 1024 $1 | awk '{if($7 != "=") print $0}' | wc -l)
	#print out reads on same chromosome within 1000bp (cis)
	near=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))<=1000) print $0}' | wc -l)
	echo $1,$near,$distal,$trans
}
export -f readtype_count

cd $dir/cells
bam_in=`ls *bam`
echo "cellid,near_cis,distal_cis,trans" > read_count.csv; parallel --jobs 20 readtype_count ::: $bam_in >> read_count.csv
```

Plotting GCC read types
```R
library(ggplot2)
library(patchwork)
setwd("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest")
dat<-read.table("./cells/read_count.csv",header=T,sep=",")
dat$total_reads<-dat$near_cis+dat$distal_cis+dat$trans

plt1<-ggplot(dat,aes(y=near_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis Reads")+theme_minimal()
plt2<-ggplot(dat,aes(y=distal_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis Reads")+theme_minimal()
plt3<-ggplot(dat,aes(y=trans,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans Reads")+theme_minimal()

plt4<-ggplot(dat,aes(y=(near_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis % Reads")+theme_minimal()+ylim(c(0,100))
plt5<-ggplot(dat,aes(y=(distal_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis % Reads")+theme_minimal()+ylim(c(0,10))
plt6<-ggplot(dat,aes(y=(trans/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans % Reads")+theme_minimal()+ylim(c(0,10))

plt<-(plt1|plt2|plt3)/(plt4|plt5|plt6)
ggsave(plt,file="./cells/read_counts.pdf")
```


## Project Library Complexity
Using Picard Tools
```bash
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
project_count() {
java -jar /volumes/seq/code/3rd_party/picard/picard-2.20.4/picard.jar EstimateLibraryComplexity I=$1 O=${1::-4}.complex_metrics.txt
}
export -f project_count

cd $dir/cells
bam_in=`ls *bam`
parallel --jobs 20 project_count ::: $bam_in
```

<!--
### Count of WGS and GCC Reads Ligation Signatures
```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
#count reads based on paired alignments
readtype_count() {
	#print out reads on same chromosome that are >= 1000 bp apart (distal)
	distal=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))>=1000) print $0}' | grep "REmotif" | wc -l)
	#print out reads on different chromosomes
	trans=$(samtools view -F 1024 $1 | awk '{if($7 != "=") print $0}' | grep "REmotif" | wc -l)
	#print out reads on same chromosome within 1000bp (cis)
	near=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))<=1000) print $0}' |grep "REmotif" |  wc -l)
	echo $1,$near,$distal,$trans
}
export -f readtype_count

cd $dir/cells
bam_in=`ls *bam`
echo "cellid,near_cis,distal_cis,trans" > read_count.csv; parallel --jobs 10 readtype_count ::: $bam_in >> read_count.csv
```
-->

### Generation of HiC Contact Matrices
Merge bam files based on CopyKit output. Then using bam2pairs from pairix to generate contacts
```bash
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

#filter reads to HiC contacts, then perform bam2pairs 
mkdir $dir/cells/contacts

#### For generation of pairix per cell
#bam_to_pairs() {
#	samtools view -F 1024 $3 | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $2 - > $1/cells/contacts/${3::-4}.contacts.bam && wait;
#	~/tools/pairix/util/bam2pairs/bam2pairs $1/cells/contacts/${3::-4}.contacts.bam $1/cells/contacts/${3::-4}
#}
#export -f bam_to_pairs

#cd $dir/cells
#bam_in=`ls *bam`
#parallel --jobs 50 bam_to_pairs $dir $ref {} ::: $bam_in &
# first argument is directory
# second is reference fasta
# third is bam file input

bamlist_merge_to_pairs() {
	dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
	ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
	samtools merge -b $1 -O SAM -@ 20 - | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $ref - > $dir/cells/contacts/${1::-13}.contacts.bam && wait;
	~/tools/pairix/util/bam2pairs/bam2pairs $dir/cells/contacts/${1::-13}.contacts.bam $dir/cells/contacts/${1::-13}
}
export -f bamlist_merge_to_pairs

cd $dir/cells
bamlist_merge_to_pairs clone_c1.bam_list.txt 
bamlist_merge_to_pairs clone_c2.bam_list.txt 
bamlist_merge_to_pairs clone_c3.bam_list.txt #skbr3
bamlist_merge_to_pairs clone_c4.bam_list.txt 

```


## Using cooler, cooltools, and EagleC to detect Structural Variants

Cooler is both python line and command line, using command line for this
https://github.com/open2c/cooler
https://github.com/open2c/cooltools

*EagleC* is a deep learning method of SV detection in bulk and single cells at multiple resolutions
https://github.com/XiaoTaoWang/EagleC

EagleC also makes use of the toolkit form the same authors called *NeoLoopFinder*, used to CNV correction in HiC Data
https://github.com/XiaoTaoWang/NeoLoopFinder

Change to cooler environment ::sunglasses::
```bash
conda deactivate #get out of r3.4 env
conda activate cooler_env #use cooler env (lower python version)
```

### Set up reference and bins
Prepare chrom.sizes file from reference fasta
```bash
#ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
#faidx $ref -i chromsizes > hg38.chrom.sizes
#previously done
```

Prepare bin of genome and the GC bins
```bash
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
FASTA_PATH=$ref
CHROMSIZES_FILE="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/hg38.chrom.sizes"

#1MB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.gc.bins"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 1000000 > $BINS_PATH #1mb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

#5KB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.gc.bins"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 5000 > $BINS_PATH  #5kb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

#10KB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.gc.bins"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 10000 > $BINS_PATH #10kb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

#50KB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.gc.bins"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 50000 > $BINS_PATH #50kb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

#500KB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.gc.bins"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 500000 > $BINS_PATH #50kb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

```

### Generate Cooler matrices from pairix data
```bash
# Note that the input pairs file happens to be space-delimited, so we convert to tab-delimited with `tr`.
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"

pairix_to_cooler_5kb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins.bed"
out_name="5kb"
cooler cload pairs -0 -c1 2 -c2 4 -p1 3 -p2 5 --temp-dir . --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_5kb

pairix_to_cooler_10kb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins.bed"
out_name="10kb"
cooler cload pairs -0 -c1 2 -c2 4 -p1 3 -p2 5 --temp-dir . --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_10kb

pairix_to_cooler_50kb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins.bed"
out_name="50kb"
cooler cload pairs -0 -c1 2 -c2 4 -p1 3 -p2 5 --temp-dir . --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_50kb


pairix_to_cooler_500kb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
out_name="500kb"
cooler cload pairs -0 -c1 2 -c2 4 -p1 3 -p2 5 --temp-dir . --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_500kb


cd $dir/cells/contacts
pairix_in=`ls clone_*pairs.gz`

#5kb
parallel --jobs 5 pairix_to_cooler_5kb ::: $pairix_in & #uses 10 cores per job

#10kb
parallel --jobs 5 pairix_to_cooler_10kb ::: $pairix_in & #uses 10 cores per job

#50kb
parallel --jobs 5 pairix_to_cooler_50kb ::: $pairix_in & #uses 10 cores per job

#500kb
parallel --jobs 5 pairix_to_cooler_500kb ::: $pairix_in & #uses 10 cores per job

# first argument is pairix gzipped file
```

<!--
### CNV normalization on HiC Data
Write out CNV segments as bedfiles at multiple resolutions for cnv correction via NeoLoopFinder
NeoLoopFinder also reports CNVs through log2 changes, making this a direct comparison.

```R
library(copykit)
library(GenomicRanges)
library(parallel)
wd_out="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"

setwd(wd_out)
tumor<-readRDS("/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/scCNA.rds")
tumor <- calcInteger(tumor, method = 'scquantum', assay = 'smoothed_bincounts') #calculate ploidy
tumor <- calcConsensus(tumor) #generate consensus

pdf("cellploidy.integer.pdf")
plotMetrics(tumor, metric = 'ploidy', label = 'ploidy_score') #plot ploidy score
dev.off()

tumor <- calcInteger(tumor, method = 'fixed', ploidy_value = median(tumor$ploidy)) #calculate integer value per consensus
tumor <- calcConsensus(tumor, consensus_by = 'subclones', assay = 'integer')


# Plot a consensus copy number heatmap 
pdf("consensus.heatmap.integer.pdf")
plotHeatmap(tumor,
            consensus = TRUE,
            label = 'subclones',
            assay = 'integer')
dev.off()
saveRDS(tumor,file="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/scCNA.consensus.rds") #save categorical consensus scCNA object


#function to define blacklist regions
make_black_list<-function(bed_in,copykit_obj){
#bins in bin bed file that are not in the copykit filtered list are considered black list,
#generate per bins.bed resolution
bins_bed<-read.table(bed_in,header=F,sep="\t")
colnames(bins_bed)<-c("chr","start","end")
bins_bed<-GRanges(bins_bed)
accepted_bed_ranges<-GRanges(copykit_obj@rowRanges) #row ranges for bedgraph
overlaps<-findOverlaps(query=bins_bed,subject=accepted_bed_ranges,minoverlap=1,ignore.strand=T)
blacklist<-bins_bed[!(range(1,nrow(bins_bed)) %in% overlaps@from),]
write.table(as.data.frame(blacklist)[1:3],col.names=F,row.names=F,quote=F,file=paste0(substr(bed_in,1,nchar(bed_in)-3),"blacklist.bed"))
print(paste("Wrote out",paste0(substr(bed_in,1,nchar(bed_in)-3),"blacklist.bed")))
}

#function to generate consensus bedGraphs for CNV correction.

#function to define blacklist regions
make_consensus_bedgraph<-function(bed_in,res,copykit_obj,clone,cores=10){
	#bins in bin bed file that are not in the copykit filtered list are considered black list,
	#generate per bins.bed resolution granges
	bins_bed<-read.table(bed_in,header=F,sep="\t")
	colnames(bins_bed)<-c("chr","start","end")
	bins_bed<-GRanges(bins_bed)

	#generate granges of copykit data
	accepted_bed_ranges<-GRanges(copykit_obj@rowRanges) #row ranges for bedgraph
	accepted_bed_ranges$cnv<-copykit_obj@consensus[clone] #add cnv data as metadata column

	#overlap the two granges
	hits<-findOverlaps(query=bins_bed,subject=accepted_bed_ranges,minoverlap=1,ignore.strand=T) #get list of intersecting
	overlaps <- pintersect(bins_bed[queryHits(hits)], accepted_bed_ranges[subjectHits(hits)]) #combine intersecting sites for overlap
	bins_bed$cnv<-mean(unlist(accepted_bed_ranges$cnv)) #set up column meta data, set to average for blacklist bed regions

	#change overlap calculation to correct for it one if bigger than the other (i think I did this right??)
	print(paste("Calculating window CNVs from CopyKit Clone consensus for",clone,"at",res))
	if(mean(width(accepted_bed_ranges))<mean(width(bins_bed))){ #if hic bins are bigger than copykit windows
		percentOverlap <- width(overlaps) / width(bins_bed[queryHits(hits)])
	} else { #if copykit windows are bigger than hic bins
		percentOverlap <- width(overlaps) / width(accepted_bed_ranges[subjectHits(hits)])}
	#lapply to 
	out<-mclapply(unique(hits@from),function(x) {
		overlaps_tmp<-hits[hits@from==x,]
		cnv_tmp<-unlist(as.data.frame(accepted_bed_ranges[overlaps_tmp@to,])[clone])
		weight_tmp<-percentOverlap[overlaps_tmp@to]
		cnv_out<-weighted.mean(cnv_tmp,w=weight_tmp) #weighted mean score by overlap percentage
		bins_bed[x,]$cnv<-cnv_out
		return(bins_bed[x,])},mc.cores=cores)
	bins_bed_cnv<-do.call("c",out)
	bins_bed_cnv<-c(bins_bed_cnv,bins_bed[-subjectHits(findOverlaps(query=bins_bed_cnv, subject=bins_bed, minoverlap=10)),] )
	bins_bed_cnv<-sortSeqlevels(bins_bed_cnv)
	bins_bed_cnv <- sort(bins_bed_cnv)
	out_cnv<-cbind(as.data.frame(bins_bed_cnv)[1:3],cnv=bins_bed_cnv$cnv)
	out_cnv$cnv<-as.integer(out_cnv$cnv)
	outname<-paste0("clone_",clone,".cnv.",res,".segmented.bedgraph")
	write.table(out_cnv,col.names=F,row.names=F,quote=F,file=outname,sep="\t")
	print(paste("Wrote out",outname))

}


#Make blacklist
lapply(bed_in_list, function(x) make_black_list(x,copykit_obj=tumor))

BINS_BED_PATH_500kb="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
clone_list<-colnames(tumor@consensus)

#Make clone specific bedgraph format for HiC windows
for(clone in clone_list){
	make_consensus_bedgraph(bed_in=BINS_BED_PATH_500kb,res="500kb",clone=clone,copykit_obj=tumor,cores=50)
}

#output values at matched bin resolutions for neoloop finder, then use segment-cnv and correct-cnv from neoloop finder afterwards
#also output a blacklist of regions that don't overlap between the bins.bed and the copykit bin filtered ranges
```
Make bigwig files for plotting

```R
library(copykit)
library(GenomicRanges)
library(parallel)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
wd_out="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive"
setwd(wd_out)
tumor<-readRDS("scCNA.rds")


#function to define blacklist regions
make_consensus_bigwig<-function(bed_in,res,copykit_obj,clone,cores=10){
	#bins in bin bed file that are not in the copykit filtered list are considered black list,
	#generate per bins.bed resolution granges
	bins_bed<-read.table(bed_in,header=F,sep="\t")
	colnames(bins_bed)<-c("chr","start","end")
	bins_bed<-GRanges(bins_bed)

	#generate granges of copykit data
	accepted_bed_ranges<-GRanges(copykit_obj@rowRanges) #row ranges for bedgraph
	accepted_bed_ranges$cnv<-copykit_obj@consensus[clone] #add cnv data as metadata column

	#overlap the two granges
	hits<-findOverlaps(query=bins_bed,subject=accepted_bed_ranges,minoverlap=1,ignore.strand=T) #get list of intersecting
	overlaps <- pintersect(bins_bed[queryHits(hits)], accepted_bed_ranges[subjectHits(hits)]) #combine intersecting sites for overlap
	bins_bed$cnv<-mean(unlist(accepted_bed_ranges$cnv),na.rm=TRUE) #set up column meta data, set to average for blacklist bed regions

	#change overlap calculation to correct for it one if bigger than the other (i think I did this right??)
	print(paste("Calculating window CNVs from CopyKit Clone consensus for",clone,"at",res))
	if(mean(width(accepted_bed_ranges))<mean(width(bins_bed))){ #if hic bins are bigger than copykit windows
		percentOverlap <- width(overlaps) / width(bins_bed[queryHits(hits)])
	} else { #if copykit windows are bigger than hic bins
		percentOverlap <- width(overlaps) / width(accepted_bed_ranges[subjectHits(hits)])}
	#lapply to 
	out<-mclapply(unique(hits@from),function(x) {
		overlaps_tmp<-hits[hits@from==x,]
		cnv_tmp<-unlist(as.data.frame(accepted_bed_ranges[overlaps_tmp@to,])[clone])
		weight_tmp<-percentOverlap[overlaps_tmp@to]
		cnv_out<-weighted.mean(cnv_tmp,w=weight_tmp,na.rm=TRUE) #weighted mean score by overlap percentage
		bins_bed[x,]$cnv<-cnv_out
		return(bins_bed[x,])},mc.cores=cores)
	bins_bed_cnv<-do.call("c",out)
	bins_bed_cnv<-c(bins_bed_cnv,bins_bed[-subjectHits(findOverlaps(query=bins_bed_cnv, subject=bins_bed, minoverlap=10)),] )
	bins_bed_cnv<-sortSeqlevels(bins_bed_cnv)
	bins_bed_cnv <- sort(bins_bed_cnv)

	out_cnv<-bins_bed_cnv
	colnames(out_cnv@elementMetadata)<-"score"
	start(out_cnv) <- start(out_cnv) + 1L
	out_name<-paste0("clone_",clone,".cnv.",res,".segmented.bigWig")
	seqinfo(out_cnv) <- seqinfo(txdb)[seqnames(seqinfo(out_cnv))]
	export(object=out_cnv, con=out_name,format="bigWig")
	print(paste("Wrote out",out_name))

}


BINS_BED_PATH_500kb="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
clone_list<-colnames(tumor@consensus)

#Make clone specific bedgraph format for HiC windows
for(clone in clone_list){
	make_consensus_bigwig(bed_in=BINS_BED_PATH_500kb,res="500kb",clone=clone,copykit_obj=tumor,cores=50)
}

```

### Detection of structural variants using Eagle C
ICE normalized output at multiresolution
```bash
conda activate EagleC
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
cd $clone_dir

eaglec_SV_detect_ICE() {
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
clone=${1::-17}
#run cooler balance on all cool matrices, using 20 cores
cooler balance -p 20 --force ${clone_dir}/${clone}.bsorted.5kb.cool
cooler balance -p 20 --force ${clone_dir}/${clone}.bsorted.10kb.cool
cooler balance -p 20 --force ${clone_dir}/${clone}.bsorted.50kb.cool

#now run predict sv

predictSV --hic-5k ${clone_dir}/${clone}.bsorted.5kb.cool \
            --hic-10k ${clone_dir}/${clone}.bsorted.10kb.cool \
            --hic-50k ${clone_dir}/${clone}.bsorted.50kb.cool \
            -O $clone -g hg38 --balance-type CNV --output-format full \
            --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999
}
export -f eaglec_SV_detect

clones=`ls clone*pairs.gz`

parallel --jobs 1 eaglec_SV_detect ::: $clones &
```

CNV normalized output at single 500kb resolution.

```bash
mamba activate EagleC
#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones" 
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
cd $clone_dir

eaglec_SV_CNV() {
	#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
	clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
	#bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive"
	bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
	in=$1
	clone=${in::-17}
	#balance through ICE first
	cooler balance -p 20 --force ${clone_dir}/${clone}.bsorted.500kb.cool
	
	#now balance by CNV
	conda activate neoloop
	correct-cnv -H ${clone_dir}/${clone}.bsorted.500kb.cool --cnv-file ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bedgraph --nproc 20 -f
	
	conda activate EagleC
	predictSV-single-resolution -H ${clone_dir}/${clone}.bsorted.500kb.cool -g hg38 --output-file ${clone_dir}/${clone}.SV.cnv.tsv --balance-type CNV --region-size 500000 --output-format full --prob-cutoff 0.5 --logFile ${clone}.eaglec.cnv.log 
	#--cache-folder ${clone_dir}/${clone}.cache
}
export -f eaglec_SV_CNV


#clones=`ls clone*pairs.gz`

#parallel --jobs 1 eaglec_SV_CNV ::: $clones & #parallel doesnt like switching between environments
eaglec_SV_CNV clone_c1.bsorted.pairs.gz
eaglec_SV_CNV clone_c2.bsorted.pairs.gz
eaglec_SV_CNV clone_c3.bsorted.pairs.gz
eaglec_SV_CNV clone_c4.bsorted.pairs.gz
eaglec_SV_CNV clone_c5.bsorted.pairs.gz
eaglec_SV_CNV clone_c6.bsorted.pairs.gz
eaglec_SV_CNV clone_c7.bsorted.pairs.gz

#Rerun with neoloopfinder output 
conda activate EagleC
#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones" 
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"

eaglec_SV_CNV_neoout() {
	#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
	clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
	#bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive"
	bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
	clone=${1::-17}
	#now balance by CNV
	conda activate EagleC
	predictSV-single-resolution -H ${clone_dir}/${clone}.bsorted.500kb.cool -g hg38 --output-file ${clone_dir}/${clone}.SV.cnv.neoloopfinder.tsv --balance-type CNV --region-size 500000 --output-format NeoLoopFinder --prob-cutoff 0.8 --logFile ${clone}.eaglec.cnv.log
	conda activate neoloop
	assemble-complexSVs -O ${clone_dir}/${clone} -B ${clone_dir}/${clone}.SV.cnv.neoloopfinder.tsv --balance-type CNV --protocol insitu --nproc 20 -H ${clone_dir}/${clone}.bsorted.500kb.cool --logFile ${clone}.neoloop.log --minimum-size 5000
}
export -f eaglec_SV_CNV_neoout

eaglec_SV_CNV_neoout clone_c1.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c2.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c3.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c4.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c5.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c6.bsorted.pairs.gz
eaglec_SV_CNV_neoout clone_c7.bsorted.pairs.gz

```

## Automate interSV plotting
Take in the SV output from EagleC, filter to significant interchr translocations, make unique and plot.

NOTE: I'll make this less specific in terms of directories, and maintain defaults in the future.
Bash script for plotting, located in /volumes/seq/projects/gccACT/src

/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots 
```bash

conda activate EagleC
#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones" 
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
cd $clone_dir

eaglec_SV_interchr_plots() {
	chrom_sizes="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/hg38.chrom.sizes"
	sv_file=$1
	clone_dir=$(dirname $sv_file)
	basename_sv=$(basename $sv_file)
	clone=${basename_sv::-11}
	cool_file=${clone_dir}"/"${clone}".bsorted.500kb.cool"
	if [ ! -d ${clone_dir}"/"${clone}"_interSVs" ] 
	then
    echo "Directory ${clone_dir}"/"${clone}"_interSVs" DOES NOT exist. Creating.."
    mkdir ${clone_dir}"/"${clone}"_interSVs"
	fi	

	cat $chrom_sizes | while read chrA end_size1 ; do cat $chrom_sizes | while read chrB end_size2 ;
			do outname="${clone_dir}/${clone}_interSVs/${clone}.interSV.${chrA}_${chrB}.png"
			plot-interSVs --cool-uri $cool_file --full-sv-file $sv_file -C $chrA $chrB --output-figure-name $outname --balance-type ICE --dpi 800
			echo "Completed file ${clone_dir}/${clone}_interSVs/${clone}.interSV.${chrA}_${chrB}.png"
		done ; done
}

export -f eaglec_SV_interchr_plots
eaglec_SV_interchr_plots $1

```

Running the bash script
```bash
#eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c1.SV.cnv.tsv

/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c2.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c3.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c4.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c5.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c6.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_interchr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c7.SV.cnv.tsv

```

## Automate intraSV plotting
Bash script for plotting, located in /volumes/seq/projects/gccACT/src
NOTE: I'll make this less specific in terms of directories, and maintain defaults in the future.

/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots 

```bash
#!/bin/bash
conda activate EagleC
#clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
#bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive"
cd $clone_dir

eaglec_SV_intrachr_plots() {
	sv_file=$1
	clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts"
	bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
	chrom_sizes="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/hg38.chrom.sizes"
	clone_dir=$(dirname $sv_file)
	basename_sv=$(basename $sv_file)
	clone=${basename_sv::-11}
	cool_file=${clone_dir}"/"${clone}".bsorted.500kb.cool"
	bedgraph_cnv=${bedgraph_dir}/${clone}.cnv.500kb.segmented.bedgraph

	if [ ! -f ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bigwig ] 
	then
	    echo "File ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bedgraph DOES NOT exist. Creating from bedgraph..."
	    sort -k1,1 -k2,2n $bedgraph_cnv > ${bedgraph_cnv}.sorted.bedgraph
	    bedGraphToBigWig ${bedgraph_cnv}.sorted.bedgraph $chrom_sizes ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bigwig
	fi	

	if [ ! -d ${clone_dir}"/"${clone}"_intraSVs" ] 
	then
	  echo "Directory ${clone_dir}"/"${clone}"_intraSVs" DOES NOT exist. Creating.."
	  mkdir ${clone_dir}"/"${clone}"_intraSVs"
	fi	

	cat $chrom_sizes | while read chr_in_size end_size ; 
	do plot-intraSVs --cool-uri $cool_file --full-sv-file $sv_file --region ${chr_in_size}:1-${end_size} --output-figure-name ${clone_dir}/${clone}_intraSVs/${clone}.intraSV.${chr_in_size}.png --cnv-file ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bigwig --coordinates-to-display 1 ${end_size} --cnv-max-value 6 --balance-type ICE --dpi 800 ; 
	echo "Completed file ${clone_dir}/${clone}_intraSVs/${clone}.intraSV.${chr_in_size}.png"
	done
}

export -f eaglec_SV_intrachr_plots

eaglec_SV_intrachr_plots $1
```

```bash
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c2.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c3.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c4.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c5.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c6.SV.cnv.tsv
/volumes/seq/projects/gccACT/src/eaglec_SV_intrachr_plots /volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clone_c7.SV.cnv.tsv

```

### Plot all by all autosome interactions
Use this for all chromosomes plots, code adapted from cooltools and inspired by:
https://github.com/bianlab-hub/zuo_ncomms_2021/blob/Hi-C_data_analysis/fig1f_plot_obs_heatmap.py


```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import os
import seaborn as sns
import cooltools.lib.plotting
import matplotlib

def trans_chr_plot(mat,chr_row,chr_col,ax_row,ax_col,axes,out,zmin,zmax,cmap,xlim,ylim):
	""" Function to plot chr by chr trans interactions (or chr by chr cis interactions if same chr given)
	Takes in c for cooler file and str for chr_row and chr_col (format "chr1","chrX", etc.)
	Takes in integers for ax_row and ax_col to add output plot as subplot.
	"""
	mat=np.log10(mat)
	ax=sns.heatmap(ax=axes[ax_row,ax_col],data=mat,cbar=False,yticklabels=False,square=False,xticklabels=False,vmin=zmin,vmax=zmax,cmap=cmap)#,
		#gridspec_kw={"truncate":False})
	ax.set_ylim(ylim,0)
	ax.set_xlim(0,xlim)
	if ax_row==0:
		axes[ax_row,ax_col].set_title(chr_col)
	if ax_col==0:
		axes[ax_row,ax_col].set_ylabel(chr_row)
	print("Completed "+chr_row+" by "+chr_col+" plot:"+out)
	return(ax)


def all_by_all_plot(infile_name,chr_count,zmin,zmax,cmap):
	"""Function to read in cooler file from given string name. And run all by all chr comparison"""
	#Read in Cooler File
	in_file=infile_name
	in_name=in_file.split(sep=".")[0]
	coolfile=in_file
	c = cooler.Cooler(coolfile)
	cooler.coarsen_cooler(coolfile,in_name+"_5mb.cool",factor=2,chunksize=10000000) #coarsen to 5mb
	c = cooler.Cooler(in_name+"_5mb.cool")
	cooler.balance_cooler(c,store=True,rescale_marginals=True)#rescale_margines=True #balance matrix ignore_diags=10,
	obs_mat = c.matrix()[:]
	chr_list=list(c.bins()[:]["chrom"].unique())
	out=''.join([in_name,'_all_by_all_log2_1Mb_obs.png'])
	#init subplot
	chr_in=len(chr_list)
	chr_sizes=pd.DataFrame(c.bins()[:]).groupby(["chrom"])["chrom"].count()
	chr_ratios=list(chr_sizes/chr_sizes[0])
	fig, axes = plt.subplots(chr_count, chr_count, figsize=(40, 40),sharex=True,sharey=True,
			gridspec_kw={'width_ratios': chr_ratios[0:chr_count],'height_ratios':chr_ratios[0:chr_count]})
		#plt.subplots_adjust(hspace=0.1,wspace=0.1)
	for i in range(0,chr_count):
		row_chrom=chr_list[i]
		ax_row=i
		xlim=chr_sizes[i]
		for j in range(0,chr_count):
			col_chrom=chr_list[j]
			ax_col=j
			ylim=chr_sizes[j]
			mat=c.matrix().fetch(row_chrom,col_chrom)
			im=trans_chr_plot(mat,row_chrom,col_chrom,ax_row,ax_col,axes,in_name,zmin,zmax,cmap,xlim,ylim)
	plt.tight_layout()
	plt.savefig(out, dpi=dpi, format='png',bbox_inches="tight")
	plt.close("all")


def all_by_all_plot_subtract(infile_name1,infile_name2,chr_count,cmap):
	"""Function to read in cooler file from given string name. And run all by all chr comparison"""
	#Read in Cooler File
	in_file1=infile_name1
	in_file2=infile_name2
	in_name1=in_file1.split(sep=".")[0]
	in_name2=in_file2.split(sep=".")[0]
	coolfile1=in_file1
	coolfile2=in_file2
	c1 = cooler.Cooler(coolfile1)
	c2 = cooler.Cooler(coolfile2)
	cooler.balance_cooler(c1,store=True) #balance matrix
	cooler.balance_cooler(c2,store=True) 
	chr_list=list(c1.bins()[:]["chrom"].unique())
	out=''.join([in_name1,"_",in_name2,'_all_by_all_log2_1Mb_obs.pdf'])
	#init subplot
	chr_in=len(chr_list)
	fig, axes = plt.subplots(chr_count, chr_count, figsize=(10, 10),sharex=True,sharey=True)
	plt.subplots_adjust(hspace=0.1,wspace=0.1)
	subplot_chr_list=[]
	for i in range(0,chr_count):
		row_chrom=chr_list[i]
		ax_row=i
		for j in range(0,chr_count):
			col_chrom=chr_list[j]
			ax_col=j
			obs_mat1 = c1.matrix().fetch(row_chrom,col_chrom)
			obs_mat2 = c2.matrix().fetch(row_chrom,col_chrom)
			mat=np.log2(obs_mat1) / np.log2(obs_mat2)
			trans_chr_plot(mat,row_chrom,col_chrom,ax_row,ax_col,axes,out,cmap) 
	print(plt)
	plt.savefig(out, dpi=dpi, format='pdf')
	plt.close("all")

#wd 
os.chdir('/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts')

#in files
in_files=["clone_c1.bsorted.cool",
	"clone_c2.bsorted.cool",
	"clone_c3.bsorted.cool",
	"clone_c4.bsorted.cool",
	"clone_c5.bsorted.cool",
	"clone_c6.bsorted.cool"]

#Set up settings for plot
dpi= 300
colormap='Reds'
zmin=-3
zmax=-1

#make all by all plot for all clones
[all_by_all_plot(x,23,zmin,zmax,cmap=colormap) for x in in_files]

#colormap="coolwarm"
#for x in in_files:
#	for y in in_files:
#		all_by_all_plot_subtract(x,y,4)

#

```

### Use NeoLoopFinder to annotate complex SVs
```bash
conda activate neoloop

assemble-complexSVs -O ${clone_dir}/${clone} -B ${clone_dir}/${clone}.SV.neoloopfinder.tsv --balance-type CNV --protocol insitu --nproc 20 -H ${clone_dir}/${clone}.bsorted.500kb.cool 



eaglec_SV_CNV() {
	clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones"
	bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive"
	clone=${1::-17}
	#balance through ICE first
	
	#now balance by CNV
	conda activate neoloop

assemble-complexSVs -O ${clone_dir}/${clone} -B ${clone_dir}/${clone}.SV.neoloopfinder.tsv --balance-type CNV --protocol insitu --nproc 20 -H ${clone_dir}/${clone}.bsorted.500kb.cool 

}
export -f eaglec_SV_CNV



```
### Generate eigengenes for compartments across chromosomes
```bash
#note this requires the cooler matrix be balanced and normalized already

cooler_eigen() {
cooltools eigs-cis -o outputs/test.eigs.100000 --view data/view_hg38.tsv --phasing-track outputs/gc.100000.tsv --n-eigs 1 $cool_file::resolutions/100000
}
export -f cooler_balance


```



### Comparison with ACT-seq MDA-MB-231
```bash
#files located here:
/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3

#list of fastq files here (only some have R1 R2):
/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3/fastq_P1_P2_P3

project_dir="/volumes/seq/projects/gccACT/mdamb231_ACTseq"
mkdir $project_dir
mkdir ${project_dir}/cells


#Map reads with BWA Mem (filter to PE reads before alignment)
bwa_align() {
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
project_dir="/volumes/seq/projects/gccACT/mdamb231_ACTseq"
fq2=$1
if [[ $fq2 == *"_R2_"* ]]; 
	then 
		fq1=$(echo $fq2 | sed 's/_R2_/_R1_/')
		base_name=$(basename $fq1)
		out_name=${base_name::-9}
		echo "Checking for files for ${fq1} and ${fq2}..."
		if [ -f $fq1 ] && [ -f $fq2 ]
		then
			echo "Running alignment for ${out_name}..."
			bwa mem $ref $fq1 $fq2 | samtools view -b - > ${project_dir}/cells/${out_name}.bam
		fi	
fi
}
export -f bwa_align

fq=$(cat /volumes/seq/projects/CNA_projects/DT_CNA/cell_line/231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3/fastq_P1_P2_P3)
parallel --jobs 30 bwa_align ::: $fq &

## Mark duplicate reads
#name sort, fix mates, sort by position, mark dup
sort_and_markdup() {
  samtools sort -T . -n -o - $1 | samtools fixmate -m - -| samtools sort -T . -o - - | samtools markdup -s - ${1::-4}.rmdup.bam 2> ${1::-4}.rmdup.stats.txt
}
export -f sort_and_markdup

bam_in=`ls *bam`
parallel --jobs 30 sort_and_markdup ::: $bam_in

#Ready to run copykit!

```

## Run CopyKit for WGS ACT-seq portion
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells")

act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
act <- findAneuploidCells(act)

# Mark low-quality cells for filtering
act <- findOutliers(act)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(act, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# Remove cells marked as low-quality and/or aneuploid from the copykit object
act <- act[,SummarizedExperiment::colData(act)$outlier == FALSE]
act <- act[,SummarizedExperiment::colData(act)$is_aneuploid == TRUE]


# kNN smooth profiles
act <- knnSmooth(act)


k_clones<-findSuggestedK(act)
# Create a umap embedding 
tumor <- runUmap(tumor)

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
tumor  <- findClusters(tumor,k_subclones=17)#output from k_clones

pdf("subclone.umap.pdf")
plotUmap(tumor, label = 'subclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
tumor <- calcConsensus(tumor)
tumor <- runConsensusPhylo(tumor)

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.pdf")
plotHeatmap(tumor, label = 'subclones',order='hclust')
dev.off()

saveRDS(tumor,file="/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
tumor<-readRDS("/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
clone_out<-data.frame(bam=paste0(row.names(tumor@colData),".bam"),clone=tumor@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("clone_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}
```

### Count of WGS and GCC Reads
```bash
dir="/volumes/seq/projects/gccACT/mdamb231_ACTseq"
#count reads based on paired alignments
readtype_count() {
	#print out reads on same chromosome that are >= 1000 bp apart (distal)
	distal=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))>=1000) print $0}' | wc -l)
	#print out reads on different chromosomes
	trans=$(samtools view -F 1024 $1 | awk '{if($7 != "=") print $0}' | wc -l)
	#print out reads on same chromosome within 1000bp (cis)
	near=$(samtools view -F 1024 $1 | awk '{if (sqrt(($9^2))<=1000) print $0}' | wc -l)
	echo $1,$near,$distal,$trans
}
export -f readtype_count

cd $dir/cells
bam_in=`ls *bam`
echo "cellid,near_cis,distal_cis,trans" > read_count.csv; parallel --jobs 10 readtype_count ::: $bam_in >> read_count.csv
```

Plotting GCC read types
```R
library(ggplot2)
library(patchwork)
setwd("/volumes/seq/projects/gccACT/mdamb231_ACTseq")
dat<-read.table("./cells/read_count.csv",header=T,sep=",")
dat$total_reads<-dat$near_cis+dat$distal_cis+dat$trans

plt1<-ggplot(dat,aes(y=near_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis Reads")+theme_minimal()
plt2<-ggplot(dat,aes(y=distal_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis Reads")+theme_minimal()
plt3<-ggplot(dat,aes(y=trans,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans Reads")+theme_minimal()

plt4<-ggplot(dat,aes(y=(near_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis % Reads")+theme_minimal()+ylim(c(0,100))
plt5<-ggplot(dat,aes(y=(distal_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis % Reads")+theme_minimal()+ylim(c(0,10))
plt6<-ggplot(dat,aes(y=(trans/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans % Reads")+theme_minimal()+ylim(c(0,10))

plt<-(plt1|plt2|plt3)/(plt4|plt5|plt6)
ggsave(plt,file="read_counts.pdf")
```

```
#correlate HiC to WGS output
#Test for number of split reads compared to HiC data
#Test for SVs?
--->

<!--
#Add compartments
#Add SV detection by contact matrix
#Use cooltools virtual4c for eccDNA interactions
#Add 
-->
<!--
```



### eccDNA Analysis
https://www.science.org/doi/10.1126/sciadv.aba2489
https://github.com/pk7zuva/Circle_finder/blob/master/circle_finder-pipeline-bwa-mem-samblaster.sh
It says read length must be at least 75bp if not enriched. But I'm going to try anyway.

```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
bwa mem <idxbase> samp.r1.fq samp.r2.fq | samblaster -e -d samp.disc.sam -s samp.split.sam | samtools view -Sb - > samp.out.bam
samtools -H $dir/230306_gccact.bam | samblaster -e --minNonOverlap 10 -d $6-$7\.disc.sam -s $6-$7\.split.sam -u $6-$7\.unmap.sam > $6-$7\.sam
```

```bash
#Use this script if your read length is >= 75
#Arg1 = Number of processors
#Arg2 = Genome or index file "/hdata1/MICRODNA-HG38/hg38.fa"
#Arg3 = fastq file 1 "1E_S1_L1-L4_R1_001.fastq"
#Arg4 = fastq file 2 "1E_S1_L1-L4_R2_001.fastq"
#Arg5 = minNonOverlap between two split reads "10"
#Arg6 = Sample name "1E"
#Arg7 = genome build "hg38"

#Usage: bash "Number of processors" "/path-of-whole-genome-file/hg38.fa" "fastq file 1" "fastq file 2" "minNonOverlap between two split reads" "Sample name" "genome build"
#bash circle_finder-pipeline-bwa-mem-samblaster.sh 16 hg38.fa 1E_S1_L1-L4_R1_001.fastq.75bp-R1.fastq 1E_S1_L1-L4_R2_001.fastq.75bp-R2.fastq 10 1E hg38

#Step 1: Mapping.
bwa mem -t $1 $2 $3 $4 | samblaster -e --minNonOverlap $5 -d $6-$7\.disc.sam -s $6-$7\.split.sam -u $6-$7\.unmap.sam > $6-$7\.sam

#Step 2: Converting (sam2bam), sorting and indexing mapped reads. Output of this step is input in step 3

samtools view -@ $1 -bS $6-$7\.sam -o $6-$7\.bam
samtools sort -@ $1 -O bam -o $6-$7\.sorted.bam $6-$7\.bam
samtools index $6-$7\.sorted.bam

samtools view -@ $1 -bS $6-$7\.disc.sam > $6-$7\.disc.bam
samtools view -@ $1 -bS $6-$7\.split.sam > $6-$7\.split.bam
samtools view -@ $1 -bS $6-$7\.unmap.sam > $6-$7\.unmap.bam

#Step 3: Extract concordant pairs with headers 

samtools view  -@ $1 -hf 0x2 $6-$7\.sorted.bam -bS > $6-$7\.concordant.bam

#Step 4: Converting bam to bed format (Remember bedtools generate 0 based co-ordinates)

bedtools bamtobed -cigar -i $6-$7\.split.bam | sed -e s/_2\\/2/\ 2/g | sed -e s/_1\\/1/\ 1/g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' | awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", $8)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", $8)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", $8)} 1' | awk 'BEGIN{FS=OFS=" "} {if (($9=="M" && $NF=="H") || ($9=="M" && $NF=="S"))  {printf ("%s\tfirst\n",$0)} else if (($9=="S" && $NF=="M") || ($9=="H" && $NF=="M")) {printf ("%s\tsecond\n",$0)} }' | awk 'BEGIN{FS=OFS="\t"} {gsub("\ ", "", $8)} 1' > $6-$7\.split.txt

#bedtools bamtobed -cigar -i $6-$7\.concordant.bam | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > $6-$7\.concordant.txt
bedtools bamtobed -cigar -i $6-$7\.sorted.bam | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > $6-$7\.concordant.txt

bedtools bamtobed -cigar -i $6-$7\.disc.bam | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > $6-$7\.disc.txt

#step 5: Calculating the Read-ID frequency. The frequency of 2 would indicate that particular Read is uniquely mapping in genome and there are only two split-mapped reads

awk '{print $4}' $6-$7\.split.txt | sort | uniq -c > $6-$7\.split.id-freq.txt
#This file "$6-$7\.split.id-freq.txt" will be used for collecting split id that have frequency equal to 4.
awk '$1=="2" {print $2}' $6-$7\.split.id-freq.txt > $6-$7\.split.id-freq2.txt
awk '$1=="4" {print $2}' $6-$7\.split.id-freq.txt > $6-$7\.split.id-freq4.txt

awk '{print $4}' $6-$7\.concordant.txt | sort | uniq -c > $6-$7\.concordant.id-freq.txt
#The following command will chose (may not be always true) one concordant and 2 split read
awk '$1=="3" {print $2}' $6-$7\.concordant.id-freq.txt > $6-$7\.concordant.id-freq3.txt
awk '$1>3 {print $2}' $6-$7\.concordant.id-freq.txt > $6-$7\.concordant.id-freqGr3.txt

#Step 6: Selecting split reads that were 1) mapped uniquely and 2) mapped on more than one loci. For normal microDNA identification no need to use the "freqGr2" file
grep -w -Ff $6-$7\.split.id-freq2.txt $6-$7\.split.txt > $6-$7\.split_freq2.txt
grep -w -Ff $6-$7\.split.id-freq4.txt $6-$7\.split.txt > $6-$7\.split_freq4.txt

#Selecting concordant pairs that were 1) mapped uniquely and 2) mapped on more than one loci (file "freqGr3.txt")
grep -w -Ff $6-$7\.concordant.id-freq3.txt $6-$7\.concordant.txt > $6-$7\.concordant_freq3.txt
grep -w -Ff $6-$7\.concordant.id-freqGr3.txt $6-$7\.concordant.txt > $6-$7\.concordant_freqGr3.txt

#Step 7: Putting split read with same id in one line
sed 'N;s/\n/\t/' $6-$7\.split_freq2.txt > $6-$7\.split_freq2.oneline.txt
sed 'N;s/\n/\t/' $6-$7\.split_freq4.txt > $6-$7\.split_freq4.oneline.txt

#Step 8: Split reads map on same chromosome and map on same strand. Finally extracting id (split read same chromosome, split read same strand), collecting all the split reads that had quality >0
awk '$1==$10 && $7==$16 && $6>0 && $15>0 {print $4} ' $6-$7\.split_freq2.oneline.txt > $6-$7\.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

#Step 9: Based on unique id I am extracting one continuously mapped reads and their partner mapped as split read (3 lines for each id) 
grep -w -Ff $6-$7\.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt $6-$7\.concordant_freq3.txt > $6-$7\.concordant_freq3.2SPLIT-1M.txt

#Step 10: Sorting based on read-id followed by length of mapped reads.
awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", $8)} 1' $6-$7\.concordant_freq3.2SPLIT-1M.txt | awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", $8)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", $8)} 1' | awk 'BEGIN{FS=OFS=" "} {if (($9=="M" && $NF=="H") || ($9=="M" && $NF=="S"))  {printf ("%s\tfirst\n",$0)} else if (($9=="S" && $NF=="M") || ($9=="H" && $NF=="M")) {printf ("%s\tsecond\n",$0)} else  {printf ("%s\tconfusing\n",$0)}}' | awk 'BEGIN{FS=OFS="\t"} {gsub("\ ", "", $8)} 1' | awk '{printf ("%s\t%d\n",$0,($3-$2)+1)}' | sort -k4,4 -k10,10n | sed 'N;N;s/\n/\t/g' | awk '{if ($5==$15) {print $0}  else if (($5=="1" && $15=="2" && $25=="1") || ($5=="2" && $15=="1" && $25=="2")) {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20)} else if (($5=="1" && $15=="2" && $25=="2") || ($5=="2" && $15=="1" && $25=="1")) {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n", $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)} }' > $6-$7\.concordant_freq3.2SPLIT-1M.inoneline.txt

#Step 11: Unique number of microDNA with number of split reads
awk '$1==$11 && $1==$21 && $7==$17 && length($8)<=12 && length($18)<=12 && length($28)<=12'  $6-$7\.concordant_freq3.2SPLIT-1M.inoneline.txt | awk '($7=="+" && $27=="-") || ($7=="-" && $27=="+")' | awk '{if ($17=="+" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\t%d\t%d\n",$1,$12,$3)} else if ($7=="+" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\t%d\t%d\n",$1,$2,$13)} else if ($17=="-" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\t%d\t%d\n",$1,$12,$3)} else if ($7=="-" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\t%d\t%d\n",$1,$2,$13)} }' | sort | uniq -c | awk '{printf ("%s\t%d\t%d\t%d\n",$2,$3,$4,$1)}' > $6-$7\.microDNA-JT.txt

rm *hg38.sam *hg38.bam
```

### HiC Data Analysis can be done with the DipC group's released hickit

https://github.com/lh3/hickit

```bash


# Map Dip-C reads and extract contacts (skip if you use your own pipeline)
seqtk mergepe read1.fq.gz read2.fq.gz | ~/tools/hickit-0.1_x64-linux/pre-dip-c - | bwa mem -5SP -p hs37d5.fa - | gzip > aln.sam.gz
~tools/hickit-0.1_x64-linux/k8 hickit.js vcf2tsv phased.vcf > phased_SNP.tsv   # extract phased SNPs from VCF
~tools/hickit-0.1_x64-linux/k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ~tools/hickit-0.1_x64-linux/k8 hickit.js chronly - | ~tools/hickit-0.1_x64-linux/k8 hickit.js bedflt par.bed - | gzip > contacts.seg.gz # for male
#./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ./k8 hickit.js chronly -y - | gzip > contacts.seg.gz # for female
~tools/hickit-0.1_x64-linux/hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz  # optional

# Impute phases (-i also works with contacts.seg.gz)
~tools/hickit-0.1_x64-linux/hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz
~tools/hickit-0.1_x64-linux/hickit -i contacts.pairs.gz --out-val=impute.val     # estimate imputation accuracy by holdout
# Infer 3D structure
~tools/hickit-0.1_x64-linux/hickit -i impute.pairs.gz -Sr1m -c1 -r10m -c5 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O imput.3dg

# 2D contact map in PNG (bin size determined by the image width)
~tools/hickit-0.1_x64-linux/hickit -i impute.pairs.gz --out-png impute.png
# Compute CpG density (optional)
~tools/hickit-0.1_x64-linux/hickit.js gfeat -r hs37d5.fa.gz imput.3dg | gzip > imput.cpg.3dg.gz
# Visualize 3D structure (requiring a graphical card)
~tools/hickit-0.1_x64-linux/hickit-gl -I imput.cpg.3dg.gz --view
```


-->