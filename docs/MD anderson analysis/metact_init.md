---
title: metACT-seq Ligation Benchmarking
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /metact_init/
category: mdanderson
---

## Initial analysis of Nextseq 2000 Run

Sequencing files located here:
```bash
/volumes/seq/flowcells/MDA/nextseq2000/2023/20230913_ACT_Meth_triple-seq_4cells/230913_VH00219_480_AACN2CGM5
```
Reads are sequenced in this format: 55-11-16-55

Set up working directory
```bash
mkdir /volumes/seq/projects/metACT/230913_metACT_benchmark
cd /volumes/seq/projects/metACT/230913_metACT_benchmark
ln -s /volumes/seq/flowcells/MDA/nextseq2000/2023/20230913_ACT_Meth_triple-seq_4cells/230913_VH00219_480_AACN2CGM5 230913_VH00219_480_AACN2CGM5
```


```bash
cd /volumes/seq/projects/metACT/230913_metACT_benchmark/230913_VH00219_480_AACN2CGM5

bcl2fastq -R /volumes/seq/projects/metACT/230913_metACT_benchmark/230913_VH00219_480_AACN2CGM5 \
-o /volumes/seq/projects/metACT/230913_metACT_benchmark/ \
-r 10 \
-p 10 \
-w 10 \
--ignore-missing-bcls --ignore-missing-filter \
--ignore-missing-positions --ignore-missing-controls \
--create-fastq-for-index-reads
```

Split out metACT and ACT reads via python script

plate_demultiplex.py

```python
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd

#set up plate
idx_list="/volumes/seq/projects/metACT/230913_metACT_benchmark/230830_metACT_sequencing_indexes.tsv"

idx_in=pd.read_csv(idx_list,sep="\t")

plate_i7=idx_in["i7_idx_seq"].unique()
plate_i5=idx_in["i5_idx_seq"].unique()

plate_i7=[i.strip() for i in plate_i7]
plate_i5=[i.strip() for i in plate_i5]

fq1=sys.argv[1] 
fq2=sys.argv[2] 
idx3=sys.argv[3] 
idx4=sys.argv[4] 
#examples
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
python ./plate_demultiplex.py \
Undetermined_S0_L001_R1_001.fastq.gz \
Undetermined_S0_L001_R2_001.fastq.gz \
Undetermined_S0_L001_I1_001.fastq.gz \
Undetermined_S0_L001_I2_001.fastq.gz

#then gzip
for i in *fastq; do gzip $i & done &
```

# Processing Samples via Copykit to start
## Alignment
```bash
#set up variables and directory
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/seq/projects/metACT/230913_metACT_benchmark"
fq1="/volumes/seq/projects/metACT/230913_metACT_benchmark/Undetermined_S0_L001_R1_001.barc.fastq.gz"
fq2="/volumes/seq/projects/metACT/230913_metACT_benchmark/Undetermined_S0_L001_R2_001.barc.fastq.gz"

#Map reads with BWA Mem
bwa mem -t 50 $ref $fq1 $fq2 | samtools view -b - > $dir/230921_metact_benchmark.bam
```


## Split out single-cells
```bash
#set up cell directory
mkdir $dir/cells
```
Split by readname bam field into cells subdir and add metadata column

```bash
samtools view $dir/230921_metact_benchmark.bam| awk -v dir=$dir 'OFS="\t" {split($1,a,":"); print $0,"XM:Z:"a[1] > "./cells/"a[1]".230921_metact_benchmark.sam"}' &
 #split out bam to cell level sam
```

## Add header to each sam and convert to bam and sort
Using parallel to save time
```bash

add_header() {
	ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
	dir="/volumes/seq/projects/metACT/230913_metACT_benchmark"
	samtools view -bT $ref $1 | samtools sort -T $dir -o ${1::-4}.bam -
}
export -f add_header
cd $dir/cells
sam_in=`ls ./*sam`
parallel --jobs 50 add_header ::: $sam_in

#remove sam files to clear up space
rm -f $dir/cells/*sam
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

## Downsample ACT-seq cells to match metACT
```bash
dir="/volumes/seq/projects/metACT/230913_metACT_benchmark"

subsample() {
	dir="/volumes/seq/projects/metACT/230913_metACT_benchmark"
	in=`echo $1 | tr -d '"'` #remove quotes
	in_bam=`ls $dir/cells/${in}*rmdup.bam`
	samtools view -bs 42.2 $in_bam > ${in_bam::-4}.subsamp2.bam
	echo "Subsampling ${in_bam}..."
}
export -f subsample

bam_in=`awk 'OFS="\t" {if($7="act") print $11}' $dir/230830_metACT_sequencing_indexes.tsv`

parallel --jobs 50 subsample ::: $bam_in

mkdir $dir/cells/subsamp
mv $dir/cells/*subsamp2.bam $dir/cells/subsamp

find $dir/cells/subsamp -size -10M -delete
```

## Run CopyKit for WGS portion (Subsampled)
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/metACT/230913_metACT_benchmark/")

tumor2 <- runVarbin("/volumes/seq/projects/metACT/230913_metACT_benchmark/cells/subsamp/",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE,min_bincount=0)

# kNN smooth profiles
tumor <- knnSmooth(tumor2)

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.subsamp.pdf")
plotHeatmap(tumor, order='hclust')
dev.off()

```
## Project Library Complexity
Using Picard Tools
```bash
dir="/volumes/seq/projects/metACT/230913_metACT_benchmark"

project_count() {
java -jar /volumes/seq/code/3rd_party/picard/picard-2.20.4/picard.jar EstimateLibraryComplexity I=$1 O=${1::-4}.complex_metrics.txt
}
export -f project_count

cd $dir/cells
bam_in=`ls *bam`
parallel --jobs 50 project_count ::: $bam_in
```


## Run Fastqc on everything and clean up

```bash
bam_in=`ls *rmdup.bam`
parallel --jobs 30 fastqc ::: $bam_in

mkdir $dir/cells/fastqc
mv *fastqc* $dir/cells/fastqc
mv *stats.txt $dir/cells/fastqc
mv *txt $dir/cells/fastqc
rm -rf *benchmark.bam #only keep duplicate marked bams

#run multiqc to aggregate
multiqc . 
#C100 had zero reads, gzipping to prevent copykit from reading in

```

```R
library(ggplot2)
library(dplyr)
metadat<-read.table("/volumes/seq/projects/metACT/230913_metACT_benchmark/230830_metACT_sequencing_indexes.tsv",sep="\t",header=T)
metadat$cellID<-paste0(metadat$i7_idx_seq,metadat$i5_idx_seq)
projected_complexity<-list.files("/volumes/seq/projects/metACT/230913_metACT_benchmark/cells/fastqc",pattern=".complex_metrics.txt$",full.names=T)

metadat$proj_compl<-NA

for(x in projected_complexity){
	print(x)
	tmp<-read.table(x,nrows=1,sep="\t",header=T)
	tmp$cellID<-lapply(strsplit(basename(x),"[.]"),"[[",1)
	if(!is.na(tmp$ESTIMATED_LIBRARY_SIZE) && tmp$cellID %in% metadat$cellID){
	metadat[metadat$cellID==tmp$cellID,]$proj_compl<-tmp$ESTIMATED_LIBRARY_SIZE
	}
}

write.table(metadat,"/volumes/seq/projects/metACT/230913_metACT_benchmark/230830_metACT_sequencing_indexes.tsv",col.names=T,row.names=F,sep="\t")
plt<-ggplot(metadat,aes(x=assay,y=log10(proj_compl)))+geom_violin()+geom_jitter()
ggsave(plt,file="ligation_benchmark_projected_complexity.pdf")

metadat %>% group_by(assay,cell_count,ligation_concentration) %>% summarize(median(proj_compl,na.rm=T))

#   assay  cell_count ligation_concentration `median(proj_compl, na.rm = T)`
#   <chr>       <int> <chr>                                            <dbl>
# 1 act             0 na                                                 NA 
# 2 act             1 na                                            1110743 
# 3 act            10 na                                            7495187 
# 4 metact          0 0.1                                                NA 
# 5 metact          1 0.05                                           210127 
# 6 metact          1 0.1                                            195268.
# 7 metact          1 0.2                                            203138 
# 8 metact          1 0.3                                            180872.
# 9 metact          1 0.4                                            268436 
#10 metact          1 0.5                                            164760.
#11 metact          1 0.6                                            219330.
#12 metact          1 0.7                                            223688.
#13 metact          1 0.8                                            253306.
#14 metact          1 0.9                                            240320.
#15 metact         10 1                                              679716.

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
setwd("/volumes/seq/projects/metACT/230913_metACT_benchmark/")

tumor2 <- runVarbin("/volumes/seq/projects/metACT/230913_metACT_benchmark/cells/",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
tumor2 <- findAneuploidCells(tumor2)

# Mark low-quality cells for filtering
tumor <- findOutliers(tumor2)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(tumor, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

####SKIPPED###
# Remove cells marked as low-quality and/or aneuploid from the copykit object
#tumor <- tumor[,SummarizedExperiment::colData(tumor)$outlier == FALSE]
#tumor <- tumor[,SummarizedExperiment::colData(tumor)$is_aneuploid == TRUE]


# kNN smooth profiles
tumor <- knnSmooth(tumor)


# Create a umap embedding 
tumor <- runUmap(tumor)
tumor<-findSuggestedK(tumor)

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

#read in metadata
metadat<-read.table("/volumes/seq/projects/metACT/230913_metACT_benchmark/230830_metACT_sequencing_indexes.tsv",sep="\t",header=T)
metadat$cellID<-paste0(metadat$i7_idx_seq,metadat$i5_idx_seq)
tumor@colData$cellID<-unlist(lapply(strsplit(row.names(tumor@colData),"[.]"),"[[",1))
#merge to rds
tumor@colData<-merge(tumor@colData,metadat,by="cellID")
row.names(tumor@colData)<-tumor@colData$sample
tumor@colData$proj_compl_log10<-log10(tumor@colData$proj_compl)
tumor@colData$reads_assigned_bins_log10<-log10(tumor@colData$reads_assigned_bins)
saveRDS(tumor,file="/volumes/seq/projects/metACT/230913_metACT_benchmark/scCNA.rds")

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.pdf")
plotHeatmap(tumor, label = c('assay',"cell_count","ligation_concentration","proj_compl_log10","reads_assigned_bins_log10"),order='hclust')
dev.off()

```



