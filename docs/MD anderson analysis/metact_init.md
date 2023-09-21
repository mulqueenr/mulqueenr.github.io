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
parallel --jobs 30 add_header ::: $sam_in

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
rm -rf *benchmark.bam #only keep duplicate marked bams

#run multiqc to aggregate
multiqc . 
#C100 had zero reads, gzipping to prevent copykit from reading in

```

```R
library(ggplot2)

metadat<-read.table("/volumes/seq/projects/metACT/230913_metACT_benchmark/230830_metACT_sequencing_indexes.csv",sep=",",header=T)
metadat$cellID<-paste0(metadat$i7_idx_seq,metadat$i5_idx_seq)
projected_complexity<-list.files("/volumes/seq/projects/metACT/230913_metACT_benchmark/cells",pattern=".complex_metrics.txt$")

metadat$proj_compl<-NA

for(x in projected_complexity){
	print(x)
	tmp<-read.table(x,nrows=1,sep="\t",header=T)
	tmp$cellID<-lapply(strsplit(basename(x),"[.]"),"[[",1)
	if(!is.na(tmp$ESTIMATED_LIBRARY_SIZE)){
	metadat[metadat$cellID==tmp$cellID,]$proj_compl<-tmp$ESTIMATED_LIBRARY_SIZE
	}
}

```
<!--

## Run CopyKit for WGS portion
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/230306_mdamb231_test/cells")

tumor2 <- runVarbin("/volumes/seq/projects/gccACT/230306_mdamb231_test/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
tumor2 <- findAneuploidCells(tumor2)

# Mark low-quality cells for filtering
tumor <- findOutliers(tumor)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(tumor, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# Remove cells marked as low-quality and/or aneuploid from the copykit object
tumor <- tumor[,SummarizedExperiment::colData(tumor)$outlier == FALSE]
tumor <- tumor[,SummarizedExperiment::colData(tumor)$is_aneuploid == TRUE]


# kNN smooth profiles
tumor <- knnSmooth(tumor)


k_clones<-findSuggestedK(tumor)
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


-->