---
title: gccACT-seq Initial Test
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /gccACT_init/
category: mdanderson
---

## Initial analysis of Nextseq 2000 Run
Ran P2 100 cycle kit with 50% of reads alloted for gccACT

Stored here:
```bash
/volumes/seq/flowcells/MDA/nextseq2000/2023/230306_VH00219_371_AACJJFWM5
```

```bash
#run transfered to /volumes/seq/tmp
cd /volumes/seq/flowcells/MDA/nextseq2000/2023/230306_VH00219_371_AACJJFWM5
bcl2fastq -R /volumes/seq/flowcells/MDA/nextseq2000/2023/230306_VH00219_371_AACJJFWM5 \
-o ~/fastq/230306_VH00219_371_AACJJFWM5/ \
-r 4 \
-p 10 \
-w 4 \
--ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls \
--create-fastq-for-index-reads
```

Split out gccACT reads via python script
403,575,855 total reads (expecting 400M)

```python
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#set up plate
plate2_i7=['AGGCAGAA', 'TCCTGAGC', 'GGACTCCT', 'TAGGCATG', 'CTCTCTAC', 'CAGAGAGG', 'GCTACGCT', 'CGAGGCTG', 'AAGAGGCA', 'GTAGAGGA', 'GCTCATGA', 'ATCTCAGG', 'ACTCGCTA', 'GGAGCTAC', 'GCGTAGTA', 'CGGAGCCT']
plate2_i5=['TCTACTCT', 'CTCCTTAC', 'TATGCAGT', 'TACTCCTT', 'AGGCTTAG', 'ATTAGACG', 'CGGAGAGA', 'CTAGTCGA', 'AGCTAGAA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA']

idx_array=[]
i=1
for i5 in plate2_i5:
	for i7 in plate2_i7:
		idx_well=[i5,i7,"C"+str(i)]
		idx_array.append(idx_well)
		i+=1

fq1=sys.argv[1] #read argument fq1="~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.fastq.gz"
fq2=sys.argv[2] #read argument fq2="~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.fastq.gz"
idx3=sys.argv[3] #idx3=~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I1_001.fastq.gz"
idx4=sys.argv[4] #idx4=~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I2_001.fastq.gz"

i=0
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1:
	with gzip.open(fq2, "rt") as handle2:
		with gzip.open(idx3, "rt") as handle3:
			with gzip.open(idx4, "rt") as handle4:
				with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
					with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
						for (title1, seq1, qual1), (title2, seq2, qual2), (title3,seq3,qual3), (title3,seq4,qual4) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2),FastqGeneralIterator(handle3),FastqGeneralIterator(handle4)):
							for j in idx_array:
								if seq3==j[1]+"AT" and seq4==j[0]+"GT":
									i+=1
									readname=j[2]+"_"+seq3+seq4
									outfile_fq1.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq1, qual1))
									outfile_fq2.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq2, qual2))

```
Running fastq splitter
This can be made a lot faster using something like a hash table, or limiting search space more.

```bash
python ~/src/plate2_fastqsplitter.py /volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I1_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I2_001.fastq.gz 

#then gzip
for i in *fastq; do gzip $i & done &
```
Got 253,517,850 reads total from assignment
Moving fastq.gz files to a project directory

Processing will be in: 
```bash
/volumes/seq/projects/gccACT/230306_mdamb231_test
```

Moving the two fastq files.
```bash
mkdir /volumes/seq/projects/gccACT/230306_mdamb231_test
cp /volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.barc.fastq.gz /volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.barc.fastq.gz /volumes/seq/projects/gccACT/230306_mdamb231_test
```

# Processing Samples via Copykit to start
## Alignment
```bash
#set up variables and directory
mkdir /volumes/seq/projects/gccACT/230306_mdamb231_test
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
fq1="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.barc.fastq.gz"
fq2="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.barc.fastq.gz"

#Map reads with BWA Mem
bwa mem -t 20 $ref $fq1 $fq2 | samtools view -b - > $dir/230306_gccact.bam
```

## Split out single-cells
```bash
#set up cell directory
mkdir $dir/cells
```
Split by readname bam field into cells subdir and add metadata column

```bash
samtools view $dir/230306_gccact.bam | awk -v dir=$dir 'OFS="\t" {split($1,a,":"); print $0,"XM:Z:"a[1] > "./cells/"a[1]".230306_gccact.sam"}'
 #split out bam to cell level sam
```

## Add header to each sam and convert to bam and sort
Using parallel to save time
```bash
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"

add_header() {
	samtools view -bT $ref $1 | samtools sort -T $dir -o ${1::-4}.bam -
}
export -f add_header
sam_in=`ls *sam`
parallel --jobs 30 add_header ::: $sam_in
```

## Mark duplicate reads
```bash
#name sort, fix mates, sort by position, mark dup
sort_and_markdup() {
  samtools sort -T . -n -o - $1 | samtools fixmate -m - -| samtools sort -T . -o - - | samtools markdup -s - ${1::-4}.rmdup.bam 2> ${1::-4}.rmdup.stats.txt
}
export -f sort_and_markdup

bam_in=`ls *bam`
parallel --jobs 30 sort_and_markdup ::: $bam_in
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
#C100 had zero reads, gzipping to prevent copykit from reading in

```

## Run CopyKit for WGS portion
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html

```R
library(copykit)
library(BiocParallel)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/230306_mdamb231_test/cells")

tumor <- runVarbin("/volumes/seq/projects/gccACT/230306_mdamb231_test/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

# Mark euploid cells if they exist
tumor <- findAneuploidCells(tumor)

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

### Count of WGS and GCC Reads
```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
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
setwd("/volumes/seq/projects/gccACT/230306_mdamb231_test")
dat<-read.table("./cells/read_count.csv",header=T,sep=",")
dat$total_reads<-dat$near_cis+dat$distal_cis+dat$trans

plt1<-ggplot(dat,aes(y=near_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis Reads")+theme_minimal()
plt2<-ggplot(dat,aes(y=distal_cis,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis Reads")+theme_minimal()
plt3<-ggplot(dat,aes(y=trans,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans Reads")+theme_minimal()

plt4<-ggplot(dat,aes(y=(near_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis % Reads")+theme_minimal()
plt5<-ggplot(dat,aes(y=(distal_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis % Reads")+theme_minimal()
plt6<-ggplot(dat,aes(y=(trans/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans % Reads")+theme_minimal()

plt<-(plt1|plt2|plt3)/(plt4|plt5|plt6)
ggsave(plt,file="read_counts.pdf")
```

### Generation of HiC Contact Matrices
Merge bam files based on CopyKit output. Then using bam2pairs from pairix to generate contacts
```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

#filter reads to HiC contacts, then perform bam2pairs 
mkdir $dir/cells/contacts

#### For generation of pairix per cell
#bam_to_pairs() {
#	samtools view -F 1024 $3 | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $2 - > $1/cells/contacts/${3::-4}.contacts.bam && wait;
#	bam2pairs $1/cells/contacts/${3::-4}.contacts.bam $1/cells/contacts/${3::-4}
#}
#export -f bam_to_pairs

#cd $dir/cells
#bam_in=`ls *bam`
#parallel --jobs 50 bam_to_pairs $dir $ref {} ::: $bam_in &
# first argument is directory
# second is reference fasta
# third is bam file input

bamlist_merge_to_pairs() {
	samtools merge -b $3 -O SAM -@ 20 - | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $2 - > $1/cells/contacts/${3::-13}.contacts.bam && wait;
	bam2pairs $1/cells/contacts/${3::-13}.contacts.bam $1/cells/contacts/${3::-13}
}
export -f bamlist_merge_to_pairs

bamlist_merge_to_pairs $dir $ref clone_c1.bam_list.txt 
bamlist_merge_to_pairs $dir $ref clone_c2.bam_list.txt 
bamlist_merge_to_pairs $dir $ref clone_c3.bam_list.txt 
bamlist_merge_to_pairs $dir $ref clone_c4.bam_list.txt 

#set variable for bam list in function


```


## Using cooler and cooltools (seems like the 4DN preferred format)
Cooler is both python line and command line, using command line for this
https://github.com/open2c/cooler
https://github.com/open2c/cooltools


Change to cooler environment ::sunglasses::
```bash
conda deactivate #get out of r3.4 env
conda activate cooler_env #use cooler env (lower python version)
```

### Set up reference and bins
Prepare chrom.sizes file from reference fasta
```bash
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
faidx $ref -i chromsizes > hg38.chrom.sizes
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
cooltools genome binnify $CHROMSIZES_FILE 1000000 > $BINS_PATH & #1mb bins, 
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

```

Generate Cooler matrices from pairix data
```bash
# Note that the input pairs file happens to be space-delimited, so we convert to tab-delimited with `tr`.
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"

pairix_to_cooler() {
cooler cload pairix -p 10 --assembly hg38 $2 $1 ${1::-9}.${3}.cool
}
export -f pairix_to_cooler

cd $dir/rm_archive/contacts
pairix_in=`ls clone_*pairs.gz`

#5kb
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins.bed"
parallel --jobs 5 pairix_to_cooler ::: $pairix_in $BINS_BED_PATH 5kb& #uses 10 cores per job
#10kb
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins.bed"
parallel --jobs 5 pairix_to_cooler ::: $pairix_in $BINS_BED_PATH 10kb & #uses 10 cores per job
#50kb
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins.bed"
parallel --jobs 5 pairix_to_cooler ::: $pairix_in $BINS_BED_PATH 50kb & #uses 10 cores per job

# first argument is pairix gzipped file
```

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
	ax=sns.heatmap(ax=axes[ax_row,ax_col],data=mat,cbar=False,yticklabels=False,square=False,xticklabels=False,vmin=zmin,vmax=zmax,cmap=cmap,
		gridspec_kw={"truncate":False})
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
	cooler.coarsen_cooler(coolfile,in_name+"_5mb.cool",factor=5,chunksize=10000000) #coarsen to 5mb
	c = cooler.Cooler(in_name+"_5mb.cool")
	cooler.balance_cooler(c,store=True) #balance matrix
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
os.chdir('/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clones')

#in files
in_files=["clone_c1.bsorted.cool",
	"clone_c2.bsorted.cool",
	"clone_c3.bsorted.cool",
	"clone_c4.bsorted.cool"]

#Set up settings for plot
dpi= 300
colormap='fall'
zmin=-3
zmax=-1

all_by_all_plot(in_files[0],4,zmin,zmax,cmap="PuBu")
all_by_all_plot(in_files[2],4,zmin,zmax,cmap="YlGn")

#make all by all plot for all clones
[all_by_all_plot(x,4,zmin,zmax) for x in in_files]

#colormap="coolwarm"
#for x in in_files:
#	for y in in_files:
#		all_by_all_plot_subtract(x,y,4)

#

#plot difference of one matrix from another
#get expected - observed?
#expected = cooltools.expected_cis(clr, view_df=hg38_arms, nproc=2, chunksize=1_000_000)


```

Detection of structural variants in HiC clone cool files. Using a kernel density estimate to locate them.


Using https://github.com/GaoLabXDU/HiSV/blob/main/HiSV_code/HiSV.py for starters

https://cooltools.readthedocs.io/en/latest/notebooks/dots.html

https://cooltools.readthedocs.io/en/latest/cooltools.html#module-cooltools.api.dotfinder

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import bioframe

import os
import seaborn as sns
import cooltools.lib.plotting
import matplotlib


#wd 
os.chdir('/volumes/seq/projects/gccACT/230306_mdamb231_test/rm_archive/contacts/clones') #note rm_archive is temporary to not interfere with mariam processing

#in files
in_files=["clone_c1.bsorted.cool",
	"clone_c2.bsorted.cool",
	"clone_c3.bsorted.cool",
	"clone_c4.bsorted.cool"]

#Set up settings for plot
dpi= 300
colormap='fall'
zmin=-3
zmax=-1
chr_count=4


# Use bioframe to fetch the genomic features from the UCSC. from https://cooltools.readthedocs.io/en/latest/notebooks/dots.html
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)


#def make_all_by_all_matrix(infile_name,chr_count,zmin,zmax,cmap):
	"""Function to read in cooler file from given string name. And run all by all chr comparison"""
	#Read in Cooler File

in_file=in_files[0]
#in_file=infile_name
in_name=in_file.split(sep=".")[0]
coolfile=in_file
c = cooler.Cooler(coolfile)
#cooler.coarsen_cooler(coolfile,in_name+"_5mb.cool",factor=5,chunksize=10000000) #coarsen to 5mb
#c = cooler.Cooler(in_name+"_5mb.cool")
cooler.balance_cooler(c,store=True) #balance matrix
obs_mat = c.matrix()[:]
obs_mat.shape

# intra-arm expected
expected = cooltools.expected_cis(
    c,
    view_df=hg38_arms,
    nproc=4,
)

binsize=10_000_000
kernels = cooltools.api.dotfinder.recommend_kernels(binsize)

# create a grid of coordinates from -5 to 5, to define round kernels
# see https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html for details
half = 5  # half width of the kernel
x, y = np.meshgrid(
    np.linspace(-half, half, 2*half + 1),
    np.linspace(-half, half, 2*half + 1),
)
# now define a donut-like mask as pixels between 2 radii: sqrt(7) and sqrt(30):
mask = (x**2+y**2 > 7) & (x**2+y**2 <= 30)
mask[:,half] = 0
mask[half,:] = 0

# lowleft mask - zero out neccessary parts
mask_ll = mask.copy()
mask_ll[:,:half] = 0
mask_ll[half:,:] = 0

# new kernels with more round donut and lowleft masks:
kernels_round = {'donut': mask,
 'vertical': kernels["vertical"].copy(),
 'horizontal': kernels["horizontal"].copy(),
 'lowleft': mask_ll}

dots_df = cooltools.dots(
    c,
    expected=expected,
    view_df=hg38_arms,
    kernels=5,
    # how far from the main diagonal to call dots:
    max_loci_separation=10_000_000,
    clustering_radius=None,
    cluster_filtering=False,
    nproc=4,
)

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


all_by_all_plot(in_files[0],4,zmin,zmax,cmap="PuBu")
all_by_all_plot(in_files[2],4,zmin,zmax,cmap="YlGn")


```

Generate eigengenes for compartments across chromosomes
```bash
cooler_eigen() {
cooler balance -p 10 -f -c 10000 $1

cooltools eigs-cis -o outputs/test.eigs.100000 --view data/view_hg38.tsv --phasing-track outputs/gc.100000.tsv --n-eigs 1 $cool_file::resolutions/100000
}
export -f cooler_balance
```


Get structural variants across merged data with HiSV
https://github.com/GaoLabXDU/HiSV

```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/cells/contacts/clones"
CHROMSIZES_FILE="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/hg38.chrom.sizes"
cd $dir

python ~/tools/HiSV/convert_type/hiccovert.py -o $dir/converted_output -r $CHROMSIZES_FILE -b 1000000 -t cool -m $dir/clone_c1.bsorted.cool  
#typo in python script name

#convert outputs with name format sample_1000kb_chr10_chr21_matrix.txt. need to be renamed and moved to directory
rename -v 'sample_1000kb_' '' $dir/converted_output/Inter_matrix/*txt
rename -v '_matrix.' '.' $dir/converted_output/Inter_matrix/*txt
mkdir $dir/hic_files
mv $dir/converted_output/Inter_matrix/*txt $dir/hic_files

rename -v 'sample_1000kb_' '' $dir/converted_output/Intra_matrix/*txt
rename -v "_matrix." "." $dir/converted_output/Intra_matrix/*txt
for file in $dir/converted_output/Intra_matrix/*txt; do mv "$file" "${file/chr*_/}"; done
mkdir $dir/hic_files
mv $dir/converted_output/Intra_matrix/*txt $dir/hic_files

python ~/tools/HiSV/HiSV_code/HiSV.py -o $dir \
               -l $CHROMSIZES_FILE \
               -f $dir/hic_files \
               -a 1000000 \
               -e 1000000 \
               -w 10 \
               -c 0.4
```


<!--
#Add compartments
#Add SV detection by contact matrix
#Use cooltools virtual4c for eccDNA interactions
#Add 
-->




<!--
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