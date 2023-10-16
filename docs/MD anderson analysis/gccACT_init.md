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

```bash
conda activate r4.2
```
```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
register(MulticoreParam(progressbar = T, workers = 50), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/230306_mdamb231_test")

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

tumor <- runUmap(tumor)
k_clones<-findSuggestedK(tumor)
# Create a umap embedding 

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
tumor  <- findClusters(tumor,k_subclones=25)#output from k_clones

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

saveRDS(tumor,file="scCNA.rds")
tumor<-readRDS("scCNA.rds")
clone_out<-data.frame(bam=paste0(getwd(),"/cells/",row.names(tumor@colData),".bam"),clone=tumor@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("clone_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}



```


## Correlating our HiC Data to ONT Data
See /gccACT_ONT/ for more details on ONT processing.

```bash
#Data stored here:
/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output

#CNV calls using QDNA algorithm:
/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/qdna_seq/20230726_1239_2D_PAO38369_output_combined.bed


#SV calls using Sniffles2:
/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/20230726_1239_2D_PAO38369_output.wf_sv.vcf.gz


```
Plotting CNV calls from ONT data together with gccACT data
```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
setwd("/volumes/seq/projects/gccACT/230306_mdamb231_test")
tumor<-readRDS(file="/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
ont<-read.table("/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/qdna_seq/20230726_1239_2D_PAO38369_output_combined.bed",sep="\t",skip=1)

colnames(ont)<-c("chrom","start","end","range","log2ratio","pass_filt","cnv")
ont$chrom<-paste0("chr",ont$chrom)
ont<-makeGRangesFromDataFrame(ont,keep.extra.columns=TRUE)

#https://rdrr.io/github/navinlabcode/copykit/src/R/plotHeatmap.R
copykit_ranges<-makeGRangesFromDataFrame(tumor@rowRanges,keep.extra.columns=TRUE)
overlap<-findOverlaps(subject=ont,query=copykit_ranges,select="first")

ont_overlaps<-cbind(as.data.frame(copykit_ranges),as.data.frame(ont)[overlap,c("log2ratio","cnv")])

plt<-plotHeatmap(tumor, label = 'subclones',order='hclust')

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.pdf")
print(plt)
dev.off()

seg_data_ont<-t(as.data.frame(ont_overlaps$log2ratio))
seg_data <- t(SummarizedExperiment::assay(tumor, "logr"))
dim(seg_data_ont)[2]==dim(seg_data)[2]

plt@matrix<-seg_data_ont
plt@row_order<-c(1)
# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.ont.pdf")
print(plt)
dev.off()
#manipulated pdf in illustrator to make the ONT data in line with other data, both pdfs are printed to match width so they can be stacked for easy viewing
```

Reading in SV VCF file

```R
library(vcfR)
library(dplyr)
vcf_file <- "/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/20230726_1239_2D_PAO38369_output.wf_sv.vcf.gz"
vcf <- read.vcfR(vcf_file, verbose = TRUE)
vcf_field_names(vcf, tag = "FORMAT")
vcf_field_names(vcf, tag = "INFO")
Z <- vcfR2tidy(vcf, format_fields = c("GT", "DP"))
dat<-Z$fix
dat %>% group_by(SVTYPE) %>% summarize(mean(SVLEN))

dat_transloc<-dat[dat$SVTYPE=="BND",]
dat_dup_inv_del<-dat[(!is.na(dat$SVLEN)),]
dat_dup_inv_del<-dat_dup_inv_del[(abs(dat_dup_inv_del$SVLEN)>1000000),]
dat_out<-rbind(dat_transloc,dat_dup_inv_del)
write.table(as.data.frame(dat_out),file="/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/20230726_1239_2D_PAO38369_output.1mb_SVs.tsv",col.names=T,row.names=T,sep="\t")

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

plt4<-ggplot(dat,aes(y=(near_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Near Cis % Reads")+theme_minimal()+ylim(c(0,100))
plt5<-ggplot(dat,aes(y=(distal_cis/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Distal Cis % Reads")+theme_minimal()+ylim(c(0,10))
plt6<-ggplot(dat,aes(y=(trans/total_reads)*100,x="Cells"))+geom_jitter()+geom_boxplot()+ylab("Trans % Reads")+theme_minimal()+ylim(c(0,10))

plt<-(plt1|plt2|plt3)/(plt4|plt5|plt6)
ggsave(plt,file="read_counts.pdf")
```

## Project Library Complexity
Using Picard Tools
```bash
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
project_count() {
java -jar /volumes/seq/code/3rd_party/picard/picard-2.20.4/picard.jar EstimateLibraryComplexity I=$1 O=${1::-4}.complex_metrics.txt
}
export -f project_count

cd $dir/cells
bam_in=`ls *bam`
parallel --jobs 20 project_count ::: $bam_in
```
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

https://github.com/mdozmorov/HiC_tools#cnv-aware-normalization

### Generation of HiC Contact Matrices
Merge bam files based on CopyKit output. Then using bam2pairs from pairix to generate contacts


```bash
conda activate EagleC
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

#filter reads to HiC contacts, then perform bam2pairs 
mkdir $dir/contacts

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
	samtools merge -b $3 -O SAM -@ 20 - | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $2 - > $1/contacts/${3::-13}.contacts.bam && wait;
	bam2pairs $1/contacts/${3::-13}.contacts.bam $1/contacts/${3::-13}
}
export -f bamlist_merge_to_pairs

cd $dir
bamlist_merge_to_pairs $dir $ref clone_c1.bam_list.txt &
bamlist_merge_to_pairs $dir $ref clone_c2.bam_list.txt &
bamlist_merge_to_pairs $dir $ref clone_c3.bam_list.txt &
bamlist_merge_to_pairs $dir $ref clone_c4.bam_list.txt &


#set variable for bam list in function


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
conda activate EagleC #use cooler env (lower python version)
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
GC_BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.gc.bed"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 1000000 > $BINS_PATH & #1mb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
bedtools nuc -fi $ref -bed $BINS_BED_PATH > $GC_BINS_PATH
awk 'OFS="\t" {print $1,$2,$3,$5}' $GC_BINS_PATH > $GC_BINS_BED_PATH

# #5KB Bins
# BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins"
# GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.gc.bins"
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins.bed"
# cooltools genome binnify $CHROMSIZES_FILE 5000 > $BINS_PATH  #5kb bins, 
# tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
# cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &
# #10KB Bins
# BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins"
# GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.gc.bins"
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins.bed"
# cooltools genome binnify $CHROMSIZES_FILE 10000 > $BINS_PATH #10kb bins, 
# tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
# cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

# #50KB Bins
# BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins"
# GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.gc.bins"
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins.bed"
# cooltools genome binnify $CHROMSIZES_FILE 50000 > $BINS_PATH #50kb bins, 
# tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
# cooltools genome gc $BINS_PATH $FASTA_PATH > $GC_BINS_PATH &

#500KB Bins
BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins"
GC_BINS_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.gc.bins"
GC_BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.gc.bed"
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
cooltools genome binnify $CHROMSIZES_FILE 500000 > $BINS_PATH #50kb bins, 
tail -n +2 $BINS_PATH |  head -n -1 > $BINS_BED_PATH #remove the header and hanging line to make it a proper bed file
bedtools nuc -fi $ref -bed $BINS_BED_PATH > $GC_BINS_PATH
awk 'OFS="\t" {print $1,$2,$3,$5}' $GC_BINS_PATH > $GC_BINS_BED_PATH

```

### Generate Cooler matrices from pairix data
```bash
# Note that the input pairs file happens to be space-delimited, so we convert to tab-delimited with `tr`.
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"

# pairix_to_cooler_5kb() {
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/5kb.bins.bed"
# out_name="5kb"
# cooler cload pairix -p 10 --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
# }
# export -f pairix_to_cooler_5kb

# pairix_to_cooler_10kb() {
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/10kb.bins.bed"
# out_name="10kb"
# cooler cload pairix -p 10 --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
# }
# export -f pairix_to_cooler_10kb

# pairix_to_cooler_50kb() {
# BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/50kb.bins.bed"
# out_name="50kb"
# cooler cload pairix -p 10 --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
# }
# export -f pairix_to_cooler_50kb


pairix_to_cooler_500kb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
out_name="500kb"
cooler cload pairix -0 -p 10 --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_500kb


pairix_to_cooler_1mb() {
BINS_BED_PATH="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"
out_name="1mb"
cooler cload pairix -0 -p 10 --assembly hg38 $BINS_BED_PATH $1 ${1::-9}.${out_name}.cool
}
export -f pairix_to_cooler_1mb



cd $dir/contacts
pairix_in=`ls clone_*pairs.gz`

#5kb
#parallel --jobs 6 pairix_to_cooler_5kb ::: $pairix_in & #uses 10 cores per job

#10kb
#parallel --jobs 6 pairix_to_cooler_10kb ::: $pairix_in & #uses 10 cores per job

#50kb
#parallel --jobs 6 pairix_to_cooler_50kb ::: $pairix_in & #uses 10 cores per job

#500kb
parallel --jobs 6 pairix_to_cooler_500kb ::: $pairix_in & #uses 10 cores per job

#1mb
parallel --jobs 6 pairix_to_cooler_1mb ::: $pairix_in & #uses 10 cores per job

# first argument is pairix gzipped file
```

### CNV normalization on HiC Data
Write out CNV segments as bedfiles at multiple resolutions for cnv correction via NeoLoopFinder
NeoLoopFinder also reports CNVs through log2 changes, making this a direct comparison.
```bash
conda activate r4.2
```

```R
library(copykit)
library(GenomicRanges)
library(parallel)
wd_out="/volumes/seq/projects/gccACT/230306_mdamb231_test/"

setwd(wd_out)
tumor<-readRDS("/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
unique(tumor$subclones)
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
saveRDS(tumor,file="scCNA.consensus.rds") #save categorical consensus scCNA object


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


BINS_BED_PATH_500kb="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/500kb.bins.bed"
BINS_BED_PATH_1mb="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"

#Make blacklist
lapply(c(BINS_BED_PATH_500kb,BINS_BED_PATH_1mb), function(x) make_black_list(x,copykit_obj=tumor))


clone_list<-colnames(tumor@consensus)

#Make clone specific bedgraph format for HiC windows
for(clone in clone_list){
	make_consensus_bedgraph(bed_in=BINS_BED_PATH_500kb,res="500kb",clone=clone,copykit_obj=tumor,cores=50)
	make_consensus_bedgraph(bed_in=BINS_BED_PATH_1mb,res="1mb",clone=clone,copykit_obj=tumor,cores=50)
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
tumor<-readRDS("scCNA.consensus.rds")


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
BINS_BED_PATH_1mb="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"

clone_list<-colnames(tumor@consensus)

#Make clone specific bedgraph format for HiC windows
for(clone in clone_list){
	make_consensus_bigwig(bed_in=BINS_BED_PATH_500kb,res="500kb",clone=clone,copykit_obj=tumor,cores=50)
	make_consensus_bigwig(bed_in=BINS_BED_PATH_1mb,res="1mb",clone=clone,copykit_obj=tumor,cores=50)
}

```

CNV normalized output at single 500kb resolution.

```bash
conda activate EagleC
```

```bash

eaglec_SV_CNV_neoout() {
clone_in=$1
cpus=$2
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/contacts"
bedgraph_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
clone=${clone_in::-17}

#balance through ICE first
cooler balance \
	--ignore-diags 3 \
	--force \
	--nproc $cpus \
	${clone_dir}/${clone}.bsorted.500kb.cool

#now balance by CNV
correct-cnv \
	--ignore-diags 3 \
	-H ${clone_dir}/${clone}.bsorted.500kb.cool \
	--cnv-file ${bedgraph_dir}/${clone}.cnv.500kb.segmented.bedgraph \
	--force \
	--nproc $cpus \
	--logFile ${clone}.cnv.norm.log

#predict SV breakpoints by cnv corrected data
predictSV-single-resolution \
	-H ${clone_dir}/${clone}.bsorted.500kb.cool \
	-g hg38 \
	--output-file ${clone_dir}/${clone}.SV.cnv.neoloopfinder.tsv \
	--balance-type CNV \
	--output-format NeoLoopFinder \
	--prob-cutoff 0.8 \
	--logFile ${clone}.eaglec.cnv.log \
	--cache-folder ${clone}.cache

#assemble cnv breakpoints
assemble-complexSVs \
	-H ${clone_dir}/${clone}.bsorted.500kb.cool \
	-O ${clone_dir}/${clone} \
	-B ${clone_dir}/${clone}.SV.cnv.neoloopfinder.tsv \
	--balance-type CNV \
	--protocol insitu \
	--nproc $cpus \
	--region-size 1000000 \
	--minimum-size 1000000 \
	--logFile ${clone}.neoloop.log 
}
export -f eaglec_SV_CNV_neoout

cd /volumes/seq/projects/gccACT/230306_mdamb231_test/contacts
clones=`ls clone*pairs.gz`
parallel --jobs 1 eaglec_SV_CNV_neoout {} 50 ::: $clones & 

```

### Plot all by all autosome interactions
Use this for all chromosomes plots, code adapted from cooltools and inspired by:
https://github.com/bianlab-hub/zuo_ncomms_2021/blob/Hi-C_data_analysis/fig1f_plot_obs_heatmap.py


```bash
conda activate EagleC
```

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
import bioframe

###Prepare Input Files Functions###
#read in eagleC called translocations, plot on top of diagonal
def prepare_eaglec_file(eaglec_in="clone_c3.SV.cnv.neoloopfinder.tsv"):
	""" Function to take in string file name of neoloopfinder output from eaglec_SV_CNV_neoout function (run above in page). Reads in files. 
	Sets transchr SVs to be above diagonal. Currently ignores cis SVs."""
	eaglec_dat=pd.read_csv(eaglec_in,sep="\t",names=["CHROM","CHR2","DIR","POS","CHR2POS","SVTYPE"])
	eaglec_dat_cis=eaglec_dat[eaglec_dat["CHROM"]==eaglec_dat["CHR2"]]
	for index, i in eaglec_dat_cis.iterrows():
		if eaglec_dat_cis.loc[index,"POS"] > eaglec_dat_cis.loc[index,"CHR2POS"]:
				print("Swapping "+eaglec_dat_cis.loc[index,"CHROM"]+" start and end cis SV.")
				tmp=eaglec_dat_cis.loc[index,:]
				eaglec_dat_cis.loc[index,"CHR2POS"]=tmp.POS
				eaglec_dat_cis.loc[index,"POS"]=tmp.CHR2POS
	eaglec_dat_trans=eaglec_dat[eaglec_dat["CHROM"]!=eaglec_dat["CHR2"]]
	for index, i in eaglec_dat_trans.iterrows():
		if (eaglec_dat_trans.loc[index,"CHROM"] in chr_number.keys()) & (eaglec_dat_trans.loc[index,"CHR2"] in chr_number.keys()):
			if chr_number[eaglec_dat_trans.loc[index,"CHROM"]] < chr_number[eaglec_dat_trans.loc[index,"CHR2"]]:
				print("Swapping "+eaglec_dat_trans.loc[index,"CHROM"]+" and "+eaglec_dat_trans.loc[index,"CHR2"])
				tmp=eaglec_dat_trans.loc[index,:]
				eaglec_dat_trans.loc[index,"CHROM"]=tmp.CHR2
				eaglec_dat_trans.loc[index,"CHR2"]=tmp.CHROM
				eaglec_dat_trans.loc[index,"POS"]=tmp.CHR2POS
				eaglec_dat_trans.loc[index,"CHR2POS"]=tmp.POS
	return([eaglec_dat_cis,eaglec_dat_trans])

#read in ont data and split for cis and trans interactions
def setup_ont_data(ont_in):
	ont_dat=pd.read_csv(ont_in,sep="\t")
	ont_dat_trans=ont_dat[(ont_dat['SVTYPE'] == "BND")]
	ont_dat_cis=ont_dat[((ont_dat['SVTYPE'].isin(['INV', 'DEL', 'DUP'])) & (ont_dat['SVLEN']>1000000))]
	for index, i in ont_dat_cis.iterrows(): #set chr and position order, so its always plotted below diagonal
		if ont_dat_cis.loc[index,"POS"] > ont_dat_cis.loc[index,"END"]:
				print("Swapping "+ont_dat_cis.loc[index,"CHROM"]+" start and end cis SV.")
				tmp=ont_dat_cis.loc[index,:]
				ont_dat_cis.loc[index,"END"]=tmp.POS
				ont_dat_cis.loc[index,"POS"]=tmp.END
	ont_dat_trans.loc[:,'CHR2POS']=[int(x.replace("[",",").replace("]",",").split(",")[1].split(":")[1]) for x in ont_dat_trans['ALT'].tolist()]
	for index, i in ont_dat_trans.iterrows(): #set chr and position order, so its always plotted below diagonal
		if (ont_dat_trans.loc[index,"CHROM"] in chr_number.keys()) & (ont_dat_trans.loc[index,"CHR2"] in chr_number.keys()):
			if chr_number[ont_dat_trans.loc[index,"CHROM"]]<chr_number[ont_dat_trans.loc[index,"CHR2"]]:
				print("Swapping "+ont_dat_trans.loc[index,"CHROM"]+" and "+ont_dat_trans.loc[index,"CHR2"])
				tmp=ont_dat_trans.loc[index,:]
				ont_dat_trans.loc[index,"CHROM"]=tmp.CHR2
				ont_dat_trans.loc[index,"CHR2"]=tmp.CHROM
				ont_dat_trans.loc[index,"POS"]=tmp.CHR2POS
				ont_dat_trans.loc[index,"CHR2POS"]=tmp.POS
	return([ont_dat_cis,ont_dat_trans])

###Contact Map Plotting Functions###
#function to plot chromosomes
def trans_chr_plot(mat,chr_row,chr_col,ax_row,ax_col,axes,zmin,zmax,cmap,xlim,ylim,mat_bins,ont_dat_trans,ont_dat_cis,eaglec_dat_cis,eaglec_dat_trans):
	""" Function to plot chr by chr trans interactions (or chr by chr cis interactions if same chr given)
	Takes in c for cooler file and str for chr_row and chr_col (format "chr1","chrX", etc.)
	Takes in integers for ax_row and ax_col to add output plot as subplot.
	"""
	mat=np.log10(mat+1e-9)#adding small value for log
	ax=sns.heatmap(ax=axes[ax_row,ax_col],data=mat,cbar=False,yticklabels=False,square=False,xticklabels=False,vmin=zmin,vmax=zmax,cmap=cmap)
	ax.set_ylim(ylim,0)
	ax.set_xlim(0,xlim)
	if ax_col==0: #set chr name if it is first column
		axes[ax_row,ax_col].set_ylabel(chr_row,fontsize=24)
	print("Completed "+chr_row+" by "+chr_col+" plot")
	if (chr_col==chr_row): #cis SVs
		if chr_row in ont_dat_cis.CHROM.tolist(): #ONT data
			print("Adding ONT data cis annotations for "+chr_row)
			add_ont_sv_cis(ax,mat_bins,
				chr_row=chr_row,chr_col=chr_col,
				ont_dat_cis=ont_dat_cis)
		if chr_row in eaglec_dat_cis.CHROM.tolist(): #Eagle C data
			print("Adding EagleC data cis annotations for "+chr_row)
			add_eaglec_sv_cis(ax,mat_bins,
				chr_row=chr_row,chr_col=chr_col,
				eaglec_dat_cis=eaglec_dat_cis)
	else: #trans SVs
		if any((ont_dat_trans["CHROM"]==chr_row) & (ont_dat_trans["CHR2"]==chr_col)): #ont data
			print("Adding ONT data trans annotations for "+chr_row+" and "+chr_col)
			add_ont_sv_trans(ax,mat_bins,
				chr_row=chr_row,chr_col=chr_col,
				ont_dat_trans=ont_dat_trans)
		if any((eaglec_dat_trans["CHROM"]==chr_row) & (eaglec_dat_trans["CHR2"]==chr_col)): #eagle C data
			print("Adding EagleC trans annotations for "+chr_row+" and "+chr_col)
			add_eaglec_sv_trans(ax,mat_bins,
				chr_row=chr_row,chr_col=chr_col,
				eaglec_dat_trans=eaglec_dat_trans)
	return(ax)

#wrapper function to make panels and loop through chromosome pairs
def all_by_all_plot(ont_dat_cis, ont_dat_trans,eaglec_dat_trans,eaglec_dat_cis,cnv_dat,infile_name="clone_c0.bsorted.1mb.cool", chr_count=4, zmin=-3, zmax=-1, cmap="Reds", dpi=300): 
	"""Function to read in cooler file from given string name. And run all by all chr comparison"""
	#Read in Cooler File
	in_file=infile_name
	in_name=in_file.split(sep=".")[0]
	coolfile=in_file
	print("Reading in "+infile_name)
	c = cooler.Cooler(coolfile)
	#cooler.coarsen_cooler(coolfile,in_name+"_5mb.cool",factor=5,chunksize=5000000) #coarsen to 5mb
	#c = cooler.Cooler(in_name+"_5mb.cool")
	print("Balancing matrix.")
	cooler.balance_cooler(c,store=True,rescale_marginals=True,ignore_diags=1)#balance matrix ignore_diags=10,
	obs_mat = c.matrix()[:]
	chr_list=list(c.bins()[:]["chrom"].unique())
	out=''.join([in_name,'_all_by_all_log2_1Mb_obs.test.png'])
	out_pdf=''.join([in_name,'_all_by_all_log2_1Mb_obs.test.pdf'])
	#init subplot
	chr_in=len(chr_list)
	chr_sizes=pd.DataFrame(c.bins()[:]).groupby(["chrom"])["chrom"].count()
	chr_ratios=list(chr_sizes/chr_sizes[0])
	height_ratio=[chr_ratios[0]/4]+[chr_ratios[0]/4]+chr_ratios[0:chr_count]#added for insulation score and cnv rows
	width_ratio=chr_ratios[0:chr_count]
	print("Calculating Insulation Values.")
	resolution=cnv_dat.loc[0,"end"]-cnv_dat.loc[0,"start"]
	windows = [1*resolution,3*resolution, 5*resolution, 10*resolution]
	insulation_table = cooltools.insulation(c, windows, verbose=True)
	fig, axes = plt.subplots(nrows=chr_count+2, ncols=chr_count, figsize=(sum(height_ratio)*5, sum(width_ratio)*5),sharex='col', sharey='row',
		gridspec_kw={'height_ratios':height_ratio,
		'width_ratios': width_ratio}) #,sharex=True,sharey=True,
	print("Adding column annotations.")
	for j in range(0,chr_count):
			col_chrom=chr_list[j]
			ax_row=0 #0 row is for cnv data
			ax_col=j
			mat=c.matrix().fetch(col_chrom)
			mat_bins=c.bins()[:]
			mat_bins=mat_bins[mat_bins["chrom"].isin([col_chrom])]
			cnv_plot(mat=mat,
				chr_col=col_chrom,
				ax_row=ax_row,
				ax_col=ax_col,
				axes=axes,
				mat_bins=mat_bins,
				cnv_dat=cnv_dat)
			add_insulation_scores(mat=mat,
				chr_col=col_chrom,
				ax_row=ax_row+1,
				ax_col=ax_col,
				axes=axes,
				mat_bins=mat_bins,
				insulation_table=insulation_table,
				resolution=resolution,
				windows=windows)
	print("Adding subpanels of chr-chr contacts.")
	for i in range(0,chr_count):
			row_chrom=chr_list[i]
			ax_row=i+2
			xlim=chr_sizes[i]
			for j in range(0,chr_count):
				col_chrom=chr_list[j]
				ax_col=j
				ylim=chr_sizes[j]
				mat=c.matrix().fetch(row_chrom,col_chrom)
				mat_bins=c.bins()[:]
				mat_bins=mat_bins[mat_bins["chrom"].isin([row_chrom,col_chrom])]
				trans_chr_plot(mat=mat,
					chr_row=row_chrom,chr_col=col_chrom,
					ax_row=ax_row,ax_col=ax_col,axes=axes,
					zmin=zmin,zmax=zmax,cmap=cmap,xlim=xlim,ylim=ylim,
					mat_bins=mat_bins,
					ont_dat_trans=ont_dat_trans,ont_dat_cis=ont_dat_cis,
					eaglec_dat_trans=eaglec_dat_trans,eaglec_dat_cis=eaglec_dat_cis)
	plt.subplots_adjust(wspace=0.02, hspace=0.02)
	plt.savefig(out, dpi=dpi, format='png',bbox_inches='tight')
	plt.savefig(out_pdf, format='pdf',bbox_inches='tight')
	plt.close("all")
	return(insulation_table)

###Add circle patches over detected SVs Functions###
#adds ont circles around cis SVs
def add_ont_sv_cis(ax,mat_bins,chr_row,chr_col,ont_dat_cis,col="green"):
		"""Function to add SV from input ont_dat_cis on same chr. Adds as a circle around the contact map. """			
		for index, i in ont_dat_cis[ont_dat_cis["CHROM"]==chr_row].iterrows():
			if ont_dat_cis.loc[index,"SVTYPE"] == "INV":
				circle_col=col #inversions
			elif ont_dat_cis.loc[index,"SVTYPE"] == "DUP":
				circle_col=col #duplications
			else:
				circle_col=col #deletions
			chr_bin_start=mat_bins[(mat_bins.chrom==chr_row)].index.min()
			bin_annot=mat_bins.loc[(mat_bins.chrom==chr_row) & (mat_bins.start>=ont_dat_cis.loc[index,"POS"]) & (mat_bins.end<ont_dat_cis.loc[index,"END"])]
			bin_annot_x=bin_annot.index.min()-chr_bin_start
			bin_annot_y=bin_annot.index.max()-chr_bin_start
			patch=plt.Circle((bin_annot_x,bin_annot_y), alpha=0.5, radius=10, clip_box=False,linestyle=":", color=circle_col,lw=3,fill=False)
			ax.add_patch(patch)

#adds ont circles around trans SVs
def add_ont_sv_trans(ax,mat_bins,chr_row,chr_col,ont_dat_trans,col="green"):
		"""Function to add SV from input ont_dat_trans on chr-chr pairs. Adds as a circle around the contact map. """			
		for index, i in ont_dat_trans[(ont_dat_trans["CHROM"]==chr_row) & (ont_dat_trans["CHR2"]==chr_col)].iterrows():
			circle_col=col #translocations
			chr_bin_start_x=mat_bins[(mat_bins.chrom==chr_row)].index.min()
			chr_bin_start_y=mat_bins[(mat_bins.chrom==chr_col)].index.min()
			bin_annot_x=mat_bins.loc[(mat_bins.chrom==chr_row) & (mat_bins.start>=ont_dat_trans.loc[index,"POS"])]
			bin_annot_y=mat_bins.loc[(mat_bins.chrom==chr_col) & (mat_bins.start>=ont_dat_trans.loc[index,"CHR2POS"])]
			bin_annot_x=bin_annot_x.index.min()-chr_bin_start_x
			bin_annot_y=bin_annot_y.index.min()-chr_bin_start_y
			patch=plt.Circle((bin_annot_y,bin_annot_x), alpha=0.5,radius=10, clip_box=False,linestyle=":", color=circle_col,lw=3,fill=False)
			ax.add_patch(patch)

#adds ont circles around cis SVs
def add_eaglec_sv_cis(ax,mat_bins,chr_row,chr_col,eaglec_dat_cis,col="blue"):
		"""Function to add SV from input eagleC_dat_cis on same chr. Adds as a circle around the contact map. """			
		for index, i in eaglec_dat_cis[eaglec_dat_cis["CHROM"]==chr_row].iterrows():
			if eaglec_dat_cis.loc[index,"SVTYPE"] == "inversion":
				circle_col=col #inversions
			elif eaglec_dat_cis.loc[index,"SVTYPE"] == "duplication":
				circle_col=col #duplications
			else:
				circle_col=col #deletions and translocations
			chr_bin_start=mat_bins[(mat_bins.chrom==chr_row)].index.min()
			bin_annot=mat_bins.loc[(mat_bins.chrom==chr_row) & (mat_bins.start>=eaglec_dat_cis.loc[index,"POS"]) & (mat_bins.end<eaglec_dat_cis.loc[index,"CHR2POS"])]
			bin_annot_x=bin_annot.index.min()-chr_bin_start
			bin_annot_y=bin_annot.index.max()-chr_bin_start
			patch=plt.Circle((bin_annot_x,bin_annot_y), alpha=0.5,radius=10, clip_box=False,linestyle='--', color=circle_col,lw=3,fill=False)
			ax.add_patch(patch)

#adds eagleC circles around trans SVs
def add_eaglec_sv_trans(ax,mat_bins,chr_row,chr_col,eaglec_dat_trans,col="blue"):
		"""Function to add SV from input eagleC on chr-chr pairs. Adds as a circle around the contact map. """					
		for index, i in eaglec_dat_trans[(eaglec_dat_trans["CHROM"]==chr_row) & (eaglec_dat_trans["CHR2"]==chr_col)].iterrows():
			circle_col=col #translocations
			chr_bin_start_x=mat_bins[(mat_bins.chrom==chr_row)].index.min()
			chr_bin_start_y=mat_bins[(mat_bins.chrom==chr_col)].index.min()
			bin_annot_x=mat_bins.loc[(mat_bins.chrom==chr_row) & (mat_bins.start>=eaglec_dat_trans.loc[index,"POS"])]
			bin_annot_y=mat_bins.loc[(mat_bins.chrom==chr_col) & (mat_bins.start>=eaglec_dat_trans.loc[index,"CHR2POS"])]
			bin_annot_x=bin_annot_x.index.min()-chr_bin_start_x
			bin_annot_y=bin_annot_y.index.min()-chr_bin_start_y
			patch=plt.Circle((bin_annot_y,bin_annot_x), alpha=0.5,radius=10, clip_box=False,linestyle='--', color=circle_col,lw=3,fill=False)
			ax.add_patch(patch)

#adds bedgraph cnv over the heatmap plot
def setup_cnv_bedgraph(cnv_in,col):
	"""Function to read in CNV bedgraph as pandas DF"""
	cnv_dat=pd.read_csv(cnv_in,sep="\t",names=["chrom","start","end","cnv"])
	cnv_col_dict=dict(zip(["0","1","2","3","4","5","6"],col))
	cnv_dat["cnv_col"]=[cnv_col_dict[str(x)] for x in cnv_dat.cnv]
	return(cnv_dat)

#adds cnv plots as annotation in top row
def cnv_plot(mat,chr_col,ax_row,ax_col,axes,mat_bins,cnv_dat):
	"""Function to add CNV profiles on top row of contact heatmaps"""
	print("Adding CNV Profile of "+chr_col)
	cnv_dat_tmp=cnv_dat[cnv_dat["chrom"].isin([chr_col])].copy()
	cnv_dat_tmp["cnv_one"]=[1 for x in cnv_dat_tmp.cnv]
	sns.barplot(ax=axes[ax_row,ax_col],x=cnv_dat_tmp.start,y=cnv_dat_tmp.cnv_one,palette=cnv_dat_tmp.cnv_col,saturation=1,lw=0)
	axes[ax_row,ax_col].set_xticks([])
	axes[ax_row,ax_col].set_yticks([])
	if ax_col==0:
		axes[ax_row,ax_col].set_ylabel("Copy \n Number")
	else:
		axes[ax_row,ax_col].set_ylabel(None)
	axes[ax_row,ax_col].set_ylim([0,1])
	axes[ax_row,ax_col].set_title(chr_col,fontsize=24)
	return(axes[ax_row,ax_col])

#adds detected TAD boundaries in second from top row
def add_insulation_scores(mat,chr_col,ax_row,ax_col,axes,mat_bins,insulation_table,resolution,windows):
	print("Adding insulation score for "+chr_col)
	region = (chr_col, min(mat_bins["start"]), max(mat_bins["end"]))
	insul_region = bioframe.select(insulation_table, region)
	insul_region=insul_region.assign(boundary_col=["black" if x is True else "white" for x in insul_region.loc[:,"is_boundary_"+str(resolution*3)]])
	insul_region=insul_region.assign(bound_one=[1 for x in range(0,insul_region.shape[0])])
	sns.barplot(ax=axes[ax_row,ax_col],x=insul_region.start,y=insul_region.bound_one,palette=insul_region.boundary_col,saturation=1,lw=0)
	axes[ax_row,ax_col].set_xticks([])
	axes[ax_row,ax_col].set_yticks([])
	if ax_col==0:
		axes[ax_row,ax_col].set_ylabel("TAD \n Boundaries")
	else:
		axes[ax_row,ax_col].set_ylabel(None)
	axes[ax_row,ax_col].set_ylim([0,1])
	print("Found "+str(sum(insul_region.loc[:,"is_boundary_"+str(resolution*3)]))+" boundaries in "+chr_col)
	return(axes[ax_row,ax_col])

os.chdir('/volumes/seq/projects/gccACT/230306_mdamb231_test/contacts')#wd 

#list chromosomes and set dictionary of order
L1=["chr"+str(x) for x in list(range(1,23))+["X","Y"]]
L2=list(range(1,25))
chr_number={k:v for k,v in zip(L1,L2)}

#read in ont data
ont_in="/volumes/seq/projects/gccACT/230808_mdamb231_ONT/20230726_1239_2D_PAO38369_output/20230726_1239_2D_PAO38369_output.1mb_SVs.tsv"
ont_dat=setup_ont_data(ont_in)

cnv_col_palette=["#053061","#4393c3","#f7f7f7","#fddbc7","#d6604d","#b2182b","#67001f"]

#make all by all plot for all clones, set up looping list
in_cool=[
"clone_c1.bsorted.1mb.cool",
"clone_c2.bsorted.1mb.cool",
"clone_c3.bsorted.1mb.cool"]

in_eaglec=[
"clone_c1.SV.cnv.neoloopfinder.tsv",
"clone_c2.SV.cnv.neoloopfinder.tsv",
"clone_c3.SV.cnv.neoloopfinder.tsv"]

in_bedgraph=[
"/volumes/seq/projects/gccACT/230306_mdamb231_test/clone_c1.cnv.1mb.segmented.bedgraph",
"/volumes/seq/projects/gccACT/230306_mdamb231_test/clone_c2.cnv.1mb.segmented.bedgraph",
"/volumes/seq/projects/gccACT/230306_mdamb231_test/clone_c3.cnv.1mb.segmented.bedgraph"]

for clone in range(0,len(in_cool)):
	eaglec_dat=prepare_eaglec_file(eaglec_in=in_eaglec[clone])
	cnv_dat=setup_cnv_bedgraph(in_bedgraph[clone],col=cnv_col_palette)
	all_by_all_plot(infile_name=in_cool[clone],
	chr_count=4,
	ont_dat_cis=ont_dat[0],
	ont_dat_trans=ont_dat[1],
	eaglec_dat_cis=eaglec_dat[0],
	eaglec_dat_trans=eaglec_dat[1],
	cnv_dat=cnv_dat,
	zmin=-3.5,
	cmap="fall")



eaglec_dat=prepare_eaglec_file(eaglec_in="clone_c1.SV.cnv.neoloopfinder.tsv")
cnv_dat=setup_cnv_bedgraph("/volumes/seq/projects/gccACT/230306_mdamb231_test/clone_c1.cnv.1mb.segmented.bedgraph",
	col=cnv_col_palette)
tab=all_by_all_plot(infile_name="clone_c1.bsorted.1mb.cool",
chr_count=4,
ont_dat_cis=ont_dat[0],
ont_dat_trans=ont_dat[1],
eaglec_dat_cis=eaglec_dat[0],
eaglec_dat_trans=eaglec_dat[1],
cnv_dat=cnv_dat,
zmin=-3.5,
cmap="Blues")


eaglec_dat=prepare_eaglec_file(eaglec_in="clone_c3.SV.cnv.neoloopfinder.tsv")
cnv_dat=setup_cnv_bedgraph("/volumes/seq/projects/gccACT/230306_mdamb231_test/clone_c4.cnv.1mb.segmented.bedgraph",
	col=cnv_col_palette)
tab=all_by_all_plot(infile_name="clone_c4.bsorted.1mb.cool",
chr_count=4,
ont_dat_cis=ont_dat[0],
ont_dat_trans=ont_dat[1],
eaglec_dat_cis=eaglec_dat[0],
eaglec_dat_trans=eaglec_dat[1],
cnv_dat=cnv_dat,
zmin=-3.5,
cmap="Reds")
#additional cmap color pallets here
#https://www.practicalpythonfordatascience.com/ap_seaborn_palette
#some good ones! magma, mako, rocket, twilight, viridis, vlag
```

### Generate eigengenes for compartments across chromosomes
```bash
conda activate EagleC
```
```python
#from here: https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html
# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools.lib.plotting
from cooltools import insulation
import itertools
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
import os 
os.chdir('/volumes/seq/projects/gccACT/230306_mdamb231_test/contacts')#wd 

resolution=500000
clone_dir="/volumes/seq/projects/gccACT/230306_mdamb231_test/contacts"
clone="clone_c3"
coolfile=clone_dir+"/"+clone+".bsorted.500kb.cool"
c = cooler.Cooler(coolfile)
windows = [1*resolution,3*resolution, 5*resolution, 10*resolution]
insulation_table = insulation(c, windows, verbose=True)
bp_formatter = EngFormatter('b')

# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


plt.rcParams['font.size'] = 12
start = 1000000
end = start+ 100*windows[0]
region = ('chr4', start, end)
norm = LogNorm(vmax=0.1, vmin=0.001)
data = c.matrix(balance=True).fetch(region)
f, ax = plt.subplots(figsize=(18, 6))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region)

ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
            insul_region['log2_insulation_score_'+str(windows[0])],
            label=f'Window {windows[0]} bp')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])

for res in windows[1:]:
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
f
plt.savefig("insulation_test.png", format='png',bbox_inches="tight")
plt.close("all")


f, ax = plt.subplots(figsize=(20, 10))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region)

ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]
ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
            weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
            strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])

plt.savefig("insulation_test.png", format='png',bbox_inches="tight")
plt.close("all")

histkwargs = dict(
    bins=10**np.linspace(-4,1,200),
    histtype='step',
    lw=2,
)

f, axs = plt.subplots(len(windows),1, sharex=True, figsize=(6,6), constrained_layout=True)
for i, (w, ax) in enumerate(zip(windows, axs)):
    ax.hist(
        insulation_table[f'boundary_strength_{w}'],
        **histkwargs
    )
    ax.text(0.02, 0.9,
             f'Window {w//1000}kb',
             ha='left',
             va='top',
             transform=ax.transAxes)
    ax.set(
        xscale='log',
        ylabel='# boundaries'
    )

axs[-1].set(xlabel='Boundary strength');
plt.savefig("insulation_boundary_test.png", format='png',bbox_inches="tight")
plt.close("all")


from skimage.filters import threshold_li, threshold_otsu

f, axs = plt.subplots(len(windows), 1, sharex=True, figsize=(6,6), constrained_layout=True)
thresholds_li = {}
thresholds_otsu = {}
for i, (w, ax) in enumerate(zip(windows, axs)):
    ax.hist(insulation_table[f'boundary_strength_{w}'],**histkwargs)
    thresholds_li[w] = threshold_li(insulation_table[f'boundary_strength_{w}'].dropna().values)
    thresholds_otsu[w] = threshold_otsu(insulation_table[f'boundary_strength_{w}'].dropna().values)
    n_boundaries_li = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_li[w]).sum()
    n_boundaries_otsu = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_otsu[w]).sum()
    ax.axvline(thresholds_li[w], c='green')
    ax.axvline(thresholds_otsu[w], c='magenta')
    ax.text(0.01, 0.9, f'Window {w//1000}kb', ha='left', va='top', transform=ax.transAxes) 
    ax.text(0.01, 0.7, f'{n_boundaries_otsu} boundaries (Otsu)', c='magenta', ha='left', va='top', transform=ax.transAxes)
    ax.text(0.01, 0.5, f'{n_boundaries_li} boundaries (Li)', c='green', ha='left', va='top', transform=ax.transAxes)
    ax.set(xscale='log', ylabel='# boundaries')
axs[-1].set(xlabel='Boundary strength')
plt.savefig("insulation_boundary_thresholding_test.png", format='png',bbox_inches="tight")
plt.close("all")

# Download CTCF bigwig file. The file size is 592 Mb, so the download might take a while:
ctcf_fc_file = cooltools.download_data("HFF_CTCF_fc", cache=True, data_dir="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta")

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
register(MulticoreParam(progressbar = T, workers = 20), default = T)
BiocParallel::bpparam()
setwd("/volumes/seq/projects/gccACT/mdamb231_ACTseq")

#act set
act <- runVarbin("/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)

colData(act)$info <- 'act'
#gccact set
gccact <- runVarbin("/volumes/seq/projects/gccACT/230306_mdamb231_test/cells",
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)
colData(gccact)$info <- 'gccact'
#merged_copykit <- cbind(act, gccact)
merged_copykit<-gccact

# Mark euploid cells if they exist
merged_copykit <- findAneuploidCells(merged_copykit)

# Mark low-quality cells for filtering
merged_copykit <- findOutliers(merged_copykit)

# Visualize cells labeled by filter and aneuploid status
pdf("outlier_qc.heatmap.pdf")
plotHeatmap(merged_copykit, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# Remove cells marked as low-quality and/or aneuploid from the copykit object
merged_copykit <- merged_copykit[,SummarizedExperiment::colData(merged_copykit)$outlier == FALSE]
merged_copykit <- merged_copykit[,SummarizedExperiment::colData(merged_copykit)$is_aneuploid == TRUE]


# kNN smooth profiles
merged_copykit <- knnSmooth(merged_copykit)
merged_copykit <- runUmap(merged_copykit)

k_clones<-findSuggestedK(merged_copykit)
# Create a umap embedding 

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
merged_copykit  <- findClusters(merged_copykit,k_subclones=21)#output from k_clones

pdf("subclone.umap.pdf")
plotUmap(merged_copykit, label = 'subclones')
dev.off()

#remove noise cluster
merged_copykit2 <- merged_copykit[,SummarizedExperiment::colData(merged_copykit)$subclones != "c0"]

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
merged_copykit2 <- calcConsensus(merged_copykit2)
merged_copykit2<- runConsensusPhylo(merged_copykit2)

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.pdf")
plotHeatmap(merged_copykit2, row_split="subclones",label = 'subclones',order='consensus_tree',n_threads=10,use_raster=T)
dev.off()

saveRDS(merged_copykit2,file="/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")
merged_copykit2<-readRDS("/volumes/seq/projects/gccACT/230306_mdamb231_test/scCNA.rds")

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

https://github.com/tanlongzhi/dip-c/#workflow

```bash
# align reads
bwa mem -5SP genome.fa R1.fq.gz R2.fq.gz | gzip > aln.sam.gz # for Nextera
#seqtk mergepe R1.fq.gz R2.fq.gz | pre-meta - | bwa mem -5SP -p genome.fa - | gzip > aln.sam.gz # for META

# extract segments
hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly -y - | gzip > contacts.seg.gz # for female
#hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly - | hickit.js bedflt par.bed - | gzip > contacts.seg.gz # for male

# resolve haplotypes via imputation
hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz
hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz

# generate 3D structures (with 3 replicates)
for rep in `seq 1 3`
do
  hickit -s${rep} -M -i impute.pairs.gz -Sr1m -c1 -r10m -c2 -b4m -b1m -O 1m.${rep}.3dg -b200k -O 200k.${rep}.3dg -D5 -b50k -O 50k.${rep}.3dg -D5 -b20k -O 20k.${rep}.3dg
done

# convert from hickit to dip-c formats, and remove repetitive regions from 3D structures
scripts/hickit_pairs_to_con.sh contacts.pairs.gz
scripts/hickit_impute_pairs_to_con.sh impute.pairs.gz
for rep in `seq 1 3`
do
  scripts/hickit_3dg_to_3dg_rescale_unit.sh 20k.${rep}.3dg
  dip-c clean3 -c impute.con.gz 20k.${rep}.dip-c.3dg > 20k.${rep}.clean.3dg # remove repetitive (contact-less) regions
done

# align replicate structures and calculate RMSD (overall value in .log file)
dip-c align -o aligned.20k. 20k.[1-3].clean.3dg 2> 20k.align.log > 20k.align.color

# convert to juicebox format for interactive viewing
# raw contacts
java -Xmx2g -jar juicer_tools.jar pre -n contacts.pairs.gz contacts.hic mm10
# haplotype-resolved contacts
scripts/con_imputed_to_juicer_pre_short.sh impute.con.gz
java -Xmx2g -jar juicer_tools.jar pre -n impute.juicer.txt.gz impute.hic color/mm10.chr.hom.len

# calculate single-cell chromatin compartment values along the genome
dip-c color2 -b1000000 -H -c color/mm10.cpg.1m.txt -s contacts.con.gz > cpg_b1m.color2 # contact-based
dip-c color -c color/mm10.cpg.20k.txt -s3 20k.1.clean.3dg > cpg_s3.color # 3D-structure-based

# calculate radial positioning
dip-c color -C 20k.1.clean.3dg > C.color

# color by chromosome number and visualize as mmCIF (viewable with pymol)
dip-c color -n color/mm10.chr.txt 20k.1.clean.3dg | dip-c vis -c /dev/stdin 20k.1.clean.3dg > 20k.1.clean.n.cif
```


To ADD:

<!--
Identification of cell cycling
https://github.com/shahcompbio/scdna_replication_tools

```bash
conda activate r4.2
```
```R
install.packages("GGally")
devtools::install_github("CL-CHEN-Lab/User_interface_for_Kronos_scRT", type = "source")

```
-->


<!--
### eccDNA Analysis
https://www.science.org/doi/10.1126/sciadv.aba2489
https://github.com/pk7zuva/Circle_finder/blob/master/circle_finder-pipeline-bwa-mem-samblaster.sh
It says read length must be at least 75bp if not enriched. But I'm going to try anyway.

```bash
#set up variables and directory
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
fq1="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.barc.fastq.gz"
fq2="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.barc.fastq.gz"

#Map reads with BWA Mem
cd $dir
bwa mem -t 20 $ref $fq1 $fq2  | samblaster -e -d samp.disc.sam -s samp.split.sam | samtools view -Sb - > samp.out.bam

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

#rm *hg38.sam *hg38.bam
```

```bash
cpu=50
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
fq1="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.barc.fastq.gz"
fq2="/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.barc.fastq.gz"
dir="/volumes/seq/projects/gccACT/230306_mdamb231_test"
minoverlap=10
sampname="230306_mdamb231_test_circfinder"

cd $dir
bash ~/src/circle_finder-pipeline-bwa-mem-samblaster.sh $cpu \
$ref \
$fq1 \
$fq2 \
$minoverlap \
$sampname hg38 &
```

-->

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