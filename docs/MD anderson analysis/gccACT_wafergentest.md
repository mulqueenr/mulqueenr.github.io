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

Use wafergen output to generate a list of cells
```bash
/volumes/lab/users/wet_lab/instruments/wafergen/Ryan/2023.06.12-137339
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
(Next time use bcl2fastq2 with proper sample sheet)

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
cd $dir
samtools view $dir/230701_wafergen_gccact.bam | awk -v dir=$dir 'OFS="\t" {split($1,a,":"); print $0,"XM:Z:"a[1] > "./cells/"a[1]".230701_wafergen_gccact.sam"}' &
 #split out bam to cell level sam
```

## Remove sam files that have less than 100000 reads

```bash
cd $dir/cells
wc -l *sam | sort -k1,1n - | head -n -1 | awk '{if($1<100000) print $2}' | xargs rm
#586 cells returned
```

## Add header to each sam and convert to bam and sort
Using parallel to save time
```bash
cd $dir/cells
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
  samtools sort -T . -n -o - $1 | samtools fixmate - -| samtools sort -T . -o - - | samtools markdup -s - ${1::-4}.rmdup.bam 2> ${1::-4}.rmdup.stats.txt
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

Generate a metadata file of wells and sample input based on wafergen output file
```python
import pandas as pd
import sys

well=sys.argv[1] 
outdir=sys.argv[2] 
#well="/volumes/lab/users/wet_lab/instruments/wafergen/Ryan/2023.06.12-137339/137339_WellList.TXT"
#well.txt is prefilted to only include information of wells chosen for dispensing.
#outdir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"

well=pd.read_csv(well,sep="\t")
well["idx_i7"]=[x.split("+")[0] for x in well["Barcode"]]
well["idx_i5"]=[x.split("+")[1] for x in well["Barcode"]]
well["idx_i5"]=[x[::-1].lower().replace("a","T").replace("t","A").replace("c","G").replace("g","C") for x in well["idx_i5"]]
well["idx"]=well["idx_i7"]+well["idx_i5"]
for i in well["Sample"].unique():
		out_df=well[well["Sample"]==i]["idx"]
		i_out=i.replace(" ","_").replace("-","_")
		out_df.to_csv(outdir+"/"+i_out+".barc.list.txt",index=False,header=False)

```

```bash
python ~/src/wafergen_metadatagenerator.py \
/volumes/lab/users/wet_lab/instruments/wafergen/Ryan/2023.06.12-137339/137339_WellList.TXT \
/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest

#arg1 is well list txt from wafergen run
#arg2 is run directory used in analysis
```

Move separate samples (wafergen defined) to different directories
This needs to be looked at again and fixed probably.
```bash
run_dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"

for i in ${run_dir}/*.barc.list.txt; do  
	out=$(basename $i)
	out_prefix=${out::-14}
	if [ ! -d "${run_dir}/cells/${out_prefix}" ]
		then 
			mkdir ${run_dir}/cells/${out_prefix}
	fi
	cat $i | find ${run_dir}/cells -type f -name {}.*.bam -exec sh -c 'mv "$0" "$1"/cells/"$2"' {} $run_dir $out_prefix
done

```

## Run CopyKit for WGS portion
Analysis from 
https://navinlabcode.github.io/CopyKit-UserGuide/quick-start.html
For this sample, cell labels were lost when I switched wafergen systems (whoops, but can still be identified by the WGS signals).
In the future I'll just use a single function per sample directory for processing.

```R
library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)
library(scquantum)
register(MulticoreParam(progressbar = T, workers = 5), default = T)
BiocParallel::bpparam()

setwd("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest")
tumor2 <- runVarbin("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/cells/sample",
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
tumor  <- findClusters(tumor,k_superclones=16, k_subclones=20)#output from k_clones

pdf("all_cells.subclone.umap.pdf")
plotUmap(tumor, label = 'subclones')
dev.off()

pdf("all_cells.superclone.umap.pdf")
plotUmap(tumor, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
tumor <- calcConsensus(tumor)
tumor <- runConsensusPhylo(tumor)
tumor <- runPhylo(tumor, metric = 'manhattan')

pdf("all_cells.subclone.phylo.pdf")
plotPhylo(tumor, label = 'subclones')
dev.off()

# Plot a copy number heatmap with clustering annotation
pdf("all_cells.subclone.heatmap.pdf")
plotHeatmap(tumor, label = c('superclones','subclones'),order='hclust')
dev.off()

saveRDS(tumor,file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/all_cells.scCNA.rds")
tumor<-readRDS("/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/all_cells.scCNA.rds")

#based on CNV profiles, superclone s1 is skbr3 and all others are mda-mb-231.
#going to split and rerun subclone analysis and heatmaping
# Remove cells marked as low-quality and/or aneuploid from the copykit object
skbr3<- tumor[,SummarizedExperiment::colData(tumor)$superclones == "s1"]
mdamb231 <- tumor[,SummarizedExperiment::colData(tumor)$superclones != "s1"]


# Create a umap embedding 
skbr3 <- runUmap(skbr3) ; mdamb231<- runUmap(mdamb231)
k_clones_skbr3<-findSuggestedK(skbr3); k_clones_mdamb231<-findSuggestedK(mdamb231); #16

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
skbr3 <- findClusters(skbr3,k_subclones=6)#output from k_clones
mdamb231 <- findClusters(mdamb231,k_subclones=13)#output from k_clones


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

saveRDS(skbr3,file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/skbr3.scCNA.rds")
saveRDS(mdamb231,file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/mdamb231.scCNA.rds")

skbr3 <- calcInteger(skbr3, method = 'scquantum', assay = 'smoothed_bincounts')
mdamb231 <- calcInteger(mdamb231, method = 'scquantum', assay = 'smoothed_bincounts')

pdf("skrb3_ploidy.score.pdf")
plotMetrics(skbr3, metric = 'ploidy', label = 'ploidy_score')
dev.off()

pdf("mdamb231_ploidy.score.pdf")
plotMetrics(mdamb231, metric = 'ploidy', label = 'ploidy_score')
dev.off()

pdf("skrb3_ploidy.heatmap.pdf")
plotHeatmap(skbr3, assay = 'integer', label = c("superclones", "subclones"), order_cells = 'consensus_tree')
dev.off()

pdf("mdamb231_ploidy.heatmap.pdf")
plotHeatmap(mdamb231, assay = 'integer', label = c("superclones", "subclones"), order_cells = 'consensus_tree')
dev.off()

saveRDS(skbr3,file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/skbr3.scCNA.rds")
saveRDS(mdamb231,file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/mdamb231.scCNA.rds")

skbr3<-readRDS(file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/skbr3.scCNA.rds")
mdamb231<-readRDS(file="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest/mdamb231.scCNA.rds")


clone_out<-data.frame(bam=paste0(row.names(skbr3@colData),".bam"),clone=skbr3@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("skbr3_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}

clone_out<-data.frame(bam=paste0(row.names(mdamb231@colData),".bam"),clone=mdamb231@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("mdamb231_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}
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


### Generation of HiC Contact Matrices
Merge bam files based on CopyKit output. Then using bam2pairs from pairix to generate contacts


```bash
conda activate EagleC
dir="/volumes/seq/projects/gccACT/230612_MDAMB231_SKBR3_Wafergentest"
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

#filter reads to HiC contacts, then perform bam2pairs 
mkdir $dir/contacts


bamlist_merge_to_pairs() {
	samtools merge -b $3 -O SAM -@ 20 - | awk '{if (sqrt(($9^2))>=1000 || $7 != "=") print $0}' | samtools view -bT $2 - > $1/contacts/${3::-13}.contacts.bam && wait;
	bam2pairs $1/contacts/${3::-13}.contacts.bam $1/contacts/${3::-13}
}
export -f bamlist_merge_to_pairs

bamlist_merge_to_pairs $dir $ref $dir/skbr3_allcells.txt


#set variable for bam list in function


```


#cooltools.virtual4c for eccDNA contacts
#https://cooltools.readthedocs.io/en/latest/cooltools.html#cooltools-api-virtual4c-module
