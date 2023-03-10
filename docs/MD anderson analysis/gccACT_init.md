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

# Processing Samples via Copykit to start
## Alignment
```bash
#set up variables and directory
mkdir /volumes/USR2/Ryan/projects/gccact/230306_mdamb231_test
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
dir="/volumes/USR2/Ryan/projects/gccact/230306_mdamb231_test"
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
dir="/volumes/USR2/Ryan/projects/gccact/230306_mdamb231_test"

add_header() {
	samtools view -bT $ref {$1::-4}.sam
	samtools sort -T $dir -o {}.bam -
}
export -f add_header
sam_in=`ls *sam`
parallel --jobs 30 sort_and_markdup ::: $sam_in
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

tumor <- runVarbin("/volumes/USR2/Ryan/projects/gccact/230306_mdamb231_test/cells",
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

saveRDS(tumor,file="/volumes/USR2/Ryan/projects/gccact/230306_mdamb231_test/scCNA.rds")
```

### eccDNA Analysis

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


```bash
bwa mem -t 10 ~/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
~/fastq/230124_hicACT_4cell/Undetermined_S0_L001_R1_001.fastq.gz | samtools view -b - > 230124_gccACT.bam &
```


