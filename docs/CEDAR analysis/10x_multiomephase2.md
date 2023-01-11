---
title: 10X Multiome Tumor All Together
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_phase2/
category: CEDAR
---

Multiome processing for 10X multiome data on Primary Tumors (Phase 1+2 and preliminary experiment combined).
Analysis was performed on Exacloud, OHSU's HPC which uses slurm as a job scheduler. So many parallelized analyses utilize slurm batch processing.

I also set up my environment to automatically paste figures to a slack channel, so you may notice many system calls like "slack -F [file] slack-channel", these are just a convience function for myself. 

*This code is a WIP and will be cleaned prior to manuscript finalization.*

## Initial Processing and QC


## Low Pass Whole Genome Sequencing Data
<!-- Done -->

## Align to hg38 with bwa-mem
<!-- Done -->

```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2
cat readme.txt 
#Run     Lane    Sample  I7 Index ID     Index1  I5 Index ID     Index2
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_A1-12        D712    AGCGATAG        D501    AGGCTATA
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_B1-12        D712    AGCGATAG        D502    GCCTCTAT
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_C1-12        D712    AGCGATAG        D503    AGGATAGG
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_D1-12        D712    AGCGATAG        D504    TCAGAGCC
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_E1-12        D712    AGCGATAG        D505    CTTCGCCT
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_F1-12        D712    AGCGATAG        D506    TAAGATTA
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_G1-12        D712    AGCGATAG        D507    ACGTCCTG
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_H1-12        D712    AGCGATAG        D508    GTCAGTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG01   N701    TAAGGCGA        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG03   N702    CGTACTAG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG04   N703    AGGCAGAA        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG05   N704    TCCTGAGC        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG06   N705    GGACTCCT        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG07   N706    TAGGCATG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG08   N707    CTCTCTAC        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG09   N710    CGAGGCTG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG10   N701    TAAGGCGA        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG11   N702    CGTACTAG        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG12   N703    AGGCAGAA        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG15   N704    TCCTGAGC        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG16   N705    GGACTCCT        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG19   N706    TAGGCATG        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG20   N707    CTCTCTAC        S506    TATGCAGT

#I'm taking the BCMM samples as the low pass whole genome for this project

cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM

for i in EXP220921HM_BCMM_WG*R1*fastq.gz; do line_count=`zcat $i | grep "^@" | wc -l`; echo $i $line_count & done & #get read count per samples


```

## Batch script for Alignment
<!-- Done -->

Using the bwa mem for alignment.

wgs_alignment.sbatch
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-15
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=30 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=5:00:00 ## ask for 3 hour on the node
#SBATCH --

fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
ref="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
file_in=`ls $fastq_dir/*WG*R1*.fastq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' `

#Align and output as bam file
R1=$file_in
R2=`echo $bam_in | awk '{ gsub("_R1_", "_R2_"); print $0}' `
output_name=${R1::-9}".batch.bam"
bwa mem -t 10 $ref $R1 $R2 | samtools sort -T . -@5 - | samtools view -b -@5 - > $output_name

```

```bash
sbatch wgs_alignment.sbatch
```

## Batch script for deduplication via Picard
<!-- Done -->

wgs_dedup.sbatch
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-15
#SBATCH --tasks-per-node=15 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 3 hour on the node
#SBATCH --

picard_dir="/home/groups/CEDAR/tools/picard-tools-1.119/"
fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
list_files=`ls $fastq_dir/*WG*R1*.bam`
bam_in=`ls $fastq_dir/*WG*R1*.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
i=$bam_in
output_name=${i::-4}".dedup.bam"
output_metrics=${i::-4}".dedup.metrics.txt"

#picard mark duplicates
java -jar ${picard_dir}/MarkDuplicates.jar \
        I=$i \
        O=$output_name \
        M=$output_metrics
```

```bash
sbatch wgs_dedup.sbatch
```

## Runing GATK4
<!-- Done -->

Following https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants

Download and set up a conda environment
Following 
```bash
#wget https://github.com/broadinstitute/gatk/archive/refs/tags/4.2.6.1.tar.gz
#tar -xvf 4.2.6.1.tar.gz
#cp gatkcondaenv.yml.template gatkcondaenv.yml
#conda env create -n gatk -f gatkcondaenv.yml
#using previously installed version
conda install -c bioconda gcnvkernel
#source activate gatk
```
### Process bam files into counts data
<!-- Done -->

Using 10Kbp windows and the cellranger genome for consistency.

```bash
source activate gatk
gatk="/home/groups/CEDAR/nishida/src/gatk-4.1.2.0/gatk"
ref_dir="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/"
ref="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
picard_dir="/home/groups/CEDAR/tools/picard-tools-1.119/"

#download mappability information
cd $ref_dir
wget https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bed.gz
gzip -d k50.umap.bed.gz
#Index file
$gatk IndexFeatureFile \
        -F k50.umap.bed

hg38_map="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/k50.umap.bed"


$gatk IndexFeatureFile \
        -L preprocessed.10000.interval_list \
        -R $ref \
        --mappability-track $hg38_map \
        -imr OVERLAPPING_ONLY \
        -O preprocessed.10000.annot.tsv

#Prepare reference file
  #idx
  cd $ref_dir
  samtools faidx genome.fa

  #dict
  cd $fastq_dir
  $gatk CreateSequenceDictionary -R $ref

#Prepare intervals
$gatk PreprocessIntervals \
    -R $ref \
    --bin-length 10000 \
    --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O preprocessed.10000.interval_list

#Annotate for downstream filtering
$gatk AnnotateIntervals \
        -L preprocessed.10000.interval_list \
        -R $ref \
        -imr OVERLAPPING_ONLY \
        --mappability-track $hg38_map \
        -O preprocessed.10000.annot.tsv


#Add read group for processing
for i in ${fastq_dir}/*WG*R1*.dedup.bam; do
  output_name=${i::-4}".RG.bam"
  RG_out=`basename $i | awk '{split($0,a,"_"); print a[4]}'`
  java -jar ${picard_dir}/AddOrReplaceReadGroups.jar \
        I=$i \
        O=$output_name\
        RGID=1 \
        RGLB=$RG_out \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=$RG_out &
  done &

#index bam files then perform bin counting
for i in ${fastq_dir}/*WG*R1*.dedup.RG.bam; do
  samtools index $i & done &

#note this can also be fed a bed file of intervals to perfectly match the single cell calling if needed
for i in ${fastq_dir}/*WG*R1*.dedup.RG.bam; do
  output_name=${i::-10}
   $gatk CollectReadCounts \
        -L preprocessed.10000.interval_list \
        -R $ref \
        -imr OVERLAPPING_ONLY \
        -I $i \
        -O ${output_name}.hdf5
    done &

#Filter intervals
$gatk FilterIntervals \
        -L preprocessed.10000.interval_list \
        --annotated-intervals preprocessed.10000.annot.tsv \
        -I EXP220921HM_BCMM_WG01_S9_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG06_S13_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG10_S17_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG16_S21_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG03_S10_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG07_S14_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG11_S18_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG19_S22_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG04_S11_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG08_S15_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG12_S19_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG20_S23_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG05_S12_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG09_S16_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG15_S20_L002_R1_001.batch.de.hdf5 \
        -imr OVERLAPPING_ONLY \
        -O preprocessed.10000.filtered.interval_list

#Set up ploidy prior (check that it is properly tab-separated)

echo """CONTIG_NAME PLOIDY_PRIOR_0  PLOIDY_PRIOR_1  PLOIDY_PRIOR_2  PLOIDY_PRIOR_3  PLOIDY_PRIOR_4  PLOIDY_PRIOR_5
chr1  0.01  0.01  0.95  0.01  0.01  0.01
chr2  0.01  0.01  0.95  0.01  0.01  0.01
chr3  0.01  0.01  0.95  0.01  0.01  0.01
chr4  0.01  0.01  0.95  0.01  0.01  0.01
chr5  0.01  0.01  0.95  0.01  0.01  0.01
chr6  0.01  0.01  0.95  0.01  0.01  0.01
chr7  0.01  0.01  0.95  0.01  0.01  0.01
chr8  0.01  0.01  0.95  0.01  0.01  0.01
chr9  0.01  0.01  0.95  0.01  0.01  0.01
chr10 0.01  0.01  0.95  0.01  0.01  0.01
chr11 0.01  0.01  0.95  0.01  0.01  0.01
chr12 0.01  0.01  0.95  0.01  0.01  0.01
chr13 0.01  0.01  0.95  0.01  0.01  0.01
chr14 0.01  0.01  0.95  0.01  0.01  0.01
chr15 0.01  0.01  0.95  0.01  0.01  0.01
chr16 0.01  0.01  0.95  0.01  0.01  0.01
chr17 0.01  0.01  0.95  0.01  0.01  0.01
chr18 0.01  0.01  0.95  0.01  0.01  0.01
chr19 0.01  0.01  0.95  0.01  0.01  0.01
chr20 0.01  0.01  0.95  0.01  0.01  0.01
chr21 0.01  0.01  0.95  0.01  0.01  0.01
chr22 0.01  0.01  0.95  0.01  0.01  0.01
chrX  0.01  0.01  0.95  0.01  0.01  0.01
chrY  1 0 0 0 0 0""" > ploidy_prior.tsv

#DetermineGermlineContigPloidy
$gatk DetermineGermlineContigPloidy \
        -L preprocessed.10000.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I EXP220921HM_BCMM_WG01_S9_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG06_S13_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG10_S17_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG16_S21_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG03_S10_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG07_S14_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG11_S18_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG19_S22_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG04_S11_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG08_S15_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG12_S19_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG20_S23_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG05_S12_L002_R1_001.batch.de.hdf5 -I EXP220921HM_BCMM_WG09_S16_L002_R1_001.batch.de.hdf5 \
        -I EXP220921HM_BCMM_WG15_S20_L002_R1_001.batch.de.hdf5 \
        --output . \
        --contig-ploidy-priors ploidy_prior.tsv \
        --output-prefix ploidy \
        --verbosity DEBUG


```

## Using HMMCopy as well
<!-- Done -->


Install HMMcopy utils
```bash
git clone https://github.com/shahcompbio/hmmcopy_utils
cd /home/groups/CEDAR/mulqueen/src/hmmcopy_utils
cmake .
make

ref_dir="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/"
ref="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
hmm_utils="/home/groups/CEDAR/mulqueen/src/hmmcopy_utils"
bowtie-build $ref ${ref::-3} #build bowtie reference index
${hmm_utils}/util/mappability/generateMap.pl -o ${ref::-3}.map.bw -i ${ref::-3} $ref #make mappability
#10kb windows
${hmm_utils}/bin/mapCounter -w 10000 ${ref::-3}.map.bw > ${ref::-3}.10000.map.wig #make windows
${hmm_utils}/bin/gcCounter -w 10000 ${ref} > ${ref::-3}.10000.gc.wig #make gc measure per window

```
Code from SCOPE and HMMCopy. Similar to s3wgs processing. 
```R

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM")

library(HMMcopy)

rfile <- system.file("/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.gc.wig") 
gfile <- system.file("/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.gc.wig") #gc content
mfile <- system.file("/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.map.bw") #mappability
normal_reads <- wigsToRangedData(rfile, gfile, mfile)


```

Using SCOPE WGSmapp and HMMcopy for analysis in R

```R
library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(parallel)
library(HMMcopy)
library(RColorBrewer)
#Initalization
bamfolder <- "/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
bamFile <- list.files(bamfolder, pattern = 'dedup.RG.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- paste("sample",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20),sep="_") #bams are ordered by sample number as well
setwd(bamfolder)


set_up_ref<-function(bins){ #modified version of SCOPE's get_bam_bed function
  genome <- BSgenome.Hsapiens.UCSC.hg38
  ref <- bins[which(as.character(seqnames(bins)) %in% paste0("chr", c(seq_len(22), "X", "Y")))] #autosomes and X Y

  #Compute mappability for each reference bin.
  mapp_gref<-mapp_hg38 #this is packaged with SCOPE, mappability across bins
  mapp <- rep(1, length(ref))
  #seqlevelsStyle(ref) <- "UCSC"
  for (chr in as.character(unique(seqnames(ref)))) {
      message("Getting mappability for ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- ref[which(as.character(seqnames(ref)) == chr)]
      mapp.chr <- rep(1, length(ref.chr))
      overlap <- as.matrix(findOverlaps(ref.chr, mapp_gref))
      for (i in unique(overlap[, 1])) {
          index.temp <- overlap[which(overlap[, 1] == i), 2]
          overlap.sub <- findOverlaps(ref.chr[i], mapp_gref[index.temp])
          overlap.intersect <- pintersect(ref.chr[i][queryHits(overlap.sub)],mapp_gref[index.temp][subjectHits(overlap.sub)])
          mapp.chr[i] <- sum((mapp_gref$score[index.temp]) * (width(overlap.intersect)))/sum(width(overlap.intersect))
      }
      mapp[chr.index] <- mapp.chr
  }

  #Compute GC for each bin, also from SCOPE
  gc <- rep(NA, length(ref))
  for (chr in unique(seqnames(ref))) {
      message("Getting GC content for chr ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- IRanges(start = start(ref)[chr.index], end = end(ref)[chr.index])
      if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
          chrtemp <- "chrX"
      }
      else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
          chrtemp <- "chrY"
      }
      else {
          chrtemp <- as.numeric(mapSeqlevels(as.character(chr), 
              "NCBI")[1])
      }
      if (length(chrtemp) == 0) 
      message("Chromosome cannot be found in NCBI database. ")
      chrm <- unmasked(genome[[chrtemp]])
      seqs <- Views(chrm, ref.chr)
      af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
      gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }

  ref@elementMetadata$gc<-gc
  ref@elementMetadata$mapp<-mapp
  return(ref)
}

get_sample_coverage<-function(bam_in="EXP220921HM_BCMM_WG01_S9_L002_R1_001.batch.dedup.RG.bam",ref,samp_name="sample_9"){
  sampname<-samp_name
    seg.dup <- read.table(system.file("extdata", "GRCh38GenomicSuperDup.tab", package = "WGSmapp"))
    gaps <- read.table(system.file("extdata", "hg38gaps.txt", package = "WGSmapp"))
    seg.dup <- seg.dup[!is.na(match(seg.dup[,1], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
    seg.dup <- GRanges(seqnames = seg.dup[,1], ranges = IRanges(start = seg.dup[,2], end = seg.dup[,3]))
    gaps <- gaps[!is.na(match(gaps[,2], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
    gaps <- GRanges(seqnames = gaps[,2], ranges = IRanges(start = gaps[,3], end = gaps[,4]))
    mask.ref <- sort(c(seg.dup, gaps))

    Y <- matrix(nrow = length(ref), ncol = length(sampname))
    rownames(Y) <- paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
    colnames(Y) <- sampname
    bamurl <- bam_in
    what <- c("rname", "pos", "mapq", "qwidth")
    flag <- scanBamFlag( isDuplicate = FALSE, isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE) # isFirstMateRead = TRUE #isPaired = TRUE,
    param <- ScanBamParam(what = what, flag = flag)
    bam <- scanBam(bamurl, param = param)[[1]]
    message("Getting coverage for sample ", ": ", sampname, "...", sep = "")
    
    bam.ref <- GRanges(seqnames = bam$rname, ranges = IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
    bam.ref <- bam.ref[bam$mapq >= 20] #Q20 threshold
    bam.ref <- suppressWarnings(bam.ref[countOverlaps(bam.ref, mask.ref) == 0])
    Y[, 1] <- countOverlaps(ref, bam.ref)
    return(Y)
}

genome <- BSgenome.Hsapiens.UCSC.hg38
bins <- tileGenome(seqinfo(genome), tilewidth = 1000 * 1000, cut.last.tile.in.chrom = TRUE) #set bins by other CNV callers
ref<-set_up_ref(bins=bins) #bins is granges of windows to use
Y<-lapply(1:length(bamFile),function(x) get_sample_coverage(bam_in=bamFile[x],ref=ref,samp_name=sampname_raw[x]))
Y<-do.call("cbind",Y)

#HMM Correction
hmmcopy_sample<-function(x){
  count<-cbind(as.data.frame(ref),Y[,x])
  samp<-sampname_raw[x]
  if(any(endsWith(samp,paste0("_",as.character(1:12))))){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/",samp,"/outs")
  } else {
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/",samp,"/outs")
  }
  colnames(count)<-c("chr","start","end","width","strand","gc","map","reads")
  count<-count[c("chr","start","end","reads","gc","map")]
  count$gc<-count$gc/100
  count<-data.table(count)
  count<-correctReadcount(count)
  count$chr<-as.character(count$chr)
  count<-count[count$chr!="chrY",]
  seg<-HMMsegment(count)
  count$state<-seg$state
  count$state<-as.character(count$state)
  count$chr<-factor(count$chr,levels=paste0("chr",c(1:22,"X")))
  count<-count[order(count$chr,count$start),]
  count$row_order<-1:nrow(count)
  cols = setNames(brewer.pal(n=6,name="RdBu"), nm = c("6","5","4","3","2","1")) # black, red, green, blue

  plt<-ggplot(count,aes(x=row_order,y=copy,color=as.character(state)))+
    scale_color_manual(values=cols)+
    geom_point(size=2.5,alpha=1)+
    ylab("")+
    xlab("")+
    ylim(-3,3)+
    facet_grid(~chr,space="free",scales="free_x")+
    theme_minimal()+
    theme(axis.text.y = element_text(size=30),
        axis.text.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.1,"lines"),
        strip.background = element_blank(), 
        legend.position="none",
        panel.border = element_rect(colour = "black", fill = NA,size=3))

  ggsave(plt,file=paste0(wd,"/",samp,"_HMMcopy.pdf"),width=2000,height=250,units="mm",limitsize=F)
  system(paste0("slack -F ",paste0(wd,"/",samp,"_HMMcopy.pdf")," ryan_todo"))
  saveRDS(count,file=paste0(wd,"/",samp,"_bulkWGS_HMMcopy.rds"))
}

hmm_y<-lapply(1:length(bamFile),function(x) hmmcopy_sample(x))


outdir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/cnv_comparison"

#Transfer hmmcopy to data location for travis
hmmcopy_save_as_tsv<-function(x){
  samp<-sampname_raw[x]
  if(any(endsWith(samp,paste0("_",as.character(1:12))))){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/",samp,"/outs")
  } else {
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/",samp,"/outs")
  }
  count<-readRDS(file=paste0(wd,"/",samp,"_bulkWGS_HMMcopy.rds"))
  write.table(count,file=paste0(outdir,"/",samp,"_Bulk_HMMcopy.tsv"),sep="\t",quote=F,col.names=T,row.names=F)
}
lapply(1:length(bamFile),function(x) hmmcopy_save_as_tsv(x))



```


### File Location
```bash
mkdir /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
sftp mulqueen@nix.ohsu.edu
#enter password
get -r /data/EXP220628HM
get -r /data/EXP220629HM

#and download WGS data
get -r /data/EXP220921HM
```


### Reference data
Using chipseq bed files held on cistrome for accessibility analysis.

```bash
/home/groups/CEDAR/mulqueen/ref/cistrome #stored reference files here, downloaded from site at .tar.gz file
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor_full_QC.txt #has information on each download peaks files
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor #has individual peak files
```

### Run cellranger-mkfastq
<!-- Done -->

Set up indexes used for mkfastq

```bash
echo """Lane,Sample,Index
*,sample_13,SI-NA-B8
*,sample_14,SI-NA-B7
*,sample_15,SI-NA-B6
*,sample_16,SI-NA-B5
*,sample_17,SI-NA-B4
*,sample_18,SI-NA-B3
*,sample_19,SI-NA-B2
*,sample_20,SI-NA-B1""" > multiome_atac_phase2.csv

echo """Lane,Sample,Index
1,sample_13,SI-TT-E3
1,sample_14,SI-TT-F3
1,sample_15,SI-TT-G3
1,sample_16,SI-TT-H3
1,sample_17,SI-TT-E4
1,sample_18,SI-TT-F4
1,sample_19,SI-TT-G4
1,sample_20,SI-TT-H4""" > multiome_rna_phase2.csv #other lane used for nextera libraries

```

Sample Sheet for RNA demultiplexing. Index 2 is workflow b

```
[Header]
EMFileVersion,4
 
[Reads]
28
90
 
[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
1,sample_13,sample_13,ACCAGACAAC,CCTAGTTCCT,phase2,SI-TT-E3
1,sample_14,sample_14,GAGAGGATAT,CCCATTTCAA,phase2,SI-TT-F3
1,sample_15,sample_15,ATGACGTCGC,ATCCTGACCT,phase2,SI-TT-G3
1,sample_16,sample_16,CCCGTTCTCG,CCAATCCGTC,phase2,SI-TT-H3
1,sample_17,sample_17,AACCACGCAT,TAACCTGAAT,phase2,SI-TT-E4
1,sample_18,sample_18,CCCACCACAA,AAGCGGAGGT,phase2,SI-TT-F4
1,sample_19,sample_19,GCGCTTATGG,CTAGCCAGGC,phase2,SI-TT-G4
1,sample_20,sample_20,AGTTTCCTGG,CTGTGTGGCA,phase2,SI-TT-H4
```
Run bcl2fastq
```bash
run_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2"
sample_sheet="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/multiome_rna_samplesheet.csv"
out_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_rna/HCVLCDRX2"

bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --output-dir=${out_dir}\
            --sample-sheet=${sample_sheet}
```
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2
Run Cellranger-arc
```bash
#conda install -c bih-cubi bcl2fastq2

cellranger-arc mkfastq --id=phase_2_atac \
                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220629HM/220713_A01058_0247_AHFJY3DRX2 \
                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/multiome_atac_phase2.csv \
                     --localcores=20 \
                     --localmem=80

#mkfastq is failing on RNA so I used bcl2fastq above
#cellranger-arc mkfastq --id=phase_2_rna \
#                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2 \
#                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/multiome_rna_phase2.csv \
#                     --localcores=20 \
#                     --barcode-mismatches=2 \
#                     --with-failed-reads \
#                     --lanes=1 \
#                     --localmem=80
```
### Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc.

### RM Libraries

```bash
for i in 13 14 15 16 17 18 19 20; do
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_atac/outs/fastq_path/HFJY3DRX2/sample_"""${i}""",sample_"""${i}""",Chromatin Accessibility
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_rna/HCVLCDRX2/phase2,sample_"""${i}""",Gene Expression""" > sample_${i}.csv ; done
```

### Run CellRanger-ARC
<!-- Done -->

```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2 #using this as analysis directory
```

Run Cellranger per sample

```bash          
for i in sample_13.csv sample_14.csv sample_15.csv sample_16.csv sample_17.csv sample_18.csv sample_19.csv sample_20.csv ; do
  outname=${i::-4};
  cellranger-arc count --id=${outname} \
   --reference=/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
   --libraries=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/${i} \
   --localcores=30 \
   --localmem=90 ; done &
   
#check web summaries
cd /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/
for i in 1 3 4 5 6 7 8 9 10 11 12; do
  cp ./sample_$i/outs/web_summary.html ./sample_$i/outs/$i.web_summary.html
  slack -F ./sample_$i/outs/$i.web_summary.html ryan_todo; done 
```

## Initial QC
<!-- Done -->

### Perform Scrublet on Data to Ensure Single-cells
<!-- Done -->

Code from tutorial here.[https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb]

Reading in h5 python matrices with code from https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb
```bash
pip install scrublet
```
Saving a python script as /home/groups/CEDAR/mulqueen/src/multiome_scrublet.py

```python
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import sys
import pandas as pd
np.random.seed(0)
x=sys.argv[1]

try:
    x = int(x)
    print(x)
except ValueError:
    # Handle the exception
    print('Passing to Preliminary Data:' + x)

if type(x) == int:
  if x < 13:
    input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_"+str(x)+"/outs"
    outname="sample_"+str(x)
  elif x >= 13:
    input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_"+str(x)+"/outs"
    outname="sample_"+str(x) 
else:
  input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/"+x+"/outs"
  outname=x


#Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.
counts_matrix = scipy.io.mmread(input_dir + '/filtered_feature_bc_matrix/matrix.mtx.gz').T.tocsc()
cellIDs=gzip.open(input_dir + '/filtered_feature_bc_matrix/barcodes.tsv.gz',"rb").read().split()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Preprocessing...
#Simulating doublets...
#Embedding transcriptomes using PCA...
#Calculating doublet scores...
#Automatically set threshold at doublet score = 0.07
#Detected doublet rate = 31.1%
#Estimated detectable doublet fraction = 60.5%
#Overall doublet rate:
#        Expected   = 5.0%
#        Estimated  = 51.4%
#Elapsed time: 785.7 seconds

df = pd.DataFrame({'cellid':cellIDs, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv(input_dir+'/'+outname+'.scrublet.tsv', index=False, sep="\t")
print("Done with sample: "+outname)
print("Saved output to: "+input_dir+'/'+outname+'.scrublet.tsv')
```

```bash
for i in 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 RM_1 RM_2 RM_3 RM_4; 
  do python /home/groups/CEDAR/mulqueen/src/multiome_scrublet.py ${i}; 
  done
```

### Use SoupX to remove ambient RNA
<!-- Done -->

```R
install.packages('SoupX')
```

```R
library(SoupX)

run_soupX_persample<-function(x,y=1){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  sc = load10X(wd)
  sc = autoEstCont(sc,tfidfMin=y) #1 is default
  out = adjustCounts(sc)
  saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
  print(paste("Finished:",outname))
}


lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,16,19,20,"RM_1","RM_2","RM_3","RM_4"),run_soupX_persample)

#14,15,17 and 18  failed to correct

#Sample 15 failed due to high homogeneity, autosupplying contamination fraction based on other samples
wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",15,"/outs")
outname<-paste0("sample_",15)
sc = load10X(wd)
sc = setContaminationFraction(sc, 0.1) #supplying contamination fraction manually based on values seen from other samples, using 10% 
out = adjustCounts(sc)
saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
print(paste("Finished:",outname))

```

## Seurat Generation and Processing
<!-- Done -->

### Seurat Object Generation for Samples
Performing seurat analysis following https://satijalab.org/signac/articles/pbmc_multiomic.html

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# set up sample loop to load the RNA and ATAC data, save to seurat object
setupseurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################

  setwd(wd)
  counts <- Read10X_h5("filtered_feature_bc_matrix.h5") #count data
  fragpath <- "atac_fragments.tsv.gz" #atac fragments
  metadata_cellranger<-read.csv("per_barcode_metrics.csv") #metadata
  row.names(metadata_cellranger)<-metadata_cellranger$barcode
  soupx_output<-readRDS("soupx_corrected_counts.rds") #load SoupX contamination corrected output
  scrublet_output<-read.table(paste0(outname,".scrublet.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
  #clean up scrublet output to add to metadata columns
  #just a hold over from a python output that I'm correcting.
  if(startsWith(scrublet_output$cellid[1],"b")){ 
  scrublet_output$cellID<-unlist(lapply(scrublet_output$cellid, function(x) substr(x,2,nchar(x))))}
  row.names(scrublet_output)<-scrublet_output$cellID
  scrublet_output<-scrublet_output[,c("doublet_scores","predicted_doublets")]

  # create a Seurat object containing the RNA data
  dat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  dat[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  #Create corrected RNA data and add to object
  dat[["SoupXRNA"]]<-CreateAssayObject(
    counts=soupx_output)

  #QC cells
  DefaultAssay(dat) <- "ATAC"
  dat <- NucleosomeSignal(dat)
  dat <- TSSEnrichment(dat)
  dat<-AddMetaData(dat,metadata=metadata_cellranger)
  dat<-AddMetaData(dat,metadata=scrublet_output)

  plt<-VlnPlot(
    object = dat,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
  ggsave(plt,file=paste0(outname,".qc.pdf"))
  system(paste0("slack -F ",outname,".qc.pdf ryan_todo"))
  saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))
}

#generate all seurat objects
lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),setupseurat)

```

Initial Merged Seurat Object from all phases

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".SeuratObject.rds"))
  dat$sample<-outname #set up sample metadata
  return(dat)}

out<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),merge_seurat)


dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"), project = "all_data")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.SeuratObject.rds")
```

### Call Peaks and Dimensionality Reduction
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")

dat<-readRDS("phase2.SeuratObject.rds")
dat
table(dat$sample)
#sample_1 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15 sample_16
#     3523     20000      5575     14071       186        30      1489       576
#sample_17 sample_18 sample_19 sample_20 sample_21 sample_22 sample_23 sample_24
#       19        37      2234       722      1851      1463       876       931
# sample_3  sample_4  sample_5  sample_6  sample_7  sample_8  sample_9
#     7698     15253      1453       666      2656       724      1451

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# call peaks using MACS2
DefaultAssay(dat)<-"ATAC"
peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")
#use this set of peaks for all samples

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

DefaultAssay(dat) <- "ATAC"
saveRDS(peaks,file="combined.peakset.rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = annotation
)

saveRDS(dat,file="phase2.SeuratObject.rds")
```

### Run Dim Reduction Per Sample
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

peaks <- readRDS(file="combined.peakset.rds")

#perform initial clustering, and remove scrublet detected doublets
single_sample_qc<-function(x){
 if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".SeuratObject.rds"))
  dat$sample<-x
  outname<-x
  }
# call peaks using MACS2
DefaultAssay(dat) <- "ATAC"

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = annotation
)

#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)
p1<-DimPlot(dat,reduction="rna_umap")+ggtitle("RNA UMAP")

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)
p2<-DimPlot(dat,reduction="atac_umap")+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted_doublets")+ggtitle("Multimodal UMAP Doublets")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-FeaturePlot(dat,reduction="multimodal_umap",features="doublet_scores")+ggtitle("Multimodal UMAP Scublet Scores")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file=paste0(wd,"/",outname,".umap.pdf"))
system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
table(dat$predicted_doublets)
if((sum(dat@meta.data$predicted_doublets=="True")/ncol(dat))<0.05){
cellids<-row.names(dat@meta.data[dat@meta.data$predicted_doublets=="False",])
}else{
cellids<-row.names(dat@meta.data[dat@meta.data$doublet_scores<quantile(dat@meta.data$doublet_scores,0.95),])
}
saveRDS(cellids,paste0(wd,"/",outname,".cellids.rds"))
dat<-subset(dat,cells=cellids)
saveRDS(dat,file=paste0(wd,"/",outname,".QC.SeuratObject.rds"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_qc)

#then rerun clustering now that they are filtered.
single_sample_qc2<-function(x){
 if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds"))
  dat$sample<-x
  outname<-x
  }
# call peaks using MACS2
DefaultAssay(dat) <- "ATAC"

#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")

p1<-DimPlot(dat,reduction="rna_umap",group.by="seurat_clusters")+ggtitle("RNA UMAP")
p2<-DimPlot(dat,reduction="atac_umap",group.by="seurat_clusters")+ggtitle("ATAC UMAP")
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP")

#Finally Plot results
plt<-(p1 | p2)/(p3)
ggsave(plt,file=paste0(wd,"/",outname,".umap2.pdf"))
system(paste0("slack -F ",paste0(wd,"/",outname,".umap2.pdf")," ryan_todo"))
saveRDS(dat,file=paste0(wd,"/",outname,".QC.SeuratObject.rds"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_qc2)
```


### Run cisTopic for ATAC Dimensionality Reduction
<!-- Done -->

Cistopic Per sample (Updated to include other directory folders)

```bash
nano /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/cistopic_per_sample.R
```
```R
#######CISTOPIC PROCESSING PER CELL LINE#################
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
library(RColorBrewer)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

single_sample_cistopic_generation<-function(x){
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds"))
  }
  cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
  saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))

  sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =atac_sub@meta.data)
  pdf(paste0(wd,"/",outname,"_model_selection.pdf"))
  par(mfrow=c(3,3))
  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,"_model_selection.pdf")," ryan_todo"))
  
  saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
  print("finshed running cistopic")

  #Add cell embeddings into seurat
  cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
  colnames(cell_embeddings)<-sub_cistopic_models@cell.names
  n_topics<-nrow(cell_embeddings)
  row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
  cell_embeddings<-as.data.frame(t(cell_embeddings))

  #Add feature loadings into seurat
  feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
  row.names(feature_loadings)<-paste0("topic_",1:n_topics)
  feature_loadings<-as.data.frame(t(feature_loadings))

  #combined cistopic results (cistopic loadings and umap with seurat object)
  cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
  print("Cistopic Loading into Seurat")
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  #plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(wd,"/",outname,".cistopic.umap.pdf"),width=10)
  print(plt1)
  #print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".cistopic.umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".QC.SeuratObject.rds"))
  }
single_sample_cistopic_generation(x=args[1])

#for i in 1 3 4 5 6 7 8 9 10 11 12 15 16 19 20 "RM_1" "RM_2" "RM_3" "RM_4"; do Rscript cistopic_per_sample.R $i; done &

```

## Public Datasets for Comparison 
<!-- Done -->

### Using Transfer Anchors for Cell identification.
<!-- Done -->

Using Swarbrick paper labels for transfer. https://pubmed.ncbi.nlm.nih.gov/34493872/

Download data
```bash
cd /home/groups/CEDAR/mulqueen/ref/swarbrick
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
tar -xvf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```
Make Seurat Object with Metadata

```R
library(Seurat)

setwd("/home/groups/CEDAR/mulqueen/ref/swarbrick")
counts<-ReadMtx(mtx="count_matrix_sparse.mtx",cells="count_matrix_barcodes.tsv",features="count_matrix_genes.tsv",feature.column=1) #sparse matrix of counts
metadata<-read.csv("metadata.csv") #metadata
row.names(metadata)<-metadata$X
# create a Seurat object containing the RNA adata
swarbrick <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)
swarbrick<-AddMetaData(swarbrick,metadata=metadata)
saveRDS(swarbrick,"/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
```

### Using EMBO paper for transfer of signatures as well. 
<!-- Done -->

https://doi.org/10.15252/embj.2020107333
Full code here: https://www.nature.com/articles/s41597-022-01236-2 data available here https://doi.org/10.6084/m9.figshare.17058077

Download data from GEO FTP server

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
wget https://figshare.com/ndownloader/articles/17058077/versions/1
unzip 1
```

Set up cell types by seurat cluster ID based on main figures.

```R
library(Seurat)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/ref/embo")
#match suerat clusters to assigned cell types in Fig 7C
##ER+ nonepi celltypes##
dat<-readRDS("SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
er_nonepi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("T cells","TAMs","CAFs","Pericytes","NA","Endothelial","TAMs_2","B cells","Myeloid","CAFs","Plasma cells","NA","NA"))
er_nonepi_cells<-setNames(names(er_nonepi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_nonepi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalSub.umap.pdf")
system("slack -F ERTotalSub.umap.pdf ryan_todo")
saveRDS(dat,file="SeuratObject_ERTotalSub.rds") #overwrite with cell types added to metadata

#match seurat clusters to assigned cell types in Fig EV4
dat<-readRDS("SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
er_nonepi_tcells<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("CD8+ effector","naive/resting","Treg","plasma","NK","NA"))
er_nonepi_tcells_cells<-setNames(names(er_nonepi_tcells[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_nonepi_tcells_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalTC.umap.pdf")
system("slack -F ERTotalTC.umap.pdf ryan_todo")
saveRDS(dat,file="SeuratObject_ERTotalTC.rds") #overwrite with cell types added to metadata


#match suerat clusters to assigned cell types in Fig 6E
dat<-readRDS("SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
er_epi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("epithelial","cycling epithelial","epithelial"))
er_epi_cells<-setNames(names(er_epi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_epi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalTum.umap.pdf")
system("slack -F ERTotalTum.umap.pdf ryan_todo")
saveRDS(dat,"SeuratObject_ERTotalTum.rds")

#ER+ All Cells
dat1<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
dat2<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
dat_tc<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
dat<-merge(dat1,dat2)
dat<-AddMetaData(dat,dat_tc$celltype,col.name="TCell_Subtype")
saveRDS(dat,"SeuratObject_ERProcessed.rds")
```

## Swarbrick Paper Label Transfer
<!-- Done -->

### Transfer Swarbrick cell types
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Using Label transfer to label cell types by Swarbrick paper
#seurat object made by AD
swarbrick<-readRDS("/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
swarbrick<-NormalizeData(swarbrick)
swarbrick<-FindVariableFeatures(swarbrick)
swarbrick<-ScaleData(swarbrick)

##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat@assays$peaks<-dat@assays$ATAC
  dat$sample<-paste0(x)
  }
  DefaultAssay(dat)<-"SoupXRNA"
  dat<-NormalizeData(dat)
  dat<-FindVariableFeatures(dat)
  dat<-ScaleData(dat)
  saveRDS(dat,file=file_in)

  transfer.anchors <- FindTransferAnchors(
    reference = swarbrick,
    reference.assay="RNA",
    query = dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = swarbrick$celltype_major,
  )

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)
  plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_label_transfer)
#
```

### Transfer EMBO Cell Types Per Sample
<!-- Done -->


```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Using Label transfer to label cell types by Embo Paper
#seurat object made by AD
embo_er<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERProcessed.rds")
DefaultAssay(embo_er)<-"RNA"
embo_er<-NormalizeData(embo_er)
embo_er<-FindVariableFeatures(embo_er)
embo_er<-ScaleData(embo_er)

##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat$sample<-paste0("sample_",x)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat$sample<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat@assays$peaks<-dat@assays$ATAC
    dat$sample<-paste0(x)
  }
  DefaultAssay(dat)<-"SoupXRNA"

  transfer.anchors <- FindTransferAnchors(
    reference = embo_er,
    reference.assay="RNA",
    query = dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = embo_er$celltype,
  )
  colnames(predictions)<-paste0("EMBO_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)
  plt1<-FeaturePlot(dat,features=c(                     
  "EMBO_prediction.score.Endothelial",       
  "EMBO_prediction.score.TAMs",              
  "EMBO_prediction.score.Pericytes",         
  "EMBO_prediction.score.CAFs",              
  "EMBO_prediction.score.T.cells",           
  "EMBO_prediction.score.Plasma.cells",      
  "EMBO_prediction.score.TAMs_2",            
  "EMBO_prediction.score.B.cells",           
  "EMBO_prediction.score.Myeloid",           
  "EMBO_prediction.score.epithelial",        
"EMBO_prediction.score.cycling.epithelial"),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by='EMBO_predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_label_transfer)


```


### Add sample metadata
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#from tumor sample information
meta_data_in<-as.data.frame(cbind("sample"=c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"),
  "diagnosis"= c("DCIS", "IDC", "IDC", "IDC", "IDC", "DCIS", "IDC", "IDC", "IDC", "IDC", "IDC", "NAT", "DCIS", "NAT", "IDC", "ILC", "IDC", "IDC", "NAT"), "molecular_type"=c(
"DCIS", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "DCIS", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "NA", "DCIS", "NA", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2+", "ER+/PR+/HER2-", "NA")))

sample_metadata_merged<-function(dat){
  dat_file_path=dat
  file_in=basename(dat)
  dir_in=dirname(dat)
  dat<-readRDS(dat) #read in as seurat object
  print(paste("Read in",dat_file_path))
  saveRDS(dat,paste0(dat_file_path,".backup")) #save a backup RDS file
  print("Made backup file")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$diagnosis,col.name="diagnosis")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$molecular_type,col.name="molecular_type")
  saveRDS(dat,dat_file_path)
  print("Finished Sample")
}

#run umap projections of merged samples
dat="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.QC.SeuratObject.rds"
sample_metadata_merged(dat)


sample_metadata_persample<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  print(paste("Read in",file_in))
  saveRDS(dat,paste0(file_in,".backup")) #save a backup RDS file
  print("Made backup file")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$diagnosis,col.name="diagnosis")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$molecular_type,col.name="molecular_type")
  saveRDS(dat,file_in)
  print("Finished Sample")
}
#add metadata to each sample
lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x) sample_metadata_persample(x))


```


## Create Seurat object with filtered cells
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".QC.SeuratObject.rds"))
  dat$sample<-outname #set up sample metadata
  return(dat)}

out<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),merge_seurat)


dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"), project = "all_data")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.QC.SeuratObject.rds")
```

### Add EMBO cell predictions to merged seurat object
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat_merged<-readRDS(file="phase2.QC.SeuratObject.rds")

# set up sample loop to load the RNA and ATAC data, save to seurat object
embo_metaextractor<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".QC.SeuratObject.rds"))
  dat_met<-dat@meta.data[startsWith(prefix="EMBO_prediction.",colnames(dat@meta.data))]
  row.names(dat_met)<-paste(outname,row.names(dat_met),sep="_")#set up sample metadata
  return(dat_met)}

out_met<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),embo_metaextractor)

met<-do.call("rbind",out_met)

dat_merged<-AddMetaData(dat_merged,met)
saveRDS(dat_merged,file="phase2.QC.SeuratObject.rds")
```


### Cistopic on merged samples
<!-- Done -->

Filter cells to the follow criteria
* scrublet reported predicted_doublets is false
* at least 500 ATAC fragments, no more than 10k 
* at least 100 gene expression exonic mapped UMIs, no more than 10k

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Read in and filter merged seurat object
atac_sub<-readRDS("phase2.QC.SeuratObject.rds")
atac_sub<-subset(atac_sub,cells= Reduce(intersect,
  list(which(atac_sub$predicted_doublets=="False"), 
  which(atac_sub$atac_fragments>= 250),
  which(atac_sub$atac_fragments <= 10000),
  which(atac_sub$gex_exonic_umis>=100),
  which(atac_sub$gex_exonic_umis <= 10000)))) #QC filters
saveRDS(atac_sub,file="phase2.QC.filt.SeuratObject.rds")
wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
outname<-paste0("phase2")

#Perform cistopic
cistopic_counts_frmt<-atac_sub$peaks@counts
row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
print("made cistopic object")
sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))

#Model selection
sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =atac_sub@meta.data)
pdf(paste0(wd,"/",outname,"_model_selection.pdf"))
par(mfrow=c(3,3))
sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
dev.off()
system(paste0("slack -F ",paste0(wd,"/",outname,"_model_selection.pdf")," ryan_todo"))

saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
print("finshed running cistopic")

#Add cell embeddings into seurat
cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
colnames(cell_embeddings)<-sub_cistopic_models@cell.names
n_topics<-nrow(cell_embeddings)
row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
cell_embeddings<-as.data.frame(t(cell_embeddings))

#Add feature loadings into seurat
feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
row.names(feature_loadings)<-paste0("topic_",1:n_topics)
feature_loadings<-as.data.frame(t(feature_loadings))

#combined cistopic results (cistopic loadings and umap with seurat object)
cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
atac_sub@reductions$cistopic<-cistopic_obj
n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("sample","seurat_clusters"))
plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
pdf(paste0(wd,"/",outname,".umap.pdf"),width=10)
print(plt1)
print(plt2)
dev.off()
system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
saveRDS(atac_sub,"phase2.QC.filt.SeuratObject.rds")
```


### Perform Merged Object Clustering
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")


#set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="#ed2024", "DCIS"="#bebebe","ILC"="#009444","NAT"="#aed8e6")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  sample_cols=c('sample_1'='#c7cfc5', 'sample_3'='#b4b2b2', 'sample_4'='#1e1f1d', 'sample_5'='#84205f', 'sample_6'='#a8d7b2', 'sample_7'='#7ecdc3', 'sample_8'='#242a26', 'sample_9'='#82717c', 'sample_10'='#146674', 'sample_11'='#88c25f', 'sample_12'='#708e3b', 'sample_15'='#2e3b80', 'sample_16'='#111114', 'sample_19'='#49823f', 'sample_20'='#b7dddb', 'RM_1'='#3e59a3', 'RM_2'='#b2433b', 'RM_3'='#12a0d6', 'RM_4'='#465243') 
  ########################################

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")


#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)


#Try multimodal with cistopic
# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(object = dat, 
  reduction.list = list("pca", "cistopic"), 
  dims.list = list(1:50, 1:ncol(dat@reductions$cistopic)), 
  modality.weight.name = "RNA.weight", 
  verbose = TRUE )
# build a joint UMAP visualization
dat <- RunUMAP(object = dat, 
  nn.name = "weighted.nn", 
  reduction.name="multimodal_umap", 
  assay = "SoupXRNA", 
  verbose = TRUE )

#plot cistopic umap too
alpha_val=0.33

p1<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="predicted.id",
  cols=alpha(type_cols,alpha_val))+
ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")

p2<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="sample",
  cols=alpha(sample_cols,alpha_val))+
ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")

p3<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="diagnosis",
  cols=alpha(diag_cols,alpha_val))+
ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")


#Finally Plot results
plt<-(p1/p2/p3)
ggsave(plt,file="phase2_filt.multimodal.umap.pdf",width=10,height=30)
system("slack -F phase2_filt.multimodal.umap.pdf ryan_todo")

#Also Plot EMBO designation of cell types
met<-as.data.frame(dat@meta.data)
embo_list<-  c("EMBO_prediction.score.Endothelial", "EMBO_prediction.score.NA","EMBO_prediction.score.TAMs","EMBO_prediction.score.Pericytes", "EMBO_prediction.score.CAFs","EMBO_prediction.score.T.cells", "EMBO_prediction.score.Plasma.cells","EMBO_prediction.score.TAMs_2","EMBO_prediction.score.B.cells", "EMBO_prediction.score.Myeloid", "EMBO_prediction.score.epithelial","EMBO_prediction.score.cycling.epithelial")
max_embo<-lapply(1:nrow(met),function(i) embo_list[which(met[i,embo_list]==max(met[i,embo_list],na.rm=T))])
max_embo<-unlist(lapply(1:length(max_embo),function(i) do.call("paste",as.list(max_embo[[i]]))))
max_embo<-unlist(lapply(max_embo,function(i) gsub("EMBO_","",i)))
max_embo<-unlist(lapply(max_embo,function(i) gsub("prediction.score.","",i)))
names(max_embo)<-row.names(met)
dat<-AddMetaData(dat,max_embo,col.name="EMBO_predicted.id")


p4<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="predicted.id")+
ggtitle("Swarbrick")


p5<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="EMBO_predicted.id")+
ggtitle("EMBO")

plt<-(p4/p5)
ggsave(plt,file="celltype_umap.pdf")
system("slack -F celltype_umap.pdf ryan_todo")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```


## Bar plots across cells
<!--Done-->
```R
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr) 
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
 

###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    
       
diag_cols<-c("IDC"="red", "DCIS"="grey")

molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
########################################



dat<-readRDS("phase2.QC.filt.SeuratObject.rds")


#Set up metadata and set up facet labels as factors for ordering
metadat<-as.data.frame(dat@meta.data)
metadat$diagnosis = factor(metadat$diagnosis, levels=c("NAT","DCIS","IDC","ILC"), labels=c("NAT","DCIS","IDC","ILC")) 
metadat$molecular_type = factor(metadat$molecular_type, levels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-"), labels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-")) 

#Cells PF
metadat$epi<-"Nonepi"
metadat[metadat$predicted.id %in% c("Cancer Epithelial","Normal Epithelial"),]$epi<-"Epi"
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,epi) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=epi,y=n))+geom_bar(stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free") #+ scale_y_continuous(trans='log10')
ggsave(plt1,file="barplot_qc_cellcount.pdf")
system("slack -F barplot_qc_cellcount.pdf ryan_todo")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_fill_manual(values=type_cols)+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="swarbrick_barplot_qc_celltype.pdf")
system("slack -F swarbrick_barplot_qc_celltype.pdf ryan_todo")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,EMBO_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=EMBO_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_fill_manual(values=embo_cell_cols)+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="Embo_barplot_qc_celltype.pdf")
system("slack -F Embo_barplot_qc_celltype.pdf ryan_todo")

```
### Vibe check on cell type prediction
<!-- Done -->

Plot out a heatmap of cell type scores per sample and prediction. I'm trying to figure out how specific they are and if results are concordant. Also plotting top 10 genes from Swarbrick gross cell type identification as a heatmap
Genes from Supplementary Table 9 (major classification).

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(seriation)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

dat_meta<-dat@meta.data


swarbrick_out<-as.data.frame(dat_meta %>% group_by(sample,predicted.id) %>% summarize(
  swarbrick_Myeloid=median(prediction.score.Myeloid,na.rm=T),
  swarbrick_B.cells=median(prediction.score.B.cells,na.rm=T),
  swarbrick_Plasmablasts=median(prediction.score.Plasmablasts,na.rm=T),
  swarbrick_T.cells=median(prediction.score.T.cells,na.rm=T),
  swarbrick_Normal.Epithelial=median(prediction.score.Normal.Epithelial,na.rm=T),
  swarbrick_Cancer.Epithelial=median(prediction.score.Cancer.Epithelial,na.rm=T),
  swarbrick_CAFs=median(prediction.score.CAFs,na.rm=T),
  swarbrick_PVL=median(prediction.score.PVL,na.rm=T),
  swarbrick_Endothelial=median(prediction.score.Endothelial,na.rm=T)
))

embo_out<-as.data.frame(dat_meta %>% group_by(sample,predicted.id) %>% summarize(
  EMBO_TAMs=median(EMBO_prediction.score.TAMs,na.rm=T),                            
  EMBO_TAMs_2=median(EMBO_prediction.score.TAMs_2,na.rm=T),            
  EMBO_B.cells=median(EMBO_prediction.score.B.cells,na.rm=T), 
  EMBO_Myeloid=median(EMBO_prediction.score.Myeloid,na.rm=T),           
  EMBO_Plasma.cells=median(EMBO_prediction.score.Plasma.cells,na.rm=T),
  EMBO_T.cells=median(EMBO_prediction.score.T.cells,na.rm=T),          
  EMBO_Epithelial=median(EMBO_prediction.score.epithelial,na.rm=T), 
  EMBO_cycling.epithelial=median(EMBO_prediction.score.cycling.epithelial,na.rm=T),       
  EMBO_CAFs=median(EMBO_prediction.score.CAFs,na.rm=T),   
  EMBO_Pericytes=median(EMBO_prediction.score.Pericytes,na.rm=T),                    
  EMBO_Endothelial=median(EMBO_prediction.score.Endothelial,na.rm=T),
))     


row.names(swarbrick_out)<-paste(swarbrick_out$sample,swarbrick_out$predicted.id)
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", 
#immune
"B-cells" ="#089099", "T-cells" ="#003147", 
#other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

side_ha<-rowAnnotation(df= data.frame(celltype=swarbrick_out$predicted.id, sample=swarbrick_out$sample),
                col=list(
                    celltype=setNames(type_cols,names(type_cols)),
                    cluster=setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique(swarbrick_out$sample))),unique(swarbrick_out$sample))
                    ))
#seriate order for nice lookin heatmap
o = seriate(swarbrick_out[,3:ncol(swarbrick_out)], method = "BEA_TSP")

SWARBRICK_SUB<-swarbrick_out[,3:ncol(swarbrick_out)]
swarbrick<-Heatmap(SWARBRICK_SUB,
  left_annotation=side_ha,
  row_order = get_order(o, 1), column_order=1:ncol(SWARBRICK_SUB),
  col=colorRamp2(c(0, max(SWARBRICK_SUB)), c("white", "blue")),
  show_row_names=T)

EMBO_SUB<-embo_out[,3:ncol(embo_out)]
EMBO<-Heatmap(EMBO_SUB,
  column_order = 1:ncol(EMBO_SUB),
  col=colorRamp2(c(0, max(EMBO_SUB,na.rm=T)), c("white", "red")),
  row_order=get_order(o,1))

pdf("predictions.heatmap.pdf",width=30)
swarbrick+EMBO
dev.off()
system("slack -F predictions.heatmap.pdf ryan_todo")


table(dat$sample)
#     RM_1      RM_2      RM_3      RM_4  sample_1 sample_10 sample_11 sample_12
#     1845      1461       869       922      3518     18789      5575     13303
#sample_15 sample_16 sample_19 sample_20  sample_3  sample_4  sample_5  sample_6
#     1486       574      2116       718      7698     15250      1444       661
# sample_7  sample_8  sample_9
#     2656       717      1445

table(dat[dat$predicted.id%in%c("Cancer Epithelial","Normal Epithelial"),]$sample)
#     RM_1      RM_2      RM_3      RM_4  sample_1 sample_10 sample_11 sample_12
#     1367       129       313       169        53      7023      4678     11791
#sample_15 sample_16 sample_19 sample_20  sample_3  sample_4  sample_5  sample_6
#      704       356      1648       284      4466      2768       792       448
# sample_7  sample_8  sample_9
#     1860       415      1096

```

## Run ChromVAR on all data
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)
library(BiocParallel)
register(SerialParam()) #using single core mode

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")
DefaultAssay(dat)<-"ATAC"
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(dat[["ATAC"]]))) %in% main.chroms)
dat[["ATAC"]] <- subset(dat[["ATAC"]], features = rownames(dat[["ATAC"]][keep.peaks]))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
peaks<-granges(dat[["ATAC"]])

motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, 
  pwm = pfm, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  use.counts = FALSE)

motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, 
  pwm = pfm)

dat <- SetAssayData(object = dat, 
  assay = 'ATAC', 
  slot = 'motifs', 
  new.data = motif.hg38)

dat <- RegionStats(object = dat, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

dat <- RunChromVAR( object = dat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```
## Using Signac Gene Activity Function
<!-- Done -->
This is to generate enhancer promoter linkages at genes (by proximity). Running on all data.

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(SeuratObjects)
library(EnsDb.Hsapiens.v86)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")
gene_activity<-GeneActivity(dat,process_n=10000)
saveRDS(gene_activity,file="phase2.QC.GeneActivity.rds")

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```


# Determine Tumor Cells and Clones via CNV Callers
<!-- Done -->

## InferCNV on RNA Profiles
<!-- Done -->

Immune and stromal cells were used to 
define the reference cell-inferred copy number profiles. Similar to analysis done in https://www.nature.com/articles/s41588-021-00911-1

infercnv_per_sample.R
```R
####Run InferCNV
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
args = commandArgs(trailingOnly=TRUE)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

####RUNNING INFERCNV#####
infercnv_per_sample<-function(x,prediction="EMBO"){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
    if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    #dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    #dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
  }
  dat<-readRDS("phase2.QC.filt.SeuratObject.rds") #use QC controlled bulk seurat object as input
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname

  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this
  dat$cnv_ref<-"FALSE"
  if(prediction=="EMBO"){
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("epithelial")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    }else{
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells","CAFs"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL","CAFs"))
  } 

  #write out gene order list
  gene_order<-annotation[!duplicated(annotation$gene_name),]
  gene_order<-as.data.frame(gene_order[gene_order$gene_name %in% row.names(dat),])
  gene_order<-gene_order[c("gene_name","seqnames","start","end")]
  chrorder<-paste0("chr",c(1:22,"X","Y","M"))
  gene_order$seqnames<-factor(gene_order$seqnames,levels=chrorder) # set chr order
  gene_order<-with(gene_order, gene_order[order(seqnames, start),]) #order by chr and start position
  write.table(gene_order,file="inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)
  gene_order<-read.table("inferCNV.gene_order.txt")

  counts=as.matrix(dat@assays$RNA@counts[,colnames(dat)])
  write.table(counts,file=paste0(wd,"/",outname,"_inferCNV.counts.txt"),sep="\t",col.names=T,row.names=T,quote=F)
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["cnv_ref"]))
  write.table(cell_annotation,file=paste0(wd,"/",outname,"_inferCNV.annotation.txt"),sep="\t",col.names=F,row.names=F,quote=F)

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(wd,"/",outname,"_inferCNV.counts.txt"),
                                      annotations_file=paste0(wd,"/",outname,"_inferCNV.annotation.txt"),
                                      delim="\t",
                                      gene_order_file="inferCNV.gene_order.txt",
                                      ref_group_names="TRUE")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0(wd,"/",outname,"_inferCNV"), 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               HMM_report_by="cell",
                               resume_mode=T,
                               HMM_type='i3',
                               num_threads=30)
  saveRDS(infercnv_obj,paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
  system(paste0("slack -F ",wd,"/",outname,"_inferCNV","/","infercnv.png"," -T ","\"",outname,"\"" ," ryan_todo") )
  system(paste0("slack -F ",wd,"/",outname,"_inferCNV","/","infercnv.19_HMM_predHMMi3.hmm_mode-samples.Pnorm_0.5.repr_intensities.png"," -T ","\"",outname,"\"" ," ryan_todo") )

}

infercnv_per_sample(x=as.character(args[1]),prediction="EMBO")

#lapply(c(4,10),infercnv_per_sample)

```


### Batch script for InferCNV Per Sample Processing
<!-- Done -->

Calling infercnv_per_sample.R script written above

infercnv_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-18
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=30 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --qos=long_jobs
#SBATCH --time=120:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=("1" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "15" "16" "19" "20" "RM_1" "RM_2" "RM_3" "RM_4")
sample_in=${array_in[$SLURM_ARRAY_TASK_ID]}
multiome_dir="/home/groups/CEDAR/mulqueen/projects/multiome"

srun Rscript ${multiome_dir}/infercnv_per_sample.R $sample_in


```

### Job submit all infercnv processing runs.
<!-- Done -->

```bash
sbatch infercnv_slurm.sh
```

## Run CaSpER on RNA profiles
<!-- Done -->

casper_per_sample.R
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
library(CaSpER) 
args = commandArgs(trailingOnly=TRUE)

casper_per_sample<-function(x,prediction="EMBO"){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  dat<-readRDS("phase2.QC.filt.SeuratObject.rds") #use QC controlled bulk seurat object as input
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname

  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  system(paste0("mkdir ",dir_in,"/casper"))
  bam_location<-paste0(dir_in,"/gex_possorted_bam.bam")
  BAFExtract_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract"
  hg38_list_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list" #downloaded from https://github.com/akdess/BAFExtract
  hg38_folder_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/"
  baf_sample_directory<-paste0(dir_in,"/casper")

  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  if(prediction=="EMBO"){
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("epithelial")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    }else{
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells","CAFs"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL","CAFs"))
  } 

  control<-names(dat$cnv_ref == "TRUE") #pulling this from the inferCNV function
  log.ge <- as.matrix(dat@assays$RNA@data)
  genes <- rownames(log.ge)
  annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
  log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
  rownames(log.ge) <- annotation$Gene
  log.ge <- log2(log.ge +1)

  system(paste0("samtools view ",bam_location," | ",BAFExtract_location," -generate_compressed_pileup_per_SAM stdin ",hg38_list_location," ",baf_sample_directory," 30 0 && wait;")) #generate BAF calls
  #example of actual call: samtools view /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs/gex_possorted_bam.bam| /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper 30 0 &
  system(paste0(BAFExtract_location," -get_SNVs_per_pileup ",hg38_list_location," ",baf_sample_directory," ",hg38_folder_location," 1 1 0.1 ",baf_sample_directory,"/test.snp")) #generage snv files from BAF
  #example of actual call: /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -get_SNVs_per_pileup /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/ 1 1 0.1 /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper/test.snp
  loh <- readBAFExtractOutput ( path=baf_sample_directory, sequencing.type="bulk") 
  names(loh) <- gsub(".snp", "", names(loh))
  load(paste0(hg38_folder_location,"/maf.rda")) ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
  loh<- list()
  loh[[1]] <- maf
  names(loh) <- sample_name
  loh.name.mapping <- data.frame (loh.name= sample_name , sample.name=colnames(log.ge))

  #analysis demonstration: https://rpubs.com/akdes/673120
  object <- CreateCasperObject(raw.data=log.ge,
    loh.name.mapping=loh.name.mapping, 
    sequencing.type="single-cell", 
    cnv.scale=3, 
    loh.scale=3, 
    expr.cutoff=0.1, 
    filter="median", 
    matrix.type="normalized",
    annotation=annotation, 
    method="iterative", 
    loh=loh, 
    control.sample.ids=control, 
    cytoband=cytoband)

  saveRDS(object,paste0(dir_in,"/casper/",sample_name,".initialobj.rds"))
  #object<-readRDS(paste0(dir_in,"/casper/",sample_name,".initialobj.rds"))
  ## runCaSpER
  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")
  saveRDS(final.objects,paste0(dir_in,"/casper/",sample_name,".finalobj.rds"))

  ## summarize large scale events 
  finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75)
  final.obj <- final.objects[[9]]
  saveRDS(final.obj,paste0(dir_in,"/casper/",sample_name,".finalobj.rds"))
  saveRDS(finalChrMat,paste0(dir_in,"/casper/",sample_name,".finalchrmat.rds"))
  #final.obj<-readRDS(paste0(dir_in,"/casper/",sample_name,".finalobj.rds"))

  #Segmentations
  gamma <- 6
  all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
  segment.summary <- extractSegmentSummary(final.objects)
  loss <- segment.summary$all.summary.loss
  gain <- segment.summary$all.summary.gain
  loh <- segment.summary$all.summary.loh
  loss.final <- loss[loss$count>gamma, ]
  gain.final <- gain[gain$count>gamma, ]
  loh.final <- loh[loh$count>gamma, ]

  #summrize segmentation across genes
  all.summary<- rbind(loss.final, gain.final)
  colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
  rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
  ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
  hits <- findOverlaps(rna, ann.gr)
  genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
  genes.ann <- lapply(genes, function(x) x[!(x=="")])
  all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
  all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
  rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann) #just need to fix genes.ann
  saveRDS(rna.matrix, paste0(dir_in,"/casper/",sample_name,".finalgenemat.rds"))

}

casper_per_sample(x=as.character(args[1]))

#lapply(c(7,9,11,15,16,19,"RM_2","RM_3"), function(x) casper_per_sample(x))

```


### Batch script for Casper Per Sample Processing
<!-- Done -->

Calling casper_per_sample.R script written above

casper_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-18
#SBATCH --tasks-per-node=5 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=5 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=("1" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "15" "16" "19" "20" "RM_1" "RM_2" "RM_3" "RM_4")
sample_in=${array_in[$SLURM_ARRAY_TASK_ID]}
multiome_dir="/home/groups/CEDAR/mulqueen/projects/multiome"

srun Rscript ${multiome_dir}/casper_per_sample.R $sample_in

```

### Job submit all casper processing runs.
<!-- Done -->

```bash
sbatch casper_slurm.sh
```

## Run CopyKat on RNA profiles
<!-- Done -->

https://github.com/navinlabcode/copykat

copykat_per_sample.sh
```R
#library(devtools)
#install_github("navinlabcode/copykat")

library(Signac)
library(Seurat)
library(copykat)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
args = commandArgs(trailingOnly=TRUE)

copykat_per_sample<-function(x,prediction="EMBO"){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  dat<-readRDS("phase2.QC.filt.SeuratObject.rds") #use QC controlled bulk seurat object as input
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname

  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  system(paste0("mkdir ",dir_in,"/copykat"))
  exp.rawdata <- as.matrix(dat@assays$RNA@counts)

  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  if(prediction=="EMBO"){
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("epithelial")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    }else{
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells","CAFs"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL","CAFs"))
  } 
  cnv_ref<-row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",])
  copykat_out <- copykat(rawmat=exp.rawdata, KS.cut=0.15,LOW.DR=0.05,UP.DR=0.2,id.type="S", ngene.chr=0, win.size=25, sam.name=sample_name, distance="euclidean", norm.cell.names=cnv_ref,output.seg="FALSE", plot.genes="FALSE", genome="hg20",n.cores=10)
  saveRDS(copykat_out,paste0(dir_in,"/copykat/",sample_name,".copykat.RDS"))
}

copykat_per_sample(x=as.character(args[1]))
lapply(c(10,11,12,19,"RM_1","RM_2"),function(x) copykat_per_sample(x))

#lapply(c(1,3,4,10,11,12,19,"RM_1","RM_2"),function(x) copykat_per_sample(x))
#to set CNV discrete changes, as per correspondence suggetions with Ruli Gao, 1.5x SD threshold, 1.5 absolute distance, or use +/-0.25 as cutoff
```


### Batch script for Copykat Per Sample Processing
<!-- Done -->

Calling copykat_per_sample.R script written above

copykat_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-18
#SBATCH --tasks-per-node=5 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=5 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=20gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=("1" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "15" "16" "19" "20" "RM_1" "RM_2" "RM_3" "RM_4")
sample_in=${array_in[$SLURM_ARRAY_TASK_ID]}
multiome_dir="/home/groups/CEDAR/mulqueen/projects/multiome"

srun Rscript ${multiome_dir}/copykat_per_sample.R $sample_in

```

## CopyscAT for ATAC CNV Calling 
<!-- Done -->

Using scATAC calling algorithm copyscAT from git repo https://github.com/spcdot/CopyscAT/

### Installation...
<!-- Done -->

```R
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install_github("spcdot/copyscat")
library(CopyscAT)
```
Modifies the copyscAT python script (https://github.com/spcdot/CopyscAT/blob/master/process_fragment_file.py) to filter based on a metadata table rather than read count (since I already QC cells) then posted to a subdirectory

```bash
mkdir /home/groups/CEDAR/mulqueen/ref/copyscat
```

### Now Running samples
<!-- Done -->

Code from https://github.com/spcdot/CopyscAT/blob/master/copyscat_tutorial.R
Initialize reference genome information for CopyscAT.

```R
library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Generate tile references
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38" ,tileWidth = 1e6,outputDir = "/home/groups/CEDAR/mulqueen/ref/copyscat")

##### REGULAR WORKFLOW #####

```

copyscat_per_sample.R
<!-- Done -->

```R
library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
args = commandArgs(trailingOnly=TRUE)


#initialize the environment
initialiseEnvironment(genomeFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
                      cytobandFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=500,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#Set up copyscAT Loop per sample
copyscAT_per_sample<-function(x,prediction="EMBO",knn_in=FALSE,cores=1){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  if (knn_in==TRUE){
  knn_list<-read.table(paste0("/home/groups/CEDAR/scATACcnv/Hisham_data/bed_files/WGS_eval/knn/",sample_name,"_knn5_neighbors.csv"),
    sep=",",header=T)
  knn_list<-as.data.frame(apply(knn_list, 2, function(y) gsub("[.]", "-", y)))
  #knn list is csv format <rowid><cell><neighbor1><neighbor2><neighbor3><neighbor4>
  }
  dat<-readRDS("phase2.QC.filt.SeuratObject.rds") #use QC controlled bulk seurat object as input
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  if (knn_in == TRUE){
  system(paste0("mkdir ",dir_in,"/copyscat_knn"))
  } else {
  system(paste0("mkdir ",dir_in,"/copyscat"))
  }
  #do python script preprocessing (basically just count fragments per window per cell)
  system(paste0("python /home/groups/CEDAR/mulqueen/ref/copyscat/process_fragment_file.py ",
  " -i ",dir_in,"/atac_fragments.tsv.gz",
  " -o ",dir_in,"/copyscat/copyscat.1mb.tsv",
  " -b ","1000000",
  " -f ","500",
  " -g ","/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
  " -c ",dir_in,"/metadata.tsv")) #modification takes in metadata table to filter cells by name, ignores -f flag
  if (knn_in==TRUE){
  setOutputFile(paste0(dir_in,"/copyscat_knn"),"copyscat_out_knn")
  } else {
  setOutputFile(paste0(dir_in,"/copyscat"),"copyscat_out")
  }

  #PART 1: INITIAL DATA NORMALIZATION
  scData<-readInputTable(paste0(dir_in,"/copyscat/copyscat.1mb.tsv"))
  #here is an if else, one python script also accounts for metacell merged cells, other is strictly single cell
  if(knn_in==TRUE){
    scData2<-as.data.frame(do.call("rbind",mclapply(1:nrow(knn_list), function(x){colSums(scData[row.names(scData) %in% unlist(knn_list[x,]),])},mc.cores=cores)))
    row.names(scData2)<-knn_list$cell
    scData<-scData2
  }

  #collapse into chromosome arm level
  summaryFunction<-cutAverage
  scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)
  scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 

  #PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
  #ALTERNATE METHOD FOR CNV CALLING (with normal cells as background)
  #Using same normal cell selection as used for CASPER and InferCNV
  dat$cnv_ref<-"FALSE"
  if(prediction=="EMBO"){
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("epithelial")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    }else{
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells","CAFs"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL","CAFs"))
  } 
  control<-names(dat$cnv_ref == "TRUE") #pulling this from the inferCNV function

  #compute central tendencies based on normal cells only
  colnames(scData_collapse)<-gsub(outname,"",colnames(scData_collapse))
  control<-gsub(paste0(outname,"_"),"",control)
  control <- control[control %in% colnames(scData_collapse)] #filter control list to control cells that survived filter
  median_iqr <- computeCenters(scData_collapse %>% select(chrom,control),summaryFunction=summaryFunction)
  #setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes in as normalCells can improve accuracy of CNV calls
  candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,
    useDummyCells = FALSE,
    propDummy=0.25,
    minMix=0.01,
    deltaMean = 0.03,
    deltaBIC2 = 0.25,
    bicMinimum = 0.1,
    subsetSize=50,
    fakeCellSD = 0.09,
    uncertaintyCutoff = 0.65,
    summaryFunction=summaryFunction,
    maxClust = 4,
    mergeCutoff = 3,
    IQRCutoff = 0.25,
    medianQuantileCutoff = -1,
    normalCells=control) 
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)

  if(knn_in==TRUE){
  saveRDS(candidate_cnvs_clean,file=paste0(dir_in,"/copyscat_knn/",sample_name,"copyscat_cnvs_matrix_knn.rds"))
  }else{saveRDS(candidate_cnvs_clean,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs_matrix.rds"))}

  #to save this data you can use annotateCNV4 as per usual, using normal barcodes
  final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, expectedNormals=control, saveOutput=TRUE,
    outputSuffix = "clean_cnv_b2",sdCNV = 0.6,filterResults=FALSE,filterRange=0.4,minAlteredCellProp = 0.5)

  if(knn_in==TRUE){
  saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat_knn/",sample_name,"copyscat_cnvs_knn.rds"))
  }else{saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs.rds"))}

  print(paste("Finished sample",sample_name))
}

copyscAT_per_sample(x=as.character(args[1]),knn=FALSE)
copyscAT_per_sample(x=as.character(args[1]),knn=TRUE)

#lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),copyscAT_per_sample)
#lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),function(x) copyscAT_per_sample(x,knn_in=TRUE,cores=5))
#copyscat_dat<-readRDS(file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs_matrix.rds"))

```

### Batch script for copyscAT Per Sample Processing
<!-- Done -->

Calling copyscat_per_sample.R script written above

copyscat_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-18
#SBATCH --tasks-per-node=5 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=5 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=("1" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "15" "16" "19" "20" "RM_1" "RM_2" "RM_3" "RM_4")
sample_in=${array_in[$SLURM_ARRAY_TASK_ID]}
multiome_dir="/home/groups/CEDAR/mulqueen/projects/multiome"


srun Rscript ${multiome_dir}/copyscat_per_sample.R $sample_in

```

### Job submit all copyscAT processing runs.
<!-- Done -->

```bash
sbatch copyscat_slurm.sh
```
## HMMcopy Bulk Comparison across single-cell CNV Callers
<!--Rerun-->
HMMcopy comparison across CNV Callers and low-pass Whole genome data

hmmcopy_comparisons.R

```R
#before running R increase slave limit 
#ulimit -s 32000 # enlarge stack limit to 32 megs
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(CaSpER) 
library(SCOPE)
library(WGSmapp)
library(doParallel)
library(reshape2)
library(parallel)
library(HMMcopy)
library(RColorBrewer)
library(philentropy)
library(dendextend)
library(ggalluvial)
args = commandArgs(trailingOnly=TRUE)

x=args[1]

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
pam50_colors<-c("Basal"="red","Her2"="pink","LumA"="blue","LumB"="cyan","Normal"="grey","NA"="black")
embo_colors<-c("Basal"="green","LP"="blue","ML"="orange","Str"="red","NA"="black")
########################################


getmode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

set_up_ref<-function(bins,ref_outname,save_rds=FALSE){ #modified version of SCOPE's get_bam_bed function
  genome <- BSgenome.Hsapiens.UCSC.hg38
  ref <- bins[which(as.character(seqnames(bins)) %in% paste0("chr", c(seq_len(22), "X", "Y")))] #autosomes and X Y

  #Compute mappability for each reference bin.
  mapp_gref<-mapp_hg38 #this is packaged with SCOPE, mappability across bins
  mapp <- rep(1, length(ref))
  #seqlevelsStyle(ref) <- "UCSC"
  for (chr in as.character(unique(seqnames(ref)))) {
      message("Getting mappability for ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- ref[which(as.character(seqnames(ref)) == chr)]
      mapp.chr <- rep(1, length(ref.chr))
      overlap <- as.matrix(findOverlaps(ref.chr, mapp_gref))
      for (i in unique(overlap[, 1])) {
          index.temp <- overlap[which(overlap[, 1] == i), 2]
          overlap.sub <- findOverlaps(ref.chr[i], mapp_gref[index.temp])
          overlap.intersect <- pintersect(ref.chr[i][queryHits(overlap.sub)],mapp_gref[index.temp][subjectHits(overlap.sub)])
          mapp.chr[i] <- sum((mapp_gref$score[index.temp]) * (width(overlap.intersect)))/sum(width(overlap.intersect))
      }
      mapp[chr.index] <- mapp.chr
  }

  #Compute GC for each bin, also from SCOPE
  gc <- rep(NA, length(ref))
  for (chr in unique(seqnames(ref))) {
      message("Getting GC content for chr ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- IRanges(start = start(ref)[chr.index], end = end(ref)[chr.index])
      if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
          chrtemp <- "chrX"
      }
      else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
          chrtemp <- "chrY"
      }
      else {
          chrtemp <- as.numeric(mapSeqlevels(as.character(chr), 
              "NCBI")[1])
      }
      if (length(chrtemp) == 0) 
      message("Chromosome cannot be found in NCBI database. ")
      chrm <- unmasked(genome[[chrtemp]])
      seqs <- Views(chrm, ref.chr)
      af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
      gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }

  ref@elementMetadata$gc<-gc
  ref@elementMetadata$mapp<-mapp
  if(save_rds){
  saveRDS(ref,file=ref_outname)}
  return(ref)
}

get_sample_coverage<-function(bam_in="EXP220921HM_BCMM_WG01_S9_L002_R1_001.batch.dedup.RG.bam",ref,samp_name="sample_9"){
  sampname<-samp_name
    seg.dup <- read.table(system.file("extdata", "GRCh38GenomicSuperDup.tab", package = "WGSmapp"))
    gaps <- read.table(system.file("extdata", "hg38gaps.txt", package = "WGSmapp"))
    seg.dup <- seg.dup[!is.na(match(seg.dup[,1], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
    seg.dup <- GRanges(seqnames = seg.dup[,1], ranges = IRanges(start = seg.dup[,2], end = seg.dup[,3]))
    gaps <- gaps[!is.na(match(gaps[,2], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
    gaps <- GRanges(seqnames = gaps[,2], ranges = IRanges(start = gaps[,3], end = gaps[,4]))
    mask.ref <- sort(c(seg.dup, gaps))

    Y <- matrix(nrow = length(ref), ncol = length(sampname))
    rownames(Y) <- paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
    colnames(Y) <- sampname
    bamurl <- bam_in
    what <- c("rname", "pos", "mapq", "qwidth")
    flag <- scanBamFlag( isDuplicate = FALSE, isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE) # isFirstMateRead = TRUE #isPaired = TRUE,
    param <- ScanBamParam(what = what, flag = flag)
    bam <- scanBam(bamurl, param = param)[[1]]
    message("Getting coverage for sample ", ": ", sampname, "...", sep = "")
    
    bam.ref <- GRanges(seqnames = bam$rname, ranges = IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
    bam.ref <- bam.ref[bam$mapq >= 20] #Q20 threshold
    bam.ref <- suppressWarnings(bam.ref[countOverlaps(bam.ref, mask.ref) == 0])
    Y[, 1] <- countOverlaps(ref, bam.ref)
    return(Y)
}

get_HMMcopy_counts<-function(in_bed,outname="sample_1",bulk_bam,ref_outname=NULL,MAKE_NEW_REF=FALSE,SMALL_REF=FALSE){
    bins<-makeGRangesFromDataFrame(in_bed)
    if(MAKE_NEW_REF|SMALL_REF){
      if(MAKE_NEW_REF){
      ref<-set_up_ref(bins=bins,ref_outname=ref_outname,save_rds=TRUE)}
      else{ref<-set_up_ref(bins=bins,ref_outname=ref_outname,save_rds=FALSE)} #bins is granges of windows to use
    } else {
    ref<-readRDS(ref_outname) #bins is granges of windows to use
    }
    y<-get_sample_coverage(bam_in=bulk_bam,ref=ref,samp_name=outname)
    count<-cbind(as.data.frame(ref),y)
    colnames(count)<-c("chr","start","end","width","strand","gc","map","reads")
    count<-count[c("chr","start","end","reads","gc","map")]
    count$gc<-count$gc/100
    count<-data.table(count)
    count<-correctReadcount(count)
    count$chr<-as.character(count$chr)
    count<-count[count$chr!="chrY",]
    seg<-HMMsegment(count)
    count$state<-seg$state
    count$state<-as.character(count$state)
    count$chr<-factor(count$chr,levels=paste0("chr",c(1:22,"X")))
    count<-count[order(count$chr,count$start),]
    count$row_order<-1:nrow(count)
    return(count)
}

plot_bulk_genome<-function(count){
  plt<-ggplot(count,aes(x=row_order,y=copy,color=as.character(state)))+
    scale_color_manual(values=cols)+
    geom_point(size=1,alpha=1)+
    ylab("")+
    xlab("")+
    ylim(-3,3)+
    facet_grid(~chr,space="free",scales="free_x")+
    theme_minimal()+
    theme(axis.text.y = element_text(size=30),
        axis.text.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0.1,"lines"),
        strip.background = element_blank(), 
        legend.position="none",
        panel.border = element_rect(colour = "black", fill = NA,size=3))
  return(plt)
}


plot_singlecell_cnvs<-function(dat=dat,cnv=t(infercnv_obj@expr.data),assay="infercnv",outname=outname,wd=wd,bulk_plot=plt_100kb,chr_split=infercnv_obj@gene_order$chr,sum_windows=hmmcopy_infercnv_win,file_in,amp_value,del_value){
  #dat is full path to seurat object
  dat_file_path=file_in
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #this is same as initial run of inferCNV, just didn't save seurat object
  cnv_ref<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",]),]
  cnv<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]
  col_fun = colorRamp2(c(min(unlist(cnv)), median(unlist(cnv)), max(unlist(cnv))), c("blue", "white", "red"))

  #discretized window calls
  cnv_discrete<-matrix(0,ncol=ncol(cnv),nrow=nrow(cnv))
  cnv_discrete[which(cnv>=amp_value,arr.ind=T)]<-1
  cnv_discrete[which(cnv<=del_value,arr.ind=T)]<--1
  row.names(cnv_discrete)<-row.names(cnv)
  colnames(cnv_discrete)<-colnames(cnv)
  # discrete_col<-setNames(c("blue","white","red"),nm=c("-1","0","1"))

  # dist_method="manhattan"
  # dist_x<-philentropy::distance(cnv_discrete,method=dist_method,as.dist.obj=T,use.row.names=T)
  # dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
  # k_search<-find_k(dend,krange=2:10) #search for optimal K from 2-10
  # k_clus_number<-k_search$nc
  # k_clus_id<-k_search$pamobject$clustering
  # dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters
  # saveRDS(dend,file=paste0(wd,"/",outname,".",assay,".dend.Rds")) #save dendrogram

  # #set up heatmap annotations
  # met<-as.data.frame(dat@meta.data)
  # met_ref<-met[row.names(met) %in% row.names(cnv_ref),]
  # met<-met[row.names(met) %in% row.names(cnv),]
  # if(any(!(unique(met$PAM50_designation) %in% names(pam50_colors)))){
  #   met[met$PAM50_designation %in% unique(met$PAM50_designation)[!(unique(met$PAM50_designation) %in% names(pam50_colors))],]$PAM50_designation<-"NA"}
  # if(any(!(unique(met$EMBO_designation) %in% names(embo_colors)))){
  #   met[met$EMBO_designation %in% unique(met$EMBO_designation)[!(unique(met$EMBO_designation) %in% names(embo_colors))],]$EMBO_designation<-"NA"}
  # read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
  #   max(met$gex_exonic_umis+met$gex_intronic_umis)), 
  #   c("white","black"))

  # ha = HeatmapAnnotation(which="row",
  #   cell_type=met$predicted.id,
  #   read_count= met$gex_exonic_umis+met$gex_intronic_umis,
  #   pam_50=met$PAM50_designation,
  #   embo=met$EMBO_designation,
  #         col = list(cell_type = type_cols,
  #           read_count=read_count_col,
  #           embo=embo_colors,
  #           pam_50=pam50_colors))

  # sum_windows<-sum_windows[colnames(cnv),]
  # state_cols = setNames(brewer.pal(n=6,name="RdBu"), nm = c("6","5","4","3","2","1")) # black, red, green, blue
  # copy_col<-colorRamp2(c(min(sum_windows$HMMcopy_mean_copymetric,na.rm=TRUE),0,
  #   max(sum_windows$HMMcopy_mean_copymetric,na.rm=TRUE)), 
  #   c("blue","white","red"))

  # hwin = HeatmapAnnotation(which="column",
  #   copy_state=sum_windows$HMMcopy_mode_copystate,
  #   copy_metric=sum_windows$HMMcopy_mean_copymetric,
  #         col = list(copy_state = state_cols,
  #           copy_metric=copy_col))
  # plt1<-Heatmap(cnv,
  #     show_row_names=F,
  #     show_column_names=F,
  #     column_order=1:ncol(cnv),
  #     col=col_fun,
  #     cluster_rows=dend,
  #     left_annotation=ha,
  #     top_annotation=hwin,
  #     column_split=chr_split)

  # ha_ref = HeatmapAnnotation(which="row",
  #   cell_type=met_ref$predicted.id,
  #   read_count= met_ref$gex_exonic_umis+met_ref$gex_intronic_umis,
  #         col = list(cell_type = type_cols,
  #           read_count=read_count_col))
  # plt1_ref<-Heatmap(cnv_ref,
  #     show_row_names=F,
  #     show_column_names=F,
  #     column_order=1:ncol(cnv),
  #     col=col_fun,
  #     left_annotation=ha_ref,
  #     column_split=chr_split)

  # plt2<-Heatmap(cnv_discrete,
  #     show_row_names=F,
  #     show_column_names=F,
  #     column_order=1:ncol(cnv),
  #     col=discrete_col,
  #     cluster_rows=dend,
  #     left_annotation=ha,
  #     top_annotation=hwin,
  #     column_split=chr_split)

  #   pdf(paste0(wd,"/",outname,".",assay,".heatmap.pdf"),width=40)
  #   print(bulk_plot)
  #   print(plt1_ref)
  #   print(plt1)
  #   print(plt2)
  #   dev.off()
  #   system(paste0("slack -F ",paste0(wd,"/",outname,".",assay,".heatmap.pdf")," ryan_todo"))
    return(cnv_discrete)
}

HMMcopy_comparison<-function(x,file_in="phase2.QC.filt.SeuratObject.rds"){

    if(x %in% 1:12){
      sample_name<-paste0("sample_",x)
      wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
      outname<-paste0("sample_",x)
    }else if(x %in% 13:20){
      sample_name<-paste0("sample_",x)
      wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
      outname<-paste0("sample_",x)
    }else{
      sample_name<-x
      wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
      outname<-x
    }
    dat<-readRDS(file_in)
    dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
    dir_in<-dirname(file_in)
    bulk_bam<-bam_vec[outname]
    bulk_bam<-paste(bamfolder,bulk_bam,sep="/")

  #100kb windows
    print(paste(outname,"100kb windows"))
    ref_outname=paste(dirname(bulk_bam),"100kb_windows.rds",sep="/")
    if(x==1){ MAKE_NEW_REF=FALSE
    } else {MAKE_NEW_REF=FALSE} #only make the ref windows for the first sample
    bins <- tileGenome(seqinfo(genome), tilewidth = 100 * 1000, cut.last.tile.in.chrom = TRUE) #set bins by other CNV callers
    counts<-get_HMMcopy_counts(in_bed=bins,outname=outname,bulk_bam=bulk_bam,ref_outname=ref_outname,MAKE_NEW_REF=MAKE_NEW_REF)
    saveRDS(counts,file=paste0(wd,"/",outname,"_bulkWGS_HMMcopy.100kb.rds"))
    counts<-readRDS(file=paste0(wd,"/",outname,"_bulkWGS_HMMcopy.100kb.rds"))
    plt_100kb<-plot_bulk_genome(counts)+ggtitle(paste(outname,"100kb Bins",mean(counts$start-counts$end)))

  #InferCNV
    assay="InferCNV"
    #3 state model is here (gene by cell name data is in i3_hmm@expr.data)
    i3_hmm<-readRDS(paste0(wd,"/",outname,"_inferCNV","/19_HMM_pred.repr_intensitiesHMMi3.hmm_mode-samples.Pnorm_0.5.infercnv_obj"))
    print(paste(outname,"InferCNV windows"))
    #infercnv_obj<-readRDS(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
    infercnv_obj<-i3_hmm
    #Format Data
    infercnv_bed<-infercnv_obj@gene_order
    cnv_in<-t(infercnv_obj@expr.data)
    chr_in<-infercnv_obj@gene_order$chr
    chr_in<-factor(chr_in,levels=unique(chr_in))
    #Summarize Data over WGS
    refGR<-makeGRangesFromDataFrame(counts)
    testGR<-makeGRangesFromDataFrame(infercnv_bed)
    hits<-findOverlaps(refGR,testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    bed_overlaps<-as.data.frame(cbind(as.data.frame(hits),percentOverlap))
    hmmcopy_infercnv_win<-cbind(infercnv_bed,
      HMMcopy_mean_copymetric=unlist(lapply(1:nrow(infercnv_bed),function(x) 
          mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE))),
      HMMcopy_weightedmean_copymetric=unlist(lapply(1:nrow(infercnv_bed),function(x) 
          weighted.mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE,w=bed_overlaps[bed_overlaps$subjectHits==x,]$percentOverlap))),
      HMMcopy_mode_copystate=unlist(lapply(1:nrow(infercnv_bed),function(x) names(sort(-table(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$state)))[1])))
    #Cluster and Plot
    disc_infercnv<-plot_singlecell_cnvs(
        dat=dat,
        cnv=cnv_in,
        assay=assay,
        outname=outname,
        wd=wd,
        chr_split=chr_in,
        bulk_plot=plt_100kb,
        sum_windows=hmmcopy_infercnv_win,
        file_in=file_in,
        amp_value=1.5,
        del_value=0.5)
    write.table(sep="\t",col.names=T,row.names=T,quote=F,cnv_in,file=paste0(out_dir,"/",outname,"_scCNV_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,disc_infercnv,file=paste0(out_dir,"/",outname,"_scCNV_discrete_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,hmmcopy_infercnv_win,file=paste0(out_dir,"/",outname,"_bulkWGS_",assay,"_bins.tsv"))

  #CASPER 
    #casper discretized matrix:
    dir_in<-wd
    casper_cnv<-readRDS(paste0(dir_in,"/casper/",outname,".finalgenemat.rds"))
    #Run different segmentation scales? https://rpubs.com/akdes/673120 (section 3)
    assay="CASPER"
    print(paste(outname,"CASPER windows"))
    casper_obj<-readRDS(paste0(dir_in,"/casper/",outname,".finalobj.rds"))
    #Format Data
    casper_bed<-casper_obj@annotation[,c("Chr","start","end")]
    row.names(casper_bed)<-casper_obj@annotation$Gene
    casper_bed$Chr<-paste0("chr",casper_bed$Chr)
    colnames(casper_bed)<-c("chr","start","end")
    #Summarize Data over WGS
    refGR<-makeGRangesFromDataFrame(counts)
    testGR<-makeGRangesFromDataFrame(casper_bed)
    hits<-findOverlaps(refGR,testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    bed_overlaps<-as.data.frame(cbind(as.data.frame(hits),percentOverlap))
    hmmcopy_casper_win<-cbind(casper_bed,
      HMMcopy_mean_copymetric=unlist(lapply(1:nrow(casper_bed),function(x) 
          mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE))),
      HMMcopy_weightedmean_copymetric=unlist(lapply(1:nrow(casper_bed),function(x) 
          weighted.mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE,w=bed_overlaps[bed_overlaps$subjectHits==x,]$percentOverlap))),
      HMMcopy_mode_copystate=unlist(lapply(1:nrow(casper_bed),function(x) names(sort(-table(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$state)))[1])))

    cnv_in<-t(casper_cnv)
    cnv_in<-cnv_in[,colnames(cnv_in)%in%casper_obj@annotation.filt$Gene]
    chr_in<-casper_obj@annotation.filt[casper_obj@annotation.filt$Gene %in% colnames(cnv_in),]
    chr_in<-paste0("chr",chr_in$Chr)
    chr_in<-factor(chr_in,levels=unique(chr_in))
    hmmcopy_casper_win<-hmmcopy_casper_win[colnames(cnv_in),]
    #Cluster and Plot
    disc_casper<-plot_singlecell_cnvs(
        dat=dat,
        cnv=cnv_in,
        assay="casper",
        outname=outname,
        wd=wd,
        chr_split=chr_in,
        bulk_plot=plt_100kb,
        sum_windows=hmmcopy_casper_win,
        file_in=file_in,
        amp_value=1,
        del_value=-1)
     write.table(sep="\t",col.names=T,row.names=T,quote=F,cnv_in,file=paste0(out_dir,"/",outname,"_scCNV_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,disc_casper,file=paste0(out_dir,"/",outname,"_scCNV_discrete_",assay,".tsv"))
     write.table(sep="\t",col.names=T,row.names=T,quote=F,hmmcopy_casper_win,file=paste0(out_dir,"/",outname,"_bulkWGS_",assay,"_bins.tsv"))

  #CopyKAT 

    assay="CopyKAT"
    #to set CNV discrete changes, as per correspondence suggetions with Ruli Gao, 1.5x SD threshold, 1.5 absolute distance, or use +/-0.25 as cutoff
    print(paste(outname,"CopyKat windows"))
    copykat_obj<-readRDS(paste0(dir_in,"/copykat/",outname,".copykat.RDS"))
    #Format Data
    copykat_bed<-copykat_obj$CNAmat[1:2]
    copykat_bed$chrom<-paste0("chr",copykat_bed$chrom)
    copykat_bed[copykat_bed$chrom=="chr23",]$chrom<-"chrX"
    bed_split<-split(x=copykat_bed,f=copykat_bed$chrom)
    copykat_bed<-do.call("rbind",lapply(bed_split,function(x) {
      print(x[1,1])
      chrend<-chr_end[chr_end$chr==x[1,1],]$length
      x$chromend<-c(x$chrompos[1:length(x$chrompos)-1]+diff(x$chrompos),chrend)
      return(x)}))
    colnames(copykat_bed)<-c("chr","start","end")

    #Summarize Data over WGS
    refGR<-makeGRangesFromDataFrame(counts)
    testGR<-makeGRangesFromDataFrame(copykat_bed)
    hits<-findOverlaps(refGR,testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    bed_overlaps<-as.data.frame(cbind(as.data.frame(hits),percentOverlap))
    hmmcopy_copykat_win<-cbind(copykat_bed,
      HMMcopy_mean_copymetric=unlist(lapply(1:nrow(copykat_bed),function(x) 
          mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE))),
      HMMcopy_weightedmean_copymetric=unlist(lapply(1:nrow(copykat_bed),function(x) 
          weighted.mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE,w=bed_overlaps[bed_overlaps$subjectHits==x,]$percentOverlap))),
      HMMcopy_mode_copystate=unlist(lapply(1:nrow(copykat_bed),function(x) names(sort(-table(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$state)))[1])))


    row.names(hmmcopy_copykat_win)<-row.names(1:nrow(hmmcopy_copykat_win))
    cnv_in<-t(copykat_obj$CNAmat[,4:ncol(copykat_obj$CNAmat)])
    row.names(cnv_in)<-gsub("\\.","-",row.names(cnv_in))
    chr_in<-paste0("chr",copykat_obj$CNAmat[,1])
    chr_in<-factor(chr_in,levels=unique(chr_in))
    sd_value<-sd(unlist(cnv_in))
    norm_value<-mean(unlist(cnv_in))
    amp_value<-norm_value+(sd_value*1.5)
    del_value<-norm_value-(sd_value*1.5)
    #Cluster and Plot
    disc_copykat<-plot_singlecell_cnvs(
        dat=dat,
        cnv=cnv_in,
        assay="CopyKAT",
        outname=outname,
        wd=wd,
        chr_split=chr_in,
        bulk_plot=plt_100kb,
        sum_windows=hmmcopy_copykat_win,
        file_in=file_in,
        amp_value=amp_value,
        del_value=del_value)
     write.table(sep="\t",col.names=T,row.names=T,quote=F,cnv_in,file=paste0(out_dir,"/",outname,"_scCNV_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,disc_copykat,file=paste0(out_dir,"/",outname,"_scCNV_discrete_",assay,".tsv"))
     write.table(sep="\t",col.names=T,row.names=T,quote=F,hmmcopy_copykat_win,file=paste0(out_dir,"/",outname,"_bulkWGS_",assay,"_bins.tsv"))

  #COPYSCAT
    assay="copyscat"
    copyscat_dat<-readRDS(file=paste0(dir_in,"/copyscat/",outname,"copyscat_cnvs_matrix.rds"))
    print(paste(outname,"Copyscat windows"))
    copyscat_obj<-readRDS(file=paste0(dir_in,"/copyscat/",outname,"copyscat_cnvs.rds"))
    #Format Data
    copyscat_dat<-copyscat_dat[[1]]
    row.names(copyscat_dat)<-copyscat_dat[,1]
    copyscat_dat<-t(copyscat_dat[,2:ncol(copyscat_dat)])
    copyscat_chr<-unique(copyscat_obj[[1]]$Chrom[!(copyscat_obj[[1]]$Chrom %in% row.names(copyscat_dat))])
    copyscat_cellid<-colnames(copyscat_dat)
    copyscat_cellid<-paste(outname,colnames(copyscat_dat),sep="_")
    copyscat_cellid[length(copyscat_cellid)]<-"medianNorm"
    copyscat_unreported <- data.frame(matrix(ncol = length(copyscat_cellid), nrow = length(copyscat_chr),data=2))
    row.names(copyscat_unreported)<-copyscat_chr
    colnames(copyscat_unreported)<-copyscat_cellid
    colnames(copyscat_dat)<-copyscat_cellid
    copyscat_dat<-rbind(copyscat_dat,copyscat_unreported)
    copyscat_dat<-copyscat_dat[match(cytoband$chr,row.names(copyscat_dat)),]
    copyscat_dat<-copyscat_dat[!startsWith(prefix="NA",row.names(copyscat_dat)),]
    copyscat_bed<-cytoband[cytoband$chr %in% row.names(copyscat_dat),]
    copyscat_bed$chr<-paste0("chr",copyscat_bed$chrom)
    copyscat_bed<-copyscat_bed[,c("chr","start","end")]

    #Summarize Data over WGS
    refGR<-makeGRangesFromDataFrame(counts)
    testGR<-makeGRangesFromDataFrame(copyscat_bed)
    hits<-findOverlaps(refGR,testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    bed_overlaps<-as.data.frame(cbind(as.data.frame(hits),percentOverlap))
    hmmcopy_copyscat_win<-cbind(copyscat_bed,
      HMMcopy_mean_copymetric=unlist(lapply(1:nrow(copyscat_bed),function(x) 
          mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE))),
      HMMcopy_weightedmean_copymetric=unlist(lapply(1:nrow(copyscat_bed),function(x) 
          weighted.mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE,w=bed_overlaps[bed_overlaps$subjectHits==x,]$percentOverlap))),
      HMMcopy_mode_copystate=unlist(lapply(1:nrow(copyscat_bed),function(x) names(sort(-table(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$state)))[1])))

    cnv_in<-t(copyscat_dat)
    chr_in<-substr(colnames(cnv_in),1,nchar(colnames(cnv_in))-1)
    chr_in<-factor(chr_in,levels=unique(chr_in))
    cnv_in<-cnv_in[1:nrow(cnv_in)-1,]#remove median norm measure
    row.names(hmmcopy_copyscat_win)<-colnames(cnv_in)
    #Cluster and Plot
    disc_copyscat<-plot_singlecell_cnvs(
        dat=dat,
        cnv=cnv_in,
        assay="CopySCAT",
        outname=outname,
        wd=wd,
        chr_split=chr_in,
        bulk_plot=plt_100kb,
        sum_windows=hmmcopy_copyscat_win,
        file_in=file_in,
        amp_value=2,
        del_value=0)
     write.table(sep="\t",col.names=T,row.names=T,quote=F,cnv_in,file=paste0(out_dir,"/",outname,"_scCNV_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,disc_copyscat,file=paste0(out_dir,"/",outname,"_scCNV_discrete_",assay,".tsv"))
     write.table(sep="\t",col.names=T,row.names=T,quote=F,hmmcopy_copyscat_win,file=paste0(out_dir,"/",outname,"_bulkWGS_",assay,"_bins.tsv"))


  #RobustCNV
    assay="RobustCNV"
    robustcnv_obj<-read.csv(paste0("/home/groups/CEDAR/scATACcnv/Hisham_data/bed_files/1MB/",outname,"_1MB_robustCNV.csv"))
    #robustcnv_obj<-as.data.frame(t(read.table(paste0("/home/groups/CEDAR/scATACcnv/Hisham_data/bed_files/1MB/","sample_4_scCNV_discrete_RobustCNV.tsv")))) for sample 4
    colnames(robustcnv_obj)<-paste(outname,colnames(robustcnv_obj),sep="_")
    #Format Data
    colnames(robustcnv_obj)<-gsub(colnames(robustcnv_obj),pattern="\\.",replacement="-")
    robustcnv_bed<-read.table(paste0("/home/groups/CEDAR/scATACcnv/Hisham_data/bed_files/1MB/","window_1MB.bed"))
    colnames(robustcnv_bed)<-c("chr","start","end","win")
    robustcnv_bed<-robustcnv_bed[!robustcnv_bed$chr %in% c("chrX","chrY"),]

    #Summarize Data over WGS
    refGR<-makeGRangesFromDataFrame(counts)
    testGR<-makeGRangesFromDataFrame(robustcnv_bed)
    hits<-findOverlaps(refGR,testGR)
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    bed_overlaps<-as.data.frame(cbind(as.data.frame(hits),percentOverlap))
    hmmcopy_robustcnv_win<-cbind(robustcnv_bed,
      HMMcopy_mean_copymetric=unlist(lapply(1:nrow(robustcnv_bed),function(x) 
          mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE))),
      HMMcopy_weightedmean_copymetric=unlist(lapply(1:nrow(robustcnv_bed),function(x) 
          weighted.mean(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$copy,na.rm=TRUE,w=bed_overlaps[bed_overlaps$subjectHits==x,]$percentOverlap))),
      HMMcopy_mode_copystate=unlist(lapply(1:nrow(robustcnv_bed),function(x) names(sort(-table(counts[bed_overlaps[bed_overlaps$subjectHits==x,]$queryHits,]$state)))[1])))

    cnv_in<-t(robustcnv_obj)
    colnames(cnv_in)<-robustcnv_bed$win
    row.names(hmmcopy_robustcnv_win)<-colnames(cnv_in)
    chr_in<-hmmcopy_robustcnv_win$chr
    sd_value<-sd(unlist(cnv_in))
    norm_value<-mean(unlist(cnv_in))
    amp_value<-norm_value+(sd_value*1.5)
    del_value<-norm_value-(sd_value*1.5)
    #Cluster and Plot
    disc_robustcnv<-plot_singlecell_cnvs(
        dat=dat,
        cnv=cnv_in,
        assay="RobustCNV",
        outname=outname,
        wd=wd,
        chr_split=chr_in,
        bulk_plot=plt_100kb,
        sum_windows=hmmcopy_robustcnv_win,
        file_in=file_in,
        amp_value=amp_value,
        del_value=del_value)
    write.table(sep="\t",col.names=T,row.names=T,quote=F,cnv_in,file=paste0(out_dir,"/",outname,"_scCNV_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,disc_robustcnv,file=paste0(out_dir,"/",outname,"_scCNV_discrete_",assay,".tsv"))
    write.table(sep="\t",col.names=T,row.names=T,quote=F,hmmcopy_robustcnv_win,file=paste0(out_dir,"/",outname,"_bulkWGS_",assay,"_bins.tsv"))
}

bamfolder <- "/home/groups/CEDAR/mulqueen/projects/multiome/221004_wgs/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
bamFile <- list.files(bamfolder, pattern = 'dedup.RG.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- paste("sample",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20),sep="_") #bams are ordered by sample number as well #3,11,12
bam_vec<-setNames(bamFile,nm=sampname_raw)
colnames(cytoband)<-c("chrom","start","end","arm")
cytoband$chr<-paste0("chr",cytoband$chrom,cytoband$arm)
chr_end<-data.frame(chr=BSgenome.Hsapiens.UCSC.hg38@seqinfo@seqnames,length=BSgenome.Hsapiens.UCSC.hg38@seqinfo@seqlengths)
cols = setNames(brewer.pal(n=6,name="RdBu"), nm = c("6","5","4","3","2","1")) # black, red, green, blue
genome <- BSgenome.Hsapiens.UCSC.hg38
out_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/cnv_comparison"
system(paste("mkdir",out_dir))
HMMcopy_comparison(x)
lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20),HMMcopy_comparison) 

```

Writing out as a batch script for slurm job submission

compare_hmm_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=0-13
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=10 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=("1" "3" "4" "5" "6" "7" "8" "9" "10" "11" "15" "16" "19" "20") 
sample_in=${array_in[$SLURM_ARRAY_TASK_ID]}
multiome_dir="/home/groups/CEDAR/mulqueen/projects/multiome"

srun Rscript ${multiome_dir}/compare_hmm_slurm.sh $sample_in

```

Job submit all HMMcopy jobs for comparison
```bash
sbatch compare_hmm_slurm.sh
```

## Epithelial Subtyping
<!-- Done -->

### Additional Cell Signatures
<!-- Done -->
From https://github.com/yunshun/HumanBreast10X/tree/main/Signatures

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
#downloaded files from
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/Human-PosSigGenes.RData
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/ImmuneMarkers2.txt
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/PAM50.txt
```

### Use EMBO and Swarbrick Paper Cell Types to Define Signatures
Using package genefu for PAM50 pseudobulk assignment.
https://www.bioconductor.org/packages/release/bioc/vignettes/genefu/inst/doc/genefu.html


```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
library(dplyr)
library(org.Hs.eg.db)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

#Using genefu per pseudobulked sample
#data: Matrix of annotations with at least one column named "EntrezGene.ID"
#'   (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust
#'   model), dimnames being properly defined.
#do.mapping TRUE if the mapping through Entrez Gene ids must be performed
#'   (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.

#CDCA1 KNTC2 ORC6L use different names in our data
#NUF2, NDC80, ORC6 resp.
pam50_genes<-c('ACTR3B', 'ANLN', 'BAG1', 'BCL2', 'BIRC5', 'BLVRA', 'CCNB1', 'CCNE1', 'CDC20', 'CDC6', 'NUF2', 'CDH3', 'CENPF', 'CEP55', 'CXXC5', 'EGFR', 'ERBB2', 'ESR1', 'EXO1', 'FGFR4', 'FOXA1', 'FOXC1', 'GPR160', 'GRB7', 'KIF2C', 'NDC80', 'KRT14', 'KRT17', 'KRT5', 'MAPT', 'MDM2', 'MELK', 'MIA', 'MKI67', 'MLPH', 'MMP11', 'MYBL2', 'MYC', 'NAT1', 'ORC6', 'PGR', 'PHGDH', 'PTTG1', 'RRM2', 'SFRP1', 'SLC39A6', 'TMEM45B', 'TYMS', 'UBE2C', 'UBE2T')

#dat<-subset(dat,EMBO_predicted.id %in% c("epithelial","cycling.epithelial")) #trying pam50 assignment with epithelial cell subset first

sample_names<-paste(unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(1))),
  unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(2))),sep="_")
counts<-as.data.frame(t(dat[["RNA"]]@counts)) 
counts<-cbind(counts,sample_names)
counts<-as.data.frame(counts %>% group_by(sample_names) %>% summarize_all(funs(sum)))
row.names(counts)<-counts$sample_name
counts<-counts[,2:ncol(counts)]
counts<-counts[,colSums(counts)>0]
#dat_in<-as.data.frame(t(counts[x,]))
dat_in<-counts
dat_in<-dat_in[!(row.names(dat_in) %in% c("RM_4","sample_15","sample_19")),] #exclude NAT samples
dat_in<-NormalizeData(dat_in,normalization.method="CLR")
dannot<-as.data.frame(cbind(Gene.Symbol=colnames(dat_in),EntrezGene.ID=mapIds(org.Hs.eg.db, colnames(dat_in), 'ENTREZID', 'SYMBOL'),probe=colnames(dat_in)))
pam50_out<-molecular.subtyping(sbt.model="pam50",data=dat_in,annot=dannot,do.mapping=TRUE,verbose=T)

#try this as well
#pam50_out_model<-intrinsic.cluster(data=dat_in,annot=dannot,do.mapping=TRUE,std="robust",intrinsicg=pam50$centroids.map[,c("probe","EntrezGene.ID")],verbose=T,mins=0)#,mapping=dannot)
#pam50_out<-intrinsic.cluster.predict(sbt.model=pam50_out_model$model, data=dat_in, annot=dannot, do.mapping=TRUE,do.prediction.strength=TRUE,verbose=TRUE)
#saveRDS(pam50_out,file="pseudobulk_pam50.rds")

pam50_meta<-setNames(nm=row.names(dat@meta.data),pam50_out$subtype[match(dat$sample, names(pam50_out$subtype))])
dat<-AddMetaData(dat,pam50_meta,col.name="pseudobulk_genefu_pam50")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

#tried just epithelial, tried both old method (intrinsic cluster) and updated method (molecular subtyping). maybe play around with normalizing first?
#limit to epithelial? or maybe read up on proper normalization? our HER2+ isn't being labelled as such

```

Running SSpbc method as well

Using https://github.com/StaafLab/sspbc/archive/refs/heads/main.zip for multiple classifications
https://www.nature.com/articles/s41523-022-00465-3#code-availability


```R
#wget https://github.com/StaafLab/sspbc/archive/refs/heads/main.zip
#file located in /home/groups/CEDAR/mulqueen/src/sspbc/sspbc-main/package
#R CMD INSTALL sspbc_1.0.tar.gz
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
library(dplyr)
library(org.Hs.eg.db)
library(sspbc)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

#dat<-subset(dat,EMBO_predicted.id %in% c("epithelial","cycling.epithelial")) #trying pam50 assignment with epithelial cell subset first

sample_names<-paste(unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(1))),
  unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(2))),sep="_")
counts<-as.data.frame(t(dat[["RNA"]]@counts)) 
counts<-cbind(counts,sample_names)
counts<-as.data.frame(counts %>% group_by(sample_names) %>% summarize_all(funs(sum)))
row.names(counts)<-counts$sample_name
counts<-counts[,2:ncol(counts)]
counts<-counts[,colSums(counts)>0]
dat_in<-counts
dat_in<-dat_in[!(row.names(dat_in) %in% c("RM_4","sample_15","sample_19")),] #exclude NAT samples
dat_in<-as.data.frame(t(dat_in))

#set up matrix by unique entrez gene names
dat_in<-dat_in[!duplicated(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
dat_in<-dat_in[!isNA(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
row.names(dat_in)<-mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')
myresults <- applySSP(gex=as.matrix(dat_in), id=row.names(dat_in), ssp.name="ssp.pam50",id.type="EntrezGene",report=TRUE)



#dat<-readRDS("phase2.QC.filt.SeuratObject.rds")
dat_pam50<-setNames(nm=row.names(dat@meta.data),myresults[match(dat@meta.data$sample,row.names(myresults)),1])
dat<-AddMetaData(dat,dat_pam50,col.name="pseudobulk_sspbc_PAM50")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```
<!-- Done -->

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

#cell lineage
  load("/home/groups/CEDAR/mulqueen/ref/embo/Human-PosSigGenes.RData")
  ls()
  #[1] "Basal" "LP"    "ML"    "Str" #lineage types
  lineage_in=list(EMBO_Basal=Basal,EMBO_LP=LP,EMBO_ML=ML,EMBO_Str=Str)

#Immune markers
  immune_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/ImmuneMarkers2.txt",header=T,sep="\t")
  immune_in<-lapply(split(immune_in,immune_in$CellType),function(x) x$Signatures)#split up data frame to a named list of genes per cell type
  names(immune_in)<-paste0("EMBO_",names(immune_in))#rename the list just so we can track the source

#using both the given PAM50 short list and the Swarbrick supplied more extensive gene list below
  PAM50_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/PAM50.txt",header=T,sep="\t")
  PAM50_in<-lapply(split(PAM50_in,PAM50_in$Subtype),function(x) x$Gene)#split up data frame to a named list of genes per cell type
  names(PAM50_in)<-paste0("PAM50_",names(PAM50_in))
  features_in=c(immune_in,PAM50_in)   

molecular.subtyping(sbt.model="pam50")
#SCSubtype Features determined by Swarbrick manuscript (Supp Table 4)
  module_feats<-list()
  module_feats[["Basal_SC"]]=c('EMP1', 'TAGLN', 'TTYH1', 'RTN4', 'TK1', 'BUB3', 'IGLV3.25', 'FAM3C', 'TMEM123', 'KDM5B', 'KRT14', 'ALG3', 'KLK6', 'EEF2', 'NSMCE4A', 'LYST', 'DEDD', 'HLA.DRA', 'PAPOLA', 'SOX4', 'ACTR3B', 'EIF3D', 'CACYBP', 'RARRES1', 'STRA13', 'MFGE8', 'FRZB', 'SDHD', 'UCHL1', 'TMEM176A', 'CAV2', 'MARCO', 'P4HB', 'CHI3L2', 'APOE', 'ATP1B1', 'C6orf15', 'KRT6B', 'TAF1D', 'ACTA2', 'LY6D', 'SAA2', 'CYP27A1', 'DLK1', 'IGKV1.5', 'CENPW', 'RAB18', 'TNFRSF11B', 'VPS28', 'HULC', 'KRT16', 'CDKN2A', 'AHNAK2', 'SEC22B', 'CDC42EP1', 'HMGA1', 'CAV1', 'BAMBI', 'TOMM22', 'ATP6V0E2', 'MTCH2', 'PRSS21', 'HDAC2', 'ZG16B', 'GAL', 'SCGB1D2', 'S100A2', 'GSPT1', 'ARPC1B', 'NIT1', 'NEAT1', 'DSC2', 'RP1.60O19.1', 'MAL2', 'TMEM176B', 'CYP1B1', 'EIF3L', 'FKBP4', 'WFDC2', 'SAA1', 'CXCL17', 'PFDN2', 'UCP2', 'RAB11B', 'FDCSP', 'HLA.DPB1', 'PCSK1N', 'C4orf48', 'CTSC')
  module_feats[["Her2E_SC"]]=c('PSMA2', 'PPP1R1B', 'SYNGR2', 'CNPY2', 'LGALS7B', 'CYBA', 'FTH1', 'MSL1', 'IGKV3.15', 'STARD3', 'HPD', 'HMGCS2', 'ID3', 'NDUFB8', 'COTL1', 'AIM1', 'MED24', 'CEACAM6', 'FABP7', 'CRABP2', 'NR4A2', 'COX14', 'ACADM', 'PKM', 'ECH1', 'C17orf89', 'NGRN', 'ATG5', 'SNHG25', 'ETFB', 'EGLN3', 'CSNK2B', 'RHOC', 'PSENEN', 'CDK12', 'ATP5I', 'ENTHD2', 'QRSL1', 'S100A7', 'TPM1', 'ATP5C1', 'HIST1H1E', 'LGALS1', 'GRB7', 'AQP3', 'ALDH2', 'EIF3E', 'ERBB2', 'LCN2', 'SLC38A10', 'TXN', 'DBI', 'RP11.206M11.7', 'TUBB', 'CRYAB', 'CD9', 'PDSS2', 'XIST', 'MED1', 'C6orf203', 'PSMD3', 'TMC5', 'UQCRQ', 'EFHD1', 'BCAM', 'GPX1', 'EPHX1', 'AREG', 'CDK2AP2', 'SPINK8', 'PGAP3', 'NFIC', 'THRSP', 'LDHB', 'MT1X', 'HIST1H4C', 'LRRC26', 'SLC16A3', 'BACE2', 'MIEN1', 'AR', 'CRIP2', 'NME1', 'DEGS2', 'CASC3', 'FOLR1', 'SIVA1', 'SLC25A39', 'IGHG1', 'ORMDL3', 'KRT81', 'SCGB2B2', 'LINC01285', 'CXCL8', 'KRT15', 'RSU1', 'ZFP36L2', 'DKK1', 'TMED10', 'IRX3', 'S100A9', 'YWHAZ')
  module_feats[["LumA_SC"]]=c('SH3BGRL', 'HSPB1', 'PHGR1', 'SOX9', 'CEBPD', 'CITED2', 'TM4SF1', 'S100P', 'KCNK6', 'AGR3', 'MPC2', 'CXCL13', 'RNASET2', 'DDIT4', 'SCUBE2', 'KRT8', 'MZT2B', 'IFI6', 'RPS26', 'TAGLN2', 'SPTSSA', 'ZFP36L1', 'MGP', 'KDELR2', 'PPDPF', 'AZGP1', 'AP000769.1', 'MYBPC1', 'S100A1', 'TFPI2', 'JUN', 'SLC25A6', 'HSP90AB1', 'ARF5', 'PMAIP1', 'TNFRSF12A', 'FXYD3', 'RASD1', 'PYCARD', 'PYDC1', 'PHLDA2', 'BZW2', 'HOXA9', 'XBP1', 'AGR2', 'HSP90AA1') 
  module_feats[["LumB_SC"]]=c('UGCG', 'ARMT1', 'ISOC1', 'GDF15', 'ZFP36', 'PSMC5', 'DDX5', 'TMEM150C', 'NBEAL1', 'CLEC3A', 'GADD45G', 'MARCKS', 'FHL2', 'CCDC117', 'LY6E', 'GJA1', 'PSAP', 'TAF7', 'PIP', 'HSPA2', 'DSCAM.AS1', 'PSMB7', 'STARD10', 'ATF3', 'WBP11', 'MALAT1', 'C6orf48', 'HLA.DRB1', 'HIST1H2BD', 'CCND1', 'STC2', 'NR4A1', 'NPY1R', 'FOS', 'ZFAND2A', 'CFL1', 'RHOB', 'LMNA', 'SLC40A1', 'CYB5A', 'SRSF5', 'SEC61G', 'CTSD', 'DNAJC12', 'IFITM1', 'MAGED2', 'RBP1', 'TFF1', 'APLP2', 'TFF3', 'TRH', 'NUPR1', 'EMC3', 'TXNIP', 'ARPC4', 'KCNE4', 'ANPEP', 'MGST1', 'TOB1', 'ADIRF', 'TUBA1B', 'MYEOV2', 'MLLT4', 'DHRS2', 'IFITM2')
  module_feats[["proliferation_score"]]<-c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")

#Swarbrick Gene Module Classification (Supp Table 5)
gene_module<-list()
  gene_module[["gene_module_1"]]<-c('ATF3', 'JUN', 'NR4A1', 'IER2', 'DUSP1', 'ZFP36', 'JUNB', 'FOS', 'FOSB', 'PPP1R15A', 'KLF6', 'DNAJB1', 'EGR1', 'BTG2', 'HSPA1B', 'HSPA1A', 'RHOB', 'CLDN4', 'MAFF', 'GADD45B', 'IRF1', 'EFNA1', 'SERTAD1', 'TSC22D1', 'CEBPD', 'CCNL1', 'TRIB1', 'MYC', 'ELF3', 'LMNA', 'NFKBIA', 'TOB1', 'HSPB1', 'BRD2', 'MCL1', 'PNRC1', 'IER3', 'KLF4', 'ZFP36L2', 'SAT1', 'ZFP36L1', 'DNAJB4', 'PHLDA2', 'NEAT1', 'MAP3K8', 'GPRC5A', 'RASD1', 'NFKBIZ', 'CTD-3252C9.4', 'BAMBI', 'RND1', 'HES1', 'PIM3', 'SQSTM1', 'HSPH1', 'ZFAND5', 'AREG', 'CD55', 'CDKN1A', 'UBC', 'CLDN3', 'DDIT3', 'BHLHE40', 'BTG1', 'ANKRD37', 'SOCS3', 'NAMPT', 'SOX4', 'LDLR', 'TIPARP', 'TM4SF1', 'CSRNP1', 'GDF15', 'ZFAND2A', 'NR4A2', 'ERRFI1', 'RAB11FIP1', 'TRAF4', 'MYADM', 'ZC3H12A', 'HERPUD1', 'CKS2', 'BAG3', 'TGIF1', 'ID3', 'JUND', 'PMAIP1', 'TACSTD2', 'ETS2', 'DNAJA1', 'PDLIM3', 'KLF10', 'CYR61', 'MXD1', 'TNFAIP3', 'NCOA7', 'OVOL1', 'TSC22D3', 'HSP90AA1', 'HSPA6', 'C15orf48', 'RHOV', 'DUSP4', 'B4GALT1', 'SDC4', 'C8orf4', 'DNAJB6', 'ICAM1', 'DNAJA4', 'MRPL18', 'GRB7', 'HNRNPA0', 'BCL3', 'DUSP10', 'EDN1', 'FHL2', 'CXCL2', 'TNFRSF12A', 'S100P', 'HSPB8', 'INSIG1', 'PLK3', 'EZR', 'IGFBP5', 'SLC38A2', 'DNAJB9', 'H3F3B', 'TPM4', 'TNFSF10', 'RSRP1', 'ARL5B', 'ATP1B1', 'HSPA8', 'IER5', 'SCGB2A1', 'YPEL2', 'TMC5', 'FBXO32', 'MAP1LC3B', 'MIDN', 'GADD45G', 'VMP1', 'HSPA5', 'SCGB2A2', 'TUBA1A', 'WEE1', 'PDK4', 'STAT3', 'PERP', 'RBBP6', 'KCNQ1OT1', 'OSER1', 'SERP1', 'UBE2B', 'HSPE1', 'SOX9', 'MLF1', 'UBB', 'MDK', 'YPEL5', 'HMGCS1', 'PTP4A1', 'WSB1', 'CEBPB', 'EIF4A2', 'S100A10', 'ELMSAN1', 'ISG15', 'CCNI', 'CLU', 'TIMP3', 'ARL4A', 'SERPINH1', 'SCGB1D2', 'UGDH', 'FUS', 'BAG1', 'IFRD1', 'TFF1', 'SERTAD3', 'IGFBP4', 'TPM1', 'PKIB', 'MALAT1', 'XBP1', 'HEBP2', 'GEM', 'EGR2', 'ID2', 'EGR3', 'HSPD1', 'GLUL', 'DDIT4', 'CDC42EP1', 'RBM39', 'MT-ND5', 'CSNK1A1', 'SLC25A25', 'PEG10', 'DEDD2')

gene_module[["gene_module_2"]]<-c('AZGP1', 'ATP5C1', 'ATP5F1', 'NHP2', 'MGP', 'RPN2', 'C14orf2', 'NQO1', 'REEP5', 'SSR2', 'NDUFA8', 'ATP5E', 'SH3BGRL', 'PIP', 'PRDX2', 'RAB25', 'EIF3L', 'PRDX1', 'USMG5', 'DAD1', 'SEC61G', 'CCT3', 'NDUFA4', 'APOD', 'CHCHD10', 'DDIT4', 'MRPL24', 'NME1', 'DCXR', 'NDUFAB1', 'ATP5A1', 'ATP5B', 'ATOX1', 'SLC50A1', 'POLR2I', 'TIMM8B', 'VPS29', 'TIMP1', 'AHCY', 'PRDX3', 'RBM3', 'GSTM3', 'ABRACL', 'RBX1', 'PAFAH1B3', 'AP1S1', 'RPL34', 'ATPIF1', 'PGD', 'CANX', 'SELENBP1', 'ATP5J', 'PSME2', 'PSME1', 'SDHC', 'AKR1A1', 'GSTP1', 'RARRES3', 'ISCU', 'NPM1', 'SPDEF', 'BLVRB', 'NDUFB3', 'RPL36A', 'MDH1', 'MYEOV2', 'MAGED2', 'CRIP2', 'SEC11C', 'CD151', 'COPE', 'PFN2', 'ALDH2', 'SNRPD2', 'TSTD1', 'RPL13A', 'HIGD2A', 'NDUFC1', 'PYCARD', 'FIS1', 'ITM2B', 'PSMB3', 'G6PD', 'CST3', 'SH3BGRL3', 'TAGLN2', 'NDUFA1', 'TMEM183A', 'S100A10', 'NGFRAP1', 'DEGS2', 'ARPC5', 'TM7SF2', 'RPS10', 'LAMTOR5', 'TMEM256', 'UQCRB', 'TMEM141', 'KRTCAP2', 'HM13', 'NDUFS6', 'PARK7', 'PSMD4', 'NDUFB11', 'TOMM7', 'EIF6', 'UQCRHL', 'ADI1', 'VDAC1', 'C9orf16', 'ETFA', 'LSM3', 'UQCRH', 'CYB5A', 'SNRPE', 'BSG', 'SSR3', 'DPM3', 'LAMTOR4', 'RPS11', 'FAM195A', 'TMEM261', 'ATP5I', 'EIF5A', 'PIN4', 'ATXN10', 'ATP5G3', 'ARPC3', 'UBA52', 'BEX4', 'ROMO1', 'SLC25A6', 'SDCBP', 'EIF4EBP1', 'PFDN6', 'PSMA3', 'RNF7', 'SPCS2', 'CYSTM1', 'CAPG', 'CD9', 'GRHPR', 'SEPP1', 'ESF1', 'TFF3', 'ARPC1B', 'ANXA5', 'WDR83OS', 'LYPLA1', 'COMT', 'MDH2', 'DNPH1', 'RAB13', 'EIF3K', 'PTGR1', 'LGALS3', 'TPI1', 'COPZ1', 'LDHA', 'PSMD8', 'EIF2S3', 'NME3', 'EIF3E', 'MRPL13', 'ZFAND6', 'FAM162A', 'ATP6V0E1', 'TMED10', 'HNRNPA3', 'PPA1', 'SNX17', 'APOA1BP', 'TUFM', 'ECHS1', 'GLTSCR2', 'RPS27L', 'NDUFB1', 'SSBP1', 'PRDX6', 'ENO1', 'PPP4C', 'COA3', 'TCEAL4', 'MRPL54', 'LAMTOR2', 'PAIP2', 'DAP', 'RPL22L1', 'C6orf203', 'TECR', 'PEBP1', 'TMED9', 'ATP6V1F', 'ESD', 'EIF3I', 'SCO2', 'ATP5D', 'UAP1', 'TMEM258', 'COX17')

gene_module[["gene_module_3"]]<-c('HLA-B', 'HLA-A', 'VIM', 'CD74', 'SRGN', 'HLA-C', 'IFI27', 'HLA-E', 'IFITM1', 'PSMB9', 'RGCC', 'S100A4', 'HLA-DRA', 'ISG15', 'IL32', 'SPARC', 'TAGLN', 'IFITM3', 'IFITM2', 'IGFBP7', 'CALD1', 'HLA-DPB1', 'HLA-DPA1', 'B2M', 'TIMP1', 'RGS1', 'FN1', 'ACTA2', 'HLA-DRB1', 'SERPING1', 'ANXA1', 'TPM2', 'TMSB4X', 'CD69', 'CCL4', 'LAPTM5', 'GSN', 'APOE', 'STAT1', 'SPARCL1', 'IFI6', 'DUSP1', 'CXCR4', 'CCL5', 'UBE2L6', 'MYL9', 'SLC2A3', 'BST2', 'CAV1', 'CD52', 'ZFP36L2', 'HLA-DQB1', 'PDLIM1', 'TNFAIP3', 'CORO1A', 'RARRES3', 'TYMP', 'C1S', 'PTRF', 'PSME2', 'CYTIP', 'COL1A1', 'PSMB8', 'NNMT', 'HLA-DQA1', 'DUSP2', 'COL1A2', 'ARHGDIB', 'COL6A2', 'FOS', 'CCL2', 'BGN', 'ID3', 'TUBA1A', 'RAC2', 'LBH', 'HLA-DRB5', 'FCER1G', 'GBP1', 'C1QA', 'COTL1', 'LUM', 'MYL6', 'GBP2', 'BTG1', 'CD37', 'HCST', 'LIMD2', 'IFIT3', 'IL7R', 'PTPRC', 'NKG7', 'FYB', 'TAP1', 'LTB', 'S100A6', 'COL3A1', 'EMP3', 'A2M', 'JUNB', 'TPM1', 'FABP4', 'TXNIP', 'SAT1', 'FXYD5', 'CD3E', 'HLA-DMA', 'CTSC', 'TSC22D3', 'MYL12A', 'CST3', 'CNN2', 'PHLDA1', 'LYZ', 'IFI44L', 'MARCKS', 'ID1', 'DCN', 'TGFBI', 'BIRC3', 'THY1', 'LGALS1', 'GPX1', 'C1QB', 'CD2', 'CST7', 'COL6A3', 'ACAP1', 'IFI16', 'ITM2B', 'POSTN', 'LDHB', 'FLNA', 'FILIP1L', 'CDKN1A', 'IRF1', 'LGALS3', 'SERPINH1', 'EFEMP1', 'PSME1', 'SH3BGRL3', 'IL2RG', 'CD3D', 'SFRP2', 'TIMP3', 'ALOX5AP', 'GMFG', 'CYBA', 'TAGLN2', 'LAP3', 'RGS2', 'CLEC2B', 'TRBC2', 'NR4A2', 'S100A8', 'PSMB10', 'OPTN', 'CTSB', 'FTL', 'KRT17', 'AREG', 'MYH9', 'MMP7', 'COL6A1', 'GZMA', 'RNASE1', 'PCOLCE', 'PTN', 'PYCARD', 'ARPC2', 'SGK1', 'COL18A1', 'GSTP1', 'NPC2', 'SOD3', 'MFGE8', 'COL4A1', 'ADIRF', 'HLA-F', 'CD7', 'APOC1', 'TYROBP', 'C1QC', 'TAPBP', 'STK4', 'RHOH', 'RNF213', 'SOD2', 'TPM4', 'CALM1', 'CTGF', 'PNRC1', 'CD27', 'CD3G', 'PRKCDBP', 'PARP14', 'IGKC', 'IGFBP5', 'IFIT1', 'LY6E')

gene_module[["gene_module_4"]]<-c('STMN1', 'H2AFZ', 'UBE2C', 'TUBA1B', 'BIRC5', 'HMGB2', 'ZWINT', 'TUBB', 'HMGB1', 'DEK', 'CDK1', 'HMGN2', 'UBE2T', 'TK1', 'RRM2', 'RANBP1', 'TYMS', 'CENPW', 'MAD2L1', 'CKS2', 'CKS1B', 'NUSAP1', 'TUBA1C', 'PTTG1', 'KPNA2', 'PCNA', 'CENPF', 'HIST1H4C', 'CDKN3', 'UBE2S', 'CCNB1', 'HMGA1', 'DTYMK', 'SNRPB', 'CDC20', 'NASP', 'MCM7', 'PLP2', 'TUBB4B', 'PLK1', 'CCNB2', 'MKI67', 'TOP2A', 'TPX2', 'PKMYT1', 'PRC1', 'SMC4', 'CENPU', 'RAN', 'DUT', 'PA2G4', 'BUB3', 'RAD21', 'SPC25', 'HN1', 'CDCA3', 'H2AFV', 'HNRNPA2B1', 'CCNA2', 'PBK', 'LSM5', 'DNAJC9', 'RPA3', 'TMPO', 'SNRPD1', 'CENPA', 'KIF20B', 'USP1', 'H2AFX', 'PPM1G', 'NUF2', 'SNRPG', 'KIF22', 'KIAA0101', 'DEPDC1', 'RNASEH2A', 'MT2A', 'STRA13', 'ANLN', 'CACYBP', 'NCL', 'NUDT1', 'ECT2', 'LSM4', 'ASF1B', 'CENPN', 'TMEM106C', 'CCT5', 'HSPA8', 'HMMR', 'SRSF3', 'AURKB', 'GGH', 'AURKA', 'TRIP13', 'CDCA8', 'HMGB3', 'HNRNPAB', 'FAM83D', 'CDC25B', 'GGCT', 'KNSTRN', 'CCT6A', 'PTGES3', 'ANP32E', 'CENPK', 'MCM3', 'DDX21', 'HSPD1', 'SKA2', 'CALM2', 'UHRF1', 'HINT1', 'ORC6', 'MZT1', 'MIS18BP1', 'WDR34', 'NAP1L1', 'TEX30', 'SFN', 'HSPE1', 'CENPM', 'TROAP', 'CDCA5', 'RACGAP1', 'SLC25A5', 'ATAD2', 'DBF4', 'KIF23', 'CEP55', 'SIVA1', 'SAC3D1', 'PSIP1', 'CLSPN', 'CCT2', 'DLGAP5', 'PSMA4', 'SMC2', 'AP2S1', 'RAD51AP1', 'MND1', 'ILF2', 'DNMT1', 'NUCKS1', 'LMNB1', 'RFC4', 'EIF5A', 'NPM3', 'ARL6IP1', 'ASPM', 'GTSE1', 'TOMM40', 'HNRNPA1', 'GMNN', 'FEN1', 'CDCA7', 'SLBP', 'TNFRSF12A', 'TM4SF1', 'CKAP2', 'CENPE', 'SRP9', 'DDX39A', 'COMMD4', 'RBM8A', 'CALM3', 'RRM1', 'ENO1', 'ANP32B', 'SRSF7', 'FAM96A', 'TPRKB', 'FABP5', 'PPIF', 'SERPINE1', 'TACC3', 'RBBP7', 'NEK2', 'CALM1', 'GMPS', 'EMP2', 'HMG20B', 'SMC3', 'HSPA9', 'NAA20', 'NUDC', 'RPL39L', 'PRKDC', 'CDCA4', 'HIST1H1A', 'HES6', 'SUPT16H', 'PTMS', 'VDAC3', 'PSMC3', 'ATP5G1', 'PSMA3', 'PGP', 'KIF2C', 'CARHSP1')

gene_module[["gene_module_5"]]<-c('GJA1', 'SCGB2A2', 'ARMT1', 'MAGED2', 'PIP', 'SCGB1D2', 'CLTC', 'MYBPC1', 'PDZK1', 'MGP', 'SLC39A6', 'CCND1', 'SLC9A3R1', 'NAT1', 'SUB1', 'CYP4X1', 'STC2', 'CROT', 'CTSD', 'FASN', 'PBX1', 'SLC4A7', 'FOXA1', 'MCCC2', 'IDH1', 'H2AFJ', 'CYP4Z1', 'IFI27', 'TBC1D9', 'ANPEP', 'DHRS2', 'TFF3', 'LGALS3BP', 'GATA3', 'LTF', 'IFITM2', 'IFITM1', 'AHNAK', 'SEPP1', 'ACADSB', 'PDCD4', 'MUCL1', 'CERS6', 'LRRC26', 'ASS1', 'SEMA3C', 'APLP2', 'AMFR', 'CDV3', 'VTCN1', 'PREX1', 'TP53INP1', 'LRIG1', 'ANK3', 'ACLY', 'CLSTN1', 'GNB1', 'C1orf64', 'STARD10', 'CA12', 'SCGB2A1', 'MGST1', 'PSAP', 'GNAS', 'MRPS30', 'MSMB', 'DDIT4', 'TTC36', 'S100A1', 'FAM208B', 'STT3B', 'SLC38A1', 'DMKN', 'SEC14L2', 'FMO5', 'DCAF10', 'WFDC2', 'GFRA1', 'LDLRAD4', 'TXNIP', 'SCGB3A1', 'APOD', 'N4BP2L2', 'TNC', 'ADIRF', 'NPY1R', 'NBPF1', 'TMEM176A', 'GLUL', 'BMP2K', 'SLC44A1', 'GFPT1', 'PSD3', 'CCNG2', 'CGNL1', 'TMED7', 'NOVA1', 'ARCN1', 'NEK10', 'GPC6', 'SCGB1B2P', 'IGHG4', 'SYT1', 'SYNGR2', 'HSPA1A', 'ATP6AP1', 'TSPAN13', 'MT-ND2', 'NIFK', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'EVL', 'GRN', 'ERH', 'CD81', 'NUPR1', 'SELENBP1', 'C1orf56', 'LMO3', 'PLK2', 'HACD3', 'RBBP8', 'CANX', 'ENAH', 'SCD', 'CREB3L2', 'SYNCRIP', 'TBL1XR1', 'DDR1', 'ERBB3', 'CHPT1', 'BANF1', 'UGDH', 'SCUBE2', 'UQCR10', 'COX6C', 'ATP5G1', 'PRSS23', 'MYEOV2', 'PITX1', 'MT-ND4L', 'TPM1', 'HMGCS2', 'ADIPOR2', 'UGCG', 'FAM129B', 'TNIP1', 'IFI6', 'CA2', 'ESR1', 'TMBIM4', 'NFIX', 'PDCD6IP', 'CRIM1', 'ARHGEF12', 'ENTPD5', 'PATZ1', 'ZBTB41', 'UCP1', 'ANO1', 'RP11-356O9.1', 'MYB', 'ZBTB44', 'SCPEP1', 'HIPK2', 'CDK2AP1', 'CYHR1', 'SPINK8', 'FKBP10', 'ISOC1', 'CD59', 'RAMP1', 'AFF3', 'MT-CYB', 'PPP1CB', 'PKM', 'ALDH2', 'PRSS8', 'NPW', 'SPR', 'PRDX3', 'SCOC', 'TMED10', 'KIAA0196', 'NDP', 'ZSWIM7', 'AP2A1', 'PLAT', 'SUSD3', 'CRABP2', 'DNAJC12', 'DHCR24', 'PPT1', 'FAM234B', 'DDX17', 'LRP2', 'ABCD3', 'CDH1', 'NFIA') 

gene_module[["gene_module_6"]]<-c('AGR2', 'TFF3', 'SELM', 'CD63', 'CTSD', 'MDK', 'CD74', 'S100A13', 'IFITM3', 'HLA-B', 'AZGP1', 'FXYD3', 'IFITM2', 'RABAC1', 'S100A14', 'CRABP2', 'LTF', 'RARRES1', 'HLA-A', 'PPIB', 'HLA-C', 'S100A10', 'S100A9', 'TIMP1', 'DDIT4', 'S100A16', 'LGALS1', 'LAPTM4A', 'SSR4', 'S100A6', 'CD59', 'BST2', 'PDIA3', 'KRT19', 'CD9', 'FXYD5', 'SCGB2A2', 'NUCB2', 'TMED3', 'LY6E', 'CFD', 'ITM2B', 'PDZK1IP1', 'LGALS3', 'NUPR1', 'SLPI', 'CLU', 'TMED9', 'HLA-DRA', 'SPTSSB', 'TMEM59', 'KRT8', 'CALR', 'HLA-DRB1', 'IFI6', 'NNMT', 'CALML5', 'S100P', 'TFF1', 'ATP1B1', 'SPINT2', 'PDIA6', 'S100A8', 'HSP90B1', 'LMAN1', 'RARRES3', 'SELENBP1', 'CEACAM6', 'TMEM176A', 'EPCAM', 'MAGED2', 'SNCG', 'DUSP4', 'CD24', 'PERP', 'WFDC2', 'HM13', 'TMBIM6', 'C12orf57', 'DKK1', 'MAGED1', 'PYCARD', 'RAMP1', 'C11orf31', 'STOM', 'TNFSF10', 'BSG', 'TMED10', 'ASS1', 'PDLIM1', 'CST3', 'PDIA4', 'NDUFA4', 'GSTP1', 'TYMP', 'SH3BGRL3', 'PRSS23', 'P4HA1', 'MUC5B', 'S100A1', 'PSAP', 'TAGLN2', 'MGST3', 'PRDX5', 'SMIM22', 'NPC2', 'MESP1', 'MYDGF', 'ASAH1', 'APP', 'NGFRAP1', 'TMEM176B', 'C8orf4', 'KRT81', 'VIMP', 'CXCL17', 'MUC1', 'COMMD6', 'TSPAN13', 'TFPI', 'C15orf48', 'CD151', 'TACSTD2', 'PSME2', 'CLDN7', 'ATP6AP2', 'CUTA', 'MT2A', 'CYB5A', 'CD164', 'TM4SF1', 'SCGB1D2', 'GSTM3', 'EGLN3', 'LMAN2', 'IFI27', 'PPP1R1B', 'B2M', 'ANXA2', 'SARAF', 'MUCL1', 'CSRP1', 'NPW', 'SLC3A2', 'PYDC1', 'QSOX1', 'TSPAN1', 'GPX1', 'TMSB4X', 'FGG', 'GUK1', 'IL32', 'ATP6V0E1', 'BCAP31', 'CHCHD10', 'TSPO', 'TNFRSF12A', 'MT1X', 'PDE4B', 'HSPA5', 'SCD', 'SERINC2', 'PSCA', 'VAMP8', 'ELF3', 'TSC22D3', 'S100A7', 'GLUL', 'ZG16B', 'TMEM45A', 'APMAP', 'RPS26', 'CALU', 'OSTC', 'NCCRP1', 'SQLE', 'RPS28', 'SSR2', 'SOX4', 'CLEC3A', 'TMEM9', 'RPL10', 'MUC5AC', 'HLA-DPA1', 'ZNHIT1', 'AQP5', 'CAPG', 'SPINT1', 'NDFIP1', 'FKBP2', 'C1S', 'LDHA', 'NEAT1', 'RPL36A', 'S100A11', 'LCN2', 'TUBA1A', 'GSTK1', 'SEPW1', 'P4HB') 

gene_module[["gene_module_7"]]<-c('KCNQ1OT1', 'AKAP9', 'RHOB', 'SOX4', 'VEGFA', 'CCNL1', 'RSRP1', 'RRBP1', 'ELF3', 'H1FX', 'FUS', 'NEAT1', 'N4BP2L2', 'SLC38A2', 'BRD2', 'PNISR', 'CLDN4', 'MALAT1', 'SOX9', 'DDIT3', 'TAF1D', 'FOSB', 'ZNF83', 'ARGLU1', 'DSC2', 'MACF1', 'GTF2I', 'SEPP1', 'ANKRD30A', 'PRLR', 'MAFB', 'NFIA', 'ZFAS1', 'MTRNR2L12', 'RNMT', 'NUPR1', 'MT-ND6', 'RBM39', 'HSPA1A', 'HSPA1B', 'RGS16', 'SUCO', 'XIST', 'PDIA6', 'VMP1', 'SUGP2', 'LPIN1', 'NDRG1', 'PRRC2C', 'CELF1', 'HSP90B1', 'JUND', 'ACADVL', 'PTPRF', 'LMAN1', 'HEBP2', 'ATF3', 'BTG1', 'GNAS', 'TSPYL2', 'ZFP36L2', 'RHOBTB3', 'TFAP2A', 'RAB6A', 'KMT2C', 'POLR2J3', 'CTNND1', 'PRRC2B', 'RNF43', 'CAV1', 'RSPO3', 'IMPA2', 'FAM84A', 'FOS', 'IGFBP5', 'NCOA3', 'WSB1', 'MBNL2', 'MMP24-AS1', 'DDX5', 'AP000769.1', 'MIA3', 'ID2', 'HNRNPH1', 'FKBP2', 'SEL1L', 'PSAT1', 'ASNS', 'SLC3A2', 'EIF4EBP1', 'HSPH1', 'SNHG19', 'RNF19A', 'GRHL1', 'WBP1', 'SRRM2', 'RUNX1', 'ASH1L', 'HIST1H4C', 'RBM25', 'ZNF292', 'RNF213', 'PRPF38B', 'DSP', 'EPC1', 'FNBP4', 'ETV6', 'SPAG9', 'SIAH2', 'RBM33', 'CAND1', 'CEBPB', 'CD44', 'NOC2L', 'LY6E', 'ANGPTL4', 'GABPB1-AS1', 'MTSS1', 'DDX42', 'PIK3C2G', 'IAH1', 'ATL2', 'ADAM17', 'PHIP', 'MPZ', 'CYP27A1', 'IER2', 'ACTR3B', 'PDCD4', 'COLCA1', 'KIAA1324', 'TFAP2C', 'CTSC', 'MYC', 'MT1X', 'VIMP', 'SERHL2', 'YPEL3', 'MKNK2', 'ZNF552', 'CDH1', 'LUC7L3', 'DDIT4', 'HNRNPR', 'IFRD1', 'RASSF7', 'SNHG8', 'EPB41L4A-AS1', 'ZC3H11A', 'SNHG15', 'CREB3L2', 'ERBB3', 'THUMPD3-AS1', 'RBBP6', 'GPBP1', 'NARF', 'SNRNP70', 'RP11-290D2.6', 'SAT1', 'GRB7', 'H1F0', 'EDEM3', 'KIAA0907', 'ATF4', 'DNAJC3', 'DKK1', 'SF1', 'NAMPT', 'SETD5', 'DYNC1H1', 'GOLGB1', 'C4orf48', 'CLIC3', 'TECR', 'HOOK3', 'WDR60', 'TMEM101', 'SYCP2', 'C6orf62', 'METTL12', 'HIST1H2BG', 'PCMTD1', 'PWWP2A', 'HIST1H3H', 'NCK1', 'CRACR2B', 'NPW', 'RAB3GAP1', 'TMEM63A', 'MGP', 'ANKRD17', 'CALD1', 'PRKAR1A', 'PBX1', 'ATXN2L', 'FAM120A', 'SAT2', 'TAF10', 'SFRP1', 'CITED2') 

#maybe add Ecotypes as well, no given gene lists from the publication??



sample_cell_signature_transfer<-function(){
  dat_epi<-subset(dat,EMBO_predicted.id=="epithelial")

  #embo lineage
  dat_epi<-AddModuleScore(dat_epi,features=lineage_in,names=names(lineage_in),assay="SoupXRNA",seed=123,search=TRUE)

  colnames(dat_epi@meta.data)[which(colnames(dat_epi@meta.data) %in% c("Cluster1","Cluster2","Cluster3","Cluster4"))]<-c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str") #Rename them

  #Immune cell features and PAM50 canonical short list of genes
  for(i in 1:length(features_in)){
    features_in[[i]]<-features_in[[i]][features_in[[i]] %in% row.names(dat_epi[["SoupXRNA"]])] #make sure gene names match
    dat_epi<-MetaFeature(dat_epi,features=c(features_in[[i]]),meta.name=names(features_in)[i],assay="SoupXRNA")}

  #SCSubype List of genes
  #run only on epithelial cells
  module_scores<-AddModuleScore(dat_epi,features=module_feats,assay="SoupXRNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
  module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-(length(module_feats)-1),ncol(module_scores@meta.data))]
  colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like
  dat_epi<-AddMetaData(dat,metadata=module_scores)

  #Swarbrick Gene Modules
  #run only on epithelial cells
  gene_module_out<-AddModuleScore(dat_epi,features=gene_module,assay="SoupXRNA",search=TRUE,name=names(gene_module)) #use add module function to add cell scores
  gene_module_out<-gene_module_out@meta.data[seq(ncol(gene_module_out@meta.data)-(length(gene_module)-1),ncol(gene_module_out@meta.data))]#get the 7 added gene modules
  colnames(gene_module_out)<-names(gene_module) 
  dat_epi<-AddMetaData(dat_epi,metadata=gene_module_out)
  out<-dat_epi@meta.data[c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str",names(module_feats),names(gene_module))]
  return(out)
}

single_sample_PAM50_assignment<-function(){
  met<-dat@meta.data
  met<-met[met$EMBO_predicted.id %in% c("epithelial"),]
  pam50_list<-  c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str" )
  max_pam50<-lapply(1:nrow(met),function(i) pam50_list[which(met[i,pam50_list]==max(met[i,pam50_list],na.rm=T))])
  max_pam50<-unlist(lapply(1:length(max_pam50),function(i) do.call("paste",as.list(max_pam50[[i]]))))
  max_pam50<-unlist(lapply(max_pam50,function(i) gsub("EMBO_","",i)))
  names(max_pam50)<-row.names(met)
  return(max_pam50)
}

single_sample_SCtype_assignment<-function(){
  met<-dat@meta.data
  met<-met[met$EMBO_predicted.id %in% c("epithelial"),]
  scsubtype_list<-  c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")
  max_scsubtype<-lapply(1:nrow(met),function(i) scsubtype_list[which(met[i,scsubtype_list]==max(met[i,scsubtype_list],na.rm=T))])
  max_scsubtype<-unlist(lapply(1:length(max_scsubtype),function(i) do.call("paste",as.list(max_scsubtype[[i]]))))
  names(max_scsubtype)<-row.names(met)
  return(max_scsubtype)
}

#Generate scores per epithelial cell
epithelial_metadata<-sample_cell_signature_transfer()
dat<-AddMetaData(dat,metadata=epithelial_metadata) #add to master data frame metadata

#assign top PAM50 designation by epithelial cells
max_pam50<-single_sample_SCtype_assignment()
dat<-AddMetaData(dat,max_pam50,col.name="PAM50_epi_designation")

#assign top scsubtype by epithelial cells
max_scsubtype<-single_sample_PAM50_assignment()
dat<-AddMetaData(dat,max_scsubtype,col.name="SCSubtype_epi_designation")

saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```
Plot Features on Cells Per Sample

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

single_sample_epcam<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  DefaultAssay(dat)<-"SoupXRNA"
  plt<-FeaturePlot(dat,features="EPCAM",reduction="multimodal_umap",order=T)
  ggsave(plt,file=out_plot,width=10,height=10)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  plt<-VlnPlot(dat,features="EPCAM",group.by="predicted.id")
  ggsave(plt,file=paste0(out_plot,"VlnPlt.pdf"),width=10,height=10)
  system(paste0("slack -F ",paste0(out_plot,"VlnPlt.pdf")," ryan_todo"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_epcam)

```

## Run IntClust on Samples
<!-- Rerun  waiting on final CNV profiles-->

Using iC10 CRAN Package https://cran.r-project.org/web/packages/iC10/iC10.pdf

```R
#install.packages("iC10")
library(iC10)
library(Seurat)
library(Signac)
#CN = ID (Sample Name) \t chromosome_name (Chr) \t loc.start (start) loc.end (end) seg.mean (log2ratio of segment)
#OR
#CN = Row (hgnc gene names) X Column (Sample)
#Exp =  Row (hgnc gene names) X Column (Sample)

#using InferCNV(gene level CNV calling) as CN matrix, and RNA data as Exp Matrix

iC10_per_sample<-function(x){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  Idents(dat)<-dat$predicted.id
  dat_ep<-subset(dat, cells = row.names(dat@meta.data[dat@meta.data$predicted.id%in% c("Normal Epithelial", "Cancer Epithelial"),]))
  infercnv_obj<-readRDS(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
  cnv<-log2(infercnv_obj@expr.data)
  cnv<-cnv[,colnames(cnv) %in% colnames(dat_ep)]
  exp<-dat_ep[["RNA"]]@counts

  out<-matchFeatures(CN=cnv,Exp=exp,
    CN.by.feat="gene",
    Exp.by.feat="gene",
    ref=NULL)
  out<-normalizeFeatures(out, method="scale")
  out<-iC10(out)
  saveRDS(out,paste0(wd,"/",outname,"_iC10.Rds"))
  dat<-AddMetaData(dat,out$class,col.name="ic10_class")
  table(dat$ic10_class)
  saveRDS(dat,file=file_in) #overwrite old file
  print(paste("Finished Sample:",sample_name))
}

lapply(c(1,3,5,6,7,8,9,16,19,20,"RM_1","RM_2","RM_3",11,4,10,12), function(x) iC10_per_sample(x))
#done 
#15,"RM_4" not done
```

## Epithelial Subtyping Per Sample
<!-- Rerun once ic50 is done -->

```R
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)

  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
  alpha_val=0.33

epithelial_class_persample<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  dat<-readRDS(file=file_in)
  atac_sub<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial"))
  if(!("ic10_class" %in% colnames(atac_sub@meta.data))){
    atac_sub$ic10_class<-"NA"
  }
  plt_cell_count<-atac_sub@meta.data[,c("sample","predicted.id","diagnosis","molecular_type","PAM50_designation","EMBO_designation","ic10_class")]
  print(outname)
  return(plt_cell_count)
}


#grab all epithelial classifications
cell_count<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x)epithelial_class_persample(x))

saveRDS(cell_count,"sample_epithelial_designations.rds") #save nested list of cell type assignment

#plot output of celltype count per sample
out<-readRDS("sample_epithelial_designations.rds")
out<-do.call("rbind",out)
colnames(out)<-c("sample","predicted.id","diagnosis","molecular_type","SCSubtype_designation","EMBO_designation","ic10_class") #rename just for simplicity
#clean up for samples with equal values
out[!(out$SCSubtype_designation %in% c("Basal","Her2","LumA","LumB","Normal")),]$SCSubtype_designation<-NA
  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
plt1<-ggplot(out,aes(x=sample,fill=SCSubtype_designation))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt2<-ggplot(out,aes(x=sample,fill=EMBO_designation))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt3<-ggplot(out,aes(x=sample,fill=ic10_class))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt<-plt1/plt2/plt3
ggsave(plt,file="sample_epithelial_type_assignment.pdf")
system("slack -F sample_epithelial_type_assignment.pdf ryan_todo")

library(dplyr)
write.table(out,file="sample_epithelial_type_assignment.tsv",col.names=T,row.names=F,sep="\t",quote=F)
system("slack -F sample_epithelial_type_assignment.tsv ryan_todo") #note this was calculated per sample as well as in the merged data set,  the assumption is that they will be the same


```

<!-- Rerun -->
## ER binding poor and good outcome from patients, overlap with ATAC data
```bash
cd /home/groups/CEDAR/mulqueen/ref
wget http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Poor%20outcome%20ER%20regions.bed.gz
wget http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Good%20outcome%20ER%20regions.bed.gz
mv "Good outcome ER regions.bed.gz" goodoutcome_ER_regions.bed.gz
mv "Poor outcome ER regions.bed.gz" pooroutcome_ER_regions.bed.gz 

```
### Files for Travis
<!-- Rerun -->

Transfering to /home/groups/CEDAR/scATACcnv/Hisham_data/final_data:
- Metadata
- Counts matrix (atac)
- peaks bed file
- inferCNV folder
- Casper folder
- atac_possorted_bam.bam
- Seurat Object

```R
library(Signac)
library(Seurat)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat_all<-readRDS("phase2.QC.filt.SeuratObject.rds")
library(parallel)

transfer_data<-function(x){
  print(paste("Running Sample:",x))
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    dat<-subset(dat_all, sample==sample_name)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    dat<-subset(dat_all, sample==sample_name)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    dat<-subset(dat_all, sample==sample_name)
  }
  metadata<-as.data.frame(dat@meta.data)
  counts_matrix<-as.data.frame(dat[["peaks"]]@counts)
  peaks_bed_file<-as.data.frame(do.call("rbind",strsplit(row.names(dat[["peaks"]]),"-")))
  write.table(metadata,file=paste0(wd,"/","metadata.tsv"),sep="\t",col.names=T,row.names=T)
  write.table(counts_matrix,file=paste0(wd,"/","counts_matrix.tsv"),sep="\t",col.names=T,row.names=T)
  write.table(peaks_bed_file,file=paste0(wd,"/","peaks.bed"),sep="\t",col.names=F,row.names=F)
  print(paste("Finished sample:",sample_name))
  saveRDS(dat,file=paste0(wd,"/",sample_name,".QC.filt.SeuratObject.rds"))
}

mclapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x) transfer_data(x),mc.cores=6)
```

- Metadata
- Counts matrix (atac)
- peaks bed file
- inferCNV folder
- Casper folder
- atac_possorted_bam.bam
- Seurat Object

```bash
out_dir="/home/groups/CEDAR/scATACcnv/Hisham_data/final_data"
mkdir $out_dir


for i in 1 3 4 5 6 7 8 9 10 11 12; do
  sample="sample_"${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_"${i}"/outs"
  mkdir ${out_dir}"/"${sample}; 
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/sample_"${i} & 
  cp ${in_dir}"/peaks.bed" ${out_dir}"/sample_"${i} & 
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/sample_"${i} & 
  cp -r ${in_dir}"/casper" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/copykat" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/copyscat" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/sample_"${i} & 
  cp ${in_dir}"/"${sample}".QC.filt.SeuratObject.rds" ${out_dir}"/sample_"${i} &
  echo "Finished ${sample}" & done &


for i in 15 16 19 20; do
  sample="sample_"${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_"${i}"/outs"
  mkdir ${out_dir}"/"${sample}; 
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/peaks.bed" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/casper" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/copykat" ${out_dir}"/sample_"${i} &
  cp -r ${in_dir}"/copyscat" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/sample_"${i} &
  cp ${in_dir}"/"${sample}".QC.filt.SeuratObject.rds" ${out_dir}"/sample_"${i} &
  echo "Finished ${sample}" & done 


for i in "RM_1" "RM_2" "RM_3" "RM_4"; do
  sample=${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/"${i}"/outs"
  mkdir ${out_dir}"/"${sample};
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/"${sample} &
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/"${sample} &
  cp ${in_dir}"/peaks.bed" ${out_dir}"/"${sample} &
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/"${sample} &
  cp -r ${in_dir}"/casper" ${out_dir}"/"${sample} &
  cp -r ${in_dir}"/copykat" ${out_dir}"/"${sample} &
  cp -r ${in_dir}"/copyscat" ${out_dir}"/"${sample} &
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/"${sample} &
  cp ${in_dir}"/"${sample}".QC.filt.SeuratObject.rds" ${out_dir}"/"${sample} &
  echo "Finished ${sample}" & done 


```

Call peaks per cell type
<!-- Rerun -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")


dat<-readRDS("phase2.QC.SeuratObject.rds")

#call peaks per predicted.id cell type
peaks <- CallPeaks(dat, 
  assay="ATAC",
  group.by="predicted.id",
  combine.peaks=FALSE,
  macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")
  #use this set of peaks for all samples

for(i in 1:length(peaks)){
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peak_name<-unique(dat$predicted.id)[i]
  peaks_out <- keepStandardChromosomes(peaks[[i]], pruning.mode = "coarse")
  peaks_out <- subsetByOverlaps(x = peaks_out, ranges = blacklist_hg38_unified, invert = TRUE)
  print(paste0("Generated peakset for ",peak_name))
  write.table(as.data.frame(peaks_out)[1:3],file=paste0(peak_name,".bed"),sep="\t",quote=F,col.names=F,row.names=F)
}

```

```bash
cp *bed /home/groups/CEDAR/scATACcnv/Hisham_data/final_data
```

## Genome tracks of celltype markers
<!-- Rerun -->

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(SeuratObjects)
library(EnsDb.Hsapiens.v86)
library(cowplot)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

x<-"ESR1"
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
Idents(dat)<-dat$predicted.id 

cov_plots<-function(dat=dat,gene_name="ESR1",outname="test",extend=2000){
  
  #get chromvar motif name
  gene_TF <- ConvertMotifID(dat, name = gene_name,assay="ATAC")

  # get gene location
  gene_loc<-annotation[annotation$gene_name==gene_name,]
  gene_loc<-as.data.frame(gene_loc,row.names=NULL)
  gene_loc<-head(gene_loc[which(max(gene_loc$end-gene_loc$start)==c(gene_loc$end-gene_loc$start)),],n=1)
  if(gene_loc$start<gene_loc$end){
  gene_range<-paste0("chr",as.character(gene_loc$seqnames),"-",gene_loc$start-extend,"-",gene_loc$end+extend)
  }else{
  gene_range<-paste0("chr",as.character(gene_loc$seqnames),"-",gene_loc$start+extend,"-",gene_loc$end-extend)
  }
  annot_plot<-AnnotationPlot(object=dat, region=gene_range)
  peak_plot<-PeakPlot(object=dat,region=gene_range)

  cov_plot <- CoveragePlot(object = dat, region = gene_range, assay="peaks", annotation=FALSE,peaks=FALSE,links=FALSE) 
  plt_rna<-VlnPlot(dat,features=gene_name,assay="RNA",slot="counts",flip=TRUE,pt.size=0)+coord_flip()+scale_x_discrete(limits = rev(levels(Idents(dat))))
  plt_tf<-VlnPlot(dat,features=gene_TF,assay="chromvar",flip=TRUE,pt.size=0)+coord_flip()+scale_x_discrete(limits = rev(levels(Idents(dat))))

  layout<-"
  AAAAAA##
  EEEEEE##
  BBBBBBCD
  BBBBBBCD
  BBBBBBCD
  "

  plt<-wrap_plots(
    A=annot_plot,
    B=cov_plot,
    C=plt_rna,
    D=plt_tf,
    E=peak_plot,
    design=layout,heights = c(1,3,1,2,3,1,2),guides="collect")+ggtitle(outname)
  return(plt)
}


#Markers with high cell type AUC determined in section Transcription Factor Expression Markers
  prefix="dat"
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
  da_tf_markers$gene

lapply(c(da_tf_markers$gene),function(y){
  plot<-cov_plots(dat=dat,gene_name=y,outname=y)
  ggsave(plot,file=paste0("merged_coverageplt_",y,"celltype.pdf"),width=15,height=10)
  system(paste0("slack -F ",paste0("merged_coverageplt_",y,"celltype.pdf")," ryan_todo"))
  })
```

## Subclustering Per Celltype
<!-- Done -->

### Epithelial Clustering
```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cisTopic)
library(SeuratWrappers)
library(patchwork)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
library(RColorBrewer)
library(ggplot2)
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.filt.SeuratObject.rds")


#set up colors for samples
###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    
molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
########################################
alpha_val=0.33


#
celltype_cistopic_generation<-function(celltype_list=c("Cancer Epithelial","Normal Epithelial"),outname="epithelial"){
  atac_sub<-subset(dat,EMBO_predicted.id %in% celltype_list)

  cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
  saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))

  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  
  saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(outname,".CisTopicObject.Rds"))
  print("finshed running cistopic")

  #Add cell embeddings into seurat
  cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
  colnames(cell_embeddings)<-sub_cistopic_models@cell.names
  n_topics<-nrow(cell_embeddings)
  row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
  cell_embeddings<-as.data.frame(t(cell_embeddings))

  #Add feature loadings into seurat
  feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
  row.names(feature_loadings)<-paste0("topic_",1:n_topics)
  feature_loadings<-as.data.frame(t(feature_loadings))

  #combined cistopic results (cistopic loadings and umap with seurat object)
  cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
  print("Cistopic Loading into Seurat")
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  pdf(paste0(outname,".cistopic.umap.pdf"),width=10)
  print(plt1)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistopic.umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(outname,".SeuratObject.rds"))
  }

#Rerun other clustering now that data is subset
celltype_clustering<-function(x,outname){
  dat<-readRDS(x)

  #set up colors for samples
  my_cols = brewer.pal(1,"Spectral")
  alpha_val=0.33

  #RNA Processing
  DefaultAssay(dat) <- "SoupXRNA"
  dat <- SCTransform(dat)
  dat <- RunPCA(dat)
  dat<- RunUMAP(object = dat, reduction.name="rna_umap", reduction="pca", assay = "SoupXRNA", verbose = TRUE, dims=1:50 ) 
  p1<-DimPlot(dat,reduction="rna_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("RNA UMAP")+theme(legend.position="none")

  #DNA Accessibility processing
  DefaultAssay(dat) <- "peaks"
  dat <- FindTopFeatures(dat, min.cutoff = 5)
  dat <- RunTFIDF(dat)
  dat <- RunSVD(dat)
  dat<- RunUMAP(object = dat, reduction.name="atac_umap", reduction="lsi", assay = "peaks", verbose = TRUE, dims=2:40 )
  p2<-DimPlot(dat,reduction="atac_umap",group.by="EMBO_predicted.id",cols=alpha(embo_cell_cols,alpha_val))+ggtitle("ATAC UMAP")+theme(legend.position="none")


  # build a joint neighbor graph using both assays (ATAC LSI)
    dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:40),modality.weight.name = "RNA.weight", verbose = TRUE )
    # build a joint UMAP visualization
    dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_umap", assay = "SoupXRNA", verbose = TRUE ) 
    p3<-DimPlot(dat,reduction="multimodal_umap",group.by="EMBO_predicted.id",cols=alpha(embo_cell_cols,alpha_val))+ggtitle("Multimodal UMAP Doublets")+theme(legend.position="none")

  #Try multimodal with cistopic
    dat <- RunUMAP(object = dat, reduction="cistopic", reduction.name="cistopic_umap", dims=1:ncol(dat@reductions$cistopic), assay = "peaks", verbose = TRUE )
    # build a joint neighbor graph using both assays
    dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("pca", "cistopic"), dims.list = list(1:50, 1:ncol(dat@reductions$cistopic)), modality.weight.name = "RNA.weight", verbose = TRUE )
    # build a joint UMAP visualization
    dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_umap", assay = "SoupXRNA", verbose = TRUE )

  #plot cistopic umap too
  p4<-DimPlot(dat,reduction="multimodal_umap",group.by="EMBO_predicted.id",cols=alpha(embo_cell_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  p5<-DimPlot(dat,reduction="cistopic_umap",group.by="EMBO_predicted.id",cols=alpha(embo_cell_cols,alpha_val))+ggtitle("Cistopic UMAP")+theme(legend.position="none")
  p6<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  #Cluster on multimodal graph

  dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")


  #Finally Plot results
  plt<-(p1 | p2)/(p3 | p4)/(p5|p6)
  ggsave(plt,file=paste0(outname,".umap.pdf"))
  system(paste0("slack -F ",paste0(outname,".umap.pdf")," ryan_todo"))
  saveRDS(dat,file=paste0(outname,"filt.SeuratObject.rds"))
}
 
#Epithelial Cells
celltype_cistopic_generation(celltype_list=c("epithelial","cycling.epithelial"),outname="epithelial")

#Normal Epithelial Cells
dat<-subset(dat,diagnosis=="NAT")
celltype_cistopic_generation(celltype_list=c("epithelial","cycling.epithelial"),outname="normalepithelial")

#Myeloid
celltype_cistopic_generation(celltype_list=c("TAMs","TAMs_2","Myeloid"),outname="myeloid")

#T Cells
celltype_cistopic_generation(celltype_list=c("T.cells"),outname="tcell")

#Plasma cells
celltype_cistopic_generation(celltype_list=c("Plasma.cells"),outname="plasmacell")

#Fibroblasts
celltype_cistopic_generation(celltype_list=c("CAFs"),outname="fibroblast")

#B cells
celltype_cistopic_generation(celltype_list=c("B.cells"),outname="bcell")

#Endothelial
celltype_cistopic_generation(celltype_list=c("Endothelial"),outname="endothelial")

#Pericytes
celltype_cistopic_generation(celltype_list=c("Pericytes"),outname="pericyte")


```

### Integration: Now Clustering together on RNA profiles using harmony to integrate
<!-- Done -->

```R
library(harmony)
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(patchwork)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")


###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

diag_cols<-c("IDC"="red", "DCIS"="grey","NAT"="lightblue","ILC"="green")
embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    

molecular_type_cols<-c("DCIS"="grey", "ER+/PR-/HER2-"="#EBC258", "ER+/PR-/HER2+"="#F7B7BB","ER+/PR+/HER2-"="#ff6699","NA"="lightblue")
########################################


harmony_sample_integration<-function(x,outname,res=0.1){
  dat<-readRDS(x)
  dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_atac",assay.use="ATAC",reduction="cistopic",project.dim=F)
  dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_rna",assay.use="RNA",reduction="pca",project.dim=F)

  dat<-RunUMAP(dat,reduction.name="harmonyumap_rna",reduction = "harmony_rna",dims=1:dim(dat@reductions$harmony_rna)[2]) 
  dat<-RunUMAP(dat,reduction.name="harmonyumap_atac",reduction = "harmony_atac",dims=1:dim(dat@reductions$harmony_atac)[2]) 

  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(
    object = dat,
    reduction.list = list("harmony_rna", "harmony_atac"), 
    dims.list = list(1:dim(dat@reductions$harmony_rna)[2], 1:dim(dat@reductions$harmony_atac)[2]), 
    modality.weight.name = "multimodal.weight",
    weighted.nn.name="multimodal_harmony.nn",
    verbose = TRUE
  )
  dat <- FindClusters(dat, graph.name="wknn",verbose = FALSE,resolution=res)
  # build a joint UMAP Harmony visualization
  dat <- RunUMAP(object = dat, nn.name = "multimodal_harmony.nn",reduction.name="multimodal_harmony_umap", assay = "SoupXRNA", verbose = TRUE ) 

  i="EMBO_predicted.id"
  plt1<-DimPlot(dat,reduction="multimodal_umap",group.by=i,cols=embo_cell_cols)+ggtitle("Unintegrated")
  plt2<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i,cols=embo_cell_cols)+ggtitle("Integrated")
  plt<-plt1+plt2
  ggsave(plt,file=paste0(outname,".",i,".pdf"),width=20,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  i="diagnosis"
  plt1<-DimPlot(dat,reduction="multimodal_umap",group.by=i,cols=diag_cols)+ggtitle("Unintegrated")
  plt2<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i,cols=diag_cols)+ggtitle("Integrated")
  plt<-plt1+plt2
  ggsave(plt,file=paste0(outname,".",i,".pdf"),width=20,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  i="sample"
  plt1<-DimPlot(dat,reduction="multimodal_umap",group.by=i)+ggtitle("Unintegrated")
  plt2<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i)+ggtitle("Integrated")
  plt<-plt1+plt2
  ggsave(plt,file=paste0(outname,".",i,".pdf"),width=20,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  i="seurat_clusters"
  plt1<-DimPlot(dat,reduction="multimodal_umap",group.by=i)+ggtitle("Unintegrated")
  plt2<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i)+ggtitle("Integrated")
  plt<-plt1+plt2
  ggsave(plt,file=paste0(outname,".",i,".pdf"),width=20,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  saveRDS(dat,file=x)
}

harmony_sample_integration(x="phase2.QC.filt.SeuratObject.rds",outname="all_cells",res=0.1) #done
harmony_sample_integration(x="normalepithelial.SeuratObject.rds",outname="normalepithelial",res=0.1) #done
harmony_sample_integration(x="tcell.SeuratObject.rds",outname="tcell",res=0.5)

harmony_sample_integration(x="myeloid.SeuratObject.rds",outname="myeloid",res=0.2)
harmony_sample_integration(x="plasmacell.SeuratObject.rds",outname="plasmacell",res=0.1)
harmony_sample_integration(x="fibroblast.SeuratObject.rds",outname="fibroblast",res=0.3)
harmony_sample_integration(x="bcell.SeuratObject.rds",outname="bcell",res=0.2)
harmony_sample_integration(x="endothelial.SeuratObject.rds",outname="endothelial",res=0.1)
harmony_sample_integration(x="pericyte.SeuratObject.rds",outname="pericyte",res=0.1)
harmony_sample_integration(x="epithelial.SeuratObject.rds",outname="epithelial",res=0.1)

```


## Transcription Factor Expression Markers
<!-- Rerun -->
Based on seurat tutorial https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1
Using average AUC to define markers that work across modalities (RNA, Gene Activity, and TF motifs). Doing this across cell types, and then within cell types across diagnoses.

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
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")




###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    
       

diag_cols<-c("IDC"="red", "DCIS"="grey","NAT"="lightblue","ILC"="green")

molecular_type_cols<-c("DCIS"="grey", "ER+/PR-/HER2-"="#EBC258", "ER+/PR-/HER2+"="#F7B7BB","ER+/PR+/HER2-"="#ff6699","NA"="lightblue")
########################################


#Grab top overlapping TFs
topTFs <- function(markers_list,celltype, padj.cutoff = 1e-2,rna=NA,ga=NA,motifs=NA) {
  ctmarkers_rna <- dplyr::filter(
    rna, RNA.group == celltype) %>% 
    arrange(-RNA.auc)

    if(is.data.frame(motifs)) {
    ctmarkers_motif <- dplyr::filter(
      motifs, chromvar.group == celltype) %>% 
      arrange(-chromvar.auc)
    }

    if(is.data.frame(ga)) {
    ctmarkers_ga<- dplyr::filter(
      ga, GeneActivity.group == celltype) %>% 
      arrange(-GeneActivity.auc)
    }

    if(is.data.frame(motifs) && is.data.frame(ga)){    
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs <- inner_join(
        x = top_tfs ,
        y = ctmarkers_ga [,c(2, 11, 6, 7)], by = "gene"
      )
    }else if(is.data.frame(motifs)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
    } else if (is.data.frame(ga)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_ga[,c(2, 11, 6, 7)], by = "gene"
      )
    } 
  auc_colnames<-grep(".auc$",colnames(top_tfs))
  top_tfs$avg_auc <-  rowMeans(top_tfs[auc_colnames])
  top_tfs <- arrange(top_tfs, -avg_auc)
  top_tfs$celltype<-celltype
  return(top_tfs)
}


#Make volcano plot per modality
plot_volcano<-function(markers.=markers,prefix,assay){
  markers.$sig<-ifelse(markers.$padj<=0.05,"sig","non_sig")
  markers.<-markers.[!duplicated(markers.$feature),]
  markers.<-markers.[is.finite(markers.$pval),]
  markers.$label<-""
  markers.<-markers.[order(markers.$padj),]
  ident_1_labels<- row.names(head(markers.[which((markers.$sig=="sig") & (markers.$logFC>0)),] ,n=10))#significant, is ident1 specific, is top 10
  ident_2_labels<- row.names(head(markers.[which((markers.$sig=="sig") & (markers.$logFC<0)),] ,n=10))#significant, is ident1 specific, is top 10
  markers.[c(ident_1_labels,ident_2_labels),]$label<-markers.[c(ident_1_labels,ident_2_labels),]$feature
  x_scale<-max(abs(markers.$logFC))*1.1 #find proper scaling for x axis
  plt<-ggplot(markers.,aes(x=logFC,y=-log10(padj),color=sig,label=label,alpha=0.1))+
    geom_point(size=0.5)+
    theme_bw()+
    scale_fill_manual(values=c("non_sig"="#999999", "sig"="#FF0000"))+
    xlim(c(-x_scale,x_scale))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=-log10(0.05))+
    geom_text_repel(size=2,segment.size=0.1,max.overlaps=Inf,min.segment.length = 0, nudge_y = 1, segment.angle = 20,color="black") +
    theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
    ggsave(plt,file=paste0(prefix,"_",assay,"_DE_volcano.pdf"),width=5,units="in",height=5)
  system(paste0("slack -F ",prefix,"_",assay,"_DE_volcano.pdf"," ryan_todo"))
}

#Identify top markers
Identify_Marker_TFs<-function(x,group_by.="predicted.id",assay.="RNA",prefix.){
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by., 
      groups_use=unname(unlist(unique(x@meta.data[group_by.]))),
      y=unname(unlist(unique(x@meta.data[group_by.]))), 
      assay = 'data', seurat_assay = assay.)
    write.table(markers,file=paste0(prefix.,"_",assay.,"_DE_table.tsv"),sep="\t",row.names=F,col.names=T,quote=F)
    plot_volcano(markers.=markers,prefix=prefix.,assay=assay.)
    colnames(markers) <- paste(assay., colnames(markers),sep=".")
    if (assay. == "chromvar") {
      motif.names <- markers[,paste0(assay.,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="ATAC")
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
plot_top_TFs<-function(x=stromal,tf_markers=da_tf_markers,prefix="stromal",group_by.="predicted.id",CHROMVAR=TRUE,GA=TRUE,height.){
    tf_rna<-average_features(x=x,features=tf_markers$gene,assay="RNA",group_by.=group_by.)
    tf_rna<-tf_rna[row.names(tf_rna) %in% tf_markers$gene,]

  if(CHROMVAR){
    tf_motif<-average_features(x=x,features=tf_markers$chromvar.feature,assay="chromvar",group_by.=group_by.)
    tf_motif<-tf_motif[row.names(tf_motif) %in% tf_markers$chromvar.feature,]
    row.names(tf_motif)<-tf_markers[tf_markers$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
  }

  if(GA){
    tf_ga<-average_features(x=x,features=tf_markers$gene,assay="GeneActivity",group_by.=group_by.)
    tf_ga<-tf_ga[row.names(tf_ga) %in% tf_markers$gene,]
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_ga<-tf_ga[markers_list,]

  }
  if(GA&&CHROMVAR){
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
  }

    #set up heatmap seriation and order by RNA
    o = seriate(max(tf_rna) - tf_rna, method = "BEA_TSP")
    saveRDS(o,file=paste0(prefix,".geneactivity.dend.rds")) 
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o,1),]$RNA.auc)
    colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.80,0.95)),plasma(3))

  if(CHROMVAR){
    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.80,0.95)),cividis(3))
    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="ATAC",motifs = motif_list[get_order(o,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

  }
  if(GA){
    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.80,0.95)),magma(3))

  }

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))


    rna_auc<-Heatmap(side_ha_rna,
        row_order = get_order(o,1),
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))
    if(!CHROMVAR){
    rna_plot<-Heatmap(tf_rna,
        row_order = get_order(o,1),
        column_order = get_order(o,2),
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90,right_annotation=gene_ha)
    } else {
    rna_plot<-Heatmap(tf_rna,
        row_order = get_order(o,1),
        column_order = get_order(o,2),
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90)
  }
  if(GA){
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
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)

  }
  if(CHROMVAR){
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
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
  }

  if(all(CHROMVAR,GA)){
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot)
  } else if(CHROMVAR){
      plt1<-draw(rna_auc+rna_plot+motif_auc+motif_plot)
  } else {
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot)
  }


    pdf(paste0(prefix,".tf.heatmap.pdf"),height=height.)
    print(plt1)
    dev.off()

    system(paste0("slack -F ",prefix,".tf.heatmap.pdf ryan_todo"))
    system(paste0("slack -F ",prefix,".tf.heatmap.motif.pdf ryan_todo"))
}

#Final wrapper function
run_top_TFs<-function(obj=stromal,prefix="stromal",i="predicted.id",n_markers=5,CHROMVAR=TRUE,plot_height=10){
  if(CHROMVAR){
  markers<-lapply(c("RNA","GeneActivity","chromvar"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay,prefix.=prefix))
  names(markers)<-c("RNA","GeneActivity","chromvar")
  markers_out<-do.call("rbind",lapply(unique(obj@meta.data[,i]),function(x) head(topTFs(markers_list=markers,celltype=x,rna=markers$RNA,ga=markers$GeneActivity,motifs=markers$chromvar),n=n_markers))) #grab top 5 TF markers per celltype
  dim(markers_out)
  markers_out<-markers_out[!duplicated(markers_out$gene),]
  dim(markers_out)
  saveRDS(markers_out,file=paste0(prefix,"_celltype_TF_markers.RDS"))
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
  plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i,CHROMVAR=TRUE,GA=TRUE)
  } else{
  markers<-lapply(c("RNA","GeneActivity"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay,prefix.=prefix))
  names(markers)<-c("RNA","GeneActivity")
  markers_out<-do.call("rbind",lapply(unique(obj@meta.data[,i]),function(x) head(topTFs(markers_list=markers,celltype=x,rna=markers$RNA,ga=markers$GeneActivity),n=n_markers))) #grab top 5 TF markers per celltype
  dim(markers_out)
  markers_out<-markers_out[!duplicated(markers_out$gene),]
  dim(markers_out)
  saveRDS(markers_out,file=paste0(prefix,"_celltype_markers.RDS"))
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_markers.RDS"))
  plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i,CHROMVAR=FALSE,GA=TRUE,height.=plot_height)
}
}


#all cells celltypes
dat<-readRDS("phase2.QC.filt.SeuratObject.rds")
dat<-subset(dat,EMBO_predicted.id!="NA")
run_top_TFs(obj=dat,prefix="dat",i="EMBO_predicted.id",n_markers=8,CHROMVAR=FALSE,plot_height=15)
run_top_TFs(obj=dat,prefix="dat",i="EMBO_predicted.id",n_markers=8,CHROMVAR=TRUE,plot_height=15)

for (k in c("epithelial.SeuratObject.rds","myeloid.SeuratObject.rds","tcell.SeuratObject.rds","plasmacell.SeuratObject.rds","fibroblast.SeuratObject.rds","bcell.SeuratObject.rds","endothelial.SeuratObject.rds","pericyte.SeuratObject.rds","normalepithelial.SeuratObject.rds")) {
  dat_sub<-readRDS(k)
  dat_prefix<-gsub("rds$","",gsub(".SeuratObject.","",basename(k)))
  run_top_TFs(obj=dat_sub,prefix=paste0(dat_prefix,"_subcluster"),i="seurat_clusters")
  run_top_TFs(obj=dat_sub,prefix=paste0(dat_prefix,"_subcluster"),i="seurat_clusters",CHROMVAR=FALSE,n_markers=15)
}

```

## Cell Subtyping

### Normal epithelial cell subtyping.
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
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")


cov_plots<-function(atac_sub=atac_sub,gene_name,idents_in){
  plt_cov <- CoveragePlot(
    object = atac_sub,
    region = gene_name,
    features = gene_name,
    assay="ATAC",
    expression.assay = "SoupXRNA",
    extend.upstream = 10000,
    extend.downstream = 10000,
    idents=idents_in)
  plt_feat <- FeaturePlot(
    object = atac_sub,
    features = gene_name,
    raster=F,
    reduction="multimodal_harmony_umap",
    order=T)
  return((plt_feat|plt_cov)+ggtitle(gene_name))
}

sample_label_transfer<-function(in_dat,ref_dat,prefix="Tcell_"){
  transfer.anchors <- FindTransferAnchors(
    reference = ref_dat,
    reference.assay="RNA",
    query = in_dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = ref_dat$celltype,
  )
  colnames(predictions)<-paste0(prefix,colnames(predictions))

  in_dat<-AddMetaData(in_dat,metadata=predictions)
  return(in_dat)
  }

```
### Epithelial Subtyping
Continued session. 
```R
#Plotting Normal Epithelial Marker genes
  normal_epithelial_marker_genes<-c("KRT5", "ACTA2", "MYLK", "SNAI2", "NOTCH4", "DKK3", "ESR1", "PGR", "FOXA1", "TNFRSF11A", "KIT", "SOX10")
  atac_sub<-readRDS("normalepithelial.SeuratObject.rds")
  DefaultAssay(atac_sub)<-"SoupXRNA"
  Idents(atac_sub)<-atac_sub$seurat_clusters
  for (i in normal_epithelial_marker_genes){
    plt<-cov_plots(atac_sub=atac_sub,gene_name=i,idents_in=c("0","1","2"))
    ggsave(plt,file=paste0("NormEpi_",i,".featureplots.pdf"),limitsize=F)
    system(paste0("slack -F ","NormEpi_",i,".featureplots.pdf ryan_todo"))
  }
#Plotting EMBO paper defined marker sets of epithelial subtypes (pretransferred on bulk data above)
  feature_set<-c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str")
  plt1<-VlnPlot(atac_sub, features = feature_set)
  plt2<-FeaturePlot(atac_sub, features = feature_set,order=T,reduction="multimodal_harmony_umap",)
  plt<-plt1/plt2
  ggsave(plt,file="NormEpi_EMBOfeaturesets.pdf")
  system("slack -F NormEpi_EMBOfeaturesets.pdf ryan_todo")

```

### T cell Subtyping
Continued session. Analyzing T cells using marker genes and label transfer of T cell subtypes.
Using 4 clusters (res 0.3)

```R
#Plotting T cell marker genes
  tcell_marker_genes<-c("CD4","CD8A","ITGAE","NCR1","IL2RA","PDCD1","CCR7","IL7R","FOXP3","CXCL13","ZFP36","GZMK","IFIT1","MKI67")
#T cell RNA reference data 
  tcell_reference_rds<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
  DefaultAssay(tcell_reference_rds)<-"RNA"
  tcell_reference_rds<-NormalizeData(tcell_reference_rds)
  tcell_reference_rds<-FindVariableFeatures(tcell_reference_rds)
  tcell_reference_rds<-ScaleData(tcell_reference_rds)

  atac_sub<-readRDS("tcell.SeuratObject.rds")
  DefaultAssay(atac_sub)<-"SoupXRNA"
  atac_sub <- FindClusters(atac_sub, graph.name="wknn",verbose = FALSE,resolution=0.3)
  Idents(atac_sub)<-atac_sub$seurat_clusters
  atac_sub<-sample_label_transfer(in_dat=atac_sub,ref_dat=tcell_reference_rds,prefix="tcell_")

#Plotting transferred cell labels
  feature_set<-colnames(atac_sub@meta.data)[startsWith(colnames(atac_sub@meta.data),prefix="tcell_prediction.score.")]
  plt1<-DimPlot(atac_sub,group.by="seurat_clusters",reduction="multimodal_harmony_umap")
  plt2<-VlnPlot(atac_sub, features = feature_set)
  plt3<-FeaturePlot(atac_sub, features = feature_set,order=T,reduction="multimodal_harmony_umap")
  plt<-plt1/plt2/plt3
  ggsave(plt,file="tcell_EMBOfeaturesets.pdf",height=20,width=20)
  system("slack -F tcell_EMBOfeaturesets.pdf ryan_todo")

#Plot marker genes also
  for (i in tcell_marker_genes){
    plt<-cov_plots(atac_sub=atac_sub,gene_name=i,idents_in=unique(Idents(atac_sub)))
    ggsave(plt,file=paste0("tcell_",i,".featureplots.pdf"),limitsize=F)
    system(paste0("slack -F ","tcell_",i,".featureplots.pdf ryan_todo"))
  }

#Save RDS now that there are label transfer metadata columns
saveRDS(atac_sub,"tcell.SeuratObject.rds")


```

### B Cell Subtyping
Continued session. Analyzying B cells using marker genes. Marker gene list from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-22300-2/MediaObjects/41467_2021_22300_MOESM5_ESM.xlsx
```R
#Plotting B cell marker genes
  bcell_marker_genes<-c("IGHM","CD27","MS4A1","CD19","CD27","CD38","CD14") #"IGHM",#,"IGHG", #IGHD

#B cell RNA reference data 
  bcell_reference_rds<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalSub.rds") #ER+ tumor non epithelial cells

  DefaultAssay(bcell_reference_rds)<-"RNA"
  bcell_reference_rds<-NormalizeData(bcell_reference_rds)
  bcell_reference_rds<-FindVariableFeatures(bcell_reference_rds)
  bcell_reference_rds<-ScaleData(bcell_reference_rds)

  atac_sub<-readRDS("bcell.SeuratObject.rds")
  DefaultAssay(atac_sub)<-"SoupXRNA"
  Idents(atac_sub)<-atac_sub$seurat_clusters
  atac_sub <- FindClusters(atac_sub, graph.name="wknn",verbose = FALSE,resolution=0.3)
  atac_sub<-sample_label_transfer(in_dat=atac_sub,ref_dat=bcell_reference_rds,prefix="bcell_")

  feature_set<-colnames(atac_sub@meta.data)[startsWith(colnames(atac_sub@meta.data),prefix="bcell_prediction.score.")]
  plt1<-DimPlot(atac_sub,group.by="seurat_clusters",reduction="multimodal_harmony_umap")
  plt2<-VlnPlot(atac_sub, features = feature_set)
  plt3<-FeaturePlot(atac_sub, features = feature_set,order=T,reduction="multimodal_harmony_umap")
  plt<-plt1/plt2/plt3
  ggsave(plt,file="bcell_EMBOfeaturesets.pdf",height=20,width=20)
  system("slack -F bcell_EMBOfeaturesets.pdf ryan_todo")

  atac_sub<-readRDS("bcell.SeuratObject.rds")
  DefaultAssay(atac_sub)<-"SoupXRNA"
  Idents(atac_sub)<-atac_sub$seurat_clusters
  for (i in bcell_marker_genes){
    plt<-cov_plots(atac_sub=atac_sub,gene_name=i,idents_in=unique(Idents(atac_sub)))
    ggsave(plt,file=paste0("bcell_",i,".featureplots.pdf"),limitsize=F)
    system(paste0("slack -F ","bcell_",i,".featureplots.pdf ryan_todo"))
  }

bcell_list<-list()
bcell_list[["bcell_memory"]]<-c(
  'CRIP1', 'TNFRSF13B', 'CD82', 'COTL1', 'B2M', 'RPS14', 'AIM2', 'RP11-731F5.2', 'CD27', 'GAPDH', 'CLECL1', 'LGALS1', 'VIM', 'MALAT1', 'HSPA8', 'ACP5', 'S100A6', 'HLA-A', 'ITGB1', 'C1orf186', 'LSP1', 'PSMB9', 'RP5-887A10.1', 'S100A4', 'TOMM7', 'PLP2', 'ACTG1', 'KLK1', 'CD99', 'LTB', 'CAPG', 'CAPN2', 'CHCHD10', 'CD24', 'CCDC50', 'CLIC1', 'CD70', 'GPR183', 'ARHGAP24', 'SSPN', 'CIB1', 'RAC2', 'SMARCB1', 'PVT1', 'RPS29', 'NLRC5', 'CD1C', 'LDHB', 'CTSH', 'MS4A1', 'CDC42EP3', 'NEK6', 'KCNN4', 'BANK1', 'TLR10', 'GSTK1', 'ARL6IP5', 'HIGD2A', 'ARPC1B', 'VOPP1', 'AC079767.4', 'PKM', 'NCR3', 'CNPY3', 'SCIMP', 'CD48', 'FAM102A', 'RAB31', 'SOD1', 'IFITM2', 'RP5-1171I10.5', 'KIAA1551', 'SAMSN1', 'IFITM1', 'RILPL2', 'FCRL2', 'DOK3', 'NCF4', 'S100A10', 'HSPA1B', 'CPNE5', 'AHNAK', 'TCF4', 'DYNLL1', 'SMIM14', 'ID3', 'KLF6', 'TNF', 'HSPB1', 'RHOB', 'DNAJB1')

bcell_list[["bcell_naive"]]<-c('TCL1A','IL4R', 'PLPP5', 'FCER2', 'BACH2', 'CXCR4', 'CD69', 'YBX3', 'IGHD', 'BTG1', 'RPL18A', 'LAPTM5', 'CD200', 'CD37', 'CD83', 'LINC00926', 'SKAP1', 'FOXP1', 'TMEM123', 'CAMK2D', 'IL21R', 'HVCN1', 'COL19A1', 'IGLL5', 'CCR7', 'APLP2', 'SARAF', 'CLEC2B', 'TSPAN13', 'CD72', 'RP11-231C14.7', 'BCL7A', 'TAGAP', 'CLEC2D', 'SATB1', 'C1orf162', 'SNX29', 'FCRL1', 'AFF3', 'ZNF318', 'RAB30', 'CD55', 'PLEKHA2', 'HLA-DQA1', 'ADK', 'FAM129C', 'PCDH9', 'RHOH', 'ZCCHC7', 'CD79A', 'MEF2C', 'NCF1', 'SNX9', 'CDCA7L', 'REL', 'EIF1B', 'EIF2AK3', 'CD22', 'SSBP2', '43527', 'LYST', 'FOS', 'C16orf74', 'FAM26F', 'RFTN1', 'DDIT3', 'KIAA0226L', 'LAIR1', 'DGKD', 'NFKB1', 'GABPB1', 'PHACTR1', 'ATF7IP', 'SMAP2', 'SLC38A1', 'SNX2', 'SESN1', 'HIF1A', 'PNRC1', 'FCMR', 'RBM38', 'BEX4', 'LYN', 'ITPR1', 'C12orf57', 'WHSC1L1', 'BZW1', 'STK17A', 'GALNT2', 'SESTD1', 'ZFP36L1', 'DUSP1', 'HNRNPA3', 'KHDRBS2', 'BTLA', 'ABCG1', 'ITM2B', 'RNASE6', 'WDR74', 'CCND3', 'BLOC1S2', 'EAF2', 'SELL', 'RCSD1', 'TCTN1', 'IRF8', 'PIK3IP1', 'IGLC3', 'VPS37B', 'AES', 'PIM3', 'AC241585.2', 'CD79B', 'ARL4A', 'DTNBP1', 'GABPB1-AS1', 'AC245100.1', 'NUP58', 'SSH2', 'ARRDC2', 'LINC01215', 'BIRC3', 'SEC62', 'H1FX', 'LTA4H', 'PLEKHG1', '43723', 'NFKBID', 'SIPA1L1', 'FAM3C', 'TUBA1A', 'RP9', 'RCN2', 'WHAMM', 'CYTIP', 'HLA-DQB1', 'SNHG8', 'DBI', 'SLC2A3', 'ST3GAL1', 'JUN', 'CHPT1', 'PELI1', 'HIST1H1C', 'MARCKSL1', 'ZNF331')

bcell_list[["bcell_plasma"]]<-c('TXNDC5', 'AQP3', 'JSRP1', 'CHPF', 'HRASLS2', 'NT5DC2', 'RP11-290F5.1', 'PRDM1', 'SDC1', 'NUCB2', 'TRIB1', 'BIK', 'SLAMF7', 'ABCB9', 'IGHGP', 'IGHJ4', 'FKBP11', 'IGHG4', 'FNDC3B', 'TNFRSF17', 'SEMA4A', 'MAN1A1', 'SDF2L1', 'CD38', 'GAS6', 'SLC1A4', 'XBP1', 'DERL3', 'RRBP1', 'PRDX4', 'MZB1', 'SPAG4', 'CD59', 'LIME1', 'LGALS3', 'RP11-1070N10.3', 'SIL1', 'ERN1', 'PDIA4', 'HDLBP', 'VIMP', 'FKBP2', 'CLPTM1L', 'LINC00152', 'SEC11C', 'WARS', 'MANF', 'SSR3', 'MYDGF', 'HYOU1', 'TXNDC11', 'SSR4', 'ITM2C', 'SEC61A1', 'HSP90B1', 'SEC24D', 'ERLEC1', 'NANS', 'HM13', 'PDK1', 'SLC44A1', 'AC093818.1', 'SPCS3', 'ELL2', 'SEL1L', 'SPATS2', 'CCDC167', 'LMAN1', 'SEC61B', 'P4HB', 'PREB', 'JCHAIN', 'MLEC', 'PPIB', 'ZBP1', 'SRPRB', 'MYO1D', 'IGLL5', 'DSTN', 'TXNDC15', 'ISOC2', 'HSPA5', 'IGLV3-1', 'FAM46C', 'ANKRD28', 'IGHA2', 'CCPG1', 'OSTC', 'GMPPB', 'MEI1', 'KDELR2', 'DNAJC1', 'ST6GALNAC4', 'SEC61G', 'CRELD2', 'LMAN2', 'DNAJC3', 'TMEM258', 'TMED9', 'SLC17A9', 'ARFGAP3', 'SPCS2', 'SUB1', 'ARSA', 'IGHG1', 'CECR1', 'DNAJB11', 'PHPT1', 'HERPUD1', 'FNDC3A', 'SEC14L1', 'TMEM208', 'IRF4', 'RPN2', 'SRGN', 'EDEM1', 'B4GALT3', 'SURF4', 'KDELR1', 'UBE2J1', 'IDH2', 'TMED10', 'DDOST', 'RPN1', 'TP53INP1', 'KLF13', 'GSTP1', 'RABAC1', 'PDIA6', 'ATF5', 'IGHG3', 'TXN', 'IFI27L2', 'DNAJB9', 'CFLAR', 'SPCS1', 'SLC9A3R1', 'GLRX', 'SRM', 'ATOX1', 'CD27', 'LGALS1', 'GORASP2', 'IGKC', 'IFI6', 'IGHA1', 'DENND1B', 'KRTCAP2', 'CANX', 'ATF4', 'SSR1', 'COX5A', 'TECR', 'COPE', 'IGKV1-12', 'MRPS24', 'SLC3A2', 'C12orf75', 'TRAM1', 'ISG20', 'ACADVL', 'CDK2AP2', 'IGHM', 'UQCRQ', 'IGHG2', 'DAD1', 'GOLGB1', 'SRPRA', 'CALR', 'PIM2', 'PABPC4', 'ROMO1', 'RPS27L', 'DNM2', 'IGKV3-20', 'MIF', 'TMEM59', 'ZNF706', 'H1FX', 'OGT', 'IGLC2', 'HCST', 'IGHV3-7', 'PRDX5', 'NDUFA1', 'IGKV3-15') 

bcell_list[["bcell_cd14"]]<-c('TYROBP', 'LYZ', 'CST3', 'APOC1', 'C1QC', 'FCER1G', 'C1QA', 'C1QB', 'AIF1', 'CCL3L3', 'CCL4L2', 'CD14', 'RNASE1', 'CXCL2', 'MS4A6A', 'CCL2', 'PLAUR', 'FCN1', 'MAFB', 'SERPINA1', 'ANXA1', 'CXCL3', 'C15orf48', 'TMEM176B', 'CTSL', 'CEBPD', 'FCGR2A', 'IL1RN', 'STAB1', 'TMEM176A', 'CFD', 'IGSF6', 'CD163', 'CLEC7A', 'GPNMB', 'TNFAIP2', 'DAB2', 'CPVL', 'TREM2', 'FYB', 'MS4A4A', 'LILRB2', 'CSTA', 'TNFSF13B', 'FCGR3A', 'C5AR1', 'LRP1', 'TGFBI', 'CSF1R', 'CSF3R', 'PILRA', 'HNMT', 'FPR1', 'LINC01272', 'ETS2', 'MSR1', 'PTGS2', 'FCGR1A', 'JAML', 'VSIG4', 'RASSF4', 'MRC1', 'CD36', 'CD4', 'GIMAP4', 'LILRB4', 'NLRP3', 'LGALS2', 'RAB32', 'IL18', 'CD302', 'DMXL2', 'AOAH', 'SERPING1', 'LILRB3', 'PLXDC2', 'FRMD4B', 'CLEC4E', 'NRP1', 'SIGLEC1', 'FPR3', 'TNFRSF1A', 'SLCO2B1', 'C3AR1', 'LILRA5', 'PLA2G7', 'LTBR', 'ADGRE2', 'DAPK1', 'ASGR1', 'HK3', 'C10orf11', 'SPHK1', 'MERTK', 'FCGR1B', 'ACTN1', 'CD33', 'OSCAR', 'SLC8A1', 'PYGL', 'BST1', 'SECTM1', 'ANPEP', 'FN1', 'S100A9', 'TREM1', 'F13A1', 'VCAN', 'MGST2', 'SIRPA', 'A2M', 'PVRL2', 'SPRED1', 'KCTD12', 'PLXND1', 'LGALS3', 'SEPP1', 'IER3', 'NPL', 'SLC31A2', 'PLTP', 'PLAU', 'IL1B', 'ID2', 'LPAR6', 'CXCL8', 'APOE', 'FAM105A', 'CAMK1', 'SPP1', 'PDK4', 'GNA15', 'PHLDA2', 'RAB20', 'TIMP2', 'PLBD1', 'IFITM3', 'TLR4', 'ADM', 'CFP', 'CCL3', 'HAVCR2', 'DBNDD2', 'MPP1', 'MCTP1', 'TTYH3', 'SAMHD1', 'CD300A', 'SLC11A1', 'TIMP1', 'CREG1', 'SULT1A1', 'ENG', 'CCR1', 'DOK2', 'CREB5', 'GIMAP1', 'LST1', 'GIMAP7', 'GLUL', 'TSPAN4', 'PLXNB2', 'UPP1', 'MAF', 'S100A8', 'PTMS', 'CD68', 'HMOX1', 'SMCO4', 'G0S2', 'HBEGF', 'BRI3', 'RP11-670E13.6', 'CCL4', 'FOSL2', 'S100A11', 'RNF130', 'CYFIP1', 'CTSD', 'ADAM9', 'FCGRT', 'GPR34', 'PTPRE', 'ZNF385A', 'NINJ1', 'CTSB', 'ALDH2', 'DUSP6', 'ITGAX', 'PTAFR', 'DUSP3', 'LGALS3BP', 'CXCL16', 'LCP2', 'BLVRB', 'ADAP2', 'LRRC25', 'BLVRA', 'MS4A7', 'STX11', 'FNDC3B', 'RNF19B', 'CMTM3', 'C1orf54', 'SLC15A3', 'KLF4', 'NCF2', 'ACP2', 'GABARAPL1', 'LY96', 'CLEC4A', 'RAB13', 'MCOLN1', 'TRIB1', 'IL6R', 'LMNA', 'FGL2', 'TYMP', 'FOLR2', 'MGLL', 'SGK1', 'PLIN2', 'SRGN', 'TNFSF13', 'HCST', 'FTL', 'TBXAS1', 'NAGA', 'ABCA1', 'CSTB', 'FBP1', 'MIR4435-2HG', 'FABP5', 'CD63', 'IDH1', 'GPR137B', 'ARHGAP18', 'PLD3', 'SLC7A7', 'MNDA', 'EFHD2', 'IFI30', 'STOM', 'GAA', 'SAT1', 'ICAM1', 'NPC2', 'EMILIN2', 'CTSC', 'RENBP', 'AAK1', 'GNS', 'LGMN', 'ATP6V1B2', 'ANXA5', 'LGALS1', 'CEBPB', 'FTH1', 'ZYX', 'RBM47', 'C1orf162', 'AGTRAP', 'S100A10', 'FUCA1', 'ODF3B', 'GPX1', 'GNPDA1', 'GSN', 'S100A4', 'S100A6', 'NAIP', 'IFI6', 'LIMS1', 'SLC37A2', 'H2AFJ', 'CCDC88A', 'RAC1', 'GRN', 'MILR1', 'GBP2', 'VIM', 'SOD2', 'DRAM1', 'SPI1', 'HSPA1A', 'ARRB2', 'GSTP1', 'MT1G', 'ANXA2', 'ITGB2', 'GLIPR2', 'LINC00152', 'RGS10', 'SERPINB6', 'SLC16A3', 'PPT1', 'MYO1F', 'THEMIS2', 'TXN', 'NQO2', 'PTGER4', 'GLRX', 'C19orf38', 'MARCKS', 'SH3BGRL3', 'PSAP', 'PRDM1', 'SCARB2', 'GSTO1', 'MYL6', 'SERF2', 'HSBP1', 'CD151', 'GK', 'ASAH1', 'TSPO', 'GNG10', 'COMT', 'ZFHX3', 'ATF3', 'TUBA1C', 'PFKFB3', 'HEBP1', 'PYCARD', 'CD9', 'ATOX1', 'VAMP8', 'GNG5', 'SDCBP', 'AP2S1', 'LGALS9', 'CYSTM1', 'DPYD', 'SOCS3', 'CTSZ', 'RNF144B', 'MGAT1', 'HCK', 'PLK3', 'ABHD12', 'LAP3', 'H2AFY', 'PLSCR1', 'TMSB10', 'CD59', 'RAP2B', 'CARD16', 'FKBP15', 'CPPED1', 'BHLHE40', 'IQGAP2', 'APLP2', 'ISG15', 'GBP1', 'NANS', 'RHOC', 'TMSB4X', 'JOSD2', 'GAPDH', 'PCBD1', 'DUSP23', 'GLMP', 'NEAT1', 'VAMP5', 'C4orf48', 'LAMP2', 'NABP1', 'GABARAP', 'SH3BP2', 'RRBP1', 'RGCC', 'ABRACL', 'PPIF', 'PGD', '43526', 'BID', 'CTSA', 'GNAI2', 'RTN3', 'ATP5E', 'QKI', 'ATF5', 'STMN1', 'YIF1B', 'OAZ2', 'IL17RA', 'ACSL1', 'C10orf54', 'TIMM8B', 'LAMTOR2', 'ATP6V1F', 'NDUFV3', 'AKR1A1', 'ITM2B', 'PAK1', 'HSPA1B', 'RAB31', 'MGST3', 'TXNDC17', 'AGPAT2', 'ERN1', 'CTSS', 'DSE', 'DYNLL1', 'PDXK', 'BCL2A1', 'TALDO1', 'RNH1', 'SDF2L1', 'PTTG1IP', 'RUNX1', 'TPM4', 'BAG3', 'HN1', 'STXBP2', 'GRINA', 'ZEB2', 'AHR', 'MAFF', 'LACTB2', 'RPPH1', 'MKNK1', 'NOP10', 'MGAT4A', 'HEXB', 'PRDX3', 'UQCR10', 'TFEC', 'CLTA', 'DBI', 'PICALM', 'VMP1', 'RPS27L', 'ARPC5', 'DSTN', 'SRA1', 'TNFAIP3', 'ATP5J2', 'GADD45B', 'CFL1', 'BST2', 'ARF3', 'PHPT1', 'EHD4', 'COX6B1', 'EIF4EBP1', 'NAGK', 'CASP1', 'CKLF', 'ARPC1B', 'SLC3A2', 'RGS1', 'MFSD1', 'HSPB1', 'TMEM167A', 'NAMPT', 'CECR1', 'MTHFD2', 'RTN4', 'CBR1', 'CAPG', 'EHBP1L1', 'SERPINB1', 'GLA', 'KLF10', 'COX8A', 'MXD1', 'TCEB2', 'WARS', 'UBXN11', 'LAIR1', 'DYNLT1', 'ATP6V0C', 'XBP1', 'PLEK', 'NUP214', 'LINC00116', 'RNF181', 'HCFC1R1', 'MPEG1', 'ETV6', 'PLEKHB2', 'RBPJ', 'ATP5H', 'ABI3', 'MRPL41', 'FLOT1', 'SAMSN1', 'CRTAP', 'YWHAH', 'IFIT3', 'COX5B', 'IFNGR1', 'ATP6AP2', 'NEU1', 'P4HB', 'GTF2H5', 'TPI1', 'USMG5', 'ANAPC11', 'UQCR11', 'NDUFB1', 'ARL8B', 'ALDOA', 'CISD2', 'SLC43A2', 'PFN1', 'FKBP1A', 'ATP6AP1', 'GPX4', 'SCO2', 'NFKBIZ', 'SOX4', 'NCOA4', 'ATG3', 'FKBP2', 'ATP5J', 'PPIB', 'FGR', 'MT1E', 'OAZ1', 'ENO1', 'ARAP1', 'GUSB', 'UBL5', 'EPSTI1', 'PLEKHO1', 'HEXA', 'PRDX1', 'TNFRSF1B', 'CTA-29F11.1', 'CLEC2B', 'LDHA', 'GNPTG', 'CYBA', 'MID1IP1', 'SEPW1', 'NDUFS5', 'NDUFB3', 'TMEM14C', 'TUBA1B', 'NDUFA1', 'ATP6V0B', 'CALM3', 'PNP', 'AP1S2', 'MYO9B', 'CD86', 'LINC00936', 'NDUFC1', 'ACTG1', 'LITAF', 'ACTB', 'TNF', 'NDUFB2', 'C14orf2', 'UQCRQ', 'NUCB1', 'CFLAR', 'NDUFA2', 'NENF', 'CITED2', 'PSMA6', 'CD81', 'HEXIM1', 'POLR2L', 'HM13', 'BNIP3L', 'POMP', 'PRELID1', 'RGS2', 'NR4A3', 'YWHAE', 'CAPZB', 'NDUFS6', 'UBE2L6', 'LCP1', 'RBMS1', 'LAMTOR4', 'CYBB', 'IFI27L2', 'TKT', 'NDUFA3', 'PKM', 'SSR3', 'GLIPR1', 'COX7B', 'ATP6V0D1', 'MINOS1', 'CALM2', 'TLN1', 'NBPF19', 'TWF2', 'MYDGF', 'CLIC1', 'GLUD1', 'SEC11A', 'RNASEK', 'FAM46A', 'GMFG', 'C7orf73', 'MYEOV2', 'RNF149', 'SEC61G', 'ATP6V0E1', 'NDUFB7', 'C20orf24', 'TMEM258', 'CAPZA2', 'RAB10', 'RBX1', 'RBKS', 'IVNS1ABP', 'FERMT3', 'PSMB3', 'TMEM219', 'H2AFV', 'FAM96B', 'COX5A', 'NXF1', 'MMP24-AS1', 'RHOG', 'KDELR1', 'PHACTR1', 'OST4', 'NAA38', 'SHFM1', 'ARPC2', 'MIR155HG', 'SNX6', 'SERTAD1', 'TNFSF10', 'LIPA', 'BANF1', 'CALR', 'ENY2', 'CD99', 'COTL1', 'CALM1', 'GNB2', 'ATP5G3', 'SLC38A2', 'ZFAND5', 'GUK1', 'ROMO1', 'RAB5C', 'NDUFB4', 'RNF213', 'MAP3K8', 'DRAP1', 'PGAM1', 'BRK1', 'LAMTOR1', 'MYL12B', 'CTSH', 'PSMA7', 'NDUFA4', 'ATP5C1', 'CDKN1A', 'HSPA5', 'PSME2', 'H3F3A', 'H2AFZ', 'HSP90AA1', 'CKS2', 'COPE', 'MIDN', 'IRF7', 'GNAS', 'FLNA')


  module_scores<-AddModuleScore(atac_sub,features=bcell_list,assay="SoupXRNA",search=TRUE,name=names(bcell_list)) #use add module function to add cell scores
  module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-(length(module_feats)-1),ncol(module_scores@meta.data))]
  colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like
  dat_epi<-AddMetaData(dat,metadata=module_scores)



```

### TBD
"myeloid.SeuratObject.rds"
"plasmacell.SeuratObject.rds"
"fibroblast.SeuratObject.rds"
"endothelial.SeuratObject.rds"
"pericyte.SeuratObject.rds"



## Plot of Differential Genes across Normal epithelial (NAT) DCIS and IDC
<!-- Rerun -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")



atac_sub<-readRDS("epithelial.SeuratObject.rds")
Idents(atac_sub)<-atac_sub$diagnosis
DefaultAssay(atac_sub)<-"SoupXRNA"
RNA_markers <- FindMarkers(atac_sub, ident.1 = "IDC", ident.2 = c("DCIS","NAT"), min.pct = 0.1)
RNA_markers$gene_name<-row.names(RNA_markers)
DefaultAssay(atac_sub)<-"ATAC"
ATAC_markers <- FindMarkers(atac_sub,ident.1 = "IDC", ident.2 = c("DCIS","NAT"), min.pct = 0.1)
ATAC_markers$da_region<-row.names(ATAC_markers)
closest_genes <- ClosestFeature(atac_sub,ATAC_markers$da_region)
ATAC_markers<-cbind(ATAC_markers,closest_genes)

rna<-RNA_markers[RNA_markers$gene_name %in% ATAC_markers$gene_name,]
atac<-ATAC_markers[ATAC_markers$gene_name %in% rna$gene_name,]



cov_plots<-function(dat=atac_sub,gene_name,idents_in){
  plt_cov <- CoveragePlot(
    object = atac_sub,
    region = gene_name,
    features = gene_name,
    assay="ATAC",
    expression.assay = "SoupXRNA",
    extend.upstream = 5000,
    extend.downstream = 5000,
    idents=idents_in)
  plt_feat <- FeaturePlot(
    object = atac_sub,
    features = gene_name,
    raster=T,
    reduction="multimodal_umap",
    order=T)
  return((plt_feat|plt_cov)+ggtitle(gene_name))
}


DefaultAssay(atac_sub)<-"SoupXRNA"
for (i in c(rna$gene_name[1:25])){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}

for (i in c("SOX10","SOX9","SOX4","SOX2","TEAD4","RUNX1")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("FOXM1","FOXA1","FOXA3","GRHL2","FOXP1","ATF3")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("HOXB13","EN1","DLX4","TBX15","SLC6A12","PAX6","FAM83A","ERICH5")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("ESR1")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}

#TFS
#IDC Fox Family (GRHL2) 
#ILC Sox Family TEAD RUNX EGR1 RPBJ HMGA1
#DCIS STAT3/BCL9
#DCIS more likely to be invasice (Methylation IDd) HOXB13 EN1 DLX4 TBX15 SLC6A12 PAX6 

#GENES
#FAM83A
#ERICH5
```
## Comparison of cell types across diagnoses and other factors.
<!-- Rerun -->

## 3D Plotting in Blender
<!-- Rerun -->

```R
#3d umap
out_3d <- RunUMAP(object=dat, n.components=3, reduction.name="harmonyumap",reduction = "harmony", dims = 1:20)
#format
#Astrocytes    TAGGTCCGACGTACTAGGGCCTCGGTCTATGGCCTA    4.24424248742567    -1.74691044949975    -6.48374510684418    #1C7D54
#Astrocytes    ATTCAGAAGCATCGCGCAGCCAGACTCTATGGCCTA    3.60301401455387    -1.96493138894082    -6.47136162049336    #1C7D54
#Astrocytes    TCAACGAGTTCGCGATGGTCAGAGCCCGCCGATATC    5.51775913941571    -1.87741656898663    -6.76243310557264    #1C7D54
out_3d_dat<-as.data.frame(cbind(out_3d@meta.data[,c("predicted.id")],row.names(out_3d@meta.data),Embeddings(out_3d,"harmonyumap")))
colnames(out_3d_dat)[1]<-"predicted.id"
col_dat<-as.data.frame(cbind(names(type_cols),unname(type_cols)))
colnames(col_dat)<-c("predicted.id","cols")
dat_out<-merge(out_3d_dat,col_dat,by="predicted.id")
write.table(dat_out,file="multiome_tumor.tsv",sep="\t",quote=F,col.names=F,row.names=F)
system("slack -F multiome_tumor.tsv ryan_todo")

```



<!-- NOT GOING TO WORK, TOO LOW COV
### Looking for subclonality in bulk WGS libraries

Install Battenberg
https://github.com/Wedge-lab/battenberg
Reference files downloaded from https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52

```bash
R -q -e 'BiocManager::install(c("igordot/copynumber"))'
R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'
R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'
mkdir /home/groups/CEDAR/mulqueen/ref/battenberg
cd /home/groups/CEDAR/mulqueen/ref/battenberg

#placed reference files here after direct download with SFTP
```

### Pysam to make pseudobulk bam files
https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/
```python
import pysam
import csv

cluster_dict = {}
with open('clusters.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())


fin = pysam.AlignmentFile("atac_v1_pbmc_5k_possorted_bam.bam", "rb")

# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb
fouts_dict = {}
for cluster in clusters:
    fout = pysam.AlignmentFile("cluster" + cluster + ".bam", "wb", template = fin)
    fouts_dict[cluster] = fout

for read in fin:
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    else: 
        continue
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)

## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
```

-->

<!--
#ER binding poor and good outcome from patients, overlap with ATAC data
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Poor%20outcome%20ER%20regions.bed.gz
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Good%20outcome%20ER%20regions.bed.gz


-->


## R Session Info with all packages loaded

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
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
 [1] splines   grid      parallel  stats4    stats     graphics  grDevices
 [8] utils     datasets  methods   base     

other attached packages:
  [1] harmony_0.1.0                           
  [2] cowplot_1.1.1                           
  [3] iC10_1.5                                
  [4] iC10TrainingData_1.3.1                  
  [5] impute_1.64.0                           
  [6] pamr_1.56.1                             
  [7] survival_3.3-1                          
  [8] ggrepel_0.9.1                           
  [9] motifmatchr_1.12.0                      
 [10] chromVAR_1.12.0                         
 [11] forcats_0.5.1                           
 [12] purrr_0.3.4                             
 [13] readr_2.1.3                             
 [14] tidyverse_1.3.2                         
 [15] ggalluvial_0.12.3                       
 [16] dendextend_1.16.0                       
 [17] philentropy_0.6.0                       
 [18] CopyscAT_0.30                           
 [19] jsonlite_1.8.0                          
 [20] gplots_3.1.3                            
 [21] tibble_3.1.8                            
 [22] tidyr_1.1.4                             
 [23] edgeR_3.32.1                            
 [24] changepoint_2.2.3                       
 [25] zoo_1.8-10                              
 [26] FNN_1.1.3                               
 [27] Rtsne_0.16                              
 [28] fastcluster_1.2.3                       
 [29] NMF_0.24.0                              
 [30] bigmemory_4.5.36                        
 [31] cluster_2.1.3                           
 [32] rngtools_1.5.2                          
 [33] pkgmaker_0.32.2                         
 [34] registry_0.5-1                          
 [35] viridis_0.6.2                           
 [36] viridisLite_0.4.1                       
 [37] copykat_1.1.0                           
 [38] CaSpER_0.2.0                            
 [39] GOstats_2.56.0                          
 [40] graph_1.68.0                
  [41] Category_2.56.0                         
 [42] GO.db_3.12.1                            
 [43] limma_3.46.0                            
 [44] biomaRt_2.46.3                          
 [45] ape_5.6-2                               
 [46] ggnetwork_0.5.10                        
 [47] intergraph_2.0-2                        
 [48] igraph_1.3.0                            
 [49] gridExtra_2.3                           
 [50] scales_1.2.1                            
 [51] ggpubr_0.4.0                            
 [52] mclust_5.4.10                           
 [53] reshape_0.8.9                           
 [54] pheatmap_1.0.12                         
 [55] signal_0.7-7                            
 [56] Rcpp_1.0.9                              
 [57] infercnv_1.6.0                          
 [58] cicero_1.8.1                            
 [59] Gviz_1.34.1                             
 [60] monocle_2.18.0                          
 [61] DDRTree_0.1.5                           
 [62] irlba_2.3.5                             
 [63] VGAM_1.1-6                              
 [64] Matrix_1.4-1                            
 [65] BiocParallel_1.24.1                     
 [66] TFBSTools_1.28.0                        
 [67] JASPAR2020_0.99.10                      
 [68] seriation_1.3.6                         
 [69] dplyr_1.0.9                             
 [70] AUCell_1.13.3                           
 [71] TxDb.Hsapiens.UCSC.hg38.knownGene_3.10.0
 [72] org.Hs.eg.db_3.12.0                     
 [73] cisTopic_0.3.0                          
 [74] SeuratWrappers_0.3.0                    
 [75] stringr_1.4.0                           
 [76] EnsDb.Hsapiens.v86_2.99.0               
 [77] ensembldb_2.14.0                        
 [78] AnnotationFilter_1.14.0                 
 [79] GenomicFeatures_1.42.2                  
 [80] AnnotationDbi_1.52.0                    
 [81] Biobase_2.50.0                          
 [82] sp_1.5-0                                
 [83] SeuratObject_4.1.0                      
 [84] Seurat_4.1.1                            
 [85] Signac_1.5.0                            
 [86] SoupX_1.6.1                             
 [87] RColorBrewer_1.1-3                      
  [88] circlize_0.4.15                         
 [89] reshape2_1.4.4                          
 [90] ComplexHeatmap_2.6.2                    
 [91] patchwork_1.1.1                         
 [92] ggplot2_3.3.6                           
 [93] doParallel_1.0.17                       
 [94] iterators_1.0.14                        
 [95] foreach_1.5.2                           
 [96] BSgenome.Hsapiens.UCSC.hg38_1.4.3       
 [97] WGSmapp_1.2.0                           
 [98] SCOPE_1.2.0                             
 [99] BSgenome.Hsapiens.UCSC.hg19_1.4.3       
[100] BSgenome_1.58.0                         
[101] rtracklayer_1.50.0                      
[102] Rsamtools_2.6.0                         
[103] Biostrings_2.58.0                       
[104] XVector_0.30.0                          
[105] GenomicRanges_1.42.0                    
[106] GenomeInfoDb_1.26.7                     
[107] IRanges_2.24.1                          
[108] S4Vectors_0.28.1                        
[109] BiocGenerics_0.36.1                     
[110] HMMcopy_1.32.0                          
[111] data.table_1.14.4                       

loaded via a namespace (and not attached):
  [1] pbapply_1.5-0               haven_2.5.0                
  [3] lattice_0.20-45             vctrs_0.5.0                
  [5] expm_0.999-6                fastICA_1.2-3              
  [7] mgcv_1.8-40                 RBGL_1.66.0                
  [9] blob_1.2.3                  spatstat.data_2.2-0        
 [11] later_1.3.0                 DBI_1.1.3                  
 [13] R.utils_2.12.0              SingleCellExperiment_1.12.0
 [15] rappdirs_0.3.3              uwot_0.1.11                
 [17] jpeg_0.1-9                  zlibbioc_1.36.0            
 [19] rgeos_0.5-9                 htmlwidgets_1.5.4          
 [21] mvtnorm_1.1-3               GlobalOptions_0.1.2        
 [23] future_1.26.1               leiden_0.4.2               
 [25] KernSmooth_2.23-20          DT_0.23                    
 [27] promises_1.2.0.1            DelayedArray_0.16.3        
 [29] Hmisc_4.7-0                 fs_1.5.2                   
 [31] fastmatch_1.1-3             RhpcBLASctl_0.21-247.1     
 [33] digest_0.6.30               png_0.1-7                  
 [35] rjags_4-13                  qlcMatrix_0.9.7            
 [37] sctransform_0.3.3           pkgconfig_2.0.3            
 [39] docopt_0.7.1                gridBase_0.4-7             
 [41] spatstat.random_2.2-0       statnet.common_4.6.0       
 [43] lgr_0.4.3                   reticulate_1.25            
 [45] SummarizedExperiment_1.20.0 network_1.17.2             
 [47] modeltools_0.2-23           GetoptLong_1.0.5           
 [49] xfun_0.31                   tidyselect_1.2.0           
 [51] DNAcopy_1.64.0              ica_1.0-3                  
 [53] snow_0.4-4                  rlang_1.0.6                
 [55] glue_1.6.2                  modelr_0.1.9               
 [57] lambda.r_1.2.4              text2vec_0.6.1             
 [59] CNEr_1.26.0                 matrixStats_0.62.0         
 [61] MatrixGenerics_1.2.1        ggseqlogo_0.1              
 [63] ggsignif_0.6.3              httpuv_1.6.5               
 [65] class_7.3-20                TH.data_1.1-1              
 [67] seqLogo_1.56.0              annotate_1.68.0            
 [69] bit_4.0.4                   mime_0.12                  
 [71] Exact_3.1                   stringi_1.7.5              
 [73] RcppRoll_0.3.0              spatstat.sparse_2.1-1      
 [75] scattermore_0.8             bitops_1.0-7               
 [77] cli_3.4.1                   RSQLite_2.2.8              
 [79] bigmemory.sri_0.1.3         libcoin_1.0-9              
 [81] rstudioapi_0.13             TSP_1.2-0                  
 [83] GenomicAlignments_1.26.0    nlme_3.1-158               
 [85] locfit_1.5-9.4              VariantAnnotation_1.36.0   
 [87] listenv_0.8.0               SnowballC_0.7.0            
 [89] miniUI_0.1.1.1              R.oo_1.25.0                
 [91] dbplyr_2.2.1                readxl_1.4.0               
 [93] lifecycle_1.0.3             munsell_0.5.0              
 [95] cellranger_1.1.0            R.methodsS3_1.8.2          
 [97] caTools_1.18.2              codetools_0.2-18           
 [99] coda_0.19-4                 lmtest_0.9-40              
[101] htmlTable_2.4.1             xtable_1.8-4               
[103] ROCR_1.0-11                 googlesheets4_1.0.1        
[105] formatR_1.12                BiocManager_1.30.18        
[107] abind_1.4-5                 farver_2.1.1               
[109] rsparse_0.5.0               parallelly_1.32.0          
[111] RANN_2.6.1                  askpass_1.1                
[113] biovizBase_1.38.0           poweRlaw_0.70.6            
[115] sparsesvd_0.2               RcppAnnoy_0.0.19           
[117] goftest_1.2-3               futile.options_1.0.1       
[119] dichromat_2.0-0.1           future.apply_1.9.0         
[121] ellipsis_0.3.2              prettyunits_1.1.1          
[123] reprex_2.0.2                lubridate_1.8.0            
[125] googledrive_2.0.0           ggridges_0.5.3             
[127] mlapi_0.1.1                 remotes_2.4.2              
[129] slam_0.1-50                 gargle_1.2.1               
[131] argparse_2.1.5              spatstat.utils_2.3-1       
[133] doSNOW_1.0.20               htmltools_0.5.2            
[135] BiocFileCache_1.14.0        utf8_1.2.2                 
[137] plotly_4.10.0               XML_3.99-0.9               
[139] e1071_1.7-11                foreign_0.8-82             
[141] withr_2.5.0                 fitdistrplus_1.1-8         
[143] bit64_4.0.5                 rootSolve_1.8.2.3          
[145] multcomp_1.4-19             ProtGenerics_1.22.0        
[147] spatstat.core_2.4-4         combinat_0.0-8             
[149] progressr_0.10.1            rsvd_1.0.5                 
[151] memoise_2.0.1               arrow_5.0.0.2              
[153] tzdb_0.3.0                  lmom_2.9                   
[155] curl_4.3.2                  fansi_1.0.3                
[157] GSEABase_1.52.1             tensor_1.5                 
[159] checkmate_2.1.0             float_0.3-0                
[161] cachem_1.0.6                deldir_1.0-6               
[163] rjson_0.2.21                rstatix_0.7.0              
[165] clue_0.3-61                 tools_4.0.3                
[167] sandwich_3.0-2              magrittr_2.0.3             
[169] RCurl_1.98-1.9              proxy_0.4-26               
[171] car_3.1-0                   TFMPvalue_0.0.8            
[173] xml2_1.3.3                  httr_1.4.3                 
[175] assertthat_0.2.1            boot_1.3-28                
[177] globals_0.15.1              R6_2.5.1                   
[179] nnet_7.3-17                 genefilter_1.72.1          
[181] DirichletMultinomial_1.32.0 progress_1.2.2             
[183] KEGGREST_1.30.1             gtools_3.9.3               
[185] shape_1.4.6                 coin_1.4-2                 
[187] lsa_0.73.3                  carData_3.0-5              
[189] colorspace_2.0-3            generics_0.1.3             
[191] base64enc_0.1-3             pracma_2.3.8               
[193] pillar_1.8.1                Rgraphviz_2.34.0           
[195] tweenr_1.0.2                HSMMSingleCell_1.10.0      
[197] GenomeInfoDbData_1.2.4      plyr_1.8.7                 
[199] gtable_0.3.1                futile.logger_1.4.3        
[201] rvest_1.0.3                 RcisTarget_1.11.10         
[203] knitr_1.37                  latticeExtra_0.6-29        
[205] fastmap_1.1.0               Cairo_1.5-12.2             
[207] broom_1.0.0                 openssl_2.0.1              
[209] backports_1.4.1             densityClust_0.3.2         
[211] feather_0.3.5               gld_2.6.5                  
[213] hms_1.1.2                   ggforce_0.3.3              
[215] shiny_1.7.1                 polyclip_1.10-0            
[217] DescTools_0.99.45           lazyeval_0.2.2             
[219] lda_1.4.2                   Formula_1.2-4              
[221] crayon_1.5.2                MASS_7.3-57                
[223] AnnotationForge_1.32.0      rpart_4.1.16               
[225] compiler_4.0.3              spatstat.geom_2.4-0        

```

{% endcapture %} {% include details.html %} 
