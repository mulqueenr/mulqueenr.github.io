---
title: CisTopic on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /cistopic_nmt/
category: CEDAR
---

### Script to test cistopic modeling for methylation regions ###
```python
#remotes::install_version("RSQLite", version = "2.2.5") #https://stackoverflow.com/questions/67279457/error-with-r-package-biomart-and-this-dependency-rsqlite
#devtools::install_github("aertslab/RcisTarget")
#devtools::install_github("aertslab/AUCell")
#devtools::install_github("aertslab/cisTopic")
```
#Set up methylation extraction
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=30 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=20gb ## request gigabyte per cpu
#SBATCH --time=5:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/dedup_bams/*deduplicated.bam`

srun bismark_methylation_extractor \
-s --gzip --parallel 10 \
-o /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov \
--yacht \
--buffer_size 20G \
--no_header \
$files_in
```
# python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py [argv1] [argv2]

argv1 is yacht output from bismark_methylation_extraction
argv2 is a bed file with features to aggregate methylation over

```python
#!/usr/bin/python
import pybedtools
import sys
import pandas as pd
import gzip
pybedtools.helpers.set_tempdir("/home/groups/CEDAR/mulqueen/temp")

in_list=[]
in_list.append(sys.argv[1])
in_list.append(sys.argv[2])

scmet=pd.read_csv(gzip.open(in_list[0], "rb"),header=None,sep="\t")[[2,3,3,4]]
scmet.columns = ["Chr", "Start", "End", "Met"]

scmet=scmet[scmet["Met"].isin(["Z","z"])]
ref = pybedtools.BedTool(in_list[1])
scmet = pybedtools.BedTool.from_dataframe(scmet)

scmet_intersect=scmet.intersect(ref,wa=True,wb=True,nonamecheck=True).to_dataframe()[["score","strand","thickStart","thickEnd","name"]]
scmet_intersect.columns = ["chr", "start", "end", "feat","met"]
met_count=scmet_intersect[scmet_intersect["met"]=="Z"]
met_count=met_count.groupby(by=["chr", "start", "end", "feat"]).size().to_frame('met_cg').reset_index() #get count of unique rows (methylation calls across groups)

total_count=scmet_intersect.groupby(by=["chr", "start", "end", "feat"]).size().to_frame('total_cg').reset_index() #get count of unique rows (methylation calls across groups)

out_dataframe=pd.merge(total_count,met_count,on=["chr","start","end","feat"])

out_dataframe.to_csv(in_list[0].split("/")[-1].split(".")[0]+"."+in_list[1].split("/")[-1].split(".")[0]+".count.txt",header=False,index=False)
```

sbatch aggregate_over_genes.slurm.sh

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov/*deduplicated.txt.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py $files_in /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed

```

```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(stats)
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov")

meth_files<-list.files(pattern=".genes.count.txt$")
region<-"/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed"

#create cistopic object from methylation files from https://github.com/aertslab/cisTopic/blob/04cecbb9d1112fcc1a6edc28b5a506bcb49f2803/R/InitializecisTopic.R
#i modified to fit file structure

createcisTopicObjectFromMeth <- function(
  methfiles,
  regions,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 0.5,
  ...
) {
  # Prepare annotation
  #tester regions<-region
  #tester methfiles<-meth_files
  regions_frame <- read.table(regions)[,c(1:3)]
  colnames(regions_frame) <- c('seqnames', 'start', 'end')
  rownames(regions_frame) <- paste(regions_frame$seqnames, ':', regions_frame$start, '-', regions_frame$end, sep='')
  regions_granges <- makeGRangesFromDataFrame(as.data.frame(regions_frame))

  # Prepare beta matrix, cell data and region data
  beta.matrix <- Matrix(0, nrow(regions_frame), length(methfiles), sparse=TRUE)
  rownames(beta.matrix) <- rownames(regions_frame)
  colnames(beta.matrix) <- methfiles

  cell.data <- matrix(0, length(methfiles), 2)
  rownames(cell.data) <- methfiles
  colnames(cell.data) <- c('Methylated reads', 'Total reads')

  region.data <- matrix(0, nrow(regions_frame), 2)
  rownames(region.data) <- rownames(regions_frame)
  colnames(region.data) <- c('Methylated reads', 'Total reads')

  for (file in methfiles){
  	#tester file file<-methfiles[1]
    # Read data
    print(paste('Reading file ', file, '...', sep=''))
    sample <- read.table(file,sep=",")
    sample <- sample[, c(1,2,3,6,5)] #4 is feature name, i also put total reads first in my data frames
    colnames(sample) <- c('seqnames', 'start', 'end', 'Methylated', 'Total_reads')
    sample_granges <- makeGRangesFromDataFrame(as.data.frame(sample), keep.extra.columns=TRUE)

    # Fill cell data
    cell.data[file, 1] <- sum(sample_granges$Methylated)
    cell.data[file, 2] <- sum(sample_granges$Total_reads)

    # Calculate aggregate reads per region (methylated and total)
    overlaps <- findOverlaps(sample_granges, regions_granges)
    sites <- sample_granges[queryHits(overlaps)]
    sumMethylated <- aggregate(sites$Methylated, list(subjectHits(overlaps)), sum)
    sumReads <- aggregate(sites$Total_reads, list(subjectHits(overlaps)), sum)

    # Fill region data
    rownames(sumMethylated) <- rownames(regions_frame)[sumMethylated[,1]]
    region.data[rownames(sumMethylated),1] <- region.data[rownames(sumMethylated),1] + sumMethylated[,2]

    rownames(sumReads) <- rownames(regions_frame)[sumReads[,1]]
    region.data[rownames(sumReads),2] <- region.data[rownames(sumReads),2] + sumReads[,2]

    # Calculate beta
    aggrBeta <- sumMethylated[,2]/sumReads[,2]
    names(aggrBeta) <- rownames(regions_frame)[sumMethylated[,1]]
    beta.matrix[names(aggrBeta), file] <- aggrBeta
  }
  region.data<-cbind(region.data,read.table(regions)[,4])#add in gene name
  object <- createcisTopicObject(count.matrix = beta.matrix, project.name = project.name, min.cells = min.cells, min.regions = min.regions, is.acc = is.acc, ...)
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  object <- addRegionMetadata(object, region.data = as.data.frame(region.data))
  return(object)
}

dat<-createcisTopicObjectFromMeth(methfiles=meth_files,regions=region)
saveRDS(dat,file="cistopic_object.Rds")

dat <- runWarpLDAModels(dat, topic=c(5:15, 20, 25, 40, 50), seed=123, nCores=14, addModels=FALSE,tmp="/home/groups/CEDAR/mulqueen/temp/")
saveRDS(dat,file="cistopic_object.models.Rds")

```