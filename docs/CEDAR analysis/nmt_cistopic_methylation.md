---
title: CisTopic on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /cistopic_nmt/
category: CEDAR
---

## Set up genomic annotations

Generate bed files for annotation regions to sum methylation data over.

{% capture summary %} Code {% endcapture %} {% capture details %}  

Gene location information stored here:

```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes
```

Generation of bed files:
```bash
awk 'OFS="\t" {split($10,a,"\;"); split(a[1],b,"\""); if($3=="gene") print $1,$4,$5,b[2]}' genes.gtf > genes.bed
```

Regulatory information stored here:
```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds
```

Some of these features have redundancies (i.e. same start and end sites but different names). Have to filter that out.

Generation of bed files:
```bash
#Regulatory information from ensembl regulatory build
wget http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz

zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="enhancer") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > enhancers.bed
zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="promoter") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > promoters.bed
zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="CTCF_binding_site") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > ctcf.bed

awk 'OFS="\t" {print $1,1,$2}' /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt | bedtools makewindows -b - -w 100000 | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}'> 100kb.bed


``` 
{% endcapture %} {% include details.html %} 

List of genomic annotation bed files:
```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/ctcf.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed
```

## Set up methylation extraction

Perform methylation extraction with bismark methylation extractor. Output cov.gz files.

Here I'm using 10 nodes, 20 cpus per node (since parallelization is 3x5), with 20gb per node. It's a lot of resources, but hopefully it will finish fast.

{% capture summary %} Code {% endcapture %} {% capture details %}  

meth_extract.slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-374
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=12:00:00
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/dedup_bams/*deduplicated.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun bismark_methylation_extractor \
-s --gzip --parallel 5 \
-o /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out \
--genome_folder /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta \
--CX \
--bedGraph \
--buffer_size 1G \
--no_header \
$files_in

```

```bash
sbatch meth_extract.slurm.sh
```

{% endcapture %} {% include details.html %} 

Generate cytosine NOME data for region summarization.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --array=1-374
#SBATCH --tasks-per-node=10 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=5gb 
#SBATCH --time=10:00:00 
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
out_name=`echo $files_in | awk '{n=split($1,a,"/");print a[n]}' | sed s/".cov.gz"//`

srun coverage2cytosine \
--gzip  \
--dir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out \
--genome_folder /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta \
--output $out_name \
--nome-seq \
$files_in
```

```bash
sbatch nome_extract.slurm.sh
```

{% endcapture %} {% include details.html %} 

## Using python pandas and bedtools wraper to summarize methylation over regions

* argv1 is a cov.gz output from coverage2cytosine report (either CG or GpC)
* argv2 is a bed file with features to aggregate methylation over
* argv3 is output directory
* argv4 is a prefix used for the bed file annotation. Example: "genebody" or "promoter"

Note: argv2 should be a bed file of format [chr]\t[start]\t[end]\t[feature_name]

{% capture summary %} Code {% endcapture %} {% capture details %}

```python
#!/usr/bin/python
import pybedtools
import sys
import pandas as pd
import gzip
import os.path
pybedtools.helpers.set_tempdir("/home/groups/CEDAR/mulqueen/temp")

if len(sys.argv) != 5:
        print("""
* argv1 is a cov.gz output from coverage2cytosine report (either CG or GpC)
* argv2 is a bed file with features to aggregate methylation over
* argv3 is output directory
* argv4 is a prefix used for the bed file annotation. Example: "genebody" or "promoter"

Note: argv2 should be a bed file of format [chr]\t[start]\t[end]\t[feature_name]""")
        sys.exit(0)


in_list=[]
in_list.append(sys.argv[1])
#in_list.append("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/T_E8_G02_merged.deduplicated.bismark.NOMe.CpG.cov.gz")
in_list.append(sys.argv[2])
#in_list.append("/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed")
in_list.append(sys.argv[3])
#in_list.append("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")
in_list.append(sys.argv[4])
#in_list.append("genes")

outname=os.path.basename(in_list[0])[:len(os.path.basename(in_list[0]))-7]
scmet=pd.read_csv(gzip.open(in_list[0], "rb"),header=None,sep="\t")[[0,1,2,4,5]]
scmet.columns = ["chr", "start", "end", "met","unmet"]

ref = pybedtools.BedTool(in_list[1])
scmet = pybedtools.BedTool.from_dataframe(scmet)

scmet_intersect=scmet.intersect(ref,wa=True,wb=True,nonamecheck=True).to_dataframe(disable_auto_names=True, header=None)[[5,6,7,8,3,4]]
scmet_intersect.columns = ["chr", "start", "end", "feat","met","unmet"]
scmet_intersect["total_sites"]=scmet_intersect["met"]+scmet_intersect["unmet"]
out_dataframe=scmet_intersect.groupby(by=["chr", "start", "end", "feat"]).agg(met_cg=pd.NamedAgg(column="met",aggfunc="sum"),total_cg=pd.NamedAgg(column="total_sites",aggfunc="sum")).reset_index()

if in_list[2].endswith("/"):
        out_dataframe.to_csv(in_list[2]+outname+"."+in_list[3]+".count.txt.gz",header=False,index=False,compression='gzip')
else:
        out_dataframe.to_csv(in_list[2]+"/"+outname+"."+in_list[3]+".count.txt.gz",header=False,index=False,compression='gzip')

```
Example running:

```bash
python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py [argv1] [argv2] [arg3] [argv4]
```

{% endcapture %} {% include details.html %}

## Run as a slurm batch jobs.

{% capture summary %} CpG Gene body {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*CpG.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
gene 

```

{% endcapture %} {% include details.html %} 

{% capture summary %} CpG Promoters {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*CpG.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
promoter

```

{% endcapture %} {% include details.html %} 

{% capture summary %} CpG Enhancers {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*CpG.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
enhancer

```

{% endcapture %} {% include details.html %} 

{% capture summary %} CpG 100kb {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*CpG.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
100kb

```

{% endcapture %} {% include details.html %} 

{% capture summary %} GpC Gene body {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*GpC.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
gene 

```

{% endcapture %} {% include details.html %} 

{% capture summary %} GpC Promoters {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*GpC.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
promoter

```

{% endcapture %} {% include details.html %} 

{% capture summary %} GpC Enhancers {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*GpC.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
enhancer

```

{% endcapture %} {% include details.html %} 

{% capture summary %} GpC 100kb {% endcapture %} {% capture details %}  

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=10 ##we want our node to do N tasks at the same time
#SBATCH --array=1-374
#SBATCH --cpus-per-task=3 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

files_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out/*GpC.cov.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun python /home/groups/CEDAR/mulqueen/src/aggregate_methylation_over_region.py \
$files_in \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
100kb

```

{% endcapture %} {% include details.html %} 

Running all of these as batch jobs

```bash
sbatch CpG_genes.slurm.sh #done
sbatch CpG_enhancers.slurm.sh  #done
sbatch CpG_promoters.slurm.sh  
sbatch CpG_100kb.slurm.sh  #done

sbatch GpC_genes.slurm.sh
sbatch GpC_enhancers.slurm.sh
sbatch GpC_promoters.slurm.sh
sbatch GpC_100kb.slurm.sh  

```

## Running cistopic on methylated regions.

Making an R script for cistopic processing of methylation files.

* argv1 is a bed file with features to aggregate methylation over. Same as aggregate_methylation_over_region.py input.
* argv2 is the methylation moeity. This is set up in argv3 of the python script. Example is \*[GpC|CpG].[prefix].count.txt.
* argv3 is prefix from aggregate_methylation_over_region.py output. This is set up in argv3 of the python script. Example is \*.[prefix].count.txt.


```bash
Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R [argv1] [argv2] [argv3]
```

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(stats)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")

#args <- commandArgs(trailingOnly=TRUE)

args<-list()
args[1]<-"/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed" 
args[2]<-"CpG"
args[3]<-"enhancer"
args<-unlist(args)

meth_files<-list.files(pattern=paste0(args[2],".",args[3],".count.txt$")) #use prefix to determine meth files to read in
region<-args[1] #region file to read through


#create cistopic object from methylation files from https://github.com/aertslab/cisTopic/blob/04cecbb9d1112fcc1a6edc28b5a506bcb49f2803/R/InitializecisTopic.R
#i modified to fit file structure

createcisTopicObjectFromMeth <- function(
  methfiles,
  regions,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 0.5, #this could probably be modified as well
  ...
) {
  # Prepare annotation
  #tester regions<-region
  #tester methfiles<-meth_files
  regions_frame <- read.table(regions)[,c(1:4)]
  colnames(regions_frame) <- c('seqnames', 'start', 'end','feat')

  regions_frame<-regions_frame[!duplicated(regions_frame[,c(1:3)]),] #ensure assayed regions are unique

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
  regions_frame <- read.table(regions)[,c(1:4)] #repopulate regions_frame for annotation field
  regions_frame<-regions_frame[!duplicated(regions_frame),] #ensure assayed regions are unique
  region.data<-cbind(region.data,regions_frame[4])#add in gene name
  object <- createcisTopicObject(count.matrix = beta.matrix, project.name = project.name, min.cells = min.cells, min.regions = min.regions, is.acc = is.acc, ...)
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  object <- addRegionMetadata(object, region.data = as.data.frame(region.data))
  return(object)
}

dat<-createcisTopicObjectFromMeth(methfiles=meth_files,regions=region)
saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.Rds"))

dat <- runWarpLDAModels(dat, topic=c(5:15, 20, 25, 40, 50), seed=123, nCores=14, addModels=FALSE,tmp="/home/groups/CEDAR/mulqueen/temp/")
saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.models.Rds"))
dat<-readRDS(file=paste0(args[2],".",args[3],".cistopic_object.models.Rds"))

#select best model
pdf(paste0(argv[1],".cistopic_model_selection.pdf"))
par(mfrow=c(3,3))
dat <- selectModel(dat, type='maximum')
dat <- selectModel(dat, type='perplexity')
dat <- selectModel(dat, type='derivative')
dev.off()

#going with 15 topics based on derivative

dat<-cisTopic::selectModel(dat,select=15,keepModels=T)

dat <- runUmap(dat, target='cell') #running umap using cistopics implementation

#add sample cell line names as metadata
dat@cell.data$cellLine<-unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",4))

#set up treatment conditions
dat@cell.data$treatment<-substr(unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",5)),1,1) #set up T47D cell treatment first
dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$cellLine,3,3)
dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$cellLine,5,5)

dat@cell.data$cellLine_treatment<-paste(dat@cell.data$cellLine,dat@cell.data$treatment,sep="_")

dat@cell.data$perc_met<-dat@cell.data$"Methylated reads"/dat@cell.data$"Total reads"

#clean up metadata a bit more
dat@cell.data[dat@cell.data$cellLine %in% c("BSM7E6","M7C1A","M7C2B","M7E4C"),]$cellLine<-"MCF7"
dat@cell.data[dat@cell.data$cellLine %in% c("T"),]$cellLine<-"T47D"
dat@cell.data[dat@cell.data$treatment %in% c("C"),]$treatment<-"control"
dat@cell.data[dat@cell.data$treatment %in% c("E"),]$treatment<-"estrogen"

pdf(paste0(args[2],".",args[3],".cistopic_clustering.pdf")
par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine','treatment','cellLine_treatment'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)

par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('Total reads','perc_met'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)


par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.Rds"))


#should probably make this Z-scored for plotting
topicmat<-as.data.frame(dat@selected.model$document_expects)
row.names(dat@cell.data) == colnames(topicmat) # order is maintained so using cell.data info for annotation
cellline_annot<-dat@cell.data$cellLine
treatment_annot<-dat@cell.data$treatment
c_count_annot<-dat@cell.data$"Total reads"
cellline_treatment_annot<-dat@cell.data$cellLine_treatment

#color function
col_fun=colorRamp2(quantile(unlist(topicmat),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))

row.names(topicmat)<-paste0("Topic_",row.names(topicmat))
plt1<-Heatmap(topicmat,
    show_column_names=F,
    bottom_annotation=columnAnnotation(cellline=cellline_annot,treatment=treatment_annot,C_covered=c_count_annot,cellline_treatment=cellline_treatment_annot,
    	    col = list(cellline = c("T47D"="red","M7"="blue"),
               			treatment = c("control" = "white", "estrogen" = "black"),
               			cellline_treatment = c("BSM7E6_E"="green","M7C1A_C"="brown","M7C2B_C"="purple","M7E4C_E"="blue","T_C"="orange","T_E"="red"))))

pdf(paste0(args[2],".",args[3],".cistopic_heatmap.pdf")
plt1
dev.off()


#testing topic significance for cell line/ treatment
dat_stat<-rbind(topicmat,cellline_annot,treatment_annot)
dat_stat<-as.data.frame(t(dat_stat))
dat_stat[1:15] <- lapply(dat_stat[1:15], as.numeric)
dat_stat<-melt(dat_stat)
colnames(dat_stat)<-c("cellline","treatment","topic","value")

#testing topic specificity by treatment and cell line
out_stats<-as.data.frame(dat_stat %>% 
  group_by(topic) %>% 
  anova_test(value ~ treatment*cellline) %>% adjust_pvalue()
  )

write.table(out_stats,file=paste0(args[2],".",args[3],".cistopic_topic_enrichment.txt"),quote=F,sep="\t",col.names=T)



```
{% endcapture %} {% include details.html %} 

```bash
Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R [argv1] [argv2] [argv3]

#CpG features
Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed \
CpG \
gene

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed \
CpG \
promoter

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed \
CpG \
enhancer

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed \
CpG \
100kb

#GpC features
Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed \
GpC \
gene

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed \
GpC \
promoter

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed \
GpC \
enhancer

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed \
GpC \
100kb

```

<!---
Additional processing for topic heatmap generation
Following this (tutorial.)[https://www.neb.com/products/n0356-dm5ctp#Product%20Information]
```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(stats)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov")

dat<-readRDS(file="cistopic_object.Rds")

#predictive measurements of topics
pred_matrix <- predictiveDistribution(dat)

#get regions from topics
dat <- getRegionsScores(dat, method='NormTop', scale=TRUE)

#plot general annotations to topics 
dat <- annotateRegions(dat, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

#binarize cistopics to get genes driving topics
dat<-binarizecisTopics(dat,thrP=0.975)
getBedFiles(dat, path='cisTopic_beds')

dat <- GREAT(dat, genome='hg38', fold_enrichment=1.5, geneHits=1, sign=0.05, request_interval=10)
#We can visualize the enrichment results:
pdf("cistopic_GO.pdf",width=30,height=30)
ontologyDotPlot(dat, top=20, topics=c(1:15), var.y='name', order.by='Binom_Adjp_BH')
dev.off()

#Save cistopic object
saveRDS(dat,file="cistopic_object.Rds")
```
{% endcapture %} {% include details.html %} 
--->

To do

* re run on genome-wide bins
* re run on defined enhancers
* re run on promoters
* re run on accessibility (GC met)
* see if cistopic needs to be binarized
* get gene lists per topic
* change methylation rate cutoffs per gene