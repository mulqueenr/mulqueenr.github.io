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
awk 'OFS="\t" {print $1,1,$2}' /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt | bedtools makewindows -b - -w 10000 | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}'> 10kb.bed


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

{% capture summary %} CpG 10kb {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
10kb

```

{% endcapture %} {% include details.html %} 

{% capture summary %} CpG Breast Cancer Enhancers {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
bcEnhance

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


{% capture summary %} GpC 10kb {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
10kb

```

{% endcapture %} {% include details.html %} 

{% capture summary %} GpC Breast Cancer Enhancers {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
bcEnhance

```

{% endcapture %} {% include details.html %} 

Running all of these as batch jobs

```bash
sbatch CpG_genes.slurm.sh 
sbatch CpG_enhancers.slurm.sh  
sbatch CpG_promoters.slurm.sh  
sbatch CpG_100kb.slurm.sh  
sbatch CpG_10kb.slurm.sh  
sbatch CpG_bcEnhancer.slurm.sh

sbatch GpC_genes.slurm.sh
sbatch GpC_enhancers.slurm.sh
sbatch GpC_promoters.slurm.sh
sbatch GpC_100kb.slurm.sh  
sbatch GpC_100kb.slurm.sh  
sbatch GpC_bcEnhancer.slurm.sh
```

## Running cistopic on methylated regions.

Making an R script for cistopic processing of methylation files.

{% capture summary %} cistopic_methylation.R {% endcapture %} {% capture details %}  

```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(mefa4)
library(stats)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(patchwork)
library(rstatix)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Ckmeans.1d.dp)

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")

args <- commandArgs(trailingOnly=TRUE)

#args<-list()
#args[1]<-"/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed" 
#args[2]<-"CpG"
#args[3]<-"gene"
#args<-unlist(args)

meth_files<-list.files(pattern=paste0(args[2],".",args[3],".count.txt.gz$")) #use prefix to determine meth files to read in
region<-args[1] #region file to read through


#create cistopic object from methylation files from https://github.com/aertslab/cisTopic/blob/04cecbb9d1112fcc1a6edc28b5a506bcb49f2803/R/InitializecisTopic.R
#i modified to fit file structure

createcisTopicObjectFromMeth <- function(
  methfiles,
  regions,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 0.5, 
  return.norm.beta=T,
  ...
) {
  # Prepare annotation
  #tester regions<-region
  #tester methfiles<-meth_files
  #tester return.norm.beta=T
  regions_frame <- read.table(regions)[,c(1:4)]
  colnames(regions_frame) <- c('seqnames', 'start', 'end','feat')

  regions_frame<-regions_frame[!duplicated(regions_frame[,c(1:3)]),] #ensure assayed regions are unique
  regions_frame<-regions_frame[regions_frame$seqnames %in% paste0('chr',c(1:22,"X","Y")),] #remove alt chrs

  rownames(regions_frame) <- paste(regions_frame$seqnames, ':', regions_frame$start, '-', regions_frame$end, sep='')
  regions_granges <- makeGRangesFromDataFrame(as.data.frame(regions_frame))

  # Prepxare beta matrix, cell data and region data
  beta.matrix <- Matrix(0, nrow(regions_frame), length(methfiles), sparse=TRUE)
  rownames(beta.matrix) <- rownames(regions_frame)
  colnames(beta.matrix) <- methfiles

  coverage.matrix<-beta.matrix #prepare a coverage matrix of the same dims as beta matrix

  cell.data <- matrix(0, length(methfiles), 2)
  rownames(cell.data) <- methfiles
  colnames(cell.data) <- c('Methylated reads', 'Total reads')

  region.data <- matrix(0, nrow(regions_frame), 2)
  rownames(region.data) <- rownames(regions_frame)
  colnames(region.data) <- c('Methylated reads', 'Total reads')

  #read in single-cell data and format for scmet
  single_cell_in<-mclapply(methfiles,mc.cores=20, FUN= function(x) {
	  	#tester file x<-methfiles[1]
	    # Read data
	    print(paste('Reading file ', x, '...', sep=''))
	    sample <- fread(x,sep=",")
	    sample <- sample[, c(1,2,3,5,6)] 
	    colnames(sample) <- c('seqnames', 'start', 'end', 'Methylated', 'Total_reads')
	    sample_granges <- makeGRangesFromDataFrame(as.data.frame(sample), keep.extra.columns=TRUE)

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

	    total_met<-sum(region.data[,1])
	    total_count<-sum(region.data[,2])

	    # Calculate beta (mean methylation) old way of processing
  		aggrBeta <- sumMethylated[,2]/sumReads[,2]

      #calculate beta-binomial distribution to account for coverage
      if(return.norm.beta){
        #beta binomial distribution correction based on https://www.biorxiv.org/content/10.1101/2019.12.11.873398v1.full
        #sample mean methylation across features
        m=mean(sumMethylated[,2]/sumReads[,2])
        v=var(sumMethylated[,2]/sumReads[,2])
        #shape parameters of beta distribution by method of moments
        alpha=m*(m*(1-m)/v-1)
        beta=(1-m)*(m*(1-m)/v-1)
        #calculate posterior
        mhat_c=(alpha+sumMethylated[,2])/(alpha+beta+sumReads[,2])
        norm_mhat_c=mhat_c/(alpha/alpha+beta)
        aggrBeta<-norm_mhat_c
      }

      coverage_count<-region.data[,2]

  		names(aggrBeta) <- rownames(regions_frame)[sumMethylated[,1]]
  		return(list(total_met=total_met,total_count=total_count,beta=aggrBeta,coverage_count=region.data[,2]))
  })



#posterior rate calculation and normalization
  cell_met_count<-unlist(lapply(single_cell_in,"[[",1))
  cell_total_count<-unlist(lapply(single_cell_in,"[[",2))
  cell.data<-cbind(cell.data,CG_met_count=cell_met_count,CG_total_count=cell_total_count)
  beta.matrix_sc<-lapply(1:length(single_cell_in), FUN= function(x) single_cell_in[[x]][3])
  print("Plotting coverage plot per feature.")
  coverage.matrix_sc<-lapply(1:length(single_cell_in), FUN= function(x) single_cell_in[[x]][4])
  plt<-ggplot()+geom_histogram(aes(x=unlist(coverage.matrix_sc)),binwidth=1)+theme_bw()+xlim(c(0,100))
  ggsave(plt,file=paste0(args[2],".",args[3],".feature_coverage.pdf"))
  system(paste0("slack -F ",args[2],".",args[3],".feature_coverage.pdf"," ryan_todo"))

  #populate the sparse matrix with a for loop. probably not the best, but whatever
  print("Populating Sparse Beta Matrix.")
	for (i in 1:length(beta.matrix_sc)){
		beta_dat<-unlist(beta.matrix_sc[[i]]$beta,use.names=T)
		beta.matrix[names(beta_dat),methfiles[i]]<-unname(beta_dat)
	}

  print("Populating Sparse Coverage Matrix.")
        for (i in 1:length(coverage.matrix_sc)){
                cov_dat<-unlist(coverage.matrix_sc[[i]]$coverage_count,use.names=T)
                coverage.matrix[names(cov_dat),methfiles[i]]<-unname(cov_dat)
        }

  #compute CG density across regions
  #https://www.biostars.org/p/478444/
  #Get the sequences and compute the GC content
  print("Calculating region specific CG density.")
  region_seq<-getSeq(BSgenome.Hsapiens.UCSC.hg38, regions_granges)
  if(args[3]=="GpC"){
    freqs <- letterFrequency(region_seq,letters="GC")
    } else {
    freqs <- letterFrequency(region_seq,letters="CG")}
  nuc_density <- freqs/region_seq@ranges@width

  region.data<-cbind(region.data,regions_frame[4],nuc_density=nuc_density,cg_count=freqs)#add in gene name and cg density

  object <- createcisTopicObject(count.matrix = beta.matrix, project.name = project.name, min.cells = min.cells, min.regions = min.regions, is.acc = is.acc, ...)
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  object <- addRegionMetadata(object, region.data = as.data.frame(region.data))
  return(c(object,coverage.matrix))
}


dat_out<-createcisTopicObjectFromMeth(methfiles=meth_files,regions=region)
dat<-dat_out[[1]]
dat_cov<-dat_out[[2]]
print("Generating an automatic methylation cutoff by univariate k means.")
#generate plot of beta values per features
feat_dat<-Melt(dat@count.matrix)
count_dat<-Melt(as.matrix(dat_cov))

feat_var<-data.frame(var=apply(dat@count.matrix,1,sd),nuc_density=as.numeric(dat@region.data[,10]),region_cg_count=as.numeric(dat@region.data[,11]))

feat_dat2<-merge(feat_dat,feat_var,by.x="rows",by.y="row.names",all.x=T)
feat_dat2<-merge(feat_dat2,count_dat,by=c("rows","cols"))
colnames(feat_dat2)<-c("rows","cols","value","var","nuc_density","region_cg_count","cov")
row.names(feat_dat2)<-row.names(feat_dat)
feat_dat<-feat_dat2
feat_dat$perc_gc_covered<-(feat_dat$cov/feat_dat$region_cg_count)*100
#perform univariate clustering to determine optimal sample segmentation https://cran.r-project.org/web/packages/Ckmeans.1d.dp/vignettes/Ckmeans.1d.dp.html
#similar to this paper https://academic.oup.com/bioinformatics/article/36/20/5027/5866975?login=true
#provide cg density as weight?
clus<-Ckmeans.1d.dp(x=feat_dat$value,k=2,y=feat_dat$nuc_density)
cutoff_value<-max(feat_dat$value[clus$cluster==1])
print(paste0("Cutoff value set to: ",cutoff_value))

ylim_feat_cov<-quantile(probs=c(0,0.95),feat_dat$cov/feat_dat$region_cg_count)
plt1<-ggplot(feat_dat,aes(x=value,y=var))+geom_density_2d(aes(colour = after_stat(level)),binwidth = 3)+xlab("Methylation")+ylab("Feature SD")+scale_color_distiller(palette = 1,direction=-1)+geom_vline(xintercept=cutoff_value,color="red") +xlim(c(0,1))
plt2<-ggplot(feat_dat,aes(x=value,y=nuc_density))+geom_density_2d(aes(colour = after_stat(level)),binwidth = 3)+xlab("Methylation")+ylab("Feature Dinuc Density")+scale_color_distiller(palette = 2,direction=-1) + geom_vline(xintercept=cutoff_value,color="red") +xlim(c(0,1))
plt3<-ggplot(feat_dat,aes(x=value,y=perc_gc_covered))+geom_density_2d(aes(colour = after_stat(level)))+xlab("Methylation")+ylab("Feature Coverage Percent")+scale_color_distiller(palette = 3,direction=-1)+ylim(ylim_feat_cov)+geom_vline(xintercept=cutoff_value,color="red")+xlim(c(0,1))

ggsave(plt1/plt2/plt3,file=paste0(args[2],".",args[3],".feature_stats.pdf"))
system(paste0("slack -F ",paste0(args[2],".",args[3],".feature_stats.pdf"," ryan_todo"))

print("Binarize Data by automatic cutoff.")
backup<-dat@binary.count.matrix
count_mat<-as.matrix(dat@count.matrix)

# #Apply coverage filter, require at 5 or more CG measurements
 dat_cov<-as.matrix(dat_cov)
 dat_cov<-dat_cov[row.names(dat_cov) %in% row.names(count_mat),colnames(dat_cov) %in% colnames(count_mat)]
 count_mat[which(dat_cov<5,arr.ind=T)]<-0


count_mat[which(count_mat>=cutoff_value,arr.ind=T)]<-1
count_mat[which(count_mat<cutoff_value,arr.ind=T)]<-0

dat@binary.count.matrix<-Matrix(count_mat)


dat <- runWarpLDAModels(dat, topic=c(5:15, 20, 25), seed=123, nCores=12, addModels=FALSE,tmp="/home/groups/CEDAR/mulqueen/temp/")

#select best model
pdf(paste0(args[2],".",args[3],".cistopic_model_selection.pdf"))
par(mfrow=c(3,3))
dat <- selectModel(dat, type='maximum')
dat <- selectModel(dat, type='perplexity')
dat <- selectModel(dat, type='derivative')
dev.off()

system(paste0("slack -F ",args[2],".",args[3],".cistopic_model_selection.pdf"," ryan_todo") )

#going with topics counts based on derivative

dat<-cisTopic::selectModel(dat,type="derivative",keepModels=T)

dat <- runUmap(dat, target='cell') #running umap using cistopics implementation

#add sample cell line names as metadata
dat@cell.data$cellLine<-unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",1))

#set up treatment conditions
dat@cell.data$treatment<-substr(unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",2)),1,1) #set up T47D cell treatment first
dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$cellLine,3,3)
dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$cellLine,5,5)

dat@cell.data$cellLine_treatment<-paste(dat@cell.data$cellLine,dat@cell.data$treatment,sep="_")

dat@cell.data$perc_met<-dat@cell.data$"CG_met_count"/dat@cell.data$"CG_total_count"

#clean up metadata a bit more
dat@cell.data[dat@cell.data$cellLine %in% c("BSM7E6","M7C1A","M7C2B","M7E4C"),]$cellLine<-"MCF7"
dat@cell.data[dat@cell.data$cellLine %in% c("T"),]$cellLine<-"T47D"
dat@cell.data[dat@cell.data$treatment %in% c("C"),]$treatment<-"control"
dat@cell.data[dat@cell.data$treatment %in% c("E"),]$treatment<-"estrogen"

write.table(dat@cell.data,file=paste0(args[2],".",args[3],".cellData.tsv"),col.names=T,row.names=T,sep="\t",quote=F)

pdf(paste0(args[2],".",args[3],".cistopic_clustering.pdf"))
par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine','treatment','cellLine_treatment'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c("CG_total_count","perc_met"), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

system(paste0("slack -F ",args[2],".",args[3],".cistopic_clustering.pdf"," ryan_todo")) 


saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.Rds"))


#should probably make this Z-scored for plotting
topicmat<-as.data.frame(dat@selected.model$document_expects)
write.table(topicmat,file=paste0(args[2],".",args[3],".cellbytopic.cistopic.tsv"),col.names=T,row.names=T,sep="\t",quote=F)

topicmat<-as.data.frame(t(scale(t(topicmat))))
cellline_annot<-dat@cell.data$cellLine
treatment_annot<-dat@cell.data$treatment
c_count_annot<-dat@cell.data$CG_total_count
cellline_treatment_annot<-dat@cell.data$cellLine_treatment

#color function
col_fun=colorRamp2(quantile(unlist(topicmat),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))

row.names(topicmat)<-paste0("Topic_",row.names(topicmat))
plt1<-Heatmap(topicmat,
	name="Z-score",
    show_column_names=F,
    bottom_annotation=columnAnnotation(cellline=cellline_annot,treatment=treatment_annot,C_covered=c_count_annot,cellline_treatment=cellline_treatment_annot,
    	    col = list(cellline = c("T47D"="red","MCF7"="blue"),
               			treatment = c("control" = "white", "estrogen" = "black"),
               			cellline_treatment = c("BSM7E6_E"="green","M7C1A_C"="brown","M7C2B_C"="purple","M7E4C_E"="blue","T_C"="orange","T_E"="red"))))

pdf(paste0(args[2],".",args[3],".cistopic_heatmap.pdf"))
plt1
dev.off()

system(paste0("slack -F ",args[2],".",args[3],".cistopic_heatmap.pdf", " ryan_todo") )


#testing topic significance for cell line/ treatment
dat_stat<-rbind(as.data.frame(dat@selected.model$document_expects),cellline_annot,treatment_annot)
dat_stat<-as.data.frame(t(dat_stat))
dat_stat[1:(ncol(dat_stat)-2)] <- lapply(dat_stat[1: (ncol(dat_stat)-2)], as.numeric)
dat_stat<-melt(dat_stat)
colnames(dat_stat)<-c("cellline","treatment","topic","value")
dat_stat$topic<-paste0("topic_",sub('.', '', dat_stat$topic))
#testing topic specificity by treatment and cell line
out_stats<-as.data.frame(dat_stat %>% 
  group_by(topic) %>% 
  anova_test(value ~ treatment*cellline) %>% adjust_pvalue()
  )

write.table(out_stats,file=paste0(args[2],".",args[3],".cistopic_topic_enrichment.txt"),quote=F,sep="\t",col.names=T)

#region annotations
dat<-readRDS(file=paste0(args[2],".",args[3],".cistopic_object.Rds"))

pred_matrix <- predictiveDistribution(dat) #predictive measurements of topics
dat <- getRegionsScores(dat, method='NormTop', scale=TRUE) #get regions from topics
dat <- annotateRegions(dat, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db') #plot general annotations to topics 
dat<-binarizecisTopics(dat,thrP=0.975) #binarize cistopics to get genes driving topics

#grab top 10 regions per topic
topic_enrich<-melt(as.matrix(dat@region.data[which(grepl(colnames(dat@region.data),pattern="Scores_Topic"))]))
region_plt<-dat@region.data[as.data.frame(topic_enrich %>% group_by(Var2) %>% slice_max(order_by = value, n = 10))[,1],]
row.names(region_plt)<-paste(row.names(region_plt),region_plt$SYMBOL)
region_plt<-region_plt[,which(grepl(pattern="Scores_Topic",colnames(region_plt)))]
col_order<-1:ncol(region_plt)
plt2<-Heatmap(region_plt,name="Topic Scores",column_order=col_order)
pdf(paste0(args[2],".",args[3],".TopicScoreRegions.pdf"),height=20)
plt2
dev.off()
system(paste0("slack -F ",args[2],".",args[3],".TopicScoreRegions.pdf"," ryan_todo"))

#getBedFiles(dat, path=paste0(args[2],".",args[3],"_cisTopic_beds")

dat <- GREAT(dat, genome='hg38', fold_enrichment=1.5, geneHits=1, sign=0.05, request_interval=10)
pdf(paste0(args[2],".",args[3],"cistopic_GO.pdf"),width=30,height=30)
ontologyDotPlot(dat, top=20, topics=c(1:ncol(dat@selected.model$document_expects)), var.y='name', order.by='Binom_Adjp_BH')
dev.off()
system(paste0("slack -F ",paste0(args[2],".",args[3],".GO.pdf"," ryan_todo"))

#Save cistopic object
saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.Rds"))


```

{% endcapture %} {% include details.html %} 



Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R [argv1] [argv2] [argv3]


+ argv1 is a bed file with features to aggregate methylation over. Same as aggregate_methylation_over_region.py input.
+ argv2 is the methylation moeity. This is set up in argv3 of the python script. Example is \*[GpC/CpG].[prefix].count.txt.
+ argv3 is prefix from aggregate_methylation_over_region.py output. This is set up in argv3 of the python script. Example is \*.[prefix].count.txt.


```bash

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

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed \
CpG \
10kb


Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed \
CpG \
bcEnhance

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

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed \
GpC \
10kb

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R \
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed \
CpG \
bcEnhance

```

Combine multiple binarized matrices prior to running cistopic

```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(mefa4)
library(stats)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(patchwork)
library(rstatix)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Ckmeans.1d.dp)

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")

args <- commandArgs(trailingOnly=TRUE)

dat_1<-readRDS(args[1])
dat_2<-readRDS(args[2])

dat<-dat_1
dat@binary.count.matrix<-rbind(dat@binary.count.matrix,dat_2@binary.count.matrix)
dat@region.data<-rbind(dat@binary.count.matrix,dat_2@region.data)
dat <- runWarpLDAModels(dat, topic=c(5:15, 20, 25, 40, 50), seed=123, nCores=14, addModels=FALSE,tmp="/home/groups/CEDAR/mulqueen/temp/")

#select best model
pdf(paste0(args[2],".",args[3],".cistopic_model_selection.pdf"))
par(mfrow=c(3,3))
dat <- selectModel(dat, type='maximum')
dat <- selectModel(dat, type='perplexity')
dat <- selectModel(dat, type='derivative')
dev.off()

system(paste0("slack ",paste0(args[2],".",args[3],".cistopic_model_selection.pdf")) )

#going with topics counts based on derivative

dat<-cisTopic::selectModel(dat,type="derivative",keepModels=T)

dat <- runUmap(dat, target='cell') #running umap using cistopics implementation

#add sample cell line names as metadata
dat@cell.data$cellLine<-unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",1))

#set up treatment conditions
dat@cell.data$treatment<-substr(unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",2)),1,1) #set up T47D cell treatment first
dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$cellLine,3,3)
dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$treatment<-substr(dat@cell.data[startsWith(dat@cell.data$cellLine,"BSM7"),]$cellLine,5,5)

dat@cell.data$cellLine_treatment<-paste(dat@cell.data$cellLine,dat@cell.data$treatment,sep="_")

dat@cell.data$perc_met<-dat@cell.data$"CG_met_count"/dat@cell.data$"CG_total_count"

#clean up metadata a bit more
dat@cell.data[dat@cell.data$cellLine %in% c("BSM7E6","M7C1A","M7C2B","M7E4C"),]$cellLine<-"MCF7"
dat@cell.data[dat@cell.data$cellLine %in% c("T"),]$cellLine<-"T47D"
dat@cell.data[dat@cell.data$treatment %in% c("C"),]$treatment<-"control"
dat@cell.data[dat@cell.data$treatment %in% c("E"),]$treatment<-"estrogen"

pdf(paste0(args[2],".",args[3],".cistopic_clustering.pdf"))
par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine','treatment','cellLine_treatment'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)

par(mfrow=c(1,3))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c("CG_total_count","perc_met"), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)


par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

system(paste0("slack ",paste0(args[2],".",args[3],".cistopic_clustering.pdf")) )


saveRDS(dat,file=paste0(args[2],".",args[3],".cistopic_object.Rds"))
```