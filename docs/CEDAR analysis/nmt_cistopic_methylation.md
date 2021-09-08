---
title: CisTopic on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /cistopic_nmt/
category: CEDAR
---

## Set up genomic annotations

Generate bed files for annotation regions to sum methylation data over.

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

Also given a list of enhancer regions that are breast cancer specific from Hisham
````bash
#Raw file stored here:
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed

#Converted to same bed format with metadata bin names
grep -v "_" /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}' > /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed

````

ENCODE ChIP Peaks
```bash
cd /home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata
#FOXA1 
#MCF7 https://www.encodeproject.org/files/ENCFF112JVK/
wget https://www.encodeproject.org/files/ENCFF112JVK/@@download/ENCFF112JVK.bed.gz
#T47D https://www.encodeproject.org/files/ENCFF758GJL/
wget https://www.encodeproject.org/files/ENCFF758GJL/@@download/ENCFF758GJL.bed.gz

zcat ENCFF112JVK.bed.gz | awk 'OFS="\t" {print $1,$2,$3,"MCF7_"NR}' > FOXA1.bed
zcat ENCFF758GJL.bed.gz | awk 'OFS="\t" {print $1,$2,$3,"T47D_"NR}' >> FOXA1.bed
```
List of genomic annotation bed files:
```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/ctcf.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed
/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed
```

## Preparing bam files from Aaron Doe's preprocessing
Need to merge resequenced bam files and combine read 1 and read 2 into a single-end bam

```bash
mkdir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams
cd /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams
ls /home/groups/CEDAR/doe/projects/my_NMT/MCF7_reseq/merge_all_MCF7/scNMT_NOMeWorkFlow/bismarkSE/dedup/*merged.bam > files.txt
ls /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scNMT_NOMeWorkFlow/bismarkSE/*bam | grep -v "M7" >> files.txt
for i in `cat files.txt`; do outname=`basename $i`; ln -s $i $outname; done
for i in *_R1.*bt2.bam; do 
r1=$i;
r2=`echo $i | sed "s/R1/R2/g"`;
r2=`echo $r2 | sed "s/val_1_/val_2_/g"`;
outname=`echo $i | awk 'OFS="_" {split($1,a,"_");print a[1],a[2],a[3]}'`;
samtools cat -o ${outname}_merged.bam $r1 $r2 ;
done &
```

## Batch script for deduplication
Using the bismark deduplication script.

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-374
#SBATCH --tasks-per-node=30 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams/*merged.bam | wc -l`

bam_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/AD_bams/*merged.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

srun deduplicate_bismark \
--single \
--output_dir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/dedup_bams \
${bam_in}

```

## Set up methylation extraction

Perform methylation extraction with bismark methylation extractor. Output cov.gz files.

Here I'm using 10 nodes, 20 cpus per node (since parallelization is 3x5), with 20gb per node. It's a lot of resources, but hopefully it will finish fast.

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

Generate cytosine NOME data for region summarization.

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


## Using python pandas and bedtools wraper to summarize methylation over regions

* argv1 is a cov.gz output from coverage2cytosine report (either CG or GpC)
* argv2 is a bed file with features to aggregate methylation over
* argv3 is output directory
* argv4 is a prefix used for the bed file annotation. Example: "genebody" or "promoter"

Note: argv2 should be a bed file of format [chr]\t[start]\t[end]\t[feature_name]


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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
bcEnhance

```

{% endcapture %} {% include details.html %} 


{% capture summary %} CpG FOXA1 {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
foxa1

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
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
bcEnhance

```

{% endcapture %} {% include details.html %} 


{% capture summary %} GpC FOXA1 {% endcapture %} {% capture details %}  

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
/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed \
/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions \
foxa1

```
{% endcapture %} {% include details.html %} 


### Running all of these as batch jobs

```bash
sbatch CpG_genes.slurm.sh 
sbatch CpG_enhancers.slurm.sh  
sbatch CpG_promoters.slurm.sh  
sbatch CpG_100kb.slurm.sh  
sbatch CpG_10kb.slurm.sh  
sbatch CpG_bcEnhancer.slurm.sh
sbatch CpG_FOXA1.slurm.sh

sbatch GpC_genes.slurm.sh
sbatch GpC_enhancers.slurm.sh
sbatch GpC_promoters.slurm.sh
sbatch GpC_100kb.slurm.sh  
sbatch GpC_100kb.slurm.sh  
sbatch GpC_bcEnhancer.slurm.sh
sbatch GpC_FOXA1.slurm.sh

```

## Running cistopic on methylated regions.

Making an R script for cistopic processing of methylation files.

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

# read in arguments
args <- commandArgs(trailingOnly=TRUE)

#test set of arguments
#args<-list()
#args[1]<-"/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed" 
#args[2]<-"GpC"
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
  min.cells = 0,
  min.regions = 0,
  is.acc = 0.5, 
  return.norm.beta=T,
  ...
) {
  # Prepare annotation
  #tester 
  #regions<-region
  #methfiles<-meth_files
  #return.norm.beta=T
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
	  	#tester
                #x<-methfiles[1]
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


dat_out<-createcisTopicObjectFromMeth(methfiles=meth_files,regions=region) #run the read in function
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
system(paste0("slack -F ",args[2],".",args[3],".feature_stats.pdf"," ryan_todo"))

print("Binarize Data by automatic cutoff.")
backup<-dat@binary.count.matrix
count_mat<-as.matrix(dat@count.matrix)

# #Apply coverage filter, require at 5 or more CG measurements
 dat_cov<-as.matrix(dat_cov)
 dat_cov<-dat_cov[row.names(dat_cov) %in% row.names(count_mat),colnames(dat_cov) %in% colnames(count_mat)]
#removed coverage filter for now
 #count_mat[which(dat_cov<5,arr.ind=T)]<-0


count_mat[which(count_mat>=cutoff_value,arr.ind=T)]<-1
count_mat[which(count_mat<cutoff_value,arr.ind=T)]<-0


#remove cells with no methylated regions (due to low coverage)
filt_keep<-which(colMeans(count_mat)>0)
dat@binary.count.matrix<-Matrix(count_mat[,filt_keep])
dat@cell.data<-dat@cell.data[filt_keep,]
dat@cell.names<-dat@cell.names[filt_keep]

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
dat <- runUmap(dat, target='region') #running umap using cistopics implementation

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

dat <- runUmap(dat, target='region') #running umap using cistopics implementation

pdf(paste0(args[2],".",args[3],".cistopic_clustering.pdf"))
par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='region',topic_contr='Probability')
dev.off()
system(paste0("slack -F ",args[2],".",args[3],".cistopic_clustering.regions.pdf"," ryan_todo")) 

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


## Run cistopic on each region and C context

Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R [argv1] [argv2] [argv3]


+ argv1 is a bed file with features to aggregate methylation over. Same as aggregate_methylation_over_region.py input.
+ argv2 is the methylation moeity. This is set up in argv3 of the python script. Example is \*[GpC/CpG].[prefix].count.txt.
+ argv3 is prefix from aggregate_methylation_over_region.py output. This is set up in argv3 of the python script. Example is \*.[prefix].count.txt.

Generate sbatch scripts for running:
```bash

reg=("/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed" "/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed" "/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed" "/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed")
c_type=("CpG" "CpG" "CpG" "CpG" "CpG" "GpC" "GpC" "GpC" "GpC" "GpC" "CpG" "GpC")
name=( "gene" "promoter" "enhancer" "100kb" "bcEnhance" "gene" "promoter" "enhancer" "100kb" "bcEnhance" "foxa1" "foxa1")

for i in "${!reg[@]}"; do
printf '#!/bin/bash 
#SBATCH --nodes=1 #request 1 node 
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time 
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs) 
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu 
#SBATCH --time=3:00:00 ## ask for 1 hour on the node 
#SBATCH -- 
srun Rscript /home/groups/CEDAR/mulqueen/src/cistopic_methylation.R %s %s %s ' "${reg[i]}" "${c_type[i]}" "${name[i]}" > "cistopic.${c_type[i]}.${name[i]}.slurm.sh"  ; done

```

### Running all of these as batch jobs

```bash 
sbatch cistopic.CpG.promoter.slurm.sh 
sbatch cistopic.CpG.100kb.slurm.sh
sbatch cistopic.CpG.enhancer.slurm.sh 
sbatch cistopic.CpG.bcEnhance.slurm.sh
sbatch cistopic.CpG.gene.slurm.sh

sbatch cistopic.GpC.bcEnhance.slurm.sh 
sbatch cistopic.GpC.gene.slurm.sh
sbatch cistopic.GpC.100kb.slurm.sh
sbatch cistopic.GpC.enhancer.slurm.sh
sbatch cistopic.GpC.promoter.slurm.sh
sbatch cistopic.CpG.foxa1.slurm.sh
sbatch cistopic.GpC.foxa1.slurm.sh

```


# RNA Processing

## Check RNA data from Aaron Doe Processing

```bash
mkdir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna
multiqc -o /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples
```

I think his filtering is way to aggressive given the quality of the data. So I'm going to realign and generate my own counts matrix.

## Raw data is located here: 
```bash
/home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples/raw
```

### Make a symbolic link to RNA file directories

```bash
ln -s /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples/raw /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw
```

### Align raw files with batch scripting
Using the genome reference and version of STAR packages with cellranger. Trimming first 15 bp and last 3 bp as per fastqc output.
Trimming first 15 bases and last 3 bases, only aligning read 1
Perform counting (I think this accounts for duplicate reads, but I'm not entirely sure)
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-672
#SBATCH --tasks-per-node=30 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/*_R1.fq.gz | wc -l`

fq_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/*_R1.fq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
name="$(basename -- $fq_in)";
name=${name%.fq.gz} 
echo $name;
/home/groups/CEDAR/mulqueen/src/cellranger/cellranger-6.0.1/lib/bin/STAR --runMode alignReads \
--genomeDir /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star \
--readFilesIn /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/${name}.fq.gz \
--outFileNamePrefix /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam/${name}. \
--readFilesCommand zcat \
--clip3pNbases 3 \
--clip5pNbases 15 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts; 
```

## Generate a Seurat Object for our RNA

Read in ReadsPerGene.out.tab files and format into a data frame. Then make a seurat object for filtering and processing cells.

Following this: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

```R
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(cisTopic)
library(biomaRt) #convert ENSG gene name to chromosome location

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam")
ff <- list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 ) #first 4 lines are mapping QC, might be good to add as a metadata table
counts_in <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "_R1[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/", "", ff )
colnames(counts_in) <- ff
row.names(counts_in) <- counts.files[[1]]$V1

#make seurat object
dat <- CreateSeuratObject(counts = counts_in, project = "cellline", min.cells = 3, min.features = 1000)

#add mito reads to feature data
mito.features= genes(EnsDb.Hsapiens.v86, filter = ~ seq_name == "MT")
#dat[["percent.mt"]] <- PercentageFeatureSet(dat, assay="RNA",features=mito.features$gene_id)

#qc plot
pdf("rna_qc.pdf")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
dev.off()
system("slack -F rna_qc.pdf ryan_todo")

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("rna_topvar.pdf")
plot1 / plot2
dev.off()
system("slack -F rna_topvar.pdf ryan_todo")

#dim reduction
all.genes <- rownames(dat)
dat <- ScaleData(dat, features = all.genes)
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

pdf("rna_pcaloadings.pdf")
VizDimLoadings(dat, dims = 1:2, reduction = "pca")
dev.off()
system("slack -F rna_pcaloadings.pdf ryan_todo")

pdf("rna_pca.pdf")
DimPlot(dat, reduction = "pca")
dev.off()
system("slack -F rna_pca.pdf ryan_todo")

pdf("rna_pcaheatmap.pdf")
DimHeatmap(dat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
system("slack -F rna_pcaheatmap.pdf ryan_todo")

#determine dimensionality
dat <- JackStraw(dat, num.replicate = 100)
dat <- ScoreJackStraw(dat, dims = 1:20)

pdf("rna_elbowplot.pdf")
JackStrawPlot(dat, dims = 1:15)
ElbowPlot(dat)
dev.off()
system("slack -F rna_elbowplot.pdf ryan_todo")

#Cluster and UMAP
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.5)
dat <- RunUMAP(dat, dims = 1:10)

#add more cell line info to data
dat$cellLine<-"MCF7"
dat@meta.data[which(dat$orig.ident %in% c("TDC1","TDE8")),]$cellLine<-"T47D"
dat$treatment <- "E"
dat@meta.data[which(dat$orig.ident %in% c("M7C1A","M7C2B","TDC1")),]$treatment<-"C"
dat$cellLine_treatment<-paste(dat$cellLine,dat$treatment)

pdf("rna_umap.pdf")
DimPlot(dat, reduction = "umap",group.by=c("seurat_clusters","orig.ident","treatment","cellLine","cellLine_treatment"))
dev.off()
system("slack -F rna_umap.pdf ryan_todo")

Idents(object = dat) <- dat$cellLine_treatment
dat.markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
var_features<-as.data.frame(dat.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC))


biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
t2g<-getBM(attributes=c('external_gene_name','ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)
markers_out <- merge(dat.markers, t2g, by.x="gene", by.y= 'ensembl_gene_id',all.x=T)

pdf("rna_markers.pdf")
plt<-DoHeatmap(dat, features = var_features$gene, size = 4,angle = 90) + NoLegend()
dev.off()
system("slack -F rna_markers.pdf ryan_todo")
markers_out[which(markers_out$gene %in% var_features$gene),]$external_gene_name
# [1] "BCAS1"    "TFRC"     "TFRC"     "EXOSC5"   "EDN1"     "EDN1"
# [7] "LXN"      "CEACAM6"  "RAB18"    "RAB18"    "GLA"      "DKK1"
#[13] "PRLR"     "PRLR"     "EFEMP1"   "EFEMP1"   "OLFML3"   "CNN3"
#[19] "CNN3"     "CCN2"     "GNG11"    "GNG11"    "PKIB"     "GABBR2"
#[25] "SLC40A1"  "SLC40A1"  "CYP1A1"   "ABHD2"    "SLC39A6"  "PRSS23"
#[31] "GFRA1"    "RCAN1"    "RCAN1"    "TFF1"     "CLDN1"    "CLDN1"
#[37] "DEGS2"    "DEGS2"    "RAB31"    "SERPINA6" "ZNF217"   "NUPR1"
#[43] "NUPR1"    "TMEM64"   "ANAPC7"   "KRT81"    "IGHM"     "NEAT1"
#[49] "MIR23AHG"

write.table(dat.markers,sep="\t",col.names=T,quote=F,file="rna_markers.tsv")
saveRDS(dat,"/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam/SeuratObject.rds")
#Clustering looks pretty good to me just based on tutorial defaults. Using this seurat object for MOFA analysis, going to recluster using cistopic for RNA

```

## Run WarpLDA on RNA data

```R
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(cisTopic)
library(biomaRt) #convert ENSG gene name to chromosome location

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")
dat<-readRDS("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam/SeuratObject.rds")

outname="RNA.gene"

cistopic_counts_frmt<-dat@assays$RNA@counts
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)

gene_loc <- merge(data.frame("ensembl_gene_id"=row.names(cistopic_counts_frmt)), t2g, by= 'ensembl_gene_id',all.x=T)

row.names(cistopic_counts_frmt)<-paste0("chr",gene_loc$chromosome_name,":",gene_loc$start_position,"-",gene_loc$end_position)
cistopic_counts_frmt<-cistopic_counts_frmt[!endsWith(row.names(cistopic_counts_frmt),"NA"),]

sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(20:40),nCores=20,addModels=FALSE)

pdf(paste(outname,"model_selection.pdf",sep="_"))
par(mfrow=c(3,3))
sub_cistopic_models<- cisTopic::selectModel(sub_cistopic_models, type='derivative')
dev.off()
system(paste0("slack -F ",paste(outname,"model_selection.pdf",sep="_")," ryan_todo"))

dat<-cisTopic::selectModel(sub_cistopic_models,type="derivative",keepModels=T)

dat <- runUmap(dat, target='cell') #running umap using cistopics implementation
dat <- getRegionsScores(dat)
dat <- binarizecisTopics(dat)
dat <- runUmap(dat, target='region') #running umap using cistopics implementation

#add sample cell line names as metadata
dat@cell.data$cellLine<-unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",1))

#set up treatment conditions
dat@cell.data$treatment<-substr(unlist(lapply(strsplit(row.names(dat@cell.data),"_"),"[",1)),3,3)
dat@cell.data$batch<-paste(dat@cell.data$cellLine,dat@cell.data$treatment,sep="_")
dat@cell.data[startsWith(dat@cell.data$cellLine,"M7"),]$cellLine<-"MCF7"
dat@cell.data[startsWith(dat@cell.data$cellLine,"TD"),]$cellLine<-"T47D"
dat@cell.data$cellLine_treatment<-paste(dat@cell.data$cellLine,dat@cell.data$treatment)

saveRDS(dat,file=paste(outname,"cistopic_object.Rds",sep="."))

write.table(dat@cell.data,file=paste0(outname,".cellData.tsv"),col.names=T,row.names=T,sep="\t",quote=F)

pdf(paste0(outname,".cistopic_clustering.pdf"))
par(mfrow=c(1,2))
plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine','treatment','cellLine_treatment'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()
system(paste0("slack -F ",outname,".cistopic_clustering.pdf"," ryan_todo")) 

pdf(paste0(outname,".cistopic_clustering.regions.pdf"))
par(mfrow=c(2,5))
plotFeatures(dat, method='Umap', target='region',topic_contr='Probability')
dev.off()
system(paste0("slack -F ",outname,".cistopic_clustering.regions.pdf"," ryan_todo")) 

```

## Combine multiple binarized matrices prior to running cistopic

Rscript /home/groups/CEDAR/mulqueen/src/merge_cistopic_methylation.R [argv1] [argv2] [argvN]

This script will merge binarized counts matrices prior to re-running cistopic and plotting.

+ argv1 through N are cistopic object .Rds files generated from cistopic_methylation.R 
This script is written to handle 2 or more cistopic objects.

```R
library(cisTopic)
library(reshape2)
library(GenomicRanges)
library(Matrix)
library(stats)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(patchwork)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)
library(Rphenograph)

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")

args <- commandArgs(trailingOnly=TRUE)

#args<-list()
#args[1]<-"CpG.bcEnhance.cistopic_object.Rds"
#args[2]<-"GpC.promoter.cistopic_object.Rds"
#args[3]<-"CpG.promoter.cistopic_object.Rds"
#args[4]<-"GpC.bcEnhance.cistopic_object.Rds"
#args[5]<-"RNA.gene.cistopic_object.Rds"
#args<-unlist(args)


dat<-lapply(args,readRDS) #read in data as list of cistopic objects

#annotate data regions for interpretability
dat<-lapply(dat,function(x){
        annotateRegions(x, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')}
        )


c_type<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][1])) #extract c type
regions<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][2])) #extract region info

#rename cell names for consistency
dat<-lapply(dat,function(x){
        colnames(x@binary.count.matrix)<-unlist(lapply(strsplit(colnames(x@binary.count.matrix),"[.]"),"[",1))
        colnames(x@count.matrix)<-unlist(lapply(strsplit(colnames(x@count.matrix),"[.]"),"[",1))
        row.names(x@cell.data)<-unlist(lapply(strsplit(row.names(x@cell.data),"[.]"),"[",1))
        return(x)
        })

#rename cell names for consistency
dat<-lapply(dat,function(x){
        colnames(x@count.matrix)<-unlist(lapply(strsplit(colnames(x@count.matrix),"[.]"),"[",1))
        colnames(x@count.matrix)<-gsub( "_merged", "", colnames(x@count.matrix) )
        colnames(x@count.matrix)<-gsub( "^T_", "TD", colnames(x@count.matrix) )
        colnames(x@count.matrix)<-gsub("BSM7E6","M7E6A",colnames(x@count.matrix))
        colnames(x@binary.count.matrix)<-unlist(lapply(strsplit(colnames(x@binary.count.matrix),"[.]"),"[",1))
        colnames(x@binary.count.matrix)<-gsub( "_merged", "", colnames(x@binary.count.matrix) )
        colnames(x@binary.count.matrix)<-gsub( "^T_", "TD", colnames(x@binary.count.matrix) )
        colnames(x@binary.count.matrix)<-gsub("BSM7E6","M7E6A",colnames(x@binary.count.matrix))
        row.names(x@cell.data)<-unlist(lapply(strsplit(row.names(x@cell.data),"[.]"),"[",1))
        row.names(x@cell.data)<-gsub( "_merged", "", row.names(x@cell.data) )
        row.names(x@cell.data)<-gsub( "^T_", "TD", row.names(x@cell.data) )
        row.names(x@cell.data)<-gsub("BSM7E6","M7E6A",row.names(x@cell.data))
        colnames(x@count.matrix)<-gsub("_S[0-9]*_L[0-9]*","",colnames(x@count.matrix))
        colnames(x@binary.count.matrix)<-gsub("_S[0-9]*_L[0-9]*","",colnames(x@binary.count.matrix))
        row.names(x@cell.data)<-gsub("_S[0-9]*_L[0-9]*","",row.names(x@cell.data))
        return(x)
        })

#i think BSM7E6 was changed to M7E6A (which is reflected in the function above)

#determine shared cells across cistopic objects
cells_to_keep<-Reduce(intersect,lapply(dat,function(x){row.names(x@cell.data)}))

dat_merged<-dat[[1]] #use first element as cistopic object for formatting

#set up merged object
dat_merged@binary.count.matrix<-do.call(rbind,
        lapply(1:length(dat),function(x){
                row.names(dat[[x]]@binary.count.matrix)<-paste(c_type[x],regions[x],dat[[x]]@region.data$SYMBOL,row.names(dat[[x]]@binary.count.matrix),sep="_")
                tmp<-dat[[x]]@binary.count.matrix[,cells_to_keep]
                return(tmp)}))

dat_merged@count.matrix<-do.call(rbind,
        lapply(1:length(dat),function(x){
                row.names(dat[[x]]@count.matrix)<-paste(c_type[x],regions[x],dat[[x]]@region.data$SYMBOL,row.names(dat[[x]]@count.matrix),sep="_")
                tmp<-dat[[x]]@count.matrix[,cells_to_keep]
                return(tmp)}))

dat_merged@region.data<-do.call(rbind,lapply(1:length(dat),function(x){
                tmp<-dat[[x]]@region.data
                row.names(tmp)<-paste(c_type[x],regions[x],dat[[x]]@region.data$SYMBOL,row.names(tmp),sep=".")
                tmp<-tmp[c("seqnames","start","end","width","nCounts","nCells","annotation","geneChr","geneEnd","geneLength","geneStrand","geneId","transcriptId","distanceToTSS","ENSEMBL","SYMBOL","GENENAME")]
                return(tmp)}))


dat_merged@cell.names<-colnames(dat_merged@binary.count.matrix)
dat_merged@cell.data<-dat_merged@cell.data[which(cells_to_keep %in% row.names(dat[[1]]@cell.data)),]
#run warp lda and overwrite dat object
dat <- runWarpLDAModels(dat_merged, topic=c(5:15, 20, 25, 40, 50), seed=123, nCores=14, addModels=FALSE,tmp="/home/groups/CEDAR/mulqueen/temp/")

#set up output name by
out_name<-paste(unlist(lapply(1:length(c_type), function(i) paste(c_type[i],regions[i],sep="_"))),collapse=".")

#select best model
pdf(paste(out_name,"cistopic_model_selection.pdf",sep="."))
par(mfrow=c(3,3))
dat <- selectModel(dat, type='maximum')
dat <- selectModel(dat, type='perplexity')
dat <- selectModel(dat, type='derivative')
dev.off()

system(paste0("slack -F ",paste(out_name,"cistopic_model_selection.pdf",sep=".")," ryan_todo") )

#going with topics counts based on derivative

dat<-cisTopic::selectModel(dat,type="derivative",keepModels=T)
dat <- runUmap(dat, target='cell') #running umap using cistopics implementation

#run phenogram based clustering on the dim reduction
Rphenograph_out <- Rphenograph(t(dat@selected.model$document_expects), k = 100)
modularity(Rphenograph_out[[2]])
membership(Rphenograph_out[[2]])
dat@cell.data$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))

#add sample cell line names as metadata
dat@cell.data$cellLine<-"MCF7"
dat@cell.data[startsWith(dat@cell.names,"T"),]$cellLine<-"T47D"

dat@cell.data$batch<-unlist(lapply(strsplit(dat@cell.names,"_"),"[",1))

#treatment
dat@cell.data$treatment<-"control"
dat@cell.data[which(dat@cell.data$batch %in% c("TDE8","M7E6A","M7E4C")),]$treatment<-"estrogen"
dat@cell.data$cellLine_treatment<-paste(dat@cell.data$cellLine,dat@cell.data$treatment)


pdf(paste(out_name,"combined.cistopic.clustering.pdf",sep="."))
        par(mfrow=c(2,2))
        plotFeatures(dat, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine_treatment','treatment','cellLine','phenograph_cluster'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
        par(mfrow=c(2,5))
        plotFeatures(dat, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()
system(paste0("slack -F ",paste(out_name,"combined.cistopic.clustering.pdf",sep="."), " ryan_todo") )


saveRDS(dat,file=paste(out_name,"cistopic_object.Rds",sep="."))


topicmat<-as.data.frame(dat@selected.model$document_expects)
write.table(topicmat,file=paste0(out_name,".cellbytopic.cistopic.tsv"),col.names=T,row.names=T,sep="\t",quote=F)

topicmat<-as.data.frame(t(scale(t(topicmat))))
cellline_annot<-dat@cell.data$cellLine
treatment_annot<-dat@cell.data$treatment
c_count_annot<-dat@cell.data$CG_total_count
cellline_treatment_annot<-dat@cell.data$batch
cluster_annot<-dat@cell.data$phenograph_cluster

#color function
col_fun=colorRamp2(quantile(unlist(topicmat),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))

#Plot heatmap topic
row.names(topicmat)<-paste0("Topic_",row.names(topicmat))
plt1<-Heatmap(topicmat,
        name="Z-score",
    show_column_names=F,
    bottom_annotation=columnAnnotation(cellline=cellline_annot,treatment=treatment_annot,C_covered=c_count_annot,batch=cellline_treatment_annot, cluster=cluster_annot,
            col = list(cellline = c("T47D"="red","MCF7"="blue"),
                                treatment = c("control" = "white", "estrogen" = "black"))))

pdf(paste0(out_name,".cistopic_heatmap.pdf"))
plt1
dev.off()

system(paste0("slack -F ",out_name,".cistopic_heatmap.pdf", " ryan_todo") )

#set up region scores to get important features per topic
dat <- getRegionsScores(dat, method='NormTop', scale=TRUE)
dat <- binarizecisTopics(dat, thrP=0.975, plot=TRUE)
dat <- runUMAP(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)

pdf(paste0(out_name,".cistopic_region_umap.pdf"))

############TO ADD###############
# add what kind of region it is to the region based UMAP for coloring

par(mfrow=c(1,2))
plotFeatures(cisTopicObject, method='UMAP', target='region', topic_contr=NULL, colorBy=c('nCounts', 'nCells'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='UMAP', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()
system(paste0("slack -F ",out_name,".cistopic_region_umap.pdf", " ryan_todo") )


saveRDS(dat,file=paste(out_name,"cistopic_object.Rds",sep="."))


```

### Running all of these as batch jobs

```bash

reg=("CpG.gene.cistopic_object.Rds CpG.promoter.cistopic_object.Rds CpG.enhancer.cistopic_object.Rds" "GpC.gene.cistopic_object.Rds GpC.promoter.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds" "CpG.bcEnhance.cistopic_object.Rds GpC.bcEnhance.cistopic_object.Rds" "CpG.bcEnhance.cistopic_object.Rds GpC.promoter.cistopic_object.Rds CpG.promoter.cistopic_object.Rds GpC.bcEnhance.cistopic_object.Rds" "CpG.gene.cistopic_object.Rds GpC.gene.cistopic_object.Rds" "CpG.enhancer.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds" "CpG.promoter.cistopic_object.Rds GpC.promoter.cistopic_object.Rds" "CpG.enhancer.cistopic_object.Rds GpC.promoter.cistopic_object.Rds CpG.promoter.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds""CpG.gene.cistopic_object.Rds CpG.promoter.cistopic_object.Rds CpG.enhancer.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "GpC.gene.cistopic_object.Rds GpC.promoter.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.bcEnhance.cistopic_object.Rds GpC.bcEnhance.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.bcEnhance.cistopic_object.Rds GpC.promoter.cistopic_object.Rds CpG.promoter.cistopic_object.Rds GpC.bcEnhance.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.gene.cistopic_object.Rds GpC.gene.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.enhancer.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.promoter.cistopic_object.Rds GpC.promoter.cistopic_object.Rds RNA.gene.cistopic_object.Rds" "CpG.enhancer.cistopic_object.Rds GpC.promoter.cistopic_object.Rds CpG.promoter.cistopic_object.Rds GpC.enhancer.cistopic_object.Rds RNA.gene.cistopic_object.Rds") 
outname=("CpG.gene.CpG.promoter.CpG.enhancer" "GpC.gene.GpC.promoter.GpC.enhancer" "CpG.bcEnhance.GpC.bcEnhance" "CpG.bcEnhance.GpC.promoter.CpG.promoter.GpC.bcEnhance" "CpG.gene.GpC.gene" "CpG.enhancer.GpC.enhancer" "CpG.promoter.GpC.promoter" "CpG.enhancer.GpC.promoter.CpG.promoter.GpC.enhancer" "CpG.gene.CpG.promoter.CpG.enhancer.RNA.gene" "GpC.gene.GpC.promoter.GpC.enhancer.RNA.gene" "CpG.bcEnhance.GpC.bcEnhance.RNA.gene" "CpG.bcEnhance.GpC.promoter.CpG.promoter.GpC.bcEnhance.RNA.gene" "CpG.gene.GpC.gene.RNA.gene" "CpG.enhancer.GpC.enhancer.RNA.gene" "CpG.promoter.GpC.promoter.RNA.gene" "CpG.enhancer.GpC.promoter.CpG.promoter.GpC.enhancer.RNA.gene")

for i in "${!reg[@]}"; do
printf '#!/bin/bash 
#SBATCH --nodes=1 #request 1 node 
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time 
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs) 
#SBATCH --mem-per-cpu=12gb ## request gigabyte per cpu 
#SBATCH --time=5:00:00 ## ask for 1 hour on the node 
#SBATCH -- 
srun Rscript /home/groups/CEDAR/mulqueen/src/merge_cistopic_methylation.R %s ' "${reg[i]}" > "cistopic.${outname[i]}.slurm.sh"  ; done

```
```bash
sbatch cistopic.CpG.gene.CpG.promoter.CpG.enhancer.slurm.sh
sbatch cistopic.GpC.gene.GpC.promoter.GpC.enhancer.slurm.sh
sbatch cistopic.CpG.bcEnhance.GpC.bcEnhance.slurm.sh
sbatch cistopic.CpG.bcEnhance.GpC.promoter.CpG.promoter.GpC.bcEnhance.slurm.sh
sbatch cistopic.CpG.gene.GpC.gene.slurm.sh
sbatch cistopic.CpG.enhancer.GpC.enhancer.slurm.sh
sbatch cistopic.CpG.promoter.GpC.promoter.slurm.sh

#rna included
sbatch cistopic.CpG.enhancer.GpC.enhancer.RNA.gene.slurm.sh
sbatch cistopic.CpG.promoter.GpC.promoter.RNA.gene.slurm.sh
sbatch cistopic.GpC.gene.GpC.promoter.GpC.enhancer.RNA.gene.slurm.sh
sbatch cistopic.CpG.bcEnhance.GpC.bcEnhance.RNA.gene.slurm.sh
sbatch cistopic.CpG.bcEnhance.GpC.promoter.CpG.promoter.GpC.bcEnhance.RNA.gene.slurm.sh
sbatch cistopic.CpG.gene.GpC.gene.RNA.gene.slurm.sh
sbatch cistopic.CpG.enhancer.GpC.enhancer.RNA.gene.slurm.sh
sbatch cistopic.CpG.promoter.GpC.promoter.RNA.gene.slurm.sh
sbatch cistopic.CpG.enhancer.GpC.promoter.CpG.promoter.GpC.enhancer.RNA.gene.slurm.sh


```
