---
title: s3GCC Analysis
layout: s3data
author: Ryan Mulqueen
permalink: /s3gcc/
category: s3processing
---

# s3-GCC Processing of PDAC Samples

*Note: This processing is exploratory and not the final processing used for the manuscript*

This notebook describes processing of s3GCC samples following the splitting of a deduplicated bwa mem aligned bam file which occurs in ["s3 Preprocessing"](https://mulqueenr.github.io/s3preprocess/)

First create an annotation file from wet lab data and a working directory.

Annotation data is listed in https://docs.google.com/spreadsheets/d/1mZ34KIwmr2vdjQlnqY7v_u0-Eca8mRc-HgI2r_WICXk/edit#gid=695371319

## Generation of cell summaries

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data")

#complexity data
crc_compl<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/complexity_data/CRC-4442CRC-4671_GCC_H.complexity.txt",header=F)
colnames(crc_compl)<-c("row_no","cellID","total_reads","uniq_reads","perc_uniq")
cellline_compl<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/complexity_data/CellLine_GCC_L.complexity.txt",header=F)
colnames(cellline_compl)<-c("row_no","cellID","total_reads","uniq_reads","perc_uniq")

compl<-rbind(crc_compl,cellline_compl)

#read in indexes master file
idx<-read.table("/home/groups/oroaklab/src/scitools/scitools-dev/SCI_stdchem_Indexes.txt",col.names=c("idx_name","idx_cycles","idx_seq"))
idx_i7<-idx[idx$idx_cycles==1,]
colnames(idx_i7)<-c("i7_idx_name","i7_idx_cycle","i7_idx_seq")
idx_i5<-idx[idx$idx_cycles==2,]
colnames(idx_i5)<-c("i5_idx_name","i5_idx_cycle","i5_idx_seq")
idx_tn5<-idx[idx$idx_cycles==3,]
colnames(idx_tn5)<-c("tn5_idx_name","tn5_idx_cycle","tn5_idx_seq")

#Per Plate breakdown of data contained as sheets in https://docs.google.com/spreadsheets/d/1mZ34KIwmr2vdjQlnqY7v_u0-Eca8mRc-HgI2r_WICXk/edit#gid=1305239528
#Not all plates are a proper barnyard mix.
crc_gcc_pcr_annot<-read.table("S3GCC_crc_pcr_annotation.txt",sep="\t",header=T)
celline_gcc_tn5_annot<-read.table("S3WGSandGCC_cellline_tn5_annotation.txt",sep="\t",header=T)

cellline_cellID<-read.table("s3gcc_cellline.bbrd.q10.filt.cellIDs.list",header=F)
colnames(cellline_cellID)<-c("cellID")
crc_cellID<-read.table("s3gcc_crc.bbrd.q10.filt.cellIDs.list",header=F)
colnames(crc_cellID)<-c("cellID")

cellID<-rbind(crc_cellID,cellline_cellID)
cellID$i7_idx_seq<-substr(cellID$cellID,0,8)
cellID$tn5_idx_seq<-substr(cellID$cellID,17,24)
cellID$i5_idx_seq<-substr(cellID$cellID,9,16)
dat<-merge(cellID,idx_i7,by="i7_idx_seq")
dat<-merge(dat,idx_tn5,by="tn5_idx_seq")
dat<-merge(dat,idx_i5,by="i5_idx_seq")

dat_cellline<-merge(dat,celline_gcc_tn5_annot,by="tn5_idx_seq")
dat_crc<-merge(dat,crc_gcc_pcr_annot,by=c("i5_idx_seq","i7_idx_seq","i5_idx_cycle","i7_idx_cycle","i5_idx_name","i7_idx_name"))

dat_cellline<-dat_cellline[c("cellID","i7_idx_seq","i5_idx_seq","tn5_idx_seq","i7_idx_name","i5_idx_name","tn5_idx_name","sample")]
dat_crc<-dat_crc[c("cellID","i7_idx_seq","i5_idx_seq","tn5_idx_seq","i7_idx_name","i5_idx_name","tn5_idx_name","sample")]

dat<-rbind(dat_cellline,dat_crc)
dat$plate<-unlist(lapply(strsplit(dat$i5_idx_name,"_"),"[",2))

dat<-merge(dat,compl,by="cellID")
write.table(dat,file="s3gcc_cellsummary.tsv",col.names=T,sep="\t")

library(ggplot2)

plt<-ggplot(dat,aes(x=sample,y=log10(uniq_reads),color=as.factor(sample)))+geom_boxplot()+geom_jitter()
ggsave(plt,file="s3gcc_uniqreads.png")
```
{% endcapture %} {% include details.html %} 

## Extraction of distal reads

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
#Extract distal reads based on >=1kb
#or transchromosomal on different chromosomes
#printing out distal connections (>=1kb or different chromosomes) using awk
#Note this is for read projection but they are not used for the actual s3GCC analysis

dir = "/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/raw_alignment"

for i in CellLine_GCC_L.nsrt.bam CRC-4442CRC-4671_GCC_H.nsrt.bam;
do \
out_name=${i::-9};
echo $out_name;
#print out reads on same chromosome that are >= 1000 bp apart (distal)
((samtools view -H $i) & (samtools view $i | awk '{if (sqrt(($9^2))>=1000) print $0}')) | samtools view -bS - > $out_name.distal.bam & 
#print out reads on different chromosomes
((samtools view -H $i) & (samtools view $i | awk '{if($7 != "=") print $0}')) | samtools view -bS - > $out_name.transchr.bam &
#print out reads on same chromosome within 1000bp (cis)
((samtools view -H $i) & (samtools view $i | awk '{if (sqrt(($9^2))<=1000) print $0}')) | samtools view -bS - > $out_name.cis.bam &
done &

#Processing of splitting bam files is similar to how wgs cells were treated in s3wgs notebook.

#moved files to s3gcc_data
mv CRC-4442CRC-4671_GCC_H.cis.bam CRC-4442CRC-4671_GCC_H.distal.bam CRC-4442CRC-4671_GCC_H.transchr.bam ../s3gcc_data/ 

#filter bam files to just cellIDs which pass QC
DIR="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis"
scitools bam-filter -L $DIR/s3gcc_data/s3gcc_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3gcc_data/s3gcc_crc_H.cis $DIR/s3gcc_data/CRC-4442CRC-4671_GCC_H.cis.bam &
scitools bam-filter -L $DIR/s3gcc_data/s3gcc_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3gcc_data/s3gcc_crc_H.distal $DIR/s3gcc_data/CRC-4442CRC-4671_GCC_H.distal.bam &
scitools bam-filter -L $DIR/s3gcc_data/s3gcc_crc.bbrd.q10.filt.cellIDs.list -O $DIR/s3gcc_data/s3gcc_crc_H.trans $DIR/s3gcc_data/CRC-4442CRC-4671_GCC_H.transchr.bam &

#Split bams by cell line
for i in s3gcc_crc_H.cis.filt.bam s3gcc_crc_H.distal.filt.bam s3gcc_crc_H.trans.filt.bam;
do scitools bam-split -A sample.annot $i & done

#ensure bam files are properly mate paired
for i in s3gcc_crc_H.cis.filt.crc4442.bam s3gcc_crc_H.distal.filt.crc4442.bam \
s3gcc_crc_H.distal.filt.crc4671.bam s3gcc_crc_H.trans.filt.crc4442.bam \
s3gcc_crc_H.cis.filt.crc4671.bam s3gcc_crc_H.trans.filt.crc4671.bam ;
do samtools sort -@ 10 -n $i | samtools fixmate - - | samtools view -@ 10 -b -f 2 -F 524 - > ${i::-4}.matefixed.bam ; done &

#perform read projection
for i in  s3gcc_crc_H.cis.filt.crc4442.matefixed.bam s3gcc_crc_H.cis.filt.crc4671.matefixed.bam s3gcc_crc_H.distal.filt.crc4442.matefixed.bam s3gcc_crc_H.distal.filt.crc4671.matefixed.bam s3gcc_crc_H.trans.filt.crc4442.matefixed.bam s3gcc_crc_H.trans.filt.crc4671.matefixed.bam;
do scitools bam-project -r 500 -n 1 -X -e $i ; done &

#Collate projections for plotting
cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data
#read projections
#distal
for i in s3gcc_crc_H.cis.filt.crc4442.matefixed.read_projections s3gcc_crc_H.cis.filt.crc4671.matefixed.read_projections s3gcc_crc_H.distal.filt.crc4442.matefixed.read_projections s3gcc_crc_H.trans.filt.crc4442.matefixed.read_projections s3gcc_crc_H.trans.filt.crc4671.matefixed.read_projections
do line=${i:21:-27};
mate_type=${i:12:-40};
name=$line"_"$mate_type
awk -v name=$name 'OFS="\t" {print $1,$2,$3,name}' $i/cell_summaries.txt;
done >> ./s3gcc_projected_reads.txt 


```
{% endcapture %} {% include details.html %} 

## Projected complexity across cells

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data")
proj<-read.table("s3gcc_projected_reads.txt",header=F)

colnames(proj)<-c("proj_perc_uniq","proj_reads","proj_unique_reads","sample_read_type")


dat<- as.data.frame(proj %>% 
                    group_by(sample_read_type) %>% 
                    summarize(mean=mean(log10(proj_unique_reads/2)),
                              sd=sd(log10(proj_unique_reads/2)),
                              median=median(log10(proj_unique_reads/2))))

dat[dat$proj_perc_uniq==0.05,]
#     sample assay proj_perc_uniq    read_type     mean        sd   median
#1   crc4442 s3gcc           0.05   cis_distal 3.199874 0.4360882 3.131619
#2   crc4442 s3gcc           0.05 cis_proximal 5.443474 0.3253194 5.345948
#3   crc4442 s3gcc           0.05        trans 3.238635 0.2888201 3.170702
#284 crc4671 s3gcc           0.05   cis_distal 3.492634 0.3608255 3.511349
#285 crc4671 s3gcc           0.05 cis_proximal 5.649266 0.3503700 5.614836
#286 crc4671 s3gcc           0.05        trans 3.280166 0.3216718 3.289254

ggplot(dat,aes(x=as.numeric(proj_perc_uniq),fill = paste(sample,assay,read_type),color=paste(sample,assay,read_type)))+
geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),alpha=0.1,color=NA) +
geom_line(aes(y=mean,linetype="solid"))+
geom_line(aes(y=median,linetype="dashed")) +
theme_bw() + scale_x_reverse()+ylim(c(0,6.5))

ggsave(file="gcc_projected_contacts.png")
ggsave(file="gcc_projected_contacts.pdf")

system("slack -F gcc_projected_contacts.png ryan_todo")
system("slack -F gcc_projected_contacts.pdf ryan_todo")
```
{% endcapture %} {% include details.html %} 

## Mate fixing bam files

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
#For splitting to single cell bam files (NOT RUN)

#add read groups to each bam file
for i in s3gcc_crc_H.trans.filt.bam s3gcc_crc_H.distal.filt.bam s3gcc_crc_H.cis.filt.bam; do scitools bam-addrg $i & done
#cis is still running

#cellID is few enough that there are no IO errors, split by RG
for i in s3gcc_crc_H.trans.filt.RG.bam s3gcc_crc_H.distal.filt.RG.bam s3gcc_crc_H.cis.filt.RG.bam; do samtools split -@10 -f './single_cell_projections/%*_%!.%.' $i ; done & #perform samtools split on RG (cellIDs)


#perform parallelized duplicate removal
find . -type f -name '*matefixed.bam' | parallel -j 20 scitools bam-rmdup {} 

#perform parallelized read projection
find . -type f -name '*matefixed.bam' | parallel -j 20 scitools bam-project -r 500 -n 1 -X -e {} 

```
{% endcapture %} {% include details.html %} 

### Collate projections to single files for plotting

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data
#read projections
#distal
for i in ./single_cell_projections/projections/*distal*matefixed.read_projections; 
do cellid=${i:38:-17};
awk -v cellid=$cellid 'OFS="\t" {print $1,$2,$3,cellid,"cis_distal"}' $i/cell_summaries.txt;
done > ./s3gcc_distal_projected_reads.txt 
#cis
for i in ./single_cell_projections/projections/*cis*matefixed.read_projections; 
do cellid=${i:38:-17};
awk -v cellid=$cellid 'OFS="\t" {print $1,$2,$3,cellid,"cis_proximal"}' $i/cell_summaries.txt;
done > ./s3gcc_cis_projected_reads.txt 
#trans
for i in ./single_cell_projections/projections/*trans*matefixed.read_projections; 
do cellid=${i:38:-17};
awk -v cellid=$cellid 'OFS="\t" {print $1,$2,$3,cellid,"trans"}' $i/cell_summaries.txt;
done > ./s3gcc_trans_projected_reads.txt 

```
{% endcapture %} {% include details.html %} 

## Make an analysis directory

{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
import os
wd="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/"
sc_dir="gcc_cisdistal_triplesparse/"

if not os.path.exists(wd+sc_dir):
    os.makedirs(wd+sc_dir)

print("Directory created or already present.")
        
```
{% endcapture %} {% include details.html %} 

### Split bam files into single cell triple-sparse format

Triple-sparse format will be the input for scHiCluster. This script will take in a bam file, split by unique identifiers in the read name, and then move to the new working directory. Output will be a triple-sparse format file for each cell id, for each chromosome. It will also make the network file, which lists full path to each cell ID (needed for scHiCluster).

Currently works at a 1mb resolution.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
import glob
import pysam
from collections import defaultdict
from collections import Counter
import pandas as pd
import os
import time

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

bam_dir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/filtered_bam/" #working directory
bam_bulk=["CRC-4442CRC-4671_GCC_H.bbrd.q10.filt.bam"] #list of bam files
#Takes read 1 and read 2 paired by query name with the above function
#then rounds read 1 and read 2 starts to the nearest 1mb and takes that as the bin (drops the 6 0s for the bin names)
out_dict = defaultdict(dict)
chr_count=[]

i=0 #Counter just to check progress
start_time = time.time() #And set up a timer

#cis interactions first
for bulk in bam_bulk:
    cellline=bulk.split(".")[0]
    bam = pysam.AlignmentFile(bam_dir+bulk, 'rb')
    for read1, read2 in read_pair_generator(bam): #set paired reads together
        if read1 is not None and read2 is not None: #remove reads which dont pair
            if read1.reference_id==read2.reference_id: #ensure chromosomes are same
                if abs(read1.template_length)>=50000: #ignore reads <50kbp in length
                    i+=1
                    if i % 100000 == 0:
                        print("Processed "+str(i)+" distal cis reads in "+str(time.time()-start_time)+" seconds.")
                    if read1.reference_start < 1000000: #bin for 0 (less than 1mbp into chr)
                        r1_bin="0"
                    else:
                        r1_bin=str(read1.reference_start)[:-6] #bin for anything greater than 1mb, just cut last 6 digits off
                    if read2.reference_start < 1000000:
                        r2_bin="0"
                    else:
                        r2_bin=str(read2.reference_start)[:-6]
                    if int(r2_bin) < int(r1_bin): #quick sorting of bins
                      tmp=r2_bin
                      r2_bin=r1_bin
                      r1_bin=tmp
                    chr_out=bam.get_reference_name(read1.reference_id) #set chr name
                    chr_count.append(chr_out) #make sure chr name was seen before or added
                    cellid_out=read1.query_name.split(":")[0] #assign read to cellID
                    cellid_out=cellline+"_"+cellid_out
                    if cellid_out not in out_dict: #add any new cellIDs
                      out_dict[cellid_out][chr_out]=[[r1_bin,r2_bin]] #build out dict
                    if chr_out not in out_dict[cellid_out]:
                      out_dict[cellid_out][chr_out]=[[r1_bin,r2_bin]]
                    else:
                      out_dict[cellid_out][chr_out].append([r1_bin,r2_bin])

chr_count=len(Counter(chr_count))

wd="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/" #working directory
sc_dir="gcc_cisdistal_triplesparse/"

with open(wd+"scHiCluster.network","w") as fout:
  for cellid_out in out_dict: #generate sciHiCluster.network file (list of path/cellname)
      fout.write(wd+sc_dir+cellid_out+"\n")

for cellid_out in out_dict: #generate [cellID].[chr].txt files in triple sparse format
    for chr_out in out_dict[cellid_out]:
      with open(wd+sc_dir+cellid_out+"_"+chr_out+".txt","w") as fout:
        dict = Counter([tuple(i) for i in out_dict[cellid_out][chr_out]])
        # Creating pandas dataframe
        Output = pd.DataFrame(data ={'bin1': list([item[0] for item in list(dict.keys())]),'bin2': list([item[1] for item in list(dict.keys())]),'count': list(dict.values())})
        Output.to_csv(fout,sep="\t",header=False,index=False)
        
```
{% endcapture %} {% include details.html %} 


## Plotting distal reads per cell

{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
import glob
import pysam
from collections import defaultdict
from collections import Counter
import pandas as pd
import os
import time

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

bam_dir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/" #working directory
bam_bulk=["s3gcc_crc.bbrd.q10.filt.bam","s3gcc_cellline.bbrd.q10.filt.bam"] #list of bam files
#Takes read 1 and read 2 paired by query name with the above function
#then rounds read 1 and read 2 starts to the nearest 1mb and takes that as the bin (drops the 6 0s for the bin names)

distal_count={}
trans_count={}
local_count={}

#cis distal and trans interactions counting
for bulk in bam_bulk:
    cellline=bulk.split(".")[0]
    bam = pysam.AlignmentFile(bam_dir+bulk, 'rb')
    for read1, read2 in read_pair_generator(bam): #set paired reads together
        if read1 is not None and read2 is not None: #remove reads which dont pair
            cellid_out=read1.query_name.split(":")[0]
            if read1.reference_id==read2.reference_id: #ensure chromosomes are same
                if abs(read1.template_length)>=1000: #ignore reads <1kbp in length
                    if cellid_out in distal_count:
                        distal_count[cellid_out]+=1
                    else:
                        distal_count[cellid_out]=1
                else:
                    if cellid_out in local_count:
                        local_count[cellid_out]+=1
                    else:
                        local_count[cellid_out]=1
                    
            if read1.reference_id!=read2.reference_id: #ensure chromosomes are same
                    if cellid_out in trans_count:
                        trans_count[cellid_out]+=1
                    else:
                        trans_count[cellid_out]=1

                        
distal_df=pd.DataFrame.from_dict(distal_count,orient="index",columns=["distal_cis_interactions"])
distal_df.index.name = 'cellID'
distal_df.reset_index(inplace=True)

trans_df=pd.DataFrame.from_dict(trans_count,orient="index",columns=["trans_interactions"])
trans_df.index.name = 'cellID'
trans_df.reset_index(inplace=True)

local_df=pd.DataFrame.from_dict(local_count,orient="index",columns=["local_cis"])
local_df.index.name = 'cellID'
local_df.reset_index(inplace=True)

df = pd.merge(left=distal_df, right=trans_df, left_on='cellID', right_on='cellID')
df = pd.merge(left=df, right=local_df, left_on='cellID', right_on='cellID')
df.to_csv(bam_dir+'s3gcc_hiccontacts.csv', index=False)
```
{% endcapture %} {% include details.html %} 

### Run scHiCluster

This script will run scHiCluster, with a harded length of hg38 chromosomes on the files produced from the script above. 

It is currently set to generate 3 clusters (via K Means clustering). sciHiCluster source and example processing is available [here.](https://github.com/zhoujt1994/scHiCluster)

{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
import time
import numpy as np
from schicluster import *
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
import csv
import glob
import pandas as pd

wd="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/" #working directory
sc_dir="gcc_cisdistal_triplesparse/"

network_file=wd+"scHiCluster.network"
network = open(network_file).read().splitlines()
network = np.loadtxt(network_file,dtype=np.str)

#sizes (in bp) of autosomes in numerical order taken from UCSC hg38 
hg38dim=[248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468] #this is all autosomes, but I'm just using chr1-16 for now (others too sparse so it throws errors)
chrom = [str(i+1) for i in range(22)] #only using autosomes 1-22 due to sparsity

chromsize = {chrom[i]:hg38dim[i] for i in range(len(chrom))} #generate dicitonary of chromosome names and sizes

nc=5 #5 clusters?

chr_count = []
chr_names=["chr"+s for s in chrom]
#need to filter cells with less than a certain number of reads per chr
#building a pandas data frame for count per chr
for cellID in network:
    cellID_out=cellID.split("/")[-1]
    for file in glob.glob(cellID+"*"+".txt"):
        chr_out=file.split("/")[-1].split("_")[-1].split(".")[0]
        num_lines = sum(1 for line in open(file))
        chr_count.append([cellID_out,chr_out,num_lines])

chr_count=pd.DataFrame(chr_count,columns=["cellID","chr_out","contact_count"])
filter_list=[]

for i in chr_count.cellID.unique():
    tmp=chr_count.loc[chr_count['cellID'] == i]
    tmp=tmp[tmp["chr_out"].isin(chr_names)]
    if len(tmp)==len(chr_names) and tmp.contact_count.min()>1:
        filter_list.append(i)
        
#filter network to cells with atleast 2 reads per chr
#goal here is to group sparse cells for more thorough analysis

network_fil=["/".join(network[0].split("/")[0:-1])+"/"+i for i in filter_list] #filtering the network list to those that pass the count filter
len(network_fil)
label=["_".join(i.split("/")[-1].split("_")[0:2]) for i in network_fil]

start_time = time.time()
cluster,embedding=hicluster_cpu(network_fil,chromsize,nc=nc, pad=1, rp=0.5, prct=20, ncpus=20)
print(time.time() - start_time)
[ARI(label, KMeans(n_clusters = nc, n_init = 200).fit(embedding[:, :ndim]).labels_) for ndim in [2,5,10,20,30,40,50]]

out_dir=wd

#output annoted clusters
with open(wd+"scHiCluster.annot","w") as fout:
    for i in range(0,len(cluster)-1):
        fout.write(filter_list[i]+"\t"+str(cluster[i])+"\n")
        print(filter_list[i]+"\t"+str(cluster[i]))

        
#output embeddings
embed=pd.DataFrame(embedding,index=filter_list)
with open(wd+"scHiCluster.dims","w") as fout:
    embed.to_csv(fout,sep="\t",header=True,index=True)

```
{% endcapture %} {% include details.html %} 

### Clustering and visualization of cis-distal domain scHiCluster

Next we will load in the output PCA dims file from sciHiCluster and use UMAP on a subset of PCs which define the most variance. For this, I'm just going to use an elbow plot to define the cutoff.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(gridExtra)
library(Rphenograph)

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data")
dat<-read.table("scHiCluster.dims",header=T,row.names=1)
scHiClus<-read.table("scHiCluster.annot",header=F)
colnames(scHiClus)<-c("fullname","scHiClus")

dat_var<-data.frame(pc=seq(1,ncol(dat)),perc_var=(sapply(dat, var)/sum(sapply(dat,var)))*100)
dat_var$cum_perc<-0
dat_var$cum_perc<-unlist(lapply(1:nrow(dat_var),function(x) sum(dat_var[seq(1,x),]$perc_var)))

plt<-ggplot(dat_var,aes(x=pc,y=perc_var))+geom_line()+theme_minimal()+geom_vline(aes(xintercept =7,color="red"))
ggsave(plt,file="scHiCluster_pc_varexplained.pdf")
system("slack -F scHiCluster_pc_varexplained.pdf ryan_todo")

annot<-read.table("scHiCluster.annot",header=F)
colnames(annot)<-c("cellID","kmean_clus")#can also just dbscan or pgclus it myself


annot_samp<-read.table("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/s3wgs_gcc.cellsummary.txt",header=T)

umap_dat<-function(x){
dat_dims<-as.data.frame(uwot::umap(dat[1:x]))
row.names(dat_dims)<-row.names(dat)
dat_dims$fullname<-row.names(dat_dims)
dat_dims$cellID_idx<-unlist(lapply(strsplit(dat_dims$fullname,"_"),"[",4))
dat_dims<-merge(dat_dims,annot_samp,by="cellID_idx")
plt<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=as.factor(sample),size=1))+geom_point()+theme_bw()+ggtitle(paste("UMAP of",x,"PCs"))
return(plt)
}


pc_dat<-function(x){
dat_dims<-as.data.frame(uwot::umap(dat[c(x,x+1)]))
row.names(dat_dims)<-row.names(dat)
dat_dims$fullname<-row.names(dat_dims)
dat_dims$cellID_idx<-unlist(lapply(strsplit(dat_dims$fullname,"_"),"[",4))
dat_dims<-merge(dat_dims,annot_samp,by="cellID_idx")
plt<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=as.factor(sample),size=1))+geom_point()+theme_bw()+ggtitle(paste("PC",x,"and",x+1))
return(plt)
}


plt_list<-lapply(c(3,4,5,6,7,8,9,10,20,30,50,75,100,ncol(dat)),umap_dat)

ncol=4
out_plt<-do.call("grid.arrange", c(plt_list, ncol=ncol))
ggsave(out_plt,file="scHiCluster_umap_pctest.svg",height=20,width=40,limitsize=FALSE)
ggsave(out_plt,file="scHiCluster_umap_pctest.png",height=20,width=40,limitsize=FALSE)
system("slack -F scHiCluster_umap_pctest.png ryan_todo")

plt_list<-lapply(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),pc_dat)

ncol=4
out_plt<-do.call("grid.arrange", c(plt_list, ncol=ncol))
ggsave(out_plt,file="scHiCluster_pca_pctest.svg",height=30,width=40,limitsize=FALSE)
ggsave(out_plt,file="scHiCluster_pca_pctest.png",height=30,width=40,limitsize=FALSE)
system("slack -F scHiCluster_pca_pctest.png ryan_todo")


#selecting x PCs to use
x=10
kmeans_out <- as.data.frame(kmeans(dat[1:x], center = 5)$cluster)
kmeans_out$fullname<-row.names(kmeans_out)
dat_dims<-as.data.frame(uwot::umap(dat[1:x]))
row.names(dat_dims)<-row.names(dat)
dat_dims$fullname<-row.names(dat_dims)
dat_dims$cellID_idx<-unlist(lapply(strsplit(dat_dims$fullname,"_"),"[",4))
dat_dims<-merge(dat_dims,kmeans_out,by="fullname")
colnames(dat_dims)[ncol(dat_dims)]<-"kmeanscluster"
dat_dims<-merge(dat_dims,scHiClus,by="fullname")

dat_dims<-merge(dat_dims,annot_samp,by="cellID_idx")
dat_dims$mapq20<-as.numeric(dat_dims$mapq20)

plt1<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=log10(mapq20)))+geom_point()+theme_bw()+ggtitle("Read Depth")
plt2<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=as.factor(kmeanscluster)))+geom_point()+theme_bw()+ggtitle(paste("PG Cluster on",x,"PCs"))
plt3<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=as.factor(scHiClus)))+geom_point()+theme_bw()+ggtitle(paste("scHiClusters on",x,"PCs"))
plt4<-ggplot(dat=dat_dims,aes(x=V1,y=V2,color=as.factor(sample)))+geom_point()+theme_bw()+ggtitle("Cell Line")

ncol=2
out_plt<-grid.arrange(plt1,plt2,plt3,plt4, ncol=ncol)
ggsave(out_plt,file="scHiCluster_umap.pdf",width=10)
ggsave(out_plt,file="scHiCluster_umap.png",width=10)
system("slack -F scHiCluster_umap.pdf ryan_todo")

write.table(dat_dims,file="s3gcc_fullsummary.txt",sep="\t",quote=F,col.names=T,row.names=F)

for (i in unique(dat_dims$scHiClus)){
    out_dir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse/"
    scHiclus_net_temp<-dat_dims[dat_dims$scHiClus==i,]
    scHiclus_net_temp<-as.data.frame(cbind(scHiclus_net_temp$fullname))
    scHiclus_net_temp$V1<-paste0(out_dir,scHiclus_net_temp$V1)
    write.table(scHiclus_net_temp,file=paste0("network_",i,".txt"),col.names=F,row.names=F,quote=F)
}

```
{% endcapture %} {% include details.html %} 

### Merging cells within clusters to prepare the running of TopDom to find intrachromosomal changes in topological domains

{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
import time
import numpy as np
from schicluster import *
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
import csv
import glob
import pandas as pd

network_file_list=glob.glob("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/for_gurkan/network*txt")

chromsize_file="/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/chromLengths.txt"
hg38dim=[248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468] #this is all autosomes, but I'm just using chr1-16 for now (others too sparse so it throws errors)
chrom = [str(i+1) for i in range(16)]
chromsize = {chrom[i]:hg38dim[i] for i in range(len(chrom))}

chr_count = []
chr_names=["chr"+s for s in chrom]

res=1000000

#this had to be corrected to work for CPU
#also added a loop for single cell topdom
def merge_cpu(network, c, res, pad=1, rp=0.5, prct=-1,chromsize=chromsize):
    ngene = int(chromsize[c] / res) + 1
    start_time = time.time()
    Q_sum = np.zeros(ngene * ngene)
    for cell in network:
        Q_sum = Q_sum + impute_cpu([cell, c, ngene, pad, rp])[1]
    end_time = time.time()
    print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
    Q_sum = Q_sum.reshape(ngene, ngene)
    return Q_sum


#this also had to be fixed with chromsize as variable
def output_topdom(cell, c, Q, res, chromsize=chromsize):
    ngene, _ = Q.shape
    B = [['chr' + c, i * res, (i + 1) * res] for i in range(ngene)]
    B[-1][-1] = chromsize[c]
    C = np.concatenate((B, Q), axis=1)
    np.savetxt(cell + '_chr' + c + '.topdommatrix', C, fmt = '%s', delimiter = '\t')
    return

for network_file in network_file_list:
    network = open(network_file).read().splitlines()
    network = np.loadtxt(network_file,dtype=np.str)
    cell=network_file.split(".")[0]
    print("Processing:"+cell)
    for c in chromsize: #for whole network merged file
        Q = merge_cpu(network, c, res=1000000)
        output_topdom(cell, c, Q, res)
        output_sparse(cell, c, Q, res)

for network_file in network_file_list:
    print(network_file)
    for cell in open(network_file).read().splitlines():
        for c in chromsize:
            ngene = int(chromsize[c] / res) + 1
            start_time = time.time()
            pad=1
            rp=0.5
            Q = np.zeros(ngene * ngene)
            print("Processing"+" "+cell+" "+c)
            Q=impute_cpu([cell, c, ngene,pad,rp])[1]
            end_time = time.time()
            Q = Q.reshape(ngene, ngene)
            output_topdom(cell,c, Q, res)
```
{% endcapture %} {% include details.html %} 

### Using R TopDom to find different topological domains from merged imputed regions.

Note: as of now this script runs through each single cell as well. This isn't particularly useful, but it might be dependent on future analysis.
{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
#R v4.0
#remotes::install_github("HenrikBengtsson/TopDom", ref="master")
library(TopDom)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse")
network_files<-list.files(path="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse",pattern="*topdommatrix")
for (i in network_files){
skip_to_next<-FALSE
print(paste("Processing file:",i))
tad=tryCatch(TopDom(data = i, window.size = 5),error=function(e){skip_to_next<<-TRUE})
if(skip_to_next){
    next
} else{
out_name=strsplit(i,"[.]")[[1]][1]
out_name=paste0(out_name,".w5.domain")
write.table(tad$bed[1:dim(tad$bed)[1],2:4], file=out_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
}
```
{% endcapture %} {% include details.html %} 

### Plot raw contact frequencies of single-cells
{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse")

files_in<-list.files(pattern="_chr16.txt$")
x<-files_in[1]
plot_singlechr<-function(x){
    name_out<-paste0(basename(x),".inter.hic.pdf")
    dat<-read.table(x,head=F,sep="\t")
    dat$V1<-dat$V1+1 #make it start counting at 1 rather than 0
    dat$V2<-dat$V2+1

    for (i in 1:nrow(dat)){
        if (dat[i,]$V2 < dat[i,]$V1) {
            temp<-dat[i,]$V1
            dat[i,]$V1 <- dat[i,]$V2
            dat[i,]$V2 <- temp
            } 
    }
    #create all the country combinations
    df <- expand.grid(1:max(dat$V1,dat$V2), 1:max(dat$V1,dat$V2))
    #change names
    colnames(df) <- c('V1', 'V2')
    #add a value of 0 for the new combinations (won't affect outcome)
    df$V3 <- 0
    #row bind with original dataset
    df <- rbind(df, dat)
    dat_cast<-as.data.frame(xtabs( V3 ~ V1 + V2, aggregate(V3~V1+V2,df,sum)))
    dat_cast<-dcast(data=dat_cast,formula=V1~V2,value.var="Freq",fun.aggregate=sum)
    dat_cast<-dat_cast[2:ncol(dat_cast)]
    #dat_cast<-log10(dat_cast)
    #dat_cast[dat_cast <= -Inf] <- 0    
    plt<-Heatmap(dat_cast,column_order=1:ncol(dat_cast),row_order=1:nrow(dat_cast),col=c("white","red"))
    pdf(name_out)
    print(plt)
    dev.off()
    #if(max(dat_cast)>1.5){
    system(paste0("slack -F ",name_out," ryan_todo"))#}
}

lapply(files_in,plot_singlechr)
```
{% endcapture %} {% include details.html %} 

### Using Python for statistical analysis on TopDom matrices
Input:

domlist: a list of all the topdom output.
cluster: the cluster assignment of each cell with the same order of domlist.
celltypelist: a list of all possible cell type label.
res: the resolution of domains.
chrom: the chromosome.


Returns:

sc_dom: domain boundaries in each single cell indicated by a binary matrix.
dom_prob: domain boundary frequency in each cluster with the same order as in celltypelist.
bins: The tested bins, since the bins with 0.0 or 1.0 domain frequency in any of the clusters were not tested.
pvalue: The p-value of each bin in bins.
{% capture summary %} Code {% endcapture %} {% capture details %}  


```python
import time
import numpy as np
from schicluster import *
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
import csv
import glob
import pandas as pd
import re
from statsmodels.sandbox.stats.multicomp import multipletests as FDR

#using this filter function for subsetting lists
def Filter(string, substr):
    return [str for str in string
    if any(sub in str for sub in substr)]


#diff domain call from scihicluster, pretty extensively modified to work properly

def diff_dom(args):
    domlist, cluster, ctlist, res, c, chromsize = args
    start_time = time.time()
    ctdict = {x:i for i,x in enumerate(ctlist)}
    ngene = int(chromsize[c] / res) + 1
    bound = np.zeros((len(cluster), ngene))
    for i,cell in enumerate(domlist):
    	domain = np.loadtxt(cell, dtype = np.str)
    	C = domain[domain[:,-1]=='domain', :2].astype(int) / res
    	bound[i, np.array(list(set(C.flatten())),dtype=int)] += 1
    	
    cellfilter = (np.sum(bound, axis=1)>0)
    count = np.array([np.sum(np.logical_and([k==x for x in cluster], cellfilter)) for k in ctlist])
    prob = np.array([np.sum(bound[[x==k for x in cluster]], axis=0) for k in ctlist])
    ptmp, btmp = [],[]
    for i in range(ngene):
    	contig = [[prob[j,i], count[j]-prob[j,i]] for j in range(len(ctlist))]
    	if np.sum(np.sum(contig, axis=0)==0)==0:
    		ptmp.append(chi2_contingency(contig)[1])
    		btmp.append('\t'.join([c, str(i*res), str(i*res+res)]))
    		
    print(c, time.time()-start_time)
    return [bound, prob/count[:,None], btmp, ptmp]


network_file_list=glob.glob("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/network*txt")

domlist=glob.glob("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse/*w5.domain")

chromsize_file="/home/groups/oroaklab/adey_lab/projects/sciWGS/Public_Data/chromLengths.txt"

hg38dim=[248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468] #this is all autosomes, but I'm just using chr1-16 for now (others too sparse so it throws errors)
chrom = [str(i+1) for i in range(16)]
chromsize = {chrom[i]:hg38dim[i] for i in range(len(chrom))}

chr_count = []
chr_names=["chr"+s for s in chrom]

res=1000000

for chrom  in chromsize:
    chrom_filt_name="_chr"+chrom+"." #set up chr name
    domlist_chr=[i for i in domlist if chrom_filt_name in i] #first filter topdom list to the chromosome
    domlist_cell_order=[] #initiate empty topdom list for each chr
    cluster_list=[] #initiate empty cluster list for each chr
    for cluster in network_file_list: #loop through each cluster network file
        clus_name=[cluster.split("/")[-1].split(".")[0]] #set cluster name
        cluster_list=cluster_list+(clus_name)*len(open(cluster).read().splitlines()) #make a list of cluster assignment per cell by multiplying cluster name by number of cells
        domlist_cell_order=domlist_cell_order+Filter(domlist_chr,open(cluster).read().splitlines()) #make list of topdoms with appropriate cluster assignment and chr
    celltype_list=list(set(cluster_list)) #set celltype_list by unique celltypes
    sc_dom, dom_prob, bins, pvalue=diff_dom([domlist_cell_order,cluster_list,celltype_list,res,chrom,chromsize])
    bin_names=["chr"+chrom+":"+str(i*res)+"-"+str((i+1)*res) for i in range(np.shape(sc_dom)[1])]
    cellid_names=[cell.split("/")[-1].split("_chr")[0] for cell in domlist_cell_order]
    sc_dom_df=pd.DataFrame(sc_dom,index=cellid_names,columns=bin_names)
    dom_prob_df=pd.DataFrame(dom_prob,index=celltype_list,columns=bin_names)
    pval_df=pd.DataFrame(pvalue,index=[col.split("\t")[0]+":"+col.split("\t")[1]+"-"+col.split("\t")[2] for col in bins],columns=["pvalue"])
    sc_dom_df.to_csv("chr"+chrom+"_scDomainBoundaries.csv")
    dom_prob_df.to_csv("chr"+chrom+"_scDomainProbabilities.csv")
    pval_df.to_csv("chr"+chrom+"_scDomainProbabilities_pval.csv")

```
{% endcapture %} {% include details.html %} 

### Using Python to merge raw matrix output

This is so we maintain transchromosomal interactions.
{% capture summary %} Code {% endcapture %} {% capture details %}  

```python
#merge raw matrices for whole chr comparisons

import glob
import pysam
import pandas as pd
from collections import defaultdict
import time
import sys, os

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

    
with open("/home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/scHiCluster.annot","r") as annot:
    annot_df=pd.DataFrame([line.split() for line in annot],columns=["cellID","cluster"])
    annot_df["cluster"]=[x[1] for x in annot_df["cellID"].str.split(pat="_")] #using cell lines



annot_df["cellname"]=[line.split("_")[2] for line in annot_df["cellID"]]

bam_dir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/" #working directory
bam_bulk=["hg38.s3GCC_4671.bbrd.q10.filt.bam","hg38.s3GCC_4442.bbrd.q10.filt.bam"] #list of bam files
#Takes read 1 and read 2 paired by query name with the above function

dir_out="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/for_gurkan/"

i=0 # start counter
start_time = time.time() # start timer
chr_set=["chr"+str(i) for i in range(1,23)]

blockPrint()
for clus in annot_df["cluster"].unique():
    cell_set=list(annot_df.loc[annot_df["cluster"]==str(clus)]["cellID"])
    with open(dir_out+"s3GCC_clus"+clus+".bulk","w") as fout:
        for bulk in bam_bulk:
            cellline=bulk.split(".")[1]
            bam = pysam.AlignmentFile(bam_dir+bulk, 'rb')
            for read1, read2 in read_pair_generator(bam): #set paired reads together
                if read1 is not None and read2 is not None: #remove reads which dont pair
                    chr1=bam.get_reference_name(read1.reference_id) #set chr1 name
                    chr2=bam.get_reference_name(read2.reference_id) #set chr2 name
                    if (chr1 in chr_set) and (chr2 in chr_set):
                        if ((chr1 == chr2) and (abs(read1.template_length)>=1000)) or (chr1 is not chr2): #so same chr >1kb or diff chr
                            i+=1
                            if i % 10000 == 0:
                                enablePrint()
                                print("Processed "+str(i)+" distal reads in "+ str(time.time()-start_time) +" seconds.")
                                blockPrint()
                            str1="+"
                            str2="+"
                            if read1.is_reverse:
                                str1="-"
                            if read2.is_reverse:
                                str2="-"
                            pos1=read1.reference_start
                            pos2=read2.reference_start
                            cellid_out=cellline+"_"+read1.query_name.split(":")[0] #assign read to cellID
                            if cellid_out in cell_set: #check to make sure it was a cellID used in clustering
                                out_list=[read1.query_name,chr1,pos1,str1,chr2,pos2,str2]
                                fout.write("\t".join(str(item) for item in out_list)+"\n")


#1. Read Name (can be blank)
#2. chromosome for read 1
#3. positions for read 1 (5' end of read, one-indexed)
#4. strand of read 1 (+ or -)
#5. chromosome for read 2
#6. positions for read 2 (5' end of read, one-indexed)
#7. strand of read 2 (+ or -)

import subprocess as sp
### Using awk for final formatting of text into HOMER input.
cmd= """for i in `ls /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/for_gurkan/*bulk`; do awk 'OFS="\t" {if ($2 > $5){ print $1,$5,$6,$7,$2,$3,$4} else {print}}' $i | awk 'OFS="\t" {if (($3 > $6) && ($2 == $5)) {print $1,$5,$6,$7,$2,$3,$4} else {print}}' | sort -T . -S 2G -k2,2d -k5,5d -k3,3n -k6,6n > $i.sorted.hicSummary ; done"""

p = sp.Popen(cmd, stdin=sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE,shell=True)
p.wait()
### Using Homer for generation of the raw formatted HiC matrix
#Resolution of 2.5MB and uwing 20CPUs
cmd="""for i in clus4442 clus4671;
do /home/groups/oroaklab/src/homer/bin/makeTagDirectory s3GCC_${i} -format HiCsummary s3GCC_${i}.bulk.sorted.hicSummary;
/home/groups/oroaklab/src/homer/bin/analyzeHiC s3GCC_${i} -res 2500000 -raw -cpu 20 > $i.hic.matrix; done"""

p = sp.Popen(cmd, stdin=sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE,shell=True)
p.wait()

```
{% endcapture %} {% include details.html %} 


### Finally using R for Knight Ruiz normalization of whole genome HiC data, and visualization

I am both looping through individual clusters and performing a cluster by cluster difference in contacts.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
library(ComplexHeatmap)
library(HiCcompare)
library(circlize)
library(EnsDb.Hsapiens.v86)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/gcc_cisdistal_triplesparse")

read_dat<-function(x){
    dat<-read.table(x)
    bincount<-seqlengths(EnsDb.Hsapiens.v86)["3"] %/% 1000000    
    bins<-paste0("bin_",0:bincount)
    dat$V1<-paste0("bin_",dat$V1)
    dat$V2<-paste0("bin_",dat$V2)
    dat_out<-data.frame(matrix(0,ncol=length(bins),nrow=length(bins)))
    colnames(dat_out)<-bins
    row.names(dat_out)<-bins
    for (row in 1:nrow(dat)) {
    dat_out[dat[row,1],dat[row,2]] <- dat[row,3]}
    dat_out[which(is.na(dat_out),arr.ind=T)]<-0
    name_out<-substr(x,1,nchar(x)-4)
    plt_dat(dat_out,name_out)
    return(dat_out)
}

plt_dat<-function(x,file_name){
    dat<-x
    
    col_fun<-viridis(10)
    col_fun<-colorRamp2(c(0,1,3,20),c("#450256","#5AC865","#21908D","#F9E721"))
    pdf(paste0(file_name,".inter.hic.png"))
    plt<-Heatmap(as.matrix(dat),
    row_order=1:nrow(dat),
    column_order=1:ncol(dat),
    show_row_names=FALSE,
    show_column_names=FALSE,
    col = col_fun)
    print(plt)
    dev.off()
    system(paste0("slack -F ",file_name,".inter.hic.png", " ryan_todo"))
}

for (i in c(as.character(1:22),"X")){
    i="3"
files_in<-list.files(pattern=paste0("*chr",i,".txt"))
for (x in files_in){
dat<-read_dat(x)
}}

#plot each cluster 
lapply(c("clus4442.hic.matrix","clus4671.hic.matrix"),FUN=plt_dat)
clus4442<-norm_dat("clus4442.hic.matrix")
clus4671<-norm_dat("clus4671.hic.matrix")

plt_chr<-function(x=clus4442,y=clus4671,chr_in){
    x<-x[grepl(paste0(chr_in,"-"),colnames(x)),grepl(paste0(chr_in,"-"),rownames(x))]
    col_fun_4442 = colorRamp2(c(0,as.numeric(quantile(unlist(x), probs = c(0.90))),as.numeric(quantile(unlist(x), probs = c(0.99))))
    , c("#fff5eb", "#f16913","#7f2704"))
    plt_4442<-Heatmap(as.matrix(x),
    cluster_rows=F,
    cluster_columns=F,
    show_row_names=FALSE,
    show_column_names=FALSE,
    col=col_fun_4442,
    column_title=paste("CRC 4442",chr_in),
    row_order=1:nrow(x),
    column_order=1:ncol(x))
    
    y<-y[grepl(paste0(chr_in,"-"),colnames(y)),grepl(paste0(chr_in,"-"),rownames(y))]
    col_fun_4671 = colorRamp2(c(0,as.numeric(quantile(unlist(y), probs = c(0.90))),as.numeric(quantile(unlist(y), probs = c(0.99))))
    , c("#fff5eb", "#f16913","#7f2704"))
    plt_4671<-Heatmap(as.matrix(y),
    cluster_rows=F,
    cluster_columns=F,
    show_row_names=FALSE,
    show_column_names=FALSE,
    col=col_fun_4671,
    column_title=paste("CRC 4671",chr_in),
    row_order=1:nrow(y),
    column_order=1:ncol(y))
    
    pdf(paste0("cellline_",chr_in,".intra.hic.png"),width=10)
    print(plt_4442+plt_4671)
    dev.off()
    system(paste0("slack -F ",paste0("cellline_",chr_in,".intra.hic.png"), " ryan_todo"))
}

for (i in paste0("chr",1:20)){plt_chr(chr_in=i,x=clus4442,y=clus4671)}



match_data_chr<-function(i,dat_mat){
    out<-unlist(lapply(strsplit(colnames(dat_mat),"-"),"[",1)) %in% strsplit(row.names(dat_mat)[i],"-")[[1]][1]
    return(out)}
    
sel_chr<-function(i,j,dat_mat){
  out_row<-unlist(lapply(strsplit(colnames(dat_mat),"-"),"[",1)) %in% i
  out_col<-unlist(lapply(strsplit(colnames(dat_mat),"-"),"[",1)) %in% j
  out_matched<-out_row & out_col
  return(out_matched)
}

diff_mat<-function(x,y){
if (x !=y){
print(x)
print(y)
dat_ref<-read.table(x,skip=1)
dat_ref<-dat_ref[2:ncol(dat_ref)]
row.names(dat_ref)<-dat_ref[,1]
dat_ref<-dat_ref[2:ncol(dat_ref)]
colnames(dat_ref)<-row.names(dat_ref)
dat_ref<-as.matrix(dat_ref)
#dat_ref[do.call("rbind",lapply(1:nrow(dat_ref),FUN=match_data_chr,dat_mat=dat_ref))]<-0
#dat_ref<-dat_ref[sel_chr("chr15","chr15",dat_ref),sel_chr("chr15","chr15",dat_ref)] #This can be used for selecting individual chromosomes for comparison
dat_ref<-KRnorm(dat_ref)

dat_test<-read.table(y,skip=1)
dat_test<-dat_test[2:ncol(dat_test)]
row.names(dat_test)<-dat_test[,1]
dat_test<-dat_test[2:ncol(dat_test)]
colnames(dat_test)<-row.names(dat_test)
dat_test<-as.matrix(dat_test)
#dat_test[do.call("rbind",lapply(1:nrow(dat_test),FUN=match_data_chr,dat_mat=dat_test))]<-0
#dat_test<-dat_test[sel_chr("chr15","chr15",dat_test),sel_chr("chr15","chr15",dat_test)]

dat_test<-KRnorm(dat_test)

dat_ref<-dat_ref[(colnames(dat_ref) %in% colnames(dat_ref)) & (colnames(dat_ref) %in% colnames(dat_test)),(row.names(dat_ref) %in% row.names(dat_ref)) & (row.names(dat_ref) %in% row.names(dat_test))]
dat_test<-dat_test[(colnames(dat_test) %in% colnames(dat_ref)) & (colnames(dat_test) %in% colnames(dat_test)),(row.names(dat_test) %in% row.names(dat_ref)) & (row.names(dat_test) %in% row.names(dat_test))]

dat_diff<-dat_test-dat_ref

x<-strsplit(x,"[.]")[[1]][1]
y<-strsplit(y,"[.]")[[1]][1]

chr_split<-unlist(lapply(strsplit(row.names(dat_diff),"-"),"[",1))
chr_split<-as.numeric(unlist(lapply(strsplit(chr_split,"r"),"[",2)))

col_func = colorRamp2(as.numeric(quantile(unlist(dat_diff), probs = c(0.001,0.1,0.5,0.9,0.999))), c("#d7191c", "#fdae61", "white", "#abdda4", "#2b83ba"))
dat_diff[which(dat_diff %in% c(NaN,-Inf,Inf),arr.ind=T)]<-NaN
pdf(paste0(x,"by",y,".difftest.heatmap.pdf"))
plt<-Heatmap(as.matrix(dat_diff),
na_col="white",
cluster_rows=F,
cluster_columns=F,
show_row_names=FALSE,
show_column_names=FALSE,
row_split=chr_split,
column_split=chr_split,
col=col_func)
print(plt)
dev.off()}
}

for (y in c("clus0.hic.matrix","clus1.hic.matrix","clus2.hic.matrix")){
lapply(c("clus0.hic.matrix","clus1.hic.matrix","clus2.hic.matrix"),FUN=diff_mat,y=y)
}

#I've been having a strange bug with KRnorm, where within this function it stalls and doesn't complete, but if clusters are run individually it is fine.

```
{% endcapture %} {% include details.html %} 
