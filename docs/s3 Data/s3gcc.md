---
title: s3GCC Analysis
layout: s3data
author: Ryan Mulqueen
permalink: /s3gcc/
category: s3processing
---

# s3-GCC Processing of PDAC Samples

This notebook describes processing of s3GCC samples following the splitting of a deduplicated bwa mem aligned bam file which occurs in s3wgs processing ipynb.

First create an annotation file from wet lab data and a working directory.

Annotation data is listed in https://docs.google.com/spreadsheets/d/1mZ34KIwmr2vdjQlnqY7v_u0-Eca8mRc-HgI2r_WICXk/edit#gid=695371319


```python
R
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


```python
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


```python
R
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


```python
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

### Collate projections to single files for plotting


```python
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


```python
import os
wd="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data/"
sc_dir="gcc_cisdistal_triplesparse/"

if not os.path.exists(wd+sc_dir):
    os.makedirs(wd+sc_dir)

print("Directory created or already present.")
        
```

    Directory created or already present.


### Split bam files into single cell triple-sparse format

Triple-sparse format will be the input for scHiCluster. This script will take in a bam file, split by unique identifiers in the read name, and then move to the new working directory. Output will be a triple-sparse format file for each cell id, for each chromosome. It will also make the network file, which lists full path to each cell ID (needed for scHiCluster).

Currently works at a 1mb resolution.


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

    Processed 100000 distal cis reads in 32.036776065826416 seconds.



    ---------------------------------------------------------------------------

    KeyboardInterrupt                         Traceback (most recent call last)

    <ipython-input-2-aa6bebf048b7> in <module>
         41     cellline=bulk.split(".")[1]
         42     bam = pysam.AlignmentFile(bam_dir+bulk, 'rb')
    ---> 43     for read1, read2 in read_pair_generator(bam): #set paired reads together
         44         if read1 is not None and read2 is not None: #remove reads which dont pair
         45             if read1.reference_id==read2.reference_id: #ensure chromosomes are same


    <ipython-input-2-aa6bebf048b7> in read_pair_generator(bam)
         14     read_dict = defaultdict(lambda: [None, None])
         15     for read in bam.fetch(until_eof=True):
    ---> 16         qname = read.query_name
         17         if qname not in read_dict:
         18             if read.is_read1:


    KeyboardInterrupt: 


## Plotting distal reads per cell


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

                    
R

setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3gcc_data")
dat<-read.table("s3gcc_cellsummary.tsv",header=T)
interactions<-read.table("s3gcc_hiccontacts.csv",header=T,sep=",")
dat<-merge(dat,interactions,by="cellID")
dat$interaction_perc<-(dat$distal_cis_interactions+dat$trans_interactions)/(dat$distal_cis_interactions+dat$trans_interactions+dat$local_cis)*100

library(dplyr)
dat %>% group_by(sample) %>% summarize(mean=mean(interaction_perc))

#  sample    mean
#  <chr>    <dbl>
#1 crc_4442 11.7
#2 crc_4671 12.2
#3 gm12878   4.70
#4 hela      5.58
#5 k562      5.23

```

### Run scHiCluster

This script will run scHiCluster, with a harded length of hg38 chromosomes on the files produced from the script above. 

It is currently set to generate 3 clusters (via K Means clustering). 


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

### Clustering and visualization of cis-distal domain scHiCluster

Next we will load in the output PCA dims file from sciHiCluster and use UMAP on a subset of PCs which define the most variance. For this, I'm just going to use an elbow plot to define the cutoff.




```python
#%load_ext rpy2.ipython
#With a working rpy2 install, this would allow for python and R scripts to run in the same notebook. However, for now I'll just leave the R scrips in place to be run via copy and paste into a terminal.
```


```python
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


      File "<ipython-input-11-1d88329b234b>", line 10
        dat_var$cum_perc<-0
               ^
    SyntaxError: invalid syntax



### Merging cells within clusters to prepare the running of TopDom to find intrachromosomal changes in topological domains




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

    Processing:/home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_2
    Load and impute chromosome 1 take 1.9365289211273193 seconds
    Load and impute chromosome 2 take 2.078774929046631 seconds
    Load and impute chromosome 3 take 1.6358287334442139 seconds
    Load and impute chromosome 4 take 1.382009744644165 seconds
    Load and impute chromosome 5 take 1.681645393371582 seconds
    Load and impute chromosome 6 take 1.6281013488769531 seconds
    Load and impute chromosome 7 take 1.1158297061920166 seconds
    Load and impute chromosome 8 take 0.9876255989074707 seconds
    Load and impute chromosome 9 take 0.8976645469665527 seconds
    Load and impute chromosome 10 take 0.7673242092132568 seconds
    Load and impute chromosome 11 take 1.057448387145996 seconds
    Load and impute chromosome 12 take 0.9621901512145996 seconds
    Load and impute chromosome 13 take 0.7142288684844971 seconds
    Load and impute chromosome 14 take 0.6337606906890869 seconds
    Load and impute chromosome 15 take 0.7955880165100098 seconds
    Load and impute chromosome 16 take 0.27492737770080566 seconds
    Processing:/home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_0
    Load and impute chromosome 1 take 0.7354927062988281 seconds
    Load and impute chromosome 2 take 0.7730979919433594 seconds
    Load and impute chromosome 3 take 0.5525975227355957 seconds
    Load and impute chromosome 4 take 0.31310319900512695 seconds
    Load and impute chromosome 5 take 0.36420202255249023 seconds
    Load and impute chromosome 6 take 0.36326146125793457 seconds
    Load and impute chromosome 7 take 0.27658605575561523 seconds
    Load and impute chromosome 8 take 0.3545567989349365 seconds
    Load and impute chromosome 9 take 0.2565169334411621 seconds
    Load and impute chromosome 10 take 0.30072760581970215 seconds
    Load and impute chromosome 11 take 0.33304333686828613 seconds
    Load and impute chromosome 12 take 0.3495612144470215 seconds
    Load and impute chromosome 13 take 0.3053710460662842 seconds
    Load and impute chromosome 14 take 0.27278971672058105 seconds
    Load and impute chromosome 15 take 0.16901493072509766 seconds
    Load and impute chromosome 16 take 0.11724734306335449 seconds
    Processing:/home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_1
    Load and impute chromosome 1 take 1.4433975219726562 seconds
    Load and impute chromosome 2 take 1.4078147411346436 seconds
    Load and impute chromosome 3 take 1.1155366897583008 seconds
    Load and impute chromosome 4 take 1.0236051082611084 seconds
    Load and impute chromosome 5 take 1.0109343528747559 seconds
    Load and impute chromosome 6 take 0.9756267070770264 seconds
    Load and impute chromosome 7 take 0.859424352645874 seconds
    Load and impute chromosome 8 take 0.7219979763031006 seconds
    Load and impute chromosome 9 take 0.8222923278808594 seconds
    Load and impute chromosome 10 take 0.9217216968536377 seconds
    Load and impute chromosome 11 take 0.7935421466827393 seconds
    Load and impute chromosome 12 take 0.6880009174346924 seconds
    Load and impute chromosome 13 take 0.5640699863433838 seconds
    Load and impute chromosome 14 take 0.5225276947021484 seconds
    Load and impute chromosome 15 take 0.4962944984436035 seconds
    Load and impute chromosome 16 take 0.3136453628540039 seconds
    /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_2.txt
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAAGTGGT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCCGCCGATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTACGACA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTATCAGCGCAGCTTGTCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGAGCTCGCT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATATGGAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATGAGGCC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGCCGTGAAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGGCTTGTCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCATATGGAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCCGGAACTG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGTGTCGGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCTAAGGTCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAATTGTGAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAGTAGGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACTCACCAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAATCTGTTGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAATATGGAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACAGTAGGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGCCACAGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCACCTTGGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGCCTCAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGTTCAGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTCCTGTT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCGTAGTG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATTGGACTC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAAGCTCGCT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAAGCTCGCT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACCTTCACC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGACTTGGTAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATACTCATA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATCGATATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCAAGCTAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGTGTCGGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTCCTGTT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAATATGGAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGCCGATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATACTCATA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGACGAAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATGTAGATA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATTGCCTAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGAACCGCG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCGCCACAGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGACCTGAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGTTCCAAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAATGGCATG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGACAAC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGAGGCTTAAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATGGATCGA 16
    /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_0.txt
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTAATACAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAAAGATCAGTCCGCTAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCAAGTCCAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTAACGATAGCTCCAACGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTAAATGCCTC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACATAGAGT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTACTCACCAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACCTGCGTATAATACAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCGCCGATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCCTCTCGTC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCGGACTTGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GAAGAGTACGGTCCGCTCGTAGTG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4442_GGTTAGTTTCAGCGCAGCACGGAC 16
    /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/network_1.txt
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGATTGTGAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAAAGATCAGTCTGTTGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAACGATAGCGATCTATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAAATCCGGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACAATTAAC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAACGGACAAC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAATGAGTAAGAACCGCG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAACACTAAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTACGCCGATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTAGTGTCGGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACCTGCGTATAAGGTCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCAGGTTATA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCCGCCGATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGACCTGAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGATCTATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCGGTACCTT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGACGGTCCGCTTGCCTAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAAAGTCCAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAAGCCACAGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGAGCTGATAATCTGTTGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATCAGCGCAGCGCAAGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGAGGAGCGTC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_ATTGAGGATGACGCGATAATACAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGCTCACCAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGATCTATC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCAATGCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCCACAGG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGGCGCAAGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTGACGAAT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTAAGATCAGTTGCCTAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAAATCCGGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACCAAGTCT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAACGGCGTGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTATGAGTAAGCAATGCA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTACGGACAAC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGCGCAAGC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTAGGTACCTT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCCTGCGTATCGTAGTG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCATTGTGAA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCCCTTCACC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTGGATCGA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGCCTAG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTCGGTCCGCTTGGACTC 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTGCTGATAAGGCATTCT 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGACGGAACTG 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATACTCATA 16
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 1
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 2
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 3
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 4
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 5
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 6
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 7
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 8
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 9
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 10
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 11
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 12
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 13
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 14
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 15
    Processing /home/groups/oroaklab/adey_lab/projects/sciWGS/200522_s3WGS_CRC/gcc_cisdistal_triplesparse/s3GCC_4671_GGTTAGTTTGACGCGATCGATATC 16


### Using R TopDom to find different topological domains from merged imputed regions.

Note: as of now this script runs through each single cell as well. This isn't particularly useful, but it might be dependent on future analysis.


```python
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


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-15-79038bc3e602> in <module>
          8 import pandas as pd
          9 import re
    ---> 10 from statsmodels.sandbox.stats.multicomp import multipletests as FDR
         11 
         12 #using this filter function for subsetting lists


    ModuleNotFoundError: No module named 'statsmodels'


### Using Python to merge raw matrix output

This is so we maintain transchromosomal interactions.


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




    0



### Finally using R for Knight Ruiz normalization of whole genome HiC data, and visualization

I am both looping through individual clusters and performing a cluster by cluster difference in contacts.


```python
#R
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
