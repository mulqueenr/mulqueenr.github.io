---
title: s3 Pre-processing
layout: s3data
author: Ryan Mulqueen
permalink: /s3preprocess/
category: s3processing
---

# Detailed run breakdown and experiment meta data is in the following google doc.

https://docs.google.com/spreadsheets/d/1mZ34KIwmr2vdjQlnqY7v_u0-Eca8mRc-HgI2r_WICXk/edit#gid=1567443037

```bash
rsync of NovaSeq files from the nix node:
    
rsync -rv mulqueen@nix:/data/LIB200619AA ./ & 
#yDVZpuYk6mh7XWZMzD4m58
rsync -rv mulqueen@nix:/data/LIB200716AA ./ &
#yDVZpuYk6mh7XWZMzD4m58


#Got raw bcl files and ran bcl2fastq
module load bcl2fastq/2.19.0
bcl2fastq -R /home/groups/oroaklab/adey_lab/projects/sciWGS/200701_NovaSeqRAW/LIB200619AA/200625_A01058_0053_BHN7HJDRXX-raw -o /home/groups/oroaklab/fastq/200630_NovaSeqSP_RM_s3 --create-fastq-for-index-reads --no-lane-splitting    
bcl2fastq --no-lane-splitting --create-fastq-for-index-reads -R /home/groups/oroaklab/adey_lab/projects/sciWGS/200729_NovaSeq/LIB200716AA/200728_A01058_0066_AHTTTTDMXX

#Files located in /home/groups/oroaklab/fastq/200630_NovaSeqSP_RM_s3

```

## FastQ Dump of reads into sci format


```bash

#Processing all nextseq runs
for i in 191118_NS500556_0362_AHVYV7AFXY 191223_NS500556_0376_AHW7GNAFXY 191203_NS500556_0369_AHVYTYAFXY 200110_NS500556_0379_AHYNJVBGXC 200130_NS500556_0382_AHW7FNAFXY;
do \
scitools fastq-dump-sci-stdchem -R $i; done &

#modified fastq-dump-sci-stdchem
#Added a -V flag for reverse complimenting index read 2 for shared processing the NovaSeq files

fq_dir="/home/groups/oroaklab/fastq/200630_NovaSeqSP_RM_s3"
fq_rundir="200630_NovaSeqSP_RM_s3"
scitools fastq-dump-sci-stdchem -V -R $fq_rundir & #output files in /home/groups/oroaklab/demultiplex/200630_NovaSeqSP_RM_s3

fq_dir="/home/groups/oroaklab/fastq/200728_NovaSeq_S2_RM_s3"
fq_rundir="200728_NovaSeq_S2_RM_s3"
scitools fastq-dump-sci-stdchem -V -R $fq_rundir & #output files in /home/groups/oroaklab/demultiplex/200728_NovaSeq_S2_RM_s3

#Setting up the output folder
mkdir /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis

```

# Do the following plate splitting for all sequencing runs


```bash
#generating annotation file based on CPT i5 PCR index
demux="/home/groups/oroaklab/demultiplex"
#Tab separated experiment annotations.

nano $demux/191118_NS500556_0362_AHVYV7AFXY/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#15	WGS	GM12878	Plate 1	D	Set 1

nano $demux/191223_NS500556_0376_AHW7GNAFXY/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#5	WGS	CRC-4442	Plate 1	E	Set 1
#6	WGS	CRC-4442	Plate 2	F	Set 1
#7	WGS	CRC-4671	Plate 3	G	Set 1

nano $demux/191203_NS500556_0369_AHVYTYAFXY/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#1	ATAC	Barnyard	Plate 1	B	Set 1
#3	ATAC	Barnyard	Plate 3	C	Set 1

nano $demux/200110_NS500556_0379_AHYNJVBGXC/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#4	ATAC	Barnyard	Plate 5	A	Set 1
#5	WGS	CRC-4442	Plate 1	E	Set 1
#6	WGS	CRC-4442	Plate 2	F	Set 1
#7	WGS	CRC-4671	Plate 3	G	Set 1
#8	GCC	CRC-4442CRC-4671	Plate 4	H	Set 1

nano $demux/200130_NS500556_0382_AHW7FNAFXY/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#3	ATAC	Barnyard	Plate 2	D	Set 1


nano $demux/200630_NovaSeqSP_RM_s3/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#2	ATAC	Barnyard	Plate 4	E	Set 1
#3	ATAC	Barnyard	Plate 5	A	Set 1
#7	GCC	CRC-4442CRC-4671	Plate 4	H	Set 1
#8	WGS	CellLine	Plate 1	J	Set 1
#9	WGS	CellLine	Plate 2	K	Set 1
#12	GCC	CellLine	Plate 3	L	Set 1

nano $demux/200728_NovaSeq_S2_RM_s3/plate_split.annot.simplified.txt
#Pool	Assay	Sample	Per Experiment Plate Name	i5 Index (Nextera)	i7 Index (truseq)
#15	WGS	GM12878	Plate 1	D	Set 1
#16	WGS	GM12878	Plate 2	B	Set 1
#5	WGS	CRC-4442	Plate 1	E	Set 1
#6	WGS	CRC-4442	Plate 2	F	Set 1
#7	WGS	CRC-4671	Plate 3	G	Set 1
#8	GCC	CRC-4442CRC-4671	Plate 4	H	Set 1
#3	ATAC	Barnyard	Plate 3	C	Set 1
#9	WGS	CellLine	Plate 1	J	Set 1
#10	WGS	CellLine	Plate 2	K	Set 1
#13	GCC	CellLine	Plate 3	L	Set 1

#python script for generating annotation files

import pandas as pd
idx=pd.read_csv("/home/groups/oroaklab/src/scitools/scitools-dev/SCI_stdchem_Indexes.txt",sep="\t",names=["idx_name","idx_pos","idx_seq"])

demux="/home/groups/oroaklab/demultiplex"
run_list=["191118_NS500556_0362_AHVYV7AFXY", "191223_NS500556_0376_AHW7GNAFXY", "191203_NS500556_0369_AHVYTYAFXY", "200110_NS500556_0379_AHYNJVBGXC", "200130_NS500556_0382_AHW7FNAFXY", "200630_NovaSeqSP_RM_s3","200728_NovaSeq_S2_RM_s3"]

def plate_split_annot_generator(x):
    plate_split=pd.read_csv(demux+"/"+x+"/plate_split.annot.simplified.txt",sep="\t",header=0)
    with open(demux+"/"+x+"/plate.annot","w") as fout:
        for i in range(0,len(plate_split.index)):
            plate_sample=plate_split["Sample"][i]
            plate_assay=plate_split["Assay"][i]
            plate_idx=plate_split["i5 Index (Nextera)"][i]
            for j in idx.loc[idx['idx_pos'] == 1]["idx_seq"]:
                for k in idx.loc[(idx['idx_pos'] == 2) & [x[1]==plate_idx for x in idx["idx_name"].str.split("_")]]["idx_seq"]:
                    for l in idx.loc[idx['idx_pos'] == 3 ]["idx_seq"]:
                        fout.write(j+k+l+"\t"+plate_sample+"_"+plate_assay+"_"+plate_idx+"\n")

                        
for x in run_list:
    plate_split_annot_generator(x)
```


```bash
#Split out fastq files:

demux="/home/groups/oroaklab/demultiplex"

for i in 191118_NS500556_0362_AHVYV7AFXY 191223_NS500556_0376_AHW7GNAFXY 191203_NS500556_0369_AHVYTYAFXY 200110_NS500556_0379_AHYNJVBGXC 200130_NS500556_0382_AHW7FNAFXY 200630_NovaSeqSP_RM_s3 200728_NovaSeq_S2_RM_s3;
do \
fq1=$demux/$i/$i.1.fq.gz;
fq2=$demux/$i/$i.2.fq.gz;
scitools fastq-split -X A $demux/$i/plate.annot $fq1 $fq2 &
done &

#Generate a list of fastq files per sequencing run
for i in 191118_NS500556_0362_AHVYV7AFXY 191223_NS500556_0376_AHW7GNAFXY 191203_NS500556_0369_AHVYTYAFXY 200110_NS500556_0379_AHYNJVBGXC 200130_NS500556_0382_AHW7FNAFXY 200630_NovaSeqSP_RM_s3 200728_NovaSeq_S2_RM_s3;
do ls $demux/$i/plate*.1.fq.gz ; done | sort | uniq > pool_list.txt

outdir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis"
demux="/home/groups/oroaklab/demultiplex"

#Now going to concatenate fastq files by plate annotation

#s3 WGS on GM12878

cat \
$demux/191118_NS500556_0362_AHVYV7AFXY/plate.GM12878_WGS_D.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.GM12878_WGS_D.1.fq.gz \
> $outdir/GM12878_WGS_D.1.fq.gz &
cat \
$demux/191118_NS500556_0362_AHVYV7AFXY/plate.GM12878_WGS_D.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.GM12878_WGS_D.2.fq.gz \
> $outdir/GM12878_WGS_D.2.fq.gz &

cat \
$demux/200728_NovaSeq_S2_RM_s3/plate.GM12878_WGS_B.1.fq.gz \
> $outdir/GM12878_WGS_B.1.fq.gz &
cat \
$demux/200728_NovaSeq_S2_RM_s3/plate.GM12878_WGS_B.2.fq.gz \
> $outdir/GM12878_WGS_B.2.fq.gz &


# s3ATAC

cat \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.Barnyard_ATAC_A.1.fq.gz \
$demux/200630_NovaSeqSP_RM_s3/plate.Barnyard_ATAC_A.1.fq.gz \
> $outdir/Barnyard_ATAC_A.1.fq.gz &
cat \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.Barnyard_ATAC_A.2.fq.gz \
$demux/200630_NovaSeqSP_RM_s3/plate.Barnyard_ATAC_A.2.fq.gz \
> $outdir/Barnyard_ATAC_A.2.fq.gz &

cat \
$demux/191203_NS500556_0369_AHVYTYAFXY/plate.Barnyard_ATAC_B.1.fq.gz \
> $outdir/Barnyard_ATAC_B.1.fq.gz &
cat \
$demux/191203_NS500556_0369_AHVYTYAFXY/plate.Barnyard_ATAC_B.2.fq.gz \
> $outdir/Barnyard_ATAC_B.2.fq.gz &

cat \
$demux/191203_NS500556_0369_AHVYTYAFXY/plate.Barnyard_ATAC_C.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.Barnyard_ATAC_C.1.fq.gz \
> $outdir/Barnyard_ATAC_C.1.fq.gz &
cat \
$demux/191203_NS500556_0369_AHVYTYAFXY/plate.Barnyard_ATAC_C.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.Barnyard_ATAC_C.2.fq.gz \
> $outdir/Barnyard_ATAC_C.2.fq.gz &

cat \
$demux/200130_NS500556_0382_AHW7FNAFXY/plate.Barnyard_ATAC_D.1.fq.gz \
> $outdir/Barnyard_ATAC_D.1.fq.gz &
cat \
$demux/200130_NS500556_0382_AHW7FNAFXY/plate.Barnyard_ATAC_D.2.fq.gz \
> $outdir/Barnyard_ATAC_D.2.fq.gz &

cat \
$demux/200630_NovaSeqSP_RM_s3/plate.Barnyard_ATAC_E.1.fq.gz \
> $outdir/Barnyard_ATAC_E.1.fq.gz &
cat \
$demux/200630_NovaSeqSP_RM_s3/plate.Barnyard_ATAC_E.2.fq.gz \
> $outdir/Barnyard_ATAC_E.2.fq.gz &

#CRC Cell line s3WGS and s3GCC

cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4442_WGS_E.1.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442_WGS_E.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442_WGS_E.1.fq.gz \
> $outdir/CRC4442_WGS_E.1.fq.gz &
cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4442_WGS_E.2.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442_WGS_E.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442_WGS_E.2.fq.gz \
> $outdir/CRC4442_WGS_E.2.fq.gz &

cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4442_WGS_F.1.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442_WGS_F.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442_WGS_F.1.fq.gz \
> $outdir/CRC4442_WGS_F.1.fq.gz &
cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4442_WGS_F.2.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442_WGS_F.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442_WGS_F.2.fq.gz \
> $outdir/CRC4442_WGS_F.2.fq.gz &

cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4671_WGS_G.1.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4671_WGS_G.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4671_WGS_G.1.fq.gz \
> $outdir/CRC4671_WGS_G.1.fq.gz &
cat \
$demux/191223_NS500556_0376_AHW7GNAFXY/plate.CRC-4671_WGS_G.2.fq.gz \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4671_WGS_G.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4671_WGS_G.2.fq.gz \
> $outdir/CRC4671_WGS_G.2.fq.gz &

cat \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442CRC-4671_GCC_H.1.fq.gz \
$demux/200630_NovaSeqSP_RM_s3/plate.CRC-4442CRC-4671_GCC_H.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442CRC4671_GCC_H.1.fq.gz \
> $outdir/CRC-4442CRC-4671_GCC_H.1.fq.gz &
cat \
$demux/200110_NS500556_0379_AHYNJVBGXC/plate.CRC-4442CRC-4671_GCC_H.2.fq.gz \
$demux/200630_NovaSeqSP_RM_s3/plate.CRC-4442CRC-4671_GCC_H.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CRC4442CRC4671_GCC_H.2.fq.gz \
> $outdir/CRC-4442CRC-4671_GCC_H.2.fq.gz &


#Cell line s3WGS and s3GCC
cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_GCC_L.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_GCC_L.1.fq.gz \
> $outdir/CellLine_GCC_L.1.fq.gz &
cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_GCC_L.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_GCC_L.2.fq.gz \
> $outdir/CellLine_GCC_L.2.fq.gz &

cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_WGS_J.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_WGS_J.1.fq.gz \
> $outdir/CellLine_WGS_J.1.fq.gz &
cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_WGS_J.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_WGS_J.2.fq.gz \
> $outdir/CellLine_WGS_J.2.fq.gz &

cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_WGS_K.1.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_WGS_K.1.fq.gz \
> $outdir/CellLine_WGS_K.1.fq.gz &
cat \
$demux/200630_NovaSeqSP_RM_s3/plate.CellLine_WGS_K.2.fq.gz \
$demux/200728_NovaSeq_S2_RM_s3/plate.CellLine_WGS_K.2.fq.gz \
> $outdir/CellLine_WGS_K.2.fq.gz &



for i in Barnyard*1.fq.gz; do
fq1=$i
fq2=${fq1::-8}.2.fq.gz
out_name=${fq1::-8}
scitools fastq-align -n -r 20 -t 20 hg38mm10 $out_name $fq1 $fq2; done ;
for i in CRC*1.fq.gz; do
fq1=$i
fq2=${fq1::-8}.2.fq.gz
out_name=${fq1::-8}
scitools fastq-align -n -r 20 -t 20 hg38 $out_name $fq1 $fq2; done ;
for i in CellLine*1.fq.gz; do
fq1=$i
fq2=${fq1::-8}.2.fq.gz
out_name=${fq1::-8}
scitools fastq-align -n -r 20 -t 20 hg38 $out_name $fq1 $fq2; done ;
for i in GM12878*1.fq.gz; do
fq1=$i
fq2=${fq1::-8}.2.fq.gz
out_name=${fq1::-8}
scitools fastq-align -n -r 20 -t 20 hg38 $out_name $fq1 $fq2; done ;

```


```bash
#Performing sorted bam-rmdup on bam files
outdir="/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis"
for i in *nsrt.bam; do
scitools bam-rmdup -n -t 30 $outdir/$i ; done &

#perform a directory clean up
mkdir complexity_data; mv *complexity* complexity_data
mkdir raw_fq; mv *fq.gz raw_fq
mkdir raw_alignment; mv *nsrt.bam raw_alignment; mv *align.log raw_alignment
#note, this should be revisited and split up by data type maybe?

#bash loop on scitools bam-filter function for each plate read cut off 
#(can differ within experiment due to changes to sequenced depth)
#s3atac barnyard filter dependent on library sequencing depth
i="Barnyard_ATAC_A.bbrd.q10.bam"; scitools bam-filter -N 15000 -c 50 -C ../complexity_data/${i::-13}.complexity.txt $i &
i="Barnyard_ATAC_B.bbrd.q10.bam"; scitools bam-filter -N 1000 -c 90 -C ../complexity_data/${i::-13}.complexity.txt $i &
i="Barnyard_ATAC_C.bbrd.q10.bam"; scitools bam-filter -N 10000 -c 80 -C ../complexity_data/${i::-13}.complexity.txt $i &
i="Barnyard_ATAC_D.bbrd.q10.bam"; scitools bam-filter -N 1000 -c 80 -C ../complexity_data/${i::-13}.complexity.txt $i &
i="Barnyard_ATAC_E.bbrd.q10.bam"; scitools bam-filter -N 10000 -c 50 -C ../complexity_data/${i::-13}.complexity.txt $i &

#s3wgs/gcc crc filter -N 100000
for i in CRC*bbrd.q10.bam; do scitools bam-filter -N 100000 $i & done &

#s3wgs gm12878 -N 100000 with differing -C
scitools bam-filter -N 1000000 -C ./complexity_data/GM12878_WGS_B.complexity.txt GM12878_WGS_B.bbrd.q10.bam &
scitools bam-filter -N 1000000 -C ./complexity_data/GM12878_WGS_D.complexity.txt -c 50 GM12878_WGS_D.bbrd.q10.bam &

#s3wgs/gcc cell lines -N 10000 -c 85
for i in CellLine*bbrd.q10.bam; do scitools bam-filter -N 10000 -c 85 -C ./complexity_data/${i::-13}.complexity.txt $i & done &

##merge bams by experiment
#barnyard atac
scitools bam-merge s3atac_barnyard.bbrd.q10.filt.bam Barnyard_ATAC_A.bbrd.q10.filt.bam Barnyard_ATAC_B.bbrd.q10.filt.bam Barnyard_ATAC_C.bbrd.q10.filt.bam Barnyard_ATAC_D.bbrd.q10.filt.bam Barnyard_ATAC_E.bbrd.q10.filt.bam &
#crc wgs
scitools bam-merge s3wgs_crc.bbrd.q10.filt.bam CRC4442_WGS_E.bbrd.q10.filt.bam CRC4442_WGS_F.bbrd.q10.filt.bam CRC4671_WGS_G.bbrd.q10.filt.bam &
#crc gcc
cp CRC-4442CRC-4671_GCC_H.bbrd.q10.filt.bam s3gcc_crc.bbrd.q10.filt.bam &
#gm12878 wgs
scitools bam-merge s3wgs_gm12878.bbrd.q10.filt.bam GM12878_WGS_D.bbrd.q10.filt.bam GM12878_WGS_B.bbrd.q10.filt.bam &
#cell line wgs
scitools bam-merge s3wgs_cellline.bbrd.q10.filt.bam CellLine_WGS_K.bbrd.q10.filt.bam CellLine_WGS_J.bbrd.q10.filt.bam &
#cell line gcc
cp CellLine_GCC_L.bbrd.q10.filt.bam s3gcc_cellline.bbrd.q10.filt.bam &

mv s3gcc* ../s3gcc_data/
mv s3atac* ../s3atac_data/
mv s3wgs* ../s3wgs_data/

```


#Files for upload (after all other processing)


```bash
#s3gcc/wgs
cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3wgs_data/singlecell_bam
ls sorted_bams*bam.sorted.bam > bamfile.list
samtools cat --threads 20 -b bamfile.list > s3wgsgcc.filtered.bam &
cp ../../files_for_upload

#s3atac
cd /home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/s3atac_data
cp mm10.bbrd.q10.bam hg38.bbrd.q10.bam mm10_SeuratObject.Rds hg38_SeuratObject.Rds ../files_for_upload
cp ./subclustering/hg38_inhibitory_neuron_SeuratObject.Rds ../../files_for_upload
#rename files for readability
mv mm10.bbrd.q10.bam s3atac.mm10.filtered.bam
mv hg38.bbrd.q10.bam s3atac.hg38.filtered.bam
mv hg38_SeuratObject.Rds s3atac.hg38_SeuratObject.Rds
mv mm10_SeuratObject.Rds s3atac.mm10_SeuratObject.Rds
mv hg38_inhibitory_neuron_SeuratObject.Rds s3atac.hg38_inhibitory_neuron_SeuratObject.Rds

```
```R
#Extracting counts matrix and meta data per cell (in R)
library(Signac)
setwd("/home/groups/oroaklab/adey_lab/projects/sciWGS/200730_s3FinalAnalysis/files_for_upload")

extract_data<-function(object,prefix){
    dat<-readRDS(object)
    write.csv(as.data.frame(dat@meta.data),quote=F,file=paste0(prefix,".metadata.csv"))
    write.csv(as.data.frame(dat[["peaks"]]@counts),quote=F,file=paste0(prefix,".counts.csv"))
    
}

file_list<-c("s3atac.hg38"="s3atac.hg38_SeuratObject.Rds",
            "s3atac.mm10"="s3atac.mm10_SeuratObject.Rds",
            "s3atac.hg38.inhibNeu"="s3atac.hg38_inhibitory_neuron_SeuratObject.Rds")

lapply(1:length(file_list),function(x) extract_data(object=file_list[x],prefix=names(file_list)[x]))

#With all files in place, gzipped and tarballed
gzip * &
tar -czf s3.tgz *gz &

```
