---
title: gccACT-seq 4 Cell Test
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /gccACT_4cell/
category: mdanderson
---

## Initial analysis of Miseq Nano run on 4 cells with gccACT

Moved sequencing run raw files to 
/volumes/USR2/Ryan/seq/230124_M01842_0058_000000000-D5YMC

Ran 4 cells on the run (nano run so ~2M reads expected). 

```bash
#run transfered to /volumes/seq/tmp
cd /volumes/USR2/Ryan/seq/230124_M01842_0058_000000000-D5YMC
bcl2fastq -R /volumes/USR2/Ryan/seq/230124_M01842_0058_000000000-D5YMC \
-o ~/fastq/230124_hicACT_4cell/ \
-r 4 \
-p 10 \
-w 4 \
--ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --with-failed-reads \
--create-fastq-for-index-reads

cd /volumes/USR2/Ryan/seq/230202_M01842_0059_000000000-G677H
bcl2fastq -R /volumes/USR2/Ryan/seq/230202_M01842_0059_000000000-G677H \
-o ~/fastq/230202_hicACT_4cell/ \
-r 4 \
-p 10 \
-w 4 \
--ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --with-failed-reads \
--create-fastq-for-index-reads
```

Looks like it only wrote reads up to 86...hmmmm so no index reads. (It was because I used a 50 cycle kit)
Still going to align to see what it looks like I guess.

```bash
bwa mem -t 10 ~/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
~/fastq/230124_hicACT_4cell/Undetermined_S0_L001_R1_001.fastq.gz | samtools view -b - > 230124_gccACT.bam &
```
