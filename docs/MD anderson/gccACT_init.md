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


