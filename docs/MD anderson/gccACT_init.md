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
--ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls \
--create-fastq-for-index-reads
```

Split out gccACT reads via python script
403,575,855 total reads (expecting 400M)

```python
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#set up plate
plate2_i7=['AGGCAGAA', 'TCCTGAGC', 'GGACTCCT', 'TAGGCATG', 'CTCTCTAC', 'CAGAGAGG', 'GCTACGCT', 'CGAGGCTG', 'AAGAGGCA', 'GTAGAGGA', 'GCTCATGA', 'ATCTCAGG', 'ACTCGCTA', 'GGAGCTAC', 'GCGTAGTA', 'CGGAGCCT']
plate2_i5=['TCTACTCT', 'CTCCTTAC', 'TATGCAGT', 'TACTCCTT', 'AGGCTTAG', 'ATTAGACG', 'CGGAGAGA', 'CTAGTCGA', 'AGCTAGAA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA']

idx_array=[]
i=1
for i5 in plate2_i5:
	for i7 in plate2_i7:
		idx_well=[i5,i7,"C"+str(i)]
		idx_array.append(idx_well)
		i+=1

fq1=sys.argv[1] #read argument fq1="~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.fastq.gz"
fq2=sys.argv[2] #read argument fq2="~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.fastq.gz"
idx3=sys.argv[3] #idx3=~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I1_001.fastq.gz"
idx4=sys.argv[4] #idx4=~/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I2_001.fastq.gz"

i=0
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1:
	with gzip.open(fq2, "rt") as handle2:
		with gzip.open(idx3, "rt") as handle3:
			with gzip.open(idx4, "rt") as handle4:
				with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
					with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
						for (title1, seq1, qual1), (title2, seq2, qual2), (title3,seq3,qual3), (title3,seq4,qual4) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2),FastqGeneralIterator(handle3),FastqGeneralIterator(handle4)):
							for j in idx_array:
								if seq3==j[1]+"AT" and seq4==j[0]+"GT":
									i+=1
									readname=j[2]+"_"+seq3+seq4
									outfile_fq1.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq1, qual1))
									outfile_fq2.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq2, qual2))

```
Running fastq splitter
This can be made a lot faster using something like a hash table, or limiting search space more.

```bash
python ~/src/plate2_fastqsplitter.py /volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R1_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_R2_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I1_001.fastq.gz \ 
/volumes/USR2/Ryan/fastq/230306_VH00219_371_AACJJFWM5/Undetermined_S0_L001_I2_001.fastq.gz 

#then gzip

gzip *fastq &
```

Got 253,517,850 reads total from assignment

### HiC Data Analysis can be done with the DipC group's released hickit
https://github.com/lh3/hickit

```bash
mkdir 
ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

#Map reads with BWA Mem
bwa mem -t 20 $ref read1.fq.gz read2.fq.gz | gzip > aln.sam.gz

# Map Dip-C reads and extract contacts (skip if you use your own pipeline)
seqtk mergepe read1.fq.gz read2.fq.gz | ~/tools/hickit-0.1_x64-linux/pre-dip-c - | bwa mem -5SP -p hs37d5.fa - | gzip > aln.sam.gz
~tools/hickit-0.1_x64-linux/k8 hickit.js vcf2tsv phased.vcf > phased_SNP.tsv   # extract phased SNPs from VCF
~tools/hickit-0.1_x64-linux/k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ~tools/hickit-0.1_x64-linux/k8 hickit.js chronly - | ~tools/hickit-0.1_x64-linux/k8 hickit.js bedflt par.bed - | gzip > contacts.seg.gz # for male
#./k8 hickit.js sam2seg -v phased_SNP.tsv aln.sam.gz | ./k8 hickit.js chronly -y - | gzip > contacts.seg.gz # for female
~tools/hickit-0.1_x64-linux/hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz  # optional

# Impute phases (-i also works with contacts.seg.gz)
~tools/hickit-0.1_x64-linux/hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz
~tools/hickit-0.1_x64-linux/hickit -i contacts.pairs.gz --out-val=impute.val     # estimate imputation accuracy by holdout
# Infer 3D structure
~tools/hickit-0.1_x64-linux/hickit -i impute.pairs.gz -Sr1m -c1 -r10m -c5 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O imput.3dg

# 2D contact map in PNG (bin size determined by the image width)
~tools/hickit-0.1_x64-linux/hickit -i impute.pairs.gz --out-png impute.png
# Compute CpG density (optional)
~tools/hickit-0.1_x64-linux/hickit.js gfeat -r hs37d5.fa.gz imput.3dg | gzip > imput.cpg.3dg.gz
# Visualize 3D structure (requiring a graphical card)
~tools/hickit-0.1_x64-linux/hickit-gl -I imput.cpg.3dg.gz --view
```


```bash
bwa mem -t 10 ~/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
~/fastq/230124_hicACT_4cell/Undetermined_S0_L001_R1_001.fastq.gz | samtools view -b - > 230124_gccACT.bam &
```


