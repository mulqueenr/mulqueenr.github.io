---
title: Tbr1 Mouse Knockout
layout: analysis
author: Ryan Mulqueen and Marissa Co
permalink: /tbr1/
category: sciATAC
---

# Processing for sciATAC portion for Tbr1 patient-specific mutation mouse models

Includes conversion of bcl to fastq files, barcode assignment, fastq splitting, alignment, removal of duplicate reads, checking sequencing depth, calling peaks, and looking at TSS enrichment.

Full breakdown of experimental setup is located [here.](https://docs.google.com/spreadsheets/d/1Px1OAE8vIi3GUXPny7OaVvYJHgGESCp4fyZKnLWW0UE/edit#gid=823628902)

## BCL File Locations

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
  #"firstplates" Plates 1,10 (10% of run pool); Plates 1,2,10
  /home/groups/oroaklab/seq/madbum/201116_NS500556_0437_AH72CMAFX2
  /home/groups/oroaklab/seq/madbum/201203_NS500556_0442_AH7FJGBGXG
  #"secondplates" Plates 3,4,5,7,6(5% spike-in); Plates 8,11,13
  /home/groups/oroaklab/seq/madbum/210202_NS500556_0456_AHK5C5BGXH
  /home/groups/oroaklab/seq/madbum/210210_NS500556_0458_AHVGCTBGXG
  #"thirdplates" Plates 9,12,14,2(5% spike-in),10(5% spike-in)
  /home/groups/oroaklab/seq/madbum/210223_NS500556_0463_AHNNH7BGXH
```
{% endcapture %} {% include details.html %} 


## BCL to FASTQ Conversion
For this we use a wrapper function NextSeq2fastq which wraps around bcl2fastq (v2.19). The wrapper is just to make it easier, since it infers where the run folder and output folders are based on our directory structure on the clusters. 

This is read in bcl files from the raw run folder in:
/home/groups/oroaklab/seq/madbum

And output fastq files in:
/home/groups/oroaklab/fastq

We repeated for all runs in firstplates, secondplates, and thirdplates.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
  run_name="201203_NS500556_0442_AH7FJGBGXG"
  NextSeq2fastq -R $run_name

  #This assumes the following variables for bcl2fastq:
  # DEFAULT VARIABLES
  #$run_path = "/home/groups/oroaklab/seq/madbum";
  #$fastq_paths = "/home/groups/oroaklab/fastq,/home/groups/oroakdata/fastq";
  #$bcl2fastq_version = "bcl2fastq/2.19.0";
  #$bcl_opts = "with-failed-reads,no-lane-splitting,fastq-compression-level=9,create-fastq-for-index-reads";
  #$bcl_ignore_opts = "with-failed-reads,no-lane-splitting,fastq-compression-level=9,create-fastq-for-index-reads,ignore-missing-bcls,ignore-missing-filter,ignore-missing-positions,ignore-missing-controls";
  #$run_processing_log_file = "/home/groups/oroaklab/fastq/run_processing.log";
  #@POSSIBLE_OUTS = ("Undetermined_S0_R1_001.fastq.gz", "Undetermined_S0_R2_001.fastq.gz", "Undetermined_S0_I1_001.fastq.gz", "Undetermined_S0_I2_001.fastq.gz");

```
{% endcapture %} {% include details.html %} 


## Demultiplexing the FASTQ files

After fastq files are generated we can then demultiplex them. By this, I mean that we are going to assign our index sequences based on the index reads from the run. 

Index cycles on the Nextseq are substantially more error prone than read cycles, so we account for an error rate. For our 8 and 10 bp indexes, we allow 2 base mismatches for each (Hamming distance of 2). This is enough to still unambiguously assign the proper original primer for any index.

We use a scitools function which is a perl script to do this. This demultiplexer also assumes the same directory structure as NextSeq2fastq, meaning only the run name must be specified on our clusters.

*The script does the following:*
  * Reads in the supplied index files. By default the index file is located here: 

  ```bash
  /home/groups/oroaklab/src/scitools/scitools-dev/SCI_Indexes.txt. 
  ```
  These are of tab-separated format: 

| IDX_NAME | IDX_NUMBER | IDX_SEQUENCE |
  
  * For each index, it then creates a hash (think python dictionary) for possible base changes for each index.
  * It then reads in the fastq data, and splits index reads to appropriate lengths for sci-chemistry.
  * It assigns the proper index sequence to all indexes based on the hash lookup table.
  * It writes out properly assigned reads (all four indexes have a proper match in the hash) in the sci-format, where the read name becomes the corrected list of indexes (referred to as a library barcode).


{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
  scitools fastq-dump -R $run_name
```

{% endcapture %} {% include details.html %} 


This will output to:
```bash
 /home/groups/oroaklab/demultiplex/201116_NS500556_0437_AH72CMAFX2
 /home/groups/oroaklab/demultiplex/201203_NS500556_0442_AH7FJGBGXG
 /home/groups/oroaklab/demultiplex/210202_NS500556_0456_AHK5C5BGXH
 /home/groups/oroaklab/demultiplex/210210_NS500556_0458_AHVGCTBGXG
 /home/groups/oroaklab/demultiplex/210223_NS500556_0463_AHNNH7BGXH

```
I then set up a working directory and moved the properly assigned reads.

{% capture summary %} Code (firstplates) {% endcapture %} {% capture details %}  

```bash
  mkdir /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates

  mv /home/groups/oroaklab/demultiplex/201116_NS500556_0437_AH72CMAFX2/201116_NS500556_0437_AH72CMAFX2.1.fq.gz \
  /home/groups/oroaklab/demultiplex/201116_NS500556_0437_AH72CMAFX2/201116_NS500556_0437_AH72CMAFX2.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates

  mv /home/groups/oroaklab/demultiplex/201203_NS500556_0442_AH7FJGBGXG/201203_NS500556_0442_AH7FJGBGXG.1.fq.gz \
  /home/groups/oroaklab/demultiplex/201203_NS500556_0442_AH7FJGBGXG/201203_NS500556_0442_AH7FJGBGXG.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates
```
{% endcapture %} {% include details.html %} 


{% capture summary %} Code (secondplates) {% endcapture %} {% capture details %}  

```bash
  mkdir /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates

  mv /home/groups/oroaklab/demultiplex/210202_NS500556_0456_AHK5C5BGXH/210202_NS500556_0456_AHK5C5BGXH.1.fq.gz \
  /home/groups/oroaklab/demultiplex/210202_NS500556_0456_AHK5C5BGXH/210202_NS500556_0456_AHK5C5BGXH.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates

  mv /home/groups/oroaklab/demultiplex/210210_NS500556_0458_AHVGCTBGXG/210210_NS500556_0458_AHVGCTBGXG.unassigned.1.fq.gz \
  /home/groups/oroaklab/demultiplex/210210_NS500556_0458_AHVGCTBGXG/210210_NS500556_0458_AHVGCTBGXG.unassigned.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (thirdplates) {% endcapture %} {% capture details %}  

```bash
  mkdir /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210224_thirdplates

  mv /home/groups/oroaklab/demultiplex/210223_NS500556_0463_AHNNH7BGXH/210223_NS500556_0463_AHNNH7BGXH.1.fq.gz \
  /home/groups/oroaklab/demultiplex/210223_NS500556_0463_AHNNH7BGXH/210223_NS500556_0463_AHNNH7BGXH.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210224_thirdplates
```
{% endcapture %} {% include details.html %}


### Generation of thorough annotation file and all meta data per cell 

Now that we have reads that assign to known scitools indexes, we have to get more specific. We are going to generate a proper annotation for our experiment based on our PCR and Tn5 primers used. We will do this for all possible index combinations as a ".annot" file.
scitools assumes annot files are in the following format: BARCODE   ANNOTATION

BARCODE structure is based on the sequencer. Because we use a Nextseq most commonly, we set them up as:

| 8bp_Tn5_i7_idx | 10bp_PCR_i7_idx | 8bp_Tn5_i5_idx | 10bp_PCR_i5_idx|

We used six tn5 plates: AB,CB,CC,BB,AC,CA (listed as Tn5_i5_idx,Tn5_i7_idx) so we will also limit barcodes to just those attainable from these tn5 combinations.

For this, I am once again looking at our experimental design from  [here.](https://docs.google.com/spreadsheets/d/1Px1OAE8vIi3GUXPny7OaVvYJHgGESCp4fyZKnLWW0UE/edit#gid=823628902)

I'm going to use a scitools helper function to do this, but a simple for loop through the index master list would work as well.

{% capture summary %} Code (firstplates) {% endcapture %} {% capture details %}  

``` bash
#Since all plates are a random assortment of all Tn5 tagmentation, we can generate a simplified annotation schematic for PCR plates.

#simplified experiment annot, I just made this in the same directory with nano text editor.
#firstplates_annot.txt
Plate   Plate_SDSBSA_Condition  PCR_Index_i5    PCR_Index_i7    PCR_Cycles
1       Fresh   A       A       17
10      Old     A       B       17
2       Fresh   C       B       17

scitools make-annot \
Plate_1+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,AA=ALL \
Plate_10+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,AB=ALL \
Plate_2+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CB=ALL \
> firstplates.annot
```
{% endcapture %} {% include details.html %} 


{% capture summary %} Code (secondplates) {% endcapture %} {% capture details %}  

``` bash
#secondplates_annot.txt
Plate   Plate_SDSBSA_Condition  PCR_Index_i5    PCR_Index_i7    PCR_Cycles
3	Fresh	C	C	17
4	Fresh	C	D	17
5	Fresh	G	F	17
6	Fresh	G	G	17
7	Fresh	H	D	17
8	Fresh	H	E	17
11	Old	C	E	17
13	Old	G	H	17

scitools make-annot \
Plate_3+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CC=ALL \
Plate_4+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CD=ALL \
Plate_5+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,GF=ALL \
Plate_6+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,GG=ALL \
Plate_7+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,HD=ALL \
Plate_8+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,HE=ALL \
Plate_11+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CE=ALL \
Plate_13+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,GH=ALL \
> secondplates.annot
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (thirdplates) {% endcapture %} {% capture details %}  

``` bash
#thirdplates_annot.txt
Plate   Plate_SDSBSA_Condition  PCR_Index_i5    PCR_Index_i7    PCR_Cycles
2	Fresh	C	B	17
9	Fresh	H	F	17
10	Old	A	B	17
12	Old	C	F	17
14	Old	H	G	17

scitools make-annot \
Plate_2+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CB=ALL \
Plate_9+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,HF=ALL \
Plate_10+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,AB=ALL \
Plate_12+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,CF=ALL \
Plate_14+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,HG=ALL \
> thirdplates.annot
```
{% endcapture %} {% include details.html %}


### Splitting out our reads from the demultiplexed fastqs

Now that we know which barcodes belong to our reads, we can split them out from the full pool.

To do this we will use a scitools function that looks at fastq read 1 and read 2 and splits it into new files based on matches to our annotation.

{% capture summary %} Code (firstplates) {% endcapture %} {% capture details %}  

```bash
#Combine fastq files across runs
cat 201116_NS500556_0437_AH72CMAFX2.1.fq.gz 201203_NS500556_0442_AH7FJGBGXG.1.fq.gz > tbr1.1.fq.gz
cat 201116_NS500556_0437_AH72CMAFX2.2.fq.gz 201203_NS500556_0442_AH7FJGBGXG.2.fq.gz > tbr1.2.fq.gz

scitools fastq-split -X -A firstplates.annot \
tbr1.1.fq.gz \
tbr1.2.fq.gz &

#The -X flag tells it to not write out barcodes which don't match. Those would be other sci formatted experiments on the same run

#Annot: Plate_1, count = 85746573
#Annot: Plate_10, count = 78801321
#Annot: Plate_2, count = 70776718
```
{% endcapture %} {% include details.html %} 


{% capture summary %} Code (secondplates) {% endcapture %} {% capture details %}  

```bash
#Combine fastq files across runs
cat 210202_NS500556_0456_AHK5C5BGXH.1.fq.gz 210210_NS500556_0458_AHVGCTBGXG.unassigned.1.fq.gz > 210202_210210_merged.1.fq.gz
cat 210202_NS500556_0456_AHK5C5BGXH.2.fq.gz 210210_NS500556_0458_AHVGCTBGXG.unassigned.2.fq.gz > 210202_210210_merged.2.fq.gz

scitools fastq-split -X -A secondplates.annot \
210202_210210_merged.1.fq.gz \
210202_210210_merged.2.fq.gz &

#Annot: Plate_7, count = 83652013
#Annot: Plate_11, count = 129213990
#Annot: Plate_8, count = 90351664
#Annot: Plate_4, count = 76915941
#Annot: Plate_13, count = 85242513
#Annot: Plate_3, count = 69656977
#Annot: Plate_5, count = 77743124
#Annot: Plate_6, count = 17561992
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (thirdplates) {% endcapture %} {% capture details %}  

```bash
scitools fastq-split -X -A thirdplates.annot \
210223_NS500556_0463_AHNNH7BGXH.1.fq.gz \
210223_NS500556_0463_AHNNH7BGXH.2.fq.gz &

#Annot: Plate_12, count = 98462906
#Annot: Plate_10, count = 7507432
#Annot: Plate_14, count = 100526129
#Annot: Plate_2, count = 14610178
#Annot: Plate_9, count = 112807352
```
{% endcapture %} {% include details.html %}


## Alignment

We have our reads, so now we can align them to the mouse reference genome. ATAC data is count based. 

We use another scitools function for convenience. It wraps bwa mem. We will use -t 10 threads for alignment and -r 10 threads for samtools sort afterwards.

{% capture summary %} Code (firstplates) {% endcapture %} {% capture details %}  

```bash
  #For plate 1
  scitools fastq-align -t 10 -r 10 mm10 plate1 firstplates.Plate_1.1.fq.gz firstplates.Plate_1.2.fq.gz &
  #For plate 10
  scitools fastq-align -t 10 -r 10 mm10 plate10 firstplates.Plate_10.1.fq.gz firstplates.Plate_10.2.fq.gz &
  #For plate 2
  scitools fastq-align -t 10 -r 10 mm10 plate2 firstplates.Plate_2.1.fq.gz firstplates.Plate_2.2.fq.gz &
```
{% endcapture %} {% include details.html %} 


{% capture summary %} Code (secondplates) {% endcapture %} {% capture details %}  

```bash
  #For plate 3
  scitools fastq-align -t 10 -r 10 mm10 plate3 secondplates.Plate_3.1.fq.gz secondplates.Plate_3.2.fq.gz &
  #For plate 4
  scitools fastq-align -t 10 -r 10 mm10 plate4 secondplates.Plate_4.1.fq.gz secondplates.Plate_4.2.fq.gz &
  #For plate 5
  scitools fastq-align -t 10 -r 10 mm10 plate5 secondplates.Plate_5.1.fq.gz secondplates.Plate_5.2.fq.gz &
  #For plate 6
  scitools fastq-align -t 10 -r 10 mm10 plate6 secondplates.Plate_6.1.fq.gz secondplates.Plate_6.2.fq.gz &
  #For plate 7
  scitools fastq-align -t 10 -r 10 mm10 plate7 secondplates.Plate_7.1.fq.gz secondplates.Plate_7.2.fq.gz &
  #For plate 8
  scitools fastq-align -t 10 -r 10 mm10 plate8 secondplates.Plate_8.1.fq.gz secondplates.Plate_8.2.fq.gz &
  #For plate 11
  scitools fastq-align -t 10 -r 10 mm10 plate11 secondplates.Plate_11.1.fq.gz secondplates.Plate_11.2.fq.gz &
  #For plate 13
  scitools fastq-align -t 10 -r 10 mm10 plate13 secondplates.Plate_13.1.fq.gz secondplates.Plate_13.2.fq.gz &
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (thirdplates) {% endcapture %} {% capture details %}  

```bash
  #For plate 2
  scitools fastq-align -t 10 -r 10 mm10 plate2 thirdplates.Plate_2.1.fq.gz thirdplates.Plate_2.2.fq.gz &
  #For plate 9
  scitools fastq-align -t 10 -r 10 mm10 plate9 thirdplates.Plate_9.1.fq.gz thirdplates.Plate_9.2.fq.gz &
  #For plate 10
  scitools fastq-align -t 10 -r 10 mm10 plate10 thirdplates.Plate_10.1.fq.gz thirdplates.Plate_10.2.fq.gz &
  #For plate 12
  scitools fastq-align -t 10 -r 10 mm10 plate12 thirdplates.Plate_12.1.fq.gz thirdplates.Plate_12.2.fq.gz &
  #For plate 14
  scitools fastq-align -t 10 -r 10 mm10 plate14 thirdplates.Plate_14.1.fq.gz thirdplates.Plate_14.2.fq.gz &
```
{% endcapture %} {% include details.html %}


## Removal of Duplicate Reads

Once we have aligned reads, we can mark PCR duplicates. Because we are sampling across the genome, it is highly unlikely that we capture the same exact start and end region twice. So we can use a combination of our barcode, and the start and end positions of a read to mark duplication rates.

{% capture summary %} Code (firstplates) {% endcapture %} {% capture details %}  

```bash
  #For plate 1
  scitools bam-rmdup plate1.bam &
  #For plate 10
  scitools bam-rmdup plate10.bam &
  #For plate 2
  scitools bam-rmdup plate2.bam &

  #Once these finish, plot the complexity per cell
  scitools plot-complexity plate10.complexity.txt &
  scitools plot-complexity plate1.complexity.txt &
  scitools plot-complexity plate2.complexity.txt &
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (secondplates) {% endcapture %} {% capture details %}  

```bash
  #For plate 3
  scitools bam-rmdup plate3.bam &
  #For plate 4
  scitools bam-rmdup plate4.bam &
  #For plate 5
  scitools bam-rmdup plate5.bam &
  #For plate 6
  scitools bam-rmdup plate6.bam &
  #For plate 7
  scitools bam-rmdup plate7.bam &
  #For plate 8
  scitools bam-rmdup plate8.bam &
  #For plate 11
  scitools bam-rmdup plate11.bam &
  #For plate 13
  scitools bam-rmdup plate13.bam &

  #Once these finish, plot the complexity per cell
  scitools plot-complexity plate3.complexity.txt &
  scitools plot-complexity plate4.complexity.txt &
  scitools plot-complexity plate5.complexity.txt &
  scitools plot-complexity plate6.complexity.txt &
  scitools plot-complexity plate7.complexity.txt &
  scitools plot-complexity plate8.complexity.txt &
  scitools plot-complexity plate11.complexity.txt &
  scitools plot-complexity plate13.complexity.txt &
```
{% endcapture %} {% include details.html %}


{% capture summary %} Code (thirdplates) {% endcapture %} {% capture details %}  

```bash

  #Combine firstplates and thirdplates bam files for plates 2 and 10 (sequenced 2x)
  mv plate2.bam plate2-third.bam
  mv plate10.bam plate10-third.bam

  cp /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates/plate2.bam .
  cp /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates/plate10.bam .

  mv plate2.bam plate2-first.bam
  mv plate10.bam plate10-first.bam

  scitools bam-merge plate2.bam plate2-first.bam plate2-third.bam &
  scitools bam-merge plate10.bam plate10-first.bam plate10-third.bam &

  #For plate 2
  scitools bam-rmdup plate2.bam &
  #For plate 9
  scitools bam-rmdup plate9.bam &
  #For plate 10
  scitools bam-rmdup plate10.bam &
  #For plate 12
  scitools bam-rmdup plate12.bam &
  #For plate 14
  scitools bam-rmdup plate14.bam &

  #Once these finish, plot the complexity per cell
  scitools plot-complexity plate2.complexity.txt &
  scitools plot-complexity plate9.complexity.txt &
  scitools plot-complexity plate10.complexity.txt &
  scitools plot-complexity plate12.complexity.txt &
  scitools plot-complexity plate14.complexity.txt &
```
{% endcapture %} {% include details.html %}


## Checking Plates and Sequencing Depth

### Reads Per Plate

{% capture summary %} Code {% endcapture %} {% capture details %}

```bash
  for j in /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210224_thirdplates \
    /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates \
    /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates;
    do for i in ${j}/*Plate*1.fq.gz;
    do echo $i `zcat $i | grep "^@" - | wc -l` ; done ; done &
```
{% endcapture %} {% include details.html %}


| Plate Prep  | Plate   | Reads Devoted |
|:--------|:--------|:--------|
| 210224_thirdplates  | Plate_10 | 7507432 |
| 210224_thirdplates  | Plate_12 | 98462906 |
| 210224_thirdplates  | Plate_14 | 100526129 |
| 210224_thirdplates  | Plate_2 | 14610178 |
| 210224_thirdplates  | Plate_9 | 112807352 |
| 210212_secondplates | Plate_11 | 129213990 |
| 210212_secondplates | Plate_13 | 85242513 |
| 210212_secondplates | Plate_3 | 69656977 |
| 210212_secondplates | Plate_4 | 76915941 |
| 210212_secondplates | Plate_5 | 77743124 |
| 210212_secondplates | Plate_6 | 17561992 |
| 210212_secondplates | Plate_7 | 83652013 |
| 210212_secondplates | Plate_8 | 90351664 |
| 201117_firstplates  | Plate_10 | 78801321 |
| 201117_firstplates  | Plate_1 | 85746573 |
| 201117_firstplates  | Plate_2 | 70776718 |

### Complexity for All Plates

Make combined complexity file for all plates.

{% capture summary %} Code {% endcapture %} {% capture details %}
```bash
  #Create new working directory for all plates
  mkdir /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates
  cd /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates
  
  #Link to complexity files for all plates
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates/plate1.complexity.txt .
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates/plate*.complexity.txt .
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210224_thirdplates/plate*.complexity.txt .
  
  #Make combined complexity file
  for j in /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates;
  do for i in ${j}/plate*.complexity.txt;
  do awk 'OFS="\t" { split(FILENAME,a,"/");split(a[9],b,"[.]"); print $2,$3,$4,$5,b[1]}' $i ; done; done > allplates.complexity.txt
```
{% endcapture %} {% include details.html %}


Make table in R with cell counts, median and mean unique reads per cell, and mean percent unique reads per cell.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")
  dat<-read.table("allplates.complexity.txt",sep="\t",header=F)
  colnames(dat)<-c("cellID","total_reads","unique_reads","percent_uniq","plate_name")
  dat<-dat[dat$unique_reads>=1000,]
  library(dplyr)
  dat %>% group_by(plate_name) %>% summarize(cell_count=n(),median_uniq_reads=median(unique_reads),mean_uniq_reads=mean(unique_reads),mean_percent_uniq=mean(percent_uniq))

```

{% endcapture %} {% include details.html %}


| plate_name | cell_count | median_uniq_reads | mean_uniq_reads | mean_percent_uniq |
|:--------|:--------|:--------|:--------|:--------|
| plate1 | 3887 | 22143 | 31404 | 90.7 |
| plate2 | 6049 | 12883 | 20923 | 93.1 |
| plate3 | 4601 | 15122 | 22041 | 92.2 |
| plate4 | 4219 | 18567 | 26743 | 91.8 |
| plate5 | 3861 | 16097 | 27445 | 88.6 |
| plate6 | 2484 | 4360 | 9001 | 89.7 |
| plate7 | 4713 | 17300 | 25658 | 91.4 |
| plate8 | 3537 | 23041 | 35603 | 89.1 |
| plate9 | 3159 | 20149 | 38844 | 73.8 |
| plate10 | 5729 | 12540 | 22063 | 92.7 |
| plate11 | 4124 | 27138 | 41738 | 86.9 |
| plate12 | 5097 | 19173 | 27789 | 91.0 |
| plate13 | 3970 | 19556 | 30665 | 90.0 |
| plate14 | 4307 | 22072 | 33249 | 89.0 |

## Merging All BAM Files and Filtering

{% capture summary %} Code {% endcapture %} {% capture details %}

```bash
  #Link to final bam files for merging
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates/plate1.bbrd.q10.bam .
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210212_secondplates/plate*.bbrd.q10.bam .
  ln -s /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210224_thirdplates/plate*.bbrd.q10.bam .
  rm plate6.bbrd.q10.bam
  
  #Combine bam files
  scitools bam-merge tbr1_ko.bam plate1.bbrd.q10.bam plate2.bbrd.q10.bam \
  plate3.bbrd.q10.bam plate4.bbrd.q10.bam plate5.bbrd.q10.bam plate7.bbrd.q10.bam \
  plate8.bbrd.q10.bam plate9.bbrd.q10.bam plate10.bbrd.q10.bam plate11.bbrd.q10.bam \
  plate12.bbrd.q10.bam plate13.bbrd.q10.bam plate14.bbrd.q10.bam &
  
  #Make complexity plots for merged dataset

  #Make a scitools compatible complexity file [rownumber][cellid][totread][uniqread][percuniq]
  #Note: this will replace the complexity file for checking sequencing depth per cell
  for j in /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates; \
  do for i in ${j}/plate*.complexity.txt; do awk 'OFS="\t" {print $1,$2,$3,$4,$5}' $i ; \
  done; done > allplates.complexity.txt

  #Make a scitools compatible annotation file [cellid][annot]
  for j in /home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates; \
  do for i in ${j}/plate*.complexity.txt; do awk 'OFS="\t" { split(FILENAME,a,"/");split(a[9],b,"[.]"); \
  print $2,b[1]}' $i ; done; done > allplates.annot

  #Make merged complexity plot colored by plate
  scitools plot-complexity -A allplates.annot allplates.complexity.txt
```
{% endcapture %} {% include details.html %}


Based on the complexity plot we filtered the merged BAM file to exclude cells with <1000 unique reads per cell.

{% capture summary %} Code {% endcapture %} {% capture details %}

```bash
  
  #Filter bam
  scitools bam-filter -N 1000 tbr1_ko.bam

```
{% endcapture %} {% include details.html %}
 
 
## Peak Calling and TSS Enrichment
 
{% capture summary %} Code {% endcapture %} {% capture details %}

```bash

  #Get insert size distribution
  scitools isize tbr1_ko.filt.bam &
  
  #Get TSS enrichment per cell
  module load bedops/2.4.36
  
  #Bulk ENCODE method
  scitools bam-tssenrich -E tbr1_ko.filt.bam mm10 &

  #Per cell method
  scitools bam-tssenrich tbr1_ko.filt.bam mm10 &
```
{% endcapture %} {% include details.html %} 

### Alternative approach to peak calling

Using macs3 for improved peak calling
```bash
# call macs3 with shift and extensions
source /home/groups/oroaklab/nishida/scitools_env/bin/activate

macs3 callpeak -f BAMPE --keep-dup all -q 0.01 -t tbr1_ko.filt.bam -g mm -n tbr1_ko.macs3.keepdup &
```

## Tabix Fragment File Generation

Tabix file format is a tab separated multicolumn data structure.

| Column Number | Name | Description |
|:--------|:-------:|:--------|
|1 |chrom |  Reference genome chromosome of fragment |
|2 |chromStart | Adjusted start position of fragment on chromosome. |
|3 |chromEnd   | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
|4 |barcode | The 10x (or sci) cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment. |
|5 |duplicateCount |The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeq instrument are excluded from this count. |

{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
  #Tabix fragment file generation
  input_bam="tbr1_ko.filt.bam"
  output_name="tbr1_ko"
  tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
  bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
  samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
  $tabix -p bed $output_name.fragments.tsv.gz &

```

{% endcapture %} {% include details.html %} 


## Testing out called peaks against other available data sets

WashU study of mouse organogenesis

{% capture summary %} Code {% endcapture %} {% capture details %}
```bash
#Using https://atlas.gs.washington.edu/mouse-atac/ data set peaks on mouse development
wget http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/matrices/atac_matrix.binary.qc_filtered.peaks.txt
awk 'OFS="\t" {split($1,a,"_"); print a[1],a[2],a[3]}' atac_matrix.binary.qc_filtered.peaks.txt > washU_mousedevel.peaks.bed
bedtools intersect -u -wa -a tbr1_ko.filt.500.bed -b washU_mousedevel.peaks.bed | wc -l
#139474
wc -l tbr1_ko.filt.500.bed
#744124 tbr1_ko.filt.500.bed
wc -l washU_mousedevel.peaks.bed
#436206 washU_mousedevel.peaks.bed



#also checking macs3 peaks
bedtools intersect -u -wa -a tbr1_ko.macs3.keepdup_summits.bed -b washU_mousedevel.peaks.bed | wc -l
#24781
wc -l tbr1_ko.macs3.keepdup_summits.bed
#277699 tbr1_ko.macs3.keepdup_summits.bed

```

{% endcapture %} {% include details.html %} 

| In Ours Only | Shared | In WashU Only |
|:--------:|:-------:|:--------:|
| 604650 (81.26%) | 139474 (18.74% Ours ; 31.97% WashU)|  296732 (68.03%)|

 Using ARSN peak set generated from all ENCODE data

{% capture summary %} Code {% endcapture %} {% capture details %}

```bash 
#Using ARSN peak set
#/home/groups/oroaklab/refs/mm10/masterlistDHS/roughmerged.DHS_masterlist.bed

bedtools intersect -u -wa -a tbr1_ko.filt.500.bed -b /home/groups/oroaklab/refs/mm10/masterlistDHS/media-6.bed | wc -l
#691937
wc -l tbr1_ko.filt.500.bed
#744124 tbr1_ko.filt.500.bed
wc -l /home/groups/oroaklab/refs/mm10/masterlistDHS/media-6.bed
#1802603 /home/groups/oroaklab/refs/mm10/masterlistDHS/media-6.bed


#checking macs3
bedtools intersect -u -wa -a tbr1_ko.macs3.keepdup_summits.bed -b /home/groups/oroaklab/refs/mm10/masterlistDHS/media-6.bed | wc -l
#65571
wc -l tbr1_ko.macs3.keepdup_summits.bed
#277699 tbr1_ko.macs3.keepdup_summits.bed

```

| In Ours Only | Shared | In Masterset |
|:--------:|:-------:|:--------:|
| 52187 (7.02%) | 691937 (92.98% Ours ; 29.10% WashU)|  1685290 (70.89%)|


### Index bam and check peaks and read pile-ups via scitools plotting function

Annotation file was generated below after seurat object meta data was made

```bash
 samtools index -@ 20 -b tbr1_ko.filt.bam tbr1_ko.filt.bai & #generate an index
scitools plot-reads -A devel_line_geno.annot -B tbr1_ko.filt.500.bed -F 0.1 -p 0.1 -G mm10 tbr1_ko.filt.bam Pax6 Neurod1 Tbr1 

scitools plot-reads -A devel_line_geno.annot -B tbr1_ko.macs3_summits.bed -F 0.1 -p 0.1 -G mm10 tbr1_ko.filt.bam Pax6 Neurod1 Tbr1

```
{% endcapture %} {% include details.html %} 

Based on this, I'm going to continue with the conservative macs3 peaks for clustering.

### Counts matrix generation on called peaks
{% capture summary %} Code {% endcapture %} {% capture details %}  

```bash
  #Generate counts matrix for macs2 peak set
  scitools atac-callpeak tbr1_ko.filt.bam &
  scitools atac-counts tbr1_ko.filt.bam tbr1_ko.filt.500.bed &

  #Generate counts matrix for macs3 peak set
  scitools atac-counts -O tbr1_ko.macs3 tbr1_ko.filt.bam tbr1_ko.macs3.keepdup_summits.bed &
```
{% endcapture %} {% include details.html %} 

# sciATAC Full Processing in R

## Generating Seurat Objects

Using R v4.0 and Signac v1.0 for processing.

{% capture summary %} Code {% endcapture %} {% capture details %}  


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  # make counts matrix from sparse matrix
  #Make counts matrix from sparse matrix
  IN<-as.matrix(read.table("tbr1_ko.filt.500.counts.sparseMatrix.values.gz")) #old peak set
  #IN<-as.matrix(read.table("tbr1_ko.macs3.counts.sparseMatrix.values.gz"))
  IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
  COLS<-read.table("tbr1_ko.filt.500.counts.sparseMatrix.cols.gz") #old peak set
  #COLS<-read.table("tbr1_ko.macs3.counts.sparseMatrix.cols.gz")
  colnames(IN)<-COLS$V1
  ROWS<-read.table("tbr1_ko.filt.500.counts.sparseMatrix.rows.gz") #old peak set
  #ROWS<-read.table("tbr1_ko.macs3.counts.sparseMatrix.rows.gz")
  row.names(IN)<-ROWS$V1

  #Read in fragment path for coverage plots
  fragment.path="./tbr1_ko.fragments.tsv.gz"

  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

  # change to UCSC style since the data was mapped to mm10
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"

  #Generate ChromatinAssay Objects
  obj.chromassay <- CreateChromatinAssay(
    counts = IN,
    genome="mm10",
    min.cells = 100,
    annotation=annotations,
    sep=c("_","_"),
    fragments=fragment.path
  )

  #Create Seurat Object
  obj <- CreateSeuratObject(
    counts = obj.chromassay,
    assay = "peaks",
  )

  #Meta.data to be updated after clustering


  #saving unprocessed SeuratObject
  saveRDS(obj,file="tbr1_ko.SeuratObject.Rds") #old peak set
  #saveRDS(obj,file="tbr1_ko.macs3.SeuratObject.Rds") #old peak set

```
{% endcapture %} {% include details.html %} 


## Performing cisTopic and UMAP

Used multiple clustering and dim reduc attempts. Commented out irrelevant ones.

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  library(cisTopic)
  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")
  obj$cellID<-row.names(obj@meta.data)

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(10,20,22,24,26,28,30,40),nCores=8,addModels=FALSE)
      #atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(29,31,32,33,34,35,36,38),nCores=7,addModels=TRUE)   
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          
  cistopic_processing(seurat_input=obj,prefix="tbr1_ko")
  #cistopic_processing(seurat_input=obj,prefix="tbr1_ko.macs3.more")

  cistopic_models<-readRDS("tbr1_ko.CisTopicObject.Rds")
  #cistopic_models_more<-readRDS("tbr1_ko.macs3.more.CisTopicObject.Rds")

  #Setting up topic count selection
  pdf("tbr1_ko.cistopic_model_selection.pdf")
  par(mfrow=c(1,3))
  cistopic_models <- selectModel(cistopic_models, type='derivative')
  dev.off()
  system("slack -F tbr1_ko.cistopic_model_selection.pdf ryan_todo")

  #Setting up topic count selection
  #pdf("tbr1_ko.cistopic_modelmore_selection.pdf")
  #par(mfrow=c(1,3))
  #cistopic_models <- selectModel(cistopic_models_more, type='derivative')
  #dev.off()
  #system("slack -F tbr1_ko.cistopic_modelmore_selection.pdf ryan_todo")


  ###############################################
  #Loop through cistopic models
  cistopic_loop<-function(topic_number,object_input,models_input){
      models_input<-selectModel(models_input,select=as.numeric(topic_number))
      #perform UMAP on topics
      topic_df<-as.data.frame(models_input@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      print("Performed UMAP.")
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,as.data.frame(object_input@meta.data),by="cellID")
     
      #combine with seurat object    
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$umap<-cistopic_umap
      
      #finally plot
      plt<-DimPlot(object_input,size=0.1)+ggtitle(as.character(topic_number))
      return(plt)
  }

  library(patchwork)
  library(parallel)

  plt_list_more<-mclapply(names(cistopic_models@models), function(x) {
    cistopic_loop(topic_number=x,object_input=obj,models_input=cistopic_models)},mc.cores=1)

  plt_list_more<-wrap_plots(plt_list_more)
  ggsave(plt_list_more,file="tbr1_ko.umap_multipleTopicModels_clustering.png",height=20,width=60,limitsize=FALSE)
  system("slack -F tbr1_ko.umap_multipleTopicModels_clustering.png ryan_todo")
  ###############################################

  #set topics based on derivative
  selected_topic=28
  cisTopicObject<-cisTopic::selectModel(cistopic_models,select=selected_topic,keepModels=T)

  #saving model selected RDS
  saveRDS(cisTopicObject,file="tbr1_ko.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   
      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by="cellID")

      #Add cell embeddings into seurat
      cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
      colnames(cell_embeddings)<-cisTopicObject@cell.names
      n_topics<-nrow(cell_embeddings)
      row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
      cell_embeddings<-as.data.frame(t(cell_embeddings))

      #Add feature loadings into seurat
      feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
      row.names(feature_loadings)<-paste0("topic_",1:n_topics)
      feature_loadings<-as.data.frame(t(feature_loadings))

      #combined cistopic results (cistopic loadings and umap with seurat object)
      cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$cistopic<-cistopic_obj
      object_input@reductions$umap<-cistopic_umap

      n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))

      object_input <- FindNeighbors(
        object = object_input,
        reduction = 'cistopic',
        dims = 1:n_topics
      )
      object_input <- FindClusters(
        object = object_input,
        verbose = TRUE,
        resolution=resolution)

  ###save Seurat file
  return(object_input)}

  obj<-cistopic_wrapper(object_input=obj,cisTopicObject=cisTopicObject,resolution=0.5)
  saveRDS(obj,file="tbr1_ko.SeuratObject.Rds")
  
  #also perform simplified dim reduction
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 'q0')
  obj <- RunSVD(obj)
  saveRDS(obj,file="tbr1_ko.SeuratObject.Rds")


```
{% endcapture %} {% include details.html %} 


### Plotting and updating metadata

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  #renaming annot for simplified annotation file making
  #rename processing_ processing. *annot
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")
  dim(obj@meta.data)

  #Read in expected sci indexes
  #Set up indexes for merging
  sci_idx<-read.table("/home/groups/oroaklab/src/scitools/scitools-dev/SCI_Indexes.txt",header=F)
  sci_idx_tn5_i5<-sci_idx[sci_idx$V2=="3",]
  sci_idx_tn5_i5$set<-unlist(lapply(strsplit(sci_idx_tn5_i5$V1,"_"),"[",2)) #split string to get set name
  sci_idx_tn5_i5$row<-unlist(lapply(strsplit(sci_idx_tn5_i5$V1,"_"),"[",4)) #split string to get row name
  colnames(sci_idx_tn5_i5)<-c("tn5_i5_idx_name","tn5_i5_idx_cycle","tn5_i5_idx_seq","tn5_i5_set","tn5_i5_row")

  sci_idx_tn5_i7<-sci_idx[sci_idx$V2=="1",]
  sci_idx_tn5_i7$set<-unlist(lapply(strsplit(sci_idx_tn5_i7$V1,"_"),"[",2)) #split string to get set name
  sci_idx_tn5_i7$column<-unlist(lapply(strsplit(sci_idx_tn5_i7$V1,"_"),"[",4)) #split string to get column name
  colnames(sci_idx_tn5_i7)<-c("tn5_i7_idx_name","tn5_i7_idx_cycle","tn5_i7_idx_seq","tn5_i7_set","tn5_i7_row")

  #merge cellID with sci_idx_tn5 data frame based on idx_sequences
  obj$tn5_i7_idx_seq<-substr(obj$cellID,1,8)
  obj$pcr_i7_idx_seq<-substr(obj$cellID,9,18)
  obj$tn5_i5_idx_seq<-substr(obj$cellID,19,26)
  obj$pcr_i5_idx_seq<-substr(obj$cellID,27,36)

  obj@meta.data<-merge(obj@meta.data,sci_idx_tn5_i5,by="tn5_i5_idx_seq")
  obj@meta.data<-merge(obj@meta.data,sci_idx_tn5_i7,by="tn5_i7_idx_seq")
  row.names(obj@meta.data)<-obj$cellID

  #Read in sample information
  samp_info<-read.table("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/sample_id.txt",header=T,sep="\t")
  #Read in Tn5-linked sample annotation
  tn5_info<-read.table("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/tn5_annotation.txt",header=T,sep="\t")
  samp_info<-merge(samp_info,tn5_info,by="Sample_ID")
  samp_info$tn5_i5_set<-substr(samp_info$Plate,1,1) #split string to get set name
  samp_info$tn5_i7_set<-substr(samp_info$Plate,2,2) #split string to get set name
  dat<-obj@meta.data
  dat<-merge(dat,samp_info,by.x=c("tn5_i5_set","tn5_i5_row","tn5_i7_set","tn5_i7_row"),by.y=c("tn5_i5_set","Row","tn5_i7_set","Column"))
  row.names(dat)<-dat$cellID
  obj<-AddMetaData(obj,dat)

  dim(obj@meta.data)

  #add in litter factor from updated samp_info

  saveRDS(obj,file="tbr1_ko.SeuratObject.Rds")
  write.table(obj@meta.data,file="summary_statistics_per_cell.tsv",col.names=T,row.names=T,sep="\t",quote=F)

  obj$line_genotype<-paste(obj$line,obj$genotype,sep="_")
  plt<-DimPlot(obj,group.by=c("peaks_snn_res.0.5","Developmental.Stage","sex","line_genotype","tube_ID"))
  ggsave(plt,file="tbr1_ko.umap.png",width=20)
  ggsave(plt,file="tbr1_ko.umap.pdf",width=20)
  system("slack -F tbr1_ko.umap.pdf ryan_todo")

  plt<-FeaturePlot(obj,features="litter_factor")
  ggsave(plt,file="tbr1_ko.umap.litter_factor.png")
  ggsave(plt,file="tbr1_ko.umap.litter_factor.pdf")
  system("slack -F tbr1_ko.umap.litter_factor.pdf ryan_todo")






```
### Make annotation files

This can be helpful with some other scitools functions.

```bash
tail -n +2 summary_statistics_per_cell.tsv| awk 'OFS="\t" {print $13,$21}' - > line.annot
tail -n +2 summary_statistics_per_cell.tsv| awk 'OFS="\t" {print $13,$21}' - > devel.annot
tail -n +2 summary_statistics_per_cell.tsv| awk 'OFS="\t" {print $13,$22}' - > line.annot
tail -n +2 summary_statistics_per_cell.tsv| awk 'OFS="\t" {print $13,$23}' - > geno.annot
tail -n +2 summary_statistics_per_cell.tsv| awk 'OFS="\t" {print $13,$21"_"$22"_"$23}' - > devel_line_geno.annot

```


{% endcapture %} {% include details.html %} 


### Statistics on cell reads

{% capture summary %} Code {% endcapture %} {% capture details %}  


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)
  library(dplyr)
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")
  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")


  #Add FRIP to meta data
  frip<-read.table("tbr1_ko.filt.500.fracOnTarget.values")
  colnames(frip)<-c("cellID","frip")
  obj$FRIP<-frip[match(obj$cellID,frip$cellID,),]$frip

  compl<-read.table("allplates.complexity.txt")
  colnames(compl)<-c("roworder","cellID","total_reads","unique_reads","percent_unique")
  obj$total_reads<-compl[match(obj$cellID,compl$cellID,),]$total_reads
  obj$unique_reads<-compl[match(obj$cellID,compl$cellID,),]$unique_reads
  obj$percent_unique<-compl[match(obj$cellID,compl$cellID,),]$percent_unique
  saveRDS(obj,"tbr1_ko.SeuratObject.Rds")

  #Cluster summaries
  dat<-obj@meta.data
  dat_sum<-as.data.frame(dat %>% 
  group_by(Developmental.Stage, line,genotype) %>% 
  summarize(mean_reads=mean(unique_reads),sd_reads=sd(unique_reads),median_reads=median(unique_reads),mean_FRIP=mean(FRIP),cell_count=n(),sample_count=length(unique(Sample_ID))))
  write.table(dat_sum,"tbr1ko_summary_statistics.tsv",col.names=T,row.names=T,quote=F,sep="\t")

  plt<-ggplot(dat,aes(x=paste(line,Developmental.Stage),y=FRIP,color=genotype))+geom_jitter(size=0.5)+geom_boxplot(outlier.shape=NA)+theme_bw()+xlab("Line, DevStage")
  ggsave(plt,file="frip_values.pdf",width=10)
  system("slack -F frip_values.pdf ryan_todo")


  plt<-ggplot(dat,aes(x=paste(line,Developmental.Stage),y=log10(unique_reads),color=genotype))+geom_jitter(size=0.5)+geom_boxplot(outlier.shape=NA)+theme_bw()+xlab("Line, DevStage")+ylab("Log10 Unique Reads")
  ggsave(plt,file="unique_reads.pdf",width=10)
  system("slack -F unique_reads.pdf ryan_todo")
```

{% endcapture %} {% include details.html %} 


## Cicero for Coaccessible Networks

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(cicero)
  library(EnsDb.Mmusculus.v79)
  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

  #Cicero processing function
  cicero_processing<-function(object_input=obj,prefix="tbr_ko"){

      #Generate CDS format from Seurat object
      atac.cds <- as.cell_data_set(object_input,group_by="seurat_clusters")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      atac.cicero<-readRDS(paste(prefix,"atac_cicero_cds.Rds",sep="_"))

      genome <- seqlengths(object_input) # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      
      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      
      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      return(object_input)
  }

  obj<-cicero_processing(object_input=obj,prefix="tbr1_ko")
  saveRDS(obj,"tbr1_ko.SeuratObject.GA.Rds")
  obj<-readRDS("tbr1_ko.SeuratObject.GA.Rds")
  
  # generate unnormalized gene activity matrix
  # gene annotation sample
  annotation_generation<-function(ensdb_obj){
      annotations <- GetGRangesFromEnsDb(ensdb = ensdb_obj)
      pos <-as.data.frame(annotations,row.names=NULL)
      pos$chromosome<-paste0("chr",pos$seqnames)
      pos$gene<-pos$gene_id
      pos <- subset(pos, strand == "+")
      pos <- pos[order(pos$start),] 
      pos <- pos[!duplicated(pos$tx_id),] # remove all but the first exons per transcript
      pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS
      neg <-as.data.frame(annotations,row.names=NULL)
      neg$chromosome<-paste0("chr",neg$seqnames)
      neg$gene<-neg$gene_id
      neg <- subset(neg, strand == "-")
      neg <- neg[order(neg$start,decreasing=TRUE),] 
      neg <- neg[!duplicated(neg$tx_id),] # remove all but the first exons per transcript
      neg$end <- neg$end + 1 # make a 1 base pair marker of the TSS
      gene_annotation<- rbind(pos, neg)
      gene_annotation <- gene_annotation[,c("chromosome","start","end","gene_name")] # Make a subset of the TSS annotation columns containing just the coordinates and the gene name
      names(gene_annotation)[4] <- "gene" # Rename the gene symbol column to "gene"
      return(gene_annotation)
    }

    mm10_annotation<-annotation_generation(ensdb_obj=EnsDb.Mmusculus.v79)

  geneactivity_processing<-function(cds_input,conns_input,prefix,gene_annotation){
      atac.cds<- annotate_cds_by_site(cds_input, gene_annotation)
      unnorm_ga <- build_gene_activity_matrix(atac.cds, conns_input)
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
  }

  #mm10
  conns<-as.data.frame(readRDS("tbr1_ko_atac_cicero_conns.Rds"))
  geneactivity_processing(cds_input=as.cell_data_set(obj,group_by="seurat_clusters"),conns_input=conns,prefix="tbr_ko",gene_annotation=mm10_annotation)

  cicero_gene_activities<-readRDS("tbr_ko.unnorm_GA.Rds")  #Read in unnormalized GA
  cicero_gene_activities<-cicero_gene_activities[2:nrow(cicero_gene_activities),] #first feature is empy
  obj[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
  obj <- NormalizeData(object = obj,assay = 'GeneActivity',normalization.method = 'LogNormalize',scale.factor = median(obj$nCount_peaks))  # normalize
  saveRDS(obj,"tbr1_ko.SeuratObject.Rds")

```

{% endcapture %} {% include details.html %} 

### Plots
```R
    library(Signac)
    library(Seurat)
    library(patchwork)
    set.seed(1234)
    library(dplyr)
    library(ggplot2)
    setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

    obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
      sep="", collapse=" ")
}

  markers<-c("Tbr1","Hes5","Olig1","Dbi","Neurod6","Gad1","Apoe","Mog","C1qb")
  markers<-unlist(lapply(markers,simpleCap)) # correct case

plt<-FeaturePlot(obj,features=markers,order=T,min.cutoff="q05",max.cutoff="q95",cols=c("lightgrey","red"))
ggsave(plt,file="tbr1_ko.markers.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.markers.test.pdf ryan_todo")

plt<-DimPlot(obj,group.by=c("peaks_snn_res.0.5","seurat_clusters","Sample_ID","Developmental.Stage","line_genotype"),order=T)
ggsave(plt,file="tbr1_ko.metadata.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.metadata.test.pdf ryan_todo")


plt<-DimPlot(obj,group.by="genotype",split.by="line",order=T)
ggsave(plt,file="tbr1_ko.metadata.line.test.pdf")
system("slack -F tbr1_ko.metadata.line.test.pdf ryan_todo")

plt<-CoveragePlot(obj,region=markers,assay="peaks",group.by="peaks_snn_res.0.5",extend.upstream=2000,extend.downstream=2000,ncol=1)
ggsave(plt,file="tbr1_ko.markers.test.cov.pdf",height=5*length(markers),limitsize=F)
system("slack -F tbr1_ko.markers.test.cov.pdf ryan_todo")

plt<-VlnPlot(obj,features=markers)
ggsave(plt,file="tbr1_ko.markers.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.markers.test.pdf ryan_todo")

dat<-obj@meta.data %>% group_by(Sample_ID,sex,Developmental.Stage,line,genotype,peaks_snn_res.0.5) %>% summarize(count=n()) 

plt<-ggplot(dat,aes(y=count,fill=peaks_snn_res.0.5,x=Sample_ID))+geom_bar(position="fill",stat="identity")+facet_wrap(Developmental.Stage+line~genotype,scales="free",strip.position="left")
ggsave(plt,file="stacked_bar_plot.png",width=10,height=10)
system("slack -F stacked_bar_plot.png ryan_todo")
```

### Add TF Motif Usage through ChromVAR

### ChromVar for Transcription Factor Motifs

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
    library(Signac)
    library(Seurat)
    library(JASPAR2020)
    library(TFBSTools)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(patchwork)
    set.seed(1234)

    #lowerign cores to be used by chromvar to 10
    library(BiocParallel)
    register(MulticoreParam(5))

    setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

    obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

    #Read in data and modify to monocle CDS file
    #read in RDS file.

    # Get a list of motif position frequency matrices from the JASPAR database
    pfm <- getMatrixSet(
      x = JASPAR2020,
      opts = list(species =9606, all_versions = FALSE))

    # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
    motif.matrix.mm10 <- CreateMotifMatrix(
      features = granges(obj[["peaks"]]),
      pwm = pfm,
      genome = 'mm10',
      use.counts = FALSE)

    # Create a new Mofif object to store the results
    motif.mm10 <- CreateMotifObject(
      data = motif.matrix.mm10,
      pwm = pfm)

    # Add the Motif object to the assays and run ChromVar
    obj <- SetAssayData(
      object = obj,
      assay = 'peaks',
      slot = 'motifs',
      new.data = motif.mm10)

    obj <- RegionStats(object = obj, genome = BSgenome.Mmusculus.UCSC.mm10)
    obj <- RunChromVAR( object = obj,genome = BSgenome.Mmusculus.UCSC.mm10)
    saveRDS(obj,"tbr1_ko.SeuratObject.Rds")

```
{% endcapture %} {% include details.html %} 


### Differential Motif Accessibility

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  ###Differential TF Accessibility by cluster###

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(parallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)

  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")
  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

  #Perform One vs. rest DA enrichment
  write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

  #define DA functions for parallelization
  #Use LR test for atac data
  da_one_v_rest<-function(i,obj,group,latent.vars.="nCount_peaks"){
      da_peaks_tmp <- FindMarkers(
          object = obj,
          ident.1 = i,
          group.by = group,
          test.use = 'LR',
          latent.vars = latent.vars.,
          only.pos=T,
          logfc.threshold=0.05
          )
      da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
      da_peaks_tmp$enriched_group<-c(i)
      da_peaks_tmp$compared_group<-c("all_other_cells")
      return(da_peaks_tmp)
    }

  # da_one_v_one<-function(i,obj,group,j_list){
  #     i<-as.character(i)
  #     da_tmp_2<-list()
  #     for (j in j_list){
  #         if ( i != j){
  #         da_peaks_tmp <- FindMarkers(
  #             object = obj,
  #             ident.1 = i,
  #             ident.2 = j,
  #             group.by = group,
  #             test.use = 'LR',
  #             latent.vars = 'nCount_peaks',
  #             only.pos=T
  #             )
  #         da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
  #         da_peaks_tmp$enriched_group<-c(i)
  #         da_peaks_tmp$compared_group<-c(j)
  #         da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
  #         }
  #     }
  #     return(da_tmp_2)
  #   }

  #Perform parallel application of DA test
  n.cores=length(unique(obj$seurat_clusters))
  #set up an empty list for looping through
  DefaultAssay(obj)<-"peaks"
  obj_peaks<-list()
  obj_peaks<-mclapply(
      unique(obj@meta.data$peaks_snn_res.0.5),
      FUN=da_one_v_rest,
      obj=obj,
      group="peaks_snn_res.0.5",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1vrest DA
  obj_peaks<-do.call("rbind",obj_peaks)
  closest_genes <- ClosestFeature(obj,obj_peaks$da_region)
  obj_peaks<-cbind(obj_peaks,closest_genes)
  write.table(obj_peaks,file="tbr1_ko.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

  #set up an empty list for looping through
  DefaultAssay(obj)<-"chromvar"
  obj_chromvar<-list()
  obj_chromvar<-mclapply(
      unique(obj@meta.data$peaks_snn_res.0.5),
      FUN=da_one_v_rest,
      obj=obj,
      group="peaks_snn_res.0.5",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1vrest DA
  obj_chromvar<-do.call("rbind",obj_chromvar)
  obj_chromvar$tf_name <- unlist(lapply(unlist(lapply(obj_chromvar$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(obj_chromvar,file="tbr1_ko.onevrest.da_chromvar.txt",sep="\t",col.names=T,row.names=T,quote=F)

  DefaultAssay(obj)<-"GeneActivity"
  obj_ga<-list()
  obj_ga<-mclapply(
      unique(obj@meta.data$peaks_snn_res.0.5),
      FUN=da_one_v_rest,
      obj=obj,
      group="peaks_snn_res.0.5",
      latent.vars.="ncount_GeneActivity",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1v1 DA
  obj_ga<-do.call("rbind",do.call("rbind",obj_ga))
  write.table(obj_ga,file="tbr1_ko.onevrest.da_ga.txt",sep="\t",col.names=T,row.names=T,quote=F)


  #Now plotting
  ########################Plot out top TF for each cluster###################
  da_tf<-read.csv(file="tbr1_ko.onevrest.da_chromvar.txt",head=T,sep="\t",row.names=NULL)

  da_tf$label<-""
  for (x in unique(da_tf$enriched_group)){
  selc_genes<-as.data.frame(da_tf %>% filter(enriched_group==x)  %>% filter(p_val < 0.05) %>% arrange(rev(desc(p_val))) %>% slice(1:10))$tf_name
  da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$label<- da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$tf_name
  }


  plt<-ggplot(da_tf,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+
  geom_point(aes(alpha=0.1))+
  geom_label_repel(aes(label=label),size=2,force=5)+
  theme_bw()+ylim(c(0,9))+xlim(c(0,30))
  ggsave(plt,file="tbr1_ko.onevrest.da_chromvar.plt.pdf")
  system(paste0("slack -F ","tbr1_ko.onevrest.da_chromvar.plt.pdf"," ryan_todo"))

    ########################Plot out top DA Peaks for each cluster###################
  da_peaks<-read.csv(file="tbr1_ko.onevrest.da_peaks.txt",head=T,sep="\t",row.names=NULL)

  da_peaks$label<-""
  for (x in unique(da_peaks$enriched_group)){
  selc_genes<-as.data.frame(da_peaks %>% filter(enriched_group==x)  %>% filter(p_val < 0.05) %>% arrange(rev(desc(p_val))) %>% slice(1:10))$gene_name
  da_peaks[da_peaks$gene_name %in% selc_genes & da_peaks$enriched_group==x,]$label<- da_peaks[da_peaks$gene_name %in% selc_genes & da_peaks$enriched_group==x,]$gene_name
  }


  plt<-ggplot(da_peaks,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+
  geom_point(aes(alpha=0.1))+
  geom_label_repel(aes(label=label),size=2,force=5)+
  theme_bw()+ylim(c(0,7))+xlim(c(0,20))
  ggsave(plt,file="tbr1_ko.onevrest.da_peaks.plt.pdf")
  system(paste0("slack -F ","tbr1_ko.onevrest.da_peaks.plt.pdf"," ryan_todo"))


```
{% endcapture %} {% include details.html %} 

### Performing GREAT on DA peaks

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  #mkdir GREAT_analysis

  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)

  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

  #To perform GREAT on peaks for enrichment per cluster
  write("Performing GREAT on all enriched sites per annotation group", stderr())
  library(rGREAT)

  #format data as bed file all seurat objects have the same peak list
  write("Preparing Background Set as all called peaks.", stderr())
  obj_bg_bed<-do.call("rbind",strsplit(unlist(obj@assays$peaks@counts@Dimnames[1]),"[-]"))
  obj_bg_bed<-as.data.frame(obj_bg_bed)
  colnames(obj_bg_bed)<-c("chr","start","end")
  obj_bg_bed$start<-as.numeric(as.character(obj_bg_bed$start))
  obj_bg_bed$end<-as.numeric(as.character(obj_bg_bed$end))

  obj_da_peaks<-read.table("tbr1_ko.onevrest.da_peaks.txt",header=T)

  write("Beginning loop through all annotation groups.", stderr())
  dir.create("GREAT_analysis")
  great_processing<-function(enriched_group_input,peak_dataframe,prefix){
      #subset bed file to peaks enriched in input group
      obj_bed<-as.data.frame(do.call("rbind",strsplit(obj_da_peaks[obj_da_peaks$enriched_group==enriched_group_input,]$da_region,"-")))
      colnames(obj_bed)<-c("chr","start","end")
      obj_bed$start<-as.numeric(as.character(obj_bed$start))
      obj_bed$end<-as.numeric(as.character(obj_bed$end))
      
      #run GREAT using all peaks as background
      write(paste("Using",nrow(obj_bed), "DA peaks from",enriched_group_input), stderr())
      job = submitGreatJob(obj_bed,obj_bg_bed,species="mm10",request_interval=30)
      tb = getEnrichmentTables(job, ontology = c("GO Molecular Function", "GO Biological Process","GO Cellular Component"))
      tb = getEnrichmentTables(job, category = c("GO","Phenotype","Genes"))
      #Plot gene association
      pdf(paste0("./GREAT_analysis/",prefix,"_DApeaks_",enriched_group_input,".GeneAssociation.pdf"))
      plotRegionGeneAssociationGraphs(job)
      dev.off()

      for (j in 1:length(names(tb))){
            write(paste("Outputting DA GREAT Analysis for", enriched_group_input, as.character(names(tb))[j]), stderr())
            tabl_name<-gsub(" ","",as.character(names(tb))[j])
            write.table(as.data.frame(tb[[j]]),file=paste0("./GREAT_analysis/",prefix,"_DApeaks_",enriched_group_input,".",tabl_name,".txt"),sep="\t",col.names=T,row.names=T,quote=F)
        }
  }

  library(parallel)
  mclapply(unique(obj_da_peaks$enriched_group), FUN=great_processing, peak_dataframe=obj_da_peaks,prefix="tbr1_ko",mc.cores=10)
```
{% endcapture %} {% include details.html %} 

```R

  setwd("/home/groups/oroaklab/adey_lab/projects/tbr1_mus/210225_allplates")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Mmusculus.v79)
  library(Matrix)

  obj<-readRDS(file="tbr1_ko.SeuratObject.Rds")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
      sep="", collapse=" ")
}

markers<-c("Rit1","Tbr1","Cux1","Neurod6","Gad2","Eomes","Mki67")
markers<-unlist(lapply(markers,simpleCap)) # correct case

plt<-FeaturePlot(obj,features=markers,order=T,min.cutoff="q05",max.cutoff="q95",cols=c("lightgrey","red"))
ggsave(plt,file="tbr1_ko.markers.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.markers.test.pdf ryan_todo")

plt<-VlnPlot(obj,features=markers,group.by="seurat_clusters")
ggsave(plt,file="tbr1_ko.markers.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.markers.test.pdf ryan_todo")

row.names(obj@assays$chromvar@data)<-unlist(lapply(unlist(lapply(row.names(obj@assays$chromvar@data), function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))

plt<-FeaturePlot(obj,features=c("TBR1","CUX1","NEUROD6","GAD2","EOMES"),order=T,min.cutoff="q01",max.cutoff="q99",cols=c("lightgrey","blue"))
ggsave(plt,file="tbr1_ko.markers.test.pdf",width=20,height=20)
system("slack -F tbr1_ko.markers.test.pdf ryan_todo")

```