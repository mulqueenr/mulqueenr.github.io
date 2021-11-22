---
title: Organoids
layout: analysis
author: Ryan Mulqueen
permalink: /organoid/
category: sciATAC
---

## Processing for sciATAC portion for organoid analysis.
I ran multiple sequecing runs for the sciATAC. For now I am just processing the most recent, but I will loop back to the original Pitstop2 experiments.

### BCL File Locations

```bash
  #First Prep
  /home/groups/oroaklab/seq/hiseq/180630_AML_Pitstop

  #Second prep
  /home/groups/oroaklab/seq/madbum/200721_NS500556_0411_AHCM3CAFX2
  /home/groups/oroaklab/seq/madbum/200804_NS500556_0413_AHCMMJBGXF

  #RNA prep
  /home/groups/oroaklab/seq/madbum/191113_NS500556_0361_AHTVFLAFXY
  /home/groups/oroaklab/seq/madbum/191118_NS500556_0362_AHVYV7AFXY
  /home/groups/oroaklab/seq/madbum/191119_NS500556_0363_AHTVL7AFXY
```

### Initial Processing of Files
Includes barcode assignment, fastq splitting, alignment, removal of duplicate reads, calling peaks and looking at TSS enrichment.


```bash
  #200722 Organoid Processing
  NextSeq2fastq -R 200721_NS500556_0411_AHCM3CAFX2
  NextSeq2fastq -R 200804_NS500556_0413_AHCMMJBGXF

  outdir="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis"

  mkdir $outdir

  scitools fastq-dump -R 200721_NS500556_0411_AHCM3CAFX2 -O $outdir
  scitools fastq-dump -R 200804_NS500556_0413_AHCMMJBGXF -O $outdir

  #quick annotation generation for initial fastq splitting from other libraries sequenced on the same run
  scitools make-annot orgo+NEX,CB=ALL+NEX,BC=ALL+NEX,CA=ALL+NEX,BB=ALL+NEX,AA=ALL+NEX,BA=ALL\
  +PCR,CB=ALL+PCR,CA=ALL+PCR,CC=ALL+PCR,CD=ALL,+PCR,CE=ALL,+PCR,CF=ALL,+PCR,AD=ALL,+PCR,AE=ALL,+PCR,AF=ALL,+PCR,AG=ALL > base_orgo.annot                 

  scitools split-fastq -X -A base_orgo.annot 200721_NS500556_0411_AHCM3CAFX2.1.fq.gz 200721_NS500556_0411_AHCM3CAFX2.2.fq.gz
  scitools split-fastq -X -A base_orgo.annot 200804_NS500556_0413_AHCMMJBGXF.1.fq.gz 200804_NS500556_0413_AHCMMJBGXF.2.fq.gz

  #move all relevant files to /home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis
      
  #Align with a wrapper for bwa mem ###RUNNING
  scitools fastq-align -t 10 -r 10 hg38 orgo_prep2_1 200721_NS500556_0411_AHCM3CAFX2.1.fq.gz 200721_NS500556_0411_AHCM3CAFX2.2.fq.gz &
  scitools fastq-align -t 10 -r 10 hg38 orgo_prep2_2 200804_NS500556_0413_AHCMMJBGXF.orgo.1.fq.gz 200804_NS500556_0413_AHCMMJBGXF.orgo.2.fq.gz 
  scitools fastq-align -t 20 -r 20 hg38 orgo_prep1_1 180630.RM.1.fq.gz 180630.RM.1.fq.gz & 

  #Barcode based remove duplicates, barcodes are contained in the read name
  scitools bam-rmdup -t 10 orgo_prep2_1.bam &
  scitools bam-rmdup -t 10 orgo_prep2_2.bam &
  scitools bam-rmdup -t 10 orgo_prep1_1.bam &

  #Filter based on barcodes with >1000 unique reads
  for i in orgo_prep1_1.bbrd.q10.bam orgo_prep2_1.bbrd.q10.bam orgo_prep2_2.bbrd.q10.bam; do scitools bam-filter -N 1000 $i ; done &

  #merge bam files
  scitools bam-merge orgo.bam orgo_prep1_1.bbrd.q10.filt.bam orgo_prep2_1.bbrd.q10.filt.bam orgo_prep2_2.bbrd.q10.filt.bam

  #Look at tss enrichment to ensure libraries are of good quality
  module load bedops/2.4.36
  scitools bam-tssenrich -X -E orgo.ID.bam hg38 & #bulk ENCODE method
  scitools bam-tssenrich -X orgo.ID.bam hg38 & #single cell method

  #Did a fresh install of macs2 for py3 environment
  #pip install macs2 #for python3 macs2
  scitools callpeaks orgo.bam &
  #Count of previously annotated human fetal corticogenesis peaks
   bedtools intersect -a orgo.500.bed -b /home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/Song_2020/celltype_peaks.bed -wa -u | wc -l

  #Modifying bam file to include prep number in cellID field (to maintain single cell identity through index collisions)
  ((samtools view -H orgo.bam)&(samtools view orgo.bam |awk 'OFS="\t" {split($1,a,":");split(a[3],b,"="); $1=a[1]"_"b[2]":"a[2]":"a[3]; print $0}')) | samtools view -bS - > orgo.ID.bam &

  #Generating sparse matrix format counts matrix
  scitools atac-counts orgo.ID.bam orgo.500.bed &
```


### Generation of thorough annotation file and all meta data per cell 


```R
  #generated organoid annotation file from google sheet https://docs.google.com/spreadsheets/d/1k93smqwxmYVUMLqq9SG8UjgjFeTinERzOflWUMRN8n8/edit#gid=1394545516

  #wrote out with nano into tsv format

  R #Using R 4.0
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(reshape2)
  ####read in files
  first_prep_annot_path="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/first_prep_annot"
  second_prep_annot_path="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/second_prep_annot"
  first_prep_tn5_plates<-list.files(path=first_prep_annot_path,pattern="^tn5_plate_...tsv")
  second_prep_tn5_plates<-list.files(path=second_prep_annot_path,pattern="^tn5_plate_...tsv")

  idx<-read.table("/home/groups/oroaklab/src/scitools/scitools-dev/SCI_Indexes.txt",col.names=c("idx_name","idx_cycles","idx_seq"))

  first_prep_orgid<-read.table(paste0(first_prep_annot_path,"/","tn5_plate_orgid.tsv"),header=T,sep="\t",col.names=c("orgID","differentiation_exp","cell_line","DIV","tag_wells","perc_tag_wells","est_cells_per_pcr_plate","cells_counted","sort_gate","freezing_protocol","treatment","organoid"))
  second_prep_orgid<-read.table(paste0(second_prep_annot_path,"/","tn5_plate_orgid.tsv"),header=T,sep="\t",col.names=c("orgID","differentiation_exp","cell_line","DIV","tag_wells","perc_tag_wells","est_cells_per_pcr_plate","est_cells_per5_pcr_plate","est_cells_per10_pcr_plate"))

  compl_1<-read.table("./source_fastq/preprocessing_files/orgo_prep1_1.complexity.txt",row.names=1,col.names=c("row.names","cellID","total_reads","uniq_reads","perc_uniq"))
  compl_2<-read.table("./source_fastq/preprocessing_files/orgo_prep2_1.complexity.txt",row.names=1,col.names=c("row.names","cellID","total_reads","uniq_reads","perc_uniq"))
  compl_3<-read.table("./source_fastq/preprocessing_files/orgo_prep2_2.complexity.txt",row.names=1,col.names=c("row.names","cellID","total_reads","uniq_reads","perc_uniq"))

  #split index names to columns
  idx$idx_plateloc<-unlist(lapply(as.character(idx$idx_name),FUN=function(x) strsplit(x,"_")[[1]][4]))
  idx$idx_moleculesrc<-unlist(lapply(as.character(idx$idx_name),FUN=function(x) strsplit(x,"_")[[1]][1]))
  idx$idx_platename<-unlist(lapply(as.character(idx$idx_name),FUN=function(x) strsplit(x,"_")[[1]][2]))
  idx$idx_moleculeloc<-unlist(lapply(as.character(idx$idx_name),FUN=function(x) strsplit(x,"_")[[1]][3]))

  #subset index sequences by cycle
  idx_tn5i5<-idx[idx$idx_moleculeloc=="i5" & idx$idx_moleculesrc=="Tn5",]
  idx_tn5i7<-idx[idx$idx_moleculeloc=="i7" & idx$idx_moleculesrc=="Tn5",]
  colnames(idx_tn5i5)<-c("tn5_i5_idx_name","idx_cycles","tn5_i5_idx_seq","row","idx_moleculesrc","tn5_i5","idx_moleculeloc")
  colnames(idx_tn5i7)<-c("tn5_i7_idx_name","idx_cycles","tn5_i7_idx_seq","column","idx_moleculesrc","tn5_i7","idx_moleculeloc")
  idx_tn5i5<-idx_tn5i5[,c(1,3,4,6)]
  idx_tn5i7<-idx_tn5i7[,c(1,3,4,6)]

  #generate cellID list and split cell id by cycle
  compl_1$cellID<-paste0(compl_1$cellID,"_1")
  compl_2$cellID<-paste0(compl_2$cellID,"_2")
  compl_3$cellID<-paste0(compl_3$cellID,"_3")
  cellid<-rbind(compl_1,compl_2,compl_3)
  cellid<-as.data.frame(cellid$cellID)
  colnames(cellid)<-"cellID"
  cellid$tn5_i7_idx_seq<-substr(cellid$cellID,1,8)
  cellid$pcr_i7_idx_seq<-substr(cellid$cellID,9,18)
  cellid$tn5_i5_idx_seq<-substr(cellid$cellID,19,26)
  cellid$pcr_i5_idx_seq<-substr(cellid$cellID,27,36)
  cellid$prep<-substr(cellid$cellID,38,38)

  #generate long format organoid tn5 data (for second plate)
  dat<-data.frame()
  for (i in second_prep_tn5_plates) { 
  tn5_plate<-strsplit(strsplit(i,"_")[[1]][-1],"[.]")[[2]][1]
  tmp<-read.table(paste0(second_prep_annot_path,"/",i),header=T,row.names=1)
  tmp$row<-row.names(tmp)
  tmp<-melt(tmp)
  tmp$tn5_plate<-tn5_plate
  ifelse(nrow(dat)==0,dat<-tmp,dat<-rbind(dat,tmp))
  }
  colnames(dat)<-c("row","column","orgID","tn5_plate")
  dat$tn5_plate<-toupper(dat$tn5_plate)
  dat$tn5_i5<-substr(dat$tn5_plate,1,1)
  dat$tn5_i7<-substr(dat$tn5_plate,2,2)
  dat$column<-substr(dat$column,2,5)

  #merge organoid tn5 data with tn5 i5 index
  dat<-merge(dat,idx_tn5i5,by=c("tn5_i5","row"),all.x=T)
  #merge organoid tn5 data with tn5 i7 index
  dat<-merge(dat,idx_tn5i7,by=c("tn5_i7","column"),all.x=T)
  #merge organoid data with organoid key
  dat<-merge(dat,second_prep_orgid,by="orgID")
  #merge organoid data with complexity file
  compl<-rbind(compl_2,compl_3)
  compl<-merge(compl,cellid[cellid$prep!=1,],by="cellID")

  dat<-merge(dat,compl,by=c("tn5_i5_idx_seq","tn5_i7_idx_seq"))
  #apply read filter used on each bam
  dat<-dat[dat$uniq_reads>=1000,]
  write.table(dat,file="second_prep_summary_statistics_per_cell.tsv",col.names=T,row.names=F,quote=F,sep="\t")

  tn5plate_annot<-dat[,c("cellID","tn5_plate")]
  orgID_annot<-dat[,c("cellID","orgID")]
  diffexp_annot<-dat[,c("cellID","differentiation_exp")]
  cellline_annot<-dat[,c("cellID","cell_line")]
  div_annot<-dat[,c("cellID","DIV")]

  write.table(tn5plate_annot,file=paste0(second_prep_annot_path,"/","second_prep_tn5plate.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(orgID_annot,file=paste0(second_prep_annot_path,"/","second_prep_orgID.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(diffexp_annot,file=paste0(second_prep_annot_path,"/","second_prep_diffexp.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(cellline_annot,file=paste0(second_prep_annot_path,"/","second_prep_cellline.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(div_annot,file=paste0(second_prep_annot_path,"/","second_prep_div.annot"),sep="\t",col.names=F,row.names=F,quote=F)


  #generate long format organoid tn5 data (for first plate)
  #merge organoid data with complexity file
  compl<-compl_1
  compl<-merge(compl,cellid[cellid$prep==1,],by="cellID")

  dat<-data.frame()
  for (i in first_prep_tn5_plates) { 
  tn5_plate<-strsplit(strsplit(i,"_")[[1]][-1],"[.]")[[2]][1]
  tmp<-read.table(paste0(first_prep_annot_path,"/",i),header=T,row.names=1)
  colnames(tmp)<-paste0("X",seq(1,12))
  tmp$row<-row.names(tmp)
  tmp<-melt(tmp,id.vars="row")
  tmp$tn5_plate<-tn5_plate
  ifelse(nrow(dat)==0,dat<-tmp,dat<-rbind(dat,tmp))
  }
  colnames(dat)<-c("row","column","orgID","tn5_plate")
  dat$tn5_plate<-toupper(dat$tn5_plate)
  dat$tn5_i5<-substr(dat$tn5_plate,1,1)
  dat$tn5_i7<-substr(dat$tn5_plate,2,2)
  dat$column<-substr(dat$column,2,5)

  #merge organoid tn5 data with tn5 i5 index
  dat<-merge(dat,idx_tn5i5,by=c("tn5_i5","row"),all.x=T)
  #merge organoid tn5 data with tn5 i7 index
  dat<-merge(dat,idx_tn5i7,by=c("tn5_i7","column"),all.x=T)
  #merge organoid data with organoid key
  dat<-merge(dat,first_prep_orgid,by="orgID")
  dat<-merge(dat,compl,by=c("tn5_i5_idx_seq","tn5_i7_idx_seq"))
  #apply read filter used on each bam
  dat<-dat[dat$uniq_reads>=1000,]

  write.table(dat,file="first_prep_summary_statistics_per_cell.tsv",col.names=T,row.names=F,quote=F,sep="\t")

  tn5plate_annot<-dat[,c("cellID","tn5_plate")]
  orgID_annot<-dat[,c("cellID","orgID")]
  diffexp_annot<-dat[,c("cellID","differentiation_exp")]
  cellline_annot<-dat[,c("cellID","cell_line")]
  div_annot<-dat[,c("cellID","DIV")]
  sortgate_annot<-dat[,c("cellID","sort_gate")]
  freezing_annot<-dat[,c("cellID","freezing_protocol")]
  treatment_annot<-dat[,c("cellID","treatment")]
  organoid_annot<-dat[,c("cellID","organoid")]

  write.table(tn5plate_annot,file=paste0(first_prep_annot_path,"/","first_prep_tn5plate.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(orgID_annot,file=paste0(first_prep_annot_path,"/","first_prep_orgID.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(diffexp_annot,file=paste0(first_prep_annot_path,"/","first_prep_diffexp.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(cellline_annot,file=paste0(first_prep_annot_path,"/","first_prep_cellline.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(div_annot,file=paste0(first_prep_annot_path,"/","first_prep_div.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(sortgate_annot,file=paste0(first_prep_annot_path,"/","first_prep_sortgate.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(freezing_annot,file=paste0(first_prep_annot_path,"/","first_prep_freezeprotocol.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(treatment_annot,file=paste0(first_prep_annot_path,"/","first_prep_treatment.annot"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(organoid_annot,file=paste0(first_prep_annot_path,"/","first_prep_organoid.annot"),sep="\t",col.names=F,row.names=F,quote=F)
```

### Tabix fragment file generation

Tabix file format is a tab separated multicolumn data structure.

| Column Number | Name | Description |
|:--------|:-------:|:--------|
|1 |chrom |  Reference genome chromosome of fragment |
|2 |chromStart | Adjusted start position of fragment on chromosome. |
|3 |chromEnd   | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
|4 |barcode | The 10x (or sci) cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment. |
|5 |duplicateCount |The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeq instrument are excluded from this count. |


```bash
  input_bam="orgo.ID.bam"
  output_name="orgo"
  tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
  bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
  samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
  $tabix -p bed $output_name.fragments.tsv.gz &
```


## sciATAC Generalized Processing in R

### Generating Seurat Objects

Using R v4.0 and Signac v1.0 for processing.


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")


  read_in_sparse<-function(x){ #x is character file prefix followed by .bbrd.q10.500.counts.sparseMatrix.values.gz
  IN<-as.matrix(read.table(paste0(x,".counts.sparseMatrix.values.gz")))
  IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
  COLS<-read.table(paste0(x,".counts.sparseMatrix.cols.gz"))
  colnames(IN)<-COLS$V1
  ROWS<-read.table(paste0(x,".counts.sparseMatrix.rows.gz"))
  row.names(IN)<-ROWS$V1
  writeMM(IN,file=paste0(x,".counts.mtx")) #this is to generate counts matrices in scrublet friendly format
  return(IN)
  }

  orgo_counts<-read_in_sparse("orgo.500") # make counts matrix from sparse matrix

  #Read in fragment path for coverage plots
  orgo_fragment.path="./orgo.fragments.tsv.gz"

  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

  # change to UCSC style since the data was mapped to hg38
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

  #Generate ChromatinAssay Objects
  orgo_chromatinassay <- CreateChromatinAssay(
    counts = orgo_counts,
    genome="hg38",
    min.cells = 1,
    annotation=annotations,
    sep=c("_","_"),
    fragments=orgo_fragment.path
  )

  #Create Seurat Object
  orgo_atac <- CreateSeuratObject(
    counts = orgo_chromatinassay,
    assay = "peaks",
  )

  #Meta.data to be updated after clustering
  #saving unprocessed SeuratObject
  saveRDS(orgo_atac,file="orgo_SeuratObject.Rds")
```



### Perform Scrublet on Data to Ensure Single-cells

Code from tutorial here.[https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb]


```bash
  #using a conda environment set up by ARSN
  source /home/groups/oroaklab/nishida/scitools_env/bin/activate
```

```python
#Installing scrublet
#pip install scrublet
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.sparse import coo_matrix
import gzip
import pandas as pd

#Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.
input_dir = '/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/'

#Perform on hg38 cells
counts_matrix = scipy.io.mmread(input_dir + 'orgo.500.counts.mtx').T.tocsc() #generated during the initialization of the Seurat Object

peaks= np.array(gzip.open(input_dir+'orgo.500.counts.sparseMatrix.rows.gz', 'rt').read().split()) #This is read in to check that our data frame is in the correct orientation
cellid= gzip.open(input_dir+'orgo.500.counts.sparseMatrix.cols.gz', 'rt').read().split() #This is read in to check that our data frame is in the correct orientation
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(peaks)))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
#Preprocessing...
#Simulating doublets...
#Embedding transcriptomes using PCA...
#Calculating doublet scores...
#Automatically set threshold at doublet score = 0.07
#Detected doublet rate = 31.1%
#Estimated detectable doublet fraction = 60.5%
#Overall doublet rate:
#        Expected   = 5.0%
#        Estimated  = 51.4%
#Elapsed time: 785.7 seconds

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

df = pd.DataFrame({'cellid':cellid, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv('orgo.scrublet.tsv', index=False, sep="\t")
```


### Plotting and updating metadata


```R
  #renaming annot for simplified annotation file making
  #rename processing_ processing. *annot
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  #Set up annotation summaries to contain the same information, in same column order
  first_annot_append<-read.table("first_prep_summary_statistics_per_cell.tsv",header=T) 
  #reorder for consistent metadata
  first_annot_append<-first_annot_append[c("cellID",
                                           "tn5_plate",
                                           "column","row",
                                           "tn5_i5","tn5_i5_idx_name","tn5_i5_idx_seq",
                                           "tn5_i7","tn5_i7_idx_name","tn5_i7_idx_seq",
                                           "pcr_i7_idx_seq","pcr_i5_idx_seq",
                                           "total_reads","uniq_reads","perc_uniq",
                                           "prep","orgID","cell_line","differentiation_exp","DIV",
                                           "freezing_protocol","sort_gate","treatment","organoid")]

  second_annot_append<-read.table("second_prep_summary_statistics_per_cell.tsv",header=T)
  second_annot_append$freezing_protocol<-"Flash_Frozen" #change this for the DIV90 cirm 43 diff exp 5 organoids
  second_annot_append[(second_annot_append$DIV=="90" & second_annot_append$differentiation_exp=="5"),]$freezing_protocol<-"Slow_Freeze"
  second_annot_append$sort_gate<-"NA"
  second_annot_append$treatment<-"No"
  second_annot_append$organoid<-second_annot_append$orgID

  second_annot_append<-second_annot_append[c("cellID",
                                           "tn5_plate",
                                           "column","row",
                                           "tn5_i5","tn5_i5_idx_name","tn5_i5_idx_seq",
                                           "tn5_i7","tn5_i7_idx_name","tn5_i7_idx_seq",
                                           "pcr_i7_idx_seq","pcr_i5_idx_seq",
                                           "total_reads","uniq_reads","perc_uniq",
                                           "prep","orgID","cell_line","differentiation_exp","DIV",
                                           "freezing_protocol","sort_gate","treatment","organoid")]

  annot_append<-rbind(first_annot_append,second_annot_append)
  #orgID and prep need to be accounted for to get unique organoids (there are duplicates in orgID)

  #Add original cluster information
  original_cluster<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/180703_sciATAC_Organoids/NatureLetterData/cellID_fullcharacterization.txt",header=T)
  original_cluster$cellID<-paste0(original_cluster$cellID,"_1")
  original_cluster<-original_cluster[c("cellID","Phenograph_Cluster")]
  colnames(original_cluster)<-c("cellID","original_cluster")

  orgo_atac<-readRDS(file="orgo_SeuratObject.Rds")
  orgo_atac@meta.data$cellID<-row.names(orgo_atac@meta.data)

  #Add scrublet info
  orgo_scrub<-read.table("orgo.scrublet.tsv",header=T) #read in scrublet

  #Add complexity info
  compl_1<-read.table("source_fastq/preprocessing_files/orgo_prep1_1.complexity.txt",head=F)
  colnames(compl_1)<-c("cellID","total_reads","unique_reads","percent_unique_reads")
  compl_1$cellID<-paste0(compl_1$cellID,"_1")
  compl_2<-read.table("source_fastq/preprocessing_files/orgo_prep2_1.complexity.txt",head=F)
  colnames(compl_2)<-c("cellID","total_reads","unique_reads","percent_unique_reads")
  compl_2$cellID<-paste0(compl_2$cellID,"_2")
  compl_3<-read.table("source_fastq/preprocessing_files/orgo_prep2_2.complexity.txt",head=F)
  colnames(compl_3)<-c("cellID","total_reads","unique_reads","percent_unique_reads")
  compl_3$cellID<-paste0(compl_3$cellID,"_3")
  compl<-rbind(compl_1,compl_2,compl_3)

  #Add TSS enrichment value
  tss_enrich<-read.table("orgo.ID.TSSenrich.value",header=F)
  colnames(tss_enrich)<-c("cellID","tss_enrichment")

  #merge all data frames
  annot<-as.data.frame(orgo_atac@meta.data)
  annot<-merge(annot,annot_append,by="cellID",all.x=T)
  annot<-merge(annot,original_cluster,by="cellID",all.x=T)
  annot<-merge(annot,orgo_scrub,by.x="cellID",by.y="cellid",all.x=T)
  #annot<-merge(annot,compl,by="cellID",all.x=T)
  annot<-merge(annot,tss_enrich,by="cellID",all.x=T)
  row.names(annot)<-annot$cellID
  orgo_atac@meta.data<-annot

  #Add FRIP to meta data
  frip<-read.table("orgo.500.fracOnTarget.values")
  colnames(frip)<-c("cellID","frip")
  row.names(frip)<-frip$cellID
  frip<-frip[frip$cellID %in% row.names(orgo_atac@meta.data),]
  frip_names<-setNames(frip$frip,nm=frip$cellID)
  orgo_atac<-AddMetaData(object=orgo_atac,col.name="FRIP",metadata=frip_names)
  #excluding differentiation experiment 4
  orgo_atac<-subset(orgo_atac, differentiation_exp %in% c("5","7"))
  orgo_atac<-subset(orgo_atac,cell_line=="CIRM43") #just cirm43 cell line and two differentiations


  saveRDS(orgo_atac,file="orgo_cirm43.SeuratObject.Rds")
```


### ChromVar for Transcription Factor Motifs


```R
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  set.seed(1234)

  #lowering cores to be used by chromvar to 10
  library(BiocParallel)
  register(MulticoreParam(10))
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  #Read in data and modify to monocle CDS file
  #read in RDS file.
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(x = JASPAR2020, opts = list(species =9606, all_versions = FALSE))

  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  motif.matrix <- CreateMotifMatrix(features = granges(orgo_cirm43), pwm = pfm, genome = 'hg38', use.counts = FALSE)

  # Create a new Mofif object to store the results
  motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)

  # Add the Motif object to the assays and run ChromVar
  orgo_cirm43 <- SetAssayData(object = orgo_cirm43, assay = 'peaks', slot = 'motifs', new.data = motif)
  orgo_cirm43 <- RegionStats(object = orgo_cirm43, genome = BSgenome.Hsapiens.UCSC.hg38)
  orgo_cirm43 <- RunChromVAR( object = orgo_cirm43,genome = BSgenome.Hsapiens.UCSC.hg38)
  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.chromvar.Rds")
```

### Performing cisTopic and UMAP


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cisTopic)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.chromvar.Rds")

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(20:30),nCores=11,addModels=FALSE)

      #Setting up topic count selection
      pdf(paste(prefix,"model_selection.pdf",sep="."))
      par(mfrow=c(1,3))
      cirm43_cistopic_models <- selectModel(atac_cistopic_models, type='derivative')
      dev.off()
      system(paste0("slack -F ",paste(prefix,"model_selection.pdf",sep=".")," ryan_todo"))
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=orgo_cirm43,prefix="orgo_cirm43")
  cirm43_cistopic_models<-readRDS("orgo_cirm43.CisTopicObject.Rds")

  ###############################################
  #Loop through cistopic models
  cistopic_loop<-function(topic_number,object_input,models_input){
      models_input<-selectModel(models_input,select=topic_number)
      #perform UMAP on topics
      topic_df<-as.data.frame(models_input@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      print("Performed UMAP.")
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")
     
      #combine with seurat object    
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$umap<-cistopic_umap
      
      #finally plot
      plt<-DimPlot(object_input,group.by=c('DIV','cell_line'),size=0.1)+ggtitle(as.character(topic_number))
      return(plt)
  }

  library(patchwork)

  plt_list<-lapply(cirm43_cistopic_models@calc.params$runWarpLDAModels$topic,
                     FUN=cistopic_loop,
                     object_input=orgo_cirm43,
                     models_input=cirm43_cistopic_models)
  plt_list<-wrap_plots(plt_list)
  ggsave(plt_list,file="cirm43.umap_multipleTopicModels_clustering.png",height=20,width=60,limitsize=FALSE)
  system("slack -F cirm43.umap_multipleTopicModels_clustering.png ryan_todo")

  ###############################################

  #set topics based on derivative
  cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,type="derivative",keepModels=T)

  #saving model selected RDS
  saveRDS(cirm43_cisTopicObject,file="orgo_cirm43.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


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
        resolution=resolution
      )

  return(object_input)}

  orgo_cirm43<-cistopic_wrapper(object_input=orgo_cirm43,cisTopicObject=cirm43_cisTopicObject,resolution=0.2)

  plt<-DimPlot(orgo_cirm43,group.by=c("DIV","peaks_snn_res.0.2"))
  ggsave(plt,file="test.umap.pdf")
  system("slack -F test.umap.pdf ryan_todo")

  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")   ###save Seurat file

```



### Cicero for Coaccessible Networks


```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  #library(monocle3)
  library(cicero)
  library(EnsDb.Hsapiens.v86)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  #Cicero processing function
  cicero_processing<-function(object_input=orgo_atac,prefix="orgo_atac"){

      #Generate CDS format from Seurat object
      atac.cds <- as.cell_data_set(object_input,group_by="seurat_clusters")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      
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


  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  orgo_cirm43<-cicero_processing(object_input=orgo_cirm43,prefix="orgo_cirm43")

  saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.unnormGA.Rds")
  
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.unnormGA.Rds")

  # generate unnormalized gene activity matrix
  # gene annotation sample
  hg38_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

  pos <-as.data.frame(hg38_annotations,row.names=NULL)
  pos$chromosome<-paste0("chr",pos$seqnames)
  pos$gene<-pos$gene_id
  pos <- subset(pos, strand == "+")
  pos <- pos[order(pos$start),] 
  pos <- pos[!duplicated(pos$tx_id),] # remove all but the first exons per transcript
  pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS

  neg <-as.data.frame(hg38_annotations,row.names=NULL)
  neg$chromosome<-paste0("chr",neg$seqnames)
  neg$gene<-neg$gene_id
  neg <- subset(neg, strand == "-")
  neg <- neg[order(neg$start,decreasing=TRUE),] 
  neg <- neg[!duplicated(neg$tx_id),] # remove all but the first exons per transcript
  neg$end <- neg$end + 1 # make a 1 base pair marker of the TSS

  gene_annotation<- rbind(pos, neg)
  gene_annotation <- gene_annotation[,c("chromosome","start","end","gene_name")] # Make a subset of the TSS annotation columns containing just the coordinates and the gene name
  names(gene_annotation)[4] <- "gene" # Rename the gene symbol column to "gene"

  geneactivity_processing<-function(cds_input,conns_input,prefix){
      atac.cds<- annotate_cds_by_site(cds_input, gene_annotation)
      unnorm_ga <- build_gene_activity_matrix(atac.cds, conns_input)
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
  }

  conns<-as.data.frame(readRDS("orgo_cirm43_atac_cicero_conns.Rds"))
  orgo_cirm43.cicero<-readRDS("orgo_cirm43_atac_cicero_cds.Rds")
  geneactivity_processing(cds_input=as.cell_data_set(orgo_cirm43,group_by="seurat_clusters"),conns_input=conns,prefix="cirm43_atac")

  #Read in unnormalized GA
  cicero_gene_activities<-readRDS("cirm43_atac.unnorm_GA.Rds")
  orgo_cirm43[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 

  # normalize
  orgo_cirm43 <- NormalizeData(
    object = orgo_cirm43,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(orgo_cirm43$nCount_GeneActivity)
  )
  saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.Rds")
```

## Filter cells in data set


### Plotting and filtering cells

```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(dplyr)
  library(patchwork)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  saveRDS(orgo_cirm43,"orgo_cirm43.preQC.SeuratObject.Rds")

  #Cluster summaries
  dat<-orgo_cirm43@meta.data
  dat_sum<-as.data.frame(dat %>% 
  group_by(orgID,seurat_clusters,differentiation_exp,DIV) %>% 
  summarize(count=n()))
  write.table(dat_sum,"cirm43_cluster_summary_statistics.tsv",col.names=T,row.names=T,quote=F,sep="\t")

  #Setting doublets as those in top 5% of doublet_score
  threshold_number<-as.numeric(quantile(x=orgo_cirm43@meta.data$doublet_scores,prob=0.95))
  orgo_cirm43$predicted_doublets<-"False"
  orgo_cirm43@meta.data[as.numeric(orgo_cirm43@meta.data$doublet_scores)>=threshold_number,]$predicted_doublets<-"True"

  #Testing organoids for proper forebrain differentiation
  plt<-VlnPlot(orgo_cirm43,feature="FOXG1",group.by="orgID",log=T)+facet_wrap(orgo_cirm43$differentiation_exp ~ orgo_cirm43$DIV,nrow=2,scale="free_x")+theme_bw()
  ggsave(plt,file="cirm43.foxg1.pdf",width=12,height=15,limitsize=F)
  system("slack -F cirm43.foxg1.pdf ryan_todo")#post to ryan_todo
  #Based on this excluding organoid 3 from diff exp 5 (low foxg1 signal so failure to differentiate)

  #Testing organoids for tss enrichment
  plt<-ggplot()+geom_histogram(aes(x=orgo_cirm43$tss_enrichment),bins=100)+theme_minimal()+geom_vline(xintercept=1.5,color="red")
  ggsave(plt,file="cirm43.tssenrich.pdf")
  system("slack -F cirm43.tssenrich.pdf ryan_todo")#post to ryan_todo
  #filtering out cells with less than 1.5 for tss enrichment

  #Testing for uneven distribution of cells by pitstop treatment
  plt<-ggplot(dat[dat$DIV=="90" & dat$differentiation_exp=="5",],aes(fill=treatment,y=1,x=seurat_clusters))+geom_bar(position="fill",stat="identity")
  ggsave(plt,file="pitstop_cluster_distribution.pdf",width=40,height=30,limitsize=F)
  system(paste0("slack -F pitstop_cluster_distribution.pdf ryan_todo")) #looks good. note pitstop treatment was only done on div90 so some bias is to be expected

  orgo_cirm43$uniq_orgID<-paste(orgo_cirm43$differentiation_exp,orgo_cirm43$DIV,orgo_cirm43$orgID,sep=" ")
  orgo_cirm43$log10_uniq_reads<-log10(orgo_cirm43$uniq_reads)

  orgo_cirm43$pass_qc<-"True"
  orgo_cirm43@meta.data[(orgo_cirm43@meta.data$orgID==3 & orgo_cirm43@meta.data$differentiation_exp==5),]$pass_qc<-"False" #low foxg1 score
  orgo_cirm43@meta.data[(orgo_cirm43@meta.data$predicted_doublets=="True"),]$pass_qc<-"False" #predicted doublets
  orgo_cirm43@meta.data[(orgo_cirm43@meta.data$tss_enrichment<1.5),]$pass_qc<-"False" #low tss enrichment

  table(orgo_cirm43$pass_qc)
  #False  True
  # 3389 32436

  plt1<-DimPlot(orgo_cirm43,group.by=c('DIV','prep','uniq_orgID',"treatment",'differentiation_exp','seurat_clusters','predicted_doublets','pass_qc'),size=0.1)
  plt2<-FeaturePlot(orgo_cirm43,feat=c("doublet_scores","tss_enrichment","log10_uniq_reads"),col=c("white","red"),pt.size=0.1,order=T,ncol=3)
  plt<-plt1/plt2
  ggsave(plt,file="cirm43.qc.umap.pdf",width=40,height=30,limitsize=F)

  system(paste0("slack -F cirm43.qc.umap.pdf ryan_todo"))#post to ryan_todo

  orgo_cirm43_pit<-subset(orgo_cirm43,differentiation_exp=="5")
  orgo_cirm43_pit<-subset(orgo_cirm43_pit,DIV=="90")
  plt1<-DimPlot(orgo_cirm43_pit,group.by=c('treatment'),size=0.1)
  ggsave(plt1,file="cirm43.exp5.div90.pitstop.pdf",limitsize=F)
  system(paste0("slack -F cirm43.exp5.div90.pitstop.pdf ryan_todo"))#post to ryan_todo

  library(dplyr)
  data.frame(orgo_cirm43_pit@meta.data) %>% group_by(treatment) %>% summarize(mean=mean(uniq_reads))
  # A tibble: 3 x 2
  #  treatment       mean
  #  <chr>          <dbl>
  #1 No            17891.
  #2 Pitstop2      39128.
  #3 Pitstop2_Digi 46146

  anova(x=orgo_cirm43_pit$uniq_reads,y=orgo_cirm43_pit$treatment)


  orgo_qc<-subset(orgo_cirm43,pass_qc=="True")
  saveRDS(orgo_qc,file="orgo_cirm43.QC.SeuratObject.Rds") #save QC passing cells seurat object


```

### Rerun cistopic clustering after filtering cells

Now to rerun cistopic clustering with cells that pass QC metrics.


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cisTopic)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds") #reading in QC passing cells

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(20:30),nCores=11,addModels=FALSE)

      #Setting up topic count selection
      pdf(paste(prefix,"model_selection.pdf",sep="."))
      par(mfrow=c(1,3))
      cirm43_cistopic_models <- selectModel(atac_cistopic_models, type='derivative')
      dev.off()
      system(paste0("slack -F ",paste(prefix,"model_selection.pdf",sep=".")," ryan_todo"))
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=orgo_cirm43,prefix="orgo_cirm43.QC")
  cirm43_cistopic_models<-readRDS("orgo_cirm43.QC.CisTopicObject.Rds")

  ###############################################
  #Loop through cistopic models
  cistopic_loop<-function(topic_number,object_input,models_input){
      models_input<-selectModel(models_input,select=topic_number)
      #perform UMAP on topics
      topic_df<-as.data.frame(models_input@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      print("Performed UMAP.")
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")
     
      #combine with seurat object    
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$umap<-cistopic_umap
      
      #finally plot
      plt<-DimPlot(object_input,group.by=c('DIV','cell_line'),size=0.1)+ggtitle(as.character(topic_number))
      return(plt)
  }

  library(patchwork)

  plt_list<-lapply(cirm43_cistopic_models@calc.params$runWarpLDAModels$topic,
                     FUN=cistopic_loop,
                     object_input=orgo_cirm43,
                     models_input=cirm43_cistopic_models)
  plt_list<-wrap_plots(plt_list)
  ggsave(plt_list,file="cirm43.qc.umap_multipleTopicModels_clustering.png",height=20,width=60,limitsize=FALSE)
  system("slack -F cirm43.qc.umap_multipleTopicModels_clustering.png ryan_todo")

  ###############################################

  #set topics based on derivative
  cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,type="derivative",keepModels=T)

  #saving model selected RDS
  saveRDS(cirm43_cisTopicObject,file="orgo_cirm43.QC.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


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
        resolution=resolution
      )

  return(object_input)}
  DefaultAssay(orgo_cirm43)<-"peaks"
  orgo_cirm43<-cistopic_wrapper(object_input=orgo_cirm43,cisTopicObject=cirm43_cisTopicObject,resolution=0.2)

  plt<-DimPlot(orgo_cirm43,group.by=c("DIV","seurat_clusters"))
  ggsave(plt,file="qc.umap.pdf")
  system("slack -F qc.umap.pdf ryan_todo")

  saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")   ###save Seurat file

```
## Organoid Cell type analysis

### Celltype Assignment of Clusters

Cell Type Assignment of Organoid Clusters
Doing this in three parts.
1. Using bulk sorted RG, IPC, eN and iN RNA markers compared to our ATAC cluster gene activity scores
2. Using bulk sorted RG, IPC, eN and iN ATAC motifs compared to our ATAC cluster motifs
3. Using single-cell Primary Cortex RG, IPC, eN and iN annotated cells to define signatures and perform CCA for label transfer


```R
  #https://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  set.seed(1234)
  library(reshape2)
  library(dplyr)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(circlize)
  library(viridis)
  # Load the pre-processed scRNA-seq and scATAC-seq data

  #Public RNA
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data")
  pubprimary<-readRDS("PublicPrimary.SeuratObject.rds")
  #perform random subsampling, since cross data integration is robust to cell count and its taking forever
  #using 1/10th the cells (~10k)
  pubprimary <- subset(pubprimary, cells = sample(x = colnames(pubprimary@assays$RNA@data), size = length(colnames(pubprimary@assays$RNA@data))/10) )
  pubprimary<-SetIdent(pubprimary,value="Type")
  #subset to cell types expected to occur in organoids
  pubprimary<-subset(pubprimary,idents=c("Excitatory Neuron","IPC","Radial Glia"))

  #Our ATAC
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS(file="orgo_cirm43.QC.SeuratObject.Rds")

  #1. Using bulk sorted RG, IPC, eN and iN RNA markers compared to our ATAC cluster gene activity scores
  #Corticogenic data on basic cell types.
  #data from http://data.nemoarchive.org/5923ca16c51011e99da31f7757ebac1c/
  #described in https://www.nature.com/articles/s41586-020-2825-4?WT.ec_id=NATURE-202010&sap-outbound-id=60313C942AEB24BFE1AF8DD74FB2E05B7385E720#data-availability
  #Bulk ATAC Peaks for marker sorted cell types located in /home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/Song_2020
  #For overlapping these narrowPeak features do:

    markers<-c("CTCF","EMX1","EMX2","LHX2","PAX6","RFX4","SOX2",
               "TBR1","EOMES","NEUROD1","NEUROD2","NEUROG1","TGIF1","TGIF2",
               "DLX1","DLX2","DLX6","GSX2","LHX6",
               "POU3F3","POU3F2","TFAP4")


    #Plot motifs alongside chromvar plot
    library(ggplot2)
    library(patchwork)


    motif_order<-names(orgo_cirm43@assays$peaks@motifs@motif.names[match(markers,unlist(orgo_cirm43@assays$peaks@motifs@motif.names),nomatch=0)])
    plt<-MotifPlot(object = orgo_cirm43@assays$peaks@motifs,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file="tf.heatmap.motif.pdf",height=10,width=2,limitsize=F)
    system("slack -F tf.heatmap.motif.pdf ryan_todo")


    #Setting up chromvar matrix from CIRM43
    tfList <- getMatrixByID(JASPAR2020, ID=row.names(orgo_cirm43@assays$chromvar@data)) 
    tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
    
    motif_plots<-MotifPlot(
      object = orgo_cirm43,
      motifs = row.names(orgo_cirm43@assays$chromvar@data)[match(markers,tfList,nomatch=0)],
      assay = 'peaks'
    )

    ggsave(motif_plots,file="marker_motifs.pdf")
    system("slack -F marker_motifs.pdf ryan_todo")

    dat_tf<-orgo_cirm43@assays$chromvar@data
    row.names(dat_tf)<-tfList
    dat_tf<-data.frame(t(dat_tf))
    dat_tf$cellID<-row.names(dat_tf)

    dat_tf$cluster_ID<-orgo_cirm43@meta.data[match(orgo_cirm43@meta.data$cellID,dat_tf$cellID),]$seurat_clusters  #append cluster ID to column
    dat_tf<-dat_tf[colnames(dat_tf) %in% c("cluster_ID",markers)] #subset to markers
    dat_tf<-melt(dat_tf) #reshape to long format
    dat_tf<-as.data.frame(dat_tf %>% group_by(cluster_ID,variable) %>% summarize(mean_chromvar=mean(value))) #group by to summarize markers by column

    #plot as heatmap
    dat_tf<-dcast(dat_tf,cluster_ID~variable)
    row.names(dat_tf)<-dat_tf$cluster_ID
    dat_tf<-dat_tf[colnames(dat_tf) %in% markers]
    dat_tf[which(is.na(dat_tf),arr.ind=T)]<-0 #set na values to 0 for clustering
    dat_tf<-as.data.frame(t(dat_tf))
    clus_order<-c("5","4","0","3","2","1","6")
    column_split<-factor(unlist(lapply(strsplit(colnames(dat_tf),"_"),"[",1)),levels=c("3","1","2","0"))
    dat_tf<-dat_tf[colnames(dat_tf) %in% clus_order]
    colfun=colorRamp2(quantile(unlist(dat_tf), probs=seq(0.1,0.9,0.1)),cividis(length(seq(0.1,0.9,0.1))))


    plt<-Heatmap(dat_tf,
                row_order=match(markers,row.names(dat_tf))[!is.na(match(markers,row.names(dat_tf)))],
                column_order=clus_order,
                col=colfun)
                #column_split=column_split)
    pdf("cirm43_celltype_tfHeatmap.pdf")  
    plt
    dev.off()
    system("slack -F cirm43_celltype_tfHeatmap.pdf ryan_todo")

  #2. Using bulk sorted RG, IPC, eN and iN ATAC motifs compared to our ATAC cluster motifs
    #Setting up gene activity matrix
    markers<-c("SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3",
               "NEUROG1","SSTR2","EOMES","PPP1R17","NEUROD4",
               "SLC17A7","NEUROD6","SATB2","TBR1","SLA",
               "DLX2","DLX1","LHX6","GAD1")
    dat_ga<-orgo_cirm43@assays$GeneActivity@data
    dat_ga<-data.frame(t(dat_ga))
    dat_ga$cellID<-row.names(dat_ga)
    #append cluster ID to column
    dat_ga$cluster_ID<-orgo_cirm43@meta.data[match(orgo_cirm43@meta.data$cellID,dat_ga$cellID),]$seurat_clusters
    #subset to markers
    dat_ga<-dat_ga[colnames(dat_ga) %in% c("cluster_ID",markers)]
    #reshape to long format
    dat_ga<-melt(dat_ga)
    #group by to summarize markers by column
    dat_ga<-as.data.frame(dat_ga %>% group_by(cluster_ID,variable) %>% summarize(mean_ga=mean(value)))
    #plot as heatmap
    dat_ga<-dcast(dat_ga,cluster_ID~variable)
    row.names(dat_ga)<-dat_ga$cluster_ID
    dat_ga<-dat_ga[colnames(dat_ga) %in% markers]
    #zscore values
    dat_ga<-scale(dat_ga)
    #set na values to 0 for clustering
    dat_ga<-data.frame(t(dat_ga))
    clus_order<-c("5","4","0","3","2","1","6")
    clus_order<-paste0("X",clus_order)
    #column_split<-factor(unlist(lapply(strsplit(colnames(dat_ga),"_"),"[",1)),levels=c("X3","X1","X2","X0"))
    dat_ga<-dat_ga[colnames(dat_ga) %in% clus_order]
    
    colfun=colorRamp2(quantile(unlist(dat_ga), probs=seq(0.1,0.9,0.1)),magma(length(seq(0.1,0.9,0.1))))

    plt<-Heatmap(dat_ga,
                row_order=match(markers,row.names(dat_ga))[!is.na(match(markers,row.names(dat_ga)))],
                column_order=clus_order,
                col=colfun)
                #column_split=column_split)

    pdf("cirm43_celltype_gaHeatmap.pdf")
    plt
    dev.off()
    system("slack -F cirm43_celltype_gaHeatmap.pdf ryan_todo")

  #3. Using single-cell Primary Cortex RG, IPC, eN and iN annotated cells to define signatures and perform CCA for label transfer
    markers<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/pubprimary.markers.txt",header=T,sep="\t")

    transfer.anchors <- FindTransferAnchors(
      reference = pubprimary,
      reference.assay="RNA",
      query = orgo_cirm43,
      query.assay="GeneActivity",
      reduction   = 'cca',
      features=markers$gene,
      verbose=T
    )
    saveRDS(transfer.anchors,"orgo_cirm43.PublicPrimary.transferanchors.rds")

    orgo_cirm43[["celltype_prediction"]]<- TransferData(
      anchorset = transfer.anchors,
      refdata = pubprimary$Type,
      weight.reduction = "cca",
      dims = 1:30,
      prediction.assay=T
    )

    saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")



    plt1<-FeaturePlot(orgo_cirm43,features=c('Radial Glia'),pt.size=0.1,min.cutoff="q75",order=T,col=c("white","#EFA421"))
    plt2<-FeaturePlot(orgo_cirm43,features=c('Excitatory Neuron'),pt.size=0.1,min.cutoff="q75",order=T,col=c("white","#6B78BA"))
    plt3<-FeaturePlot(orgo_cirm43,features=c('IPC'),pt.size=0.1,min.cutoff="q75",sort=T)

    plt<-plt1/plt2/plt3
    ggsave(plt,file="cirm43.predictedid.umap.png",width=10,height=30,limitsize=F)
    ggsave(plt,file="cirm43.predictedid.umap.pdf",width=10,height=30,limitsize=F)
    system("slack -F cirm43.predictedid.umap.pdf ryan_todo")

    predictdat<-as.data.frame(t(orgo_cirm43@assays$celltype_prediction@data))[,1:3]
    predictdat$seurat_clusters<-orgo_cirm43@meta.data[row.names(predictdat),]$seurat_clusters
  
    predictdat<-melt(predictdat)
    predictdat<-as.data.frame(predictdat %>% group_by(seurat_clusters,variable) %>% summarize(average=mean(value)))

    predictdat<-dcast(predictdat,seurat_clusters~variable)
    row.names(predictdat)<-predictdat$seurat_clusters
    predictdat<-predictdat[!(colnames(predictdat) %in% c("seurat_clusters"))]
    predictdat<-as.data.frame(t(scale(predictdat,scale=T)))
    clus_order<-c("5","4","0","3","2","1","6")
    colfun=colorRamp2(c(-2,0,2),c("#000000","#FFFFFF","#FF0000"))

    predictdat<-predictdat[colnames(predictdat) %in% clus_order]
    plt<-Heatmap(predictdat,
    row_order=c("Radial Glia","IPC","Excitatory Neuron"),
    column_order=clus_order,
    col=colfun)
    pdf("predictedid.heatmap.pdf")
    plt
    dev.off()
    system("slack -F predictedid.heatmap.pdf ryan_todo")

  saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")
```

### Cell cycle testing

Seurat has a stored set of cell cycle genes that we can use to assess cell cycle signatures.

[Following this.](https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html)
Using gene lists based on cell cycle markers listed in https://www.cell.com/neuron/pdf/S0896-6273(19)30561-6.pdf
Supplementary Table 7.

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(RColorBrewer)

  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")
  
  s.genes <- c('AK2', 'SLC25A5', 'TMEM98', 'GGCT', 'DBF4', 'PAX6', 'SPAG9', 'RPS20', 'BAZ1B', 'UQCRC1', 'ANLN', 'BRCA1', 'DDX11', 'TACC3', 'VIM', 'HMGB3', 'CENPQ', 'RFC1', 'MAT2B', 'SPDL1', 'CNTLN', 'TPR', 'RCN1', 'RFC2', 'RAD51', 'POLQ', 'MPHOSPH9', 'NNAT', 'SYNE2', 'COL11A1', 'QSER1', 'SNCAIP', 'MCM10', 'YBX1', 'ASPM', 'MPPED2', 'PKM', 'RHOA', 'PRR11', 'NUCKS1', 'RAD18', 'SMC1A', 'HMMR', 'MCM2', 'CA12', 'PTPLAD1', 'ENO1', 'GTSE1', 'ACTB', 'MCM6', 'SPAG5', 'UBE2T', 'POLD3', 'JADE1', 'FBLN1', 'SLC1A3', 'XRCC5', 'KIF22', 'RBL1', 'NDC80', 'HSPB11', 'XPO1', 'GSTP1', 'SRRT', 'SF3B2', 'TFAP2C', 'TPX2', 'RPLP0', 'FUS', 'KIF4A', 'ORC6', 'ZFHX4', 'HNRNPC', 'SUPT16H', 'WDR76', 'PHGDH', 'EZR', 'MYL6', 'CLSPN', 'CDC45', 'CDC6', 'CBX5', 'MSH2', 'CDC5L', 'HNRNPH3', 'H2AFY2', 'HNRNPM', 'RANBP1', 'SNRPD3', 'CENPM', 'HMGXB4', 'MCM5', 'RPL3', 'FKBP3', 'CEP128', 'ERH', 'VRK1', 'PNN', 'GINS1', 'PLCB4', 'NOP56', 'SMCHD1', 'RBBP8', 'POLA1', 'STAG2', 'RBBP7', 'CMC2', 'UBE2I', 'CCP110', 'CEP152', 'SFRP1', 'EEF1D', 'MCM4', 'RNASEH2A', 'HNRNPUL1', 'LIG1', 'CDK6', 'PTN', 'H2AFV', 'CHCHD2', 'HSPB1', 'NUDT1', 'PPP1R17', 'LSM5', 'RPA3', 'EZH2', 'RHEB', 'LHX2', 'ELAVL2', 'NSMCE4A', 'SMC3', 'KPNB1', 'CBX1', 'PFN1', 'TMEM97', 'WHSC1', 'NCAPG', 'CCDC34', 'MDK', 'C11orf58', 'CORO1C', 'PTGES3', 'FOXM1', 'RAD51AP1', 'TIMELESS', 'GAPDH', 'CHD4', 'TPI1', 'C12orf57', 'LDHB', 'SRSF9', 'FBXO5', 'SRSF3', 'MCM3', 'E2F3', 'GMNN', 'TTK', 'ERBB2IP', 'LMNB1', 'H2AFY', 'SMAD5', 'SMC4', 'TFDP2', 'HES1', 'ECT2', 'BBX', 'NCL', 'PPM1G', 'RPS15', 'FANCL', 'SRSF7', 'MSH6', 'SRSF4', 'IVNS1ABP', 'ACADM', 'PRDX1', 'CNN3', 'CENPF', 'RPA2', 'MESDC2', 'STAG1', 'CASP8AP2', 'HMGN3', 'RPN2', 'CCND2', 'CTNNAL1', 'WDR34', 'SET', 'CNTRL', 'FAM178A', 'HELLS', 'ENY2', 'MASTL', 'EXOSC8', 'EGR1', 'TMPO', 'NFYB', 'NCAPH', 'MND1', 'CCDC18', 'CBX3', 'HNRNPA2B1', 'WIPF3', 'NPY', 'ZWINT', 'CDKN2C', 'DDX39A', 'CENPK', 'NEUROD4', 'CDK2', 'TUBA1B', 'STIL', 'HJURP', 'BAZ2B', 'EXOSC9', 'CKS2', 'SNRPC', 'HIST1H1D', 'HIST1H1A', 'GLO1', 'DEK', 'SOX9', 'PPDPF', 'SNRPD2', 'SNRPB', 'MGME1', 'MCM8', 'HNRNPR', 'RALY', 'UBA2', 'DLGAP5', 'YEATS4', 'PIN1', 'HP1BP3', 'PKMYT1', 'PAICS', 'SPECC1', 'CALU', 'HAT1', 'DUT', 'FAM64A', 'ILF3', 'PARP2', 'MIS18BP1', 'SGOL1', 'GADD45G', 'LSM4', 'DNMT1', 'AKAP12', 'GINS2', 'PSMC3IP', 'TOP2A', 'RAN', 'PCNA', 'NES', 'NASP', 'MYH10', 'TPT1', 
    'RFC3', 'ANKRD32', 'LRRCC1', 'MEIS2', 'TMEM106C', 'RBM17', 'SYNCRIP', 'ATP5G2', 'CDK4', 'HNRNPA1', 'AHI1', 'DHX9', 'RNASEH2B', 'CKAP2', 'SCRN1', 'SRSF1', 'BRIP1', 'ACTL6A', 'TRA2B', 'SMC2', 'CDK5RAP2', 'ANP32B', 'RPL35', 'RPS6', 'GGH', 'RDX', 'CTDSPL2', 'NUSAP1', 'KIF23', 'CASC5', 'RPLP1', 'KIF11', 'KIF20B', 'DNA2', 'BARD1', 'PPIG', 'MNS1', 'ZGRF1', 'HNRNPD', '44450', 'CENPE', 'HADH', 'SCAF11', 'PHLDA1', 'SNRPF', 'NEDD1', 'ASCL1', 'BRCA2', 'DIAPH3', 'TMX1', 'SERF2', 'COMMD4', 'FANCI', 'MFGE8', 'ANAPC11', 'NFIC', 'SAE1', 'PLK4', 'ITGB3BP', 'KIF2C', 'NUF2', 'ANP32E', 'DTL', 'ILF2', 'SRP9', 'PARP1', 'LBR', 'SNRPG', 'SLC20A1', 'CDCA7', 'GULP1', 'HSPD1', 'HES6', 'FANCD2', 'CENPC', 'CCNA2', 'MYO10', 'G3BP1', 'PHIP', 'MMS22L', 'CDCA5', 'NCAPG2', 'NONO', 'RBMX', 'GINS4', 'PLIN2', 'HAUS6', 'RPL7A', 'ZEB1', 'MKI67', 'SSRP1', 'RPS3', 'INCENP', 'CHEK1', 'DSN1', 'HIRIP3', 'ITGB1', 'CCT5', 'MAGI1', 'NCAPD3', 'CENPU', 'CENPJ', 'SCHIP1', 'MZT2B', 'HAUS1', 'SPC25', 'TMEM123', 'HNRNPDL', 'CENPH', 'CARHSP1', 'SMARCA5', 'HNRNPU', 'SREK1', 'CHD1', 'BUB3', 'BTG3', 'DBI', 'TMEM237', 'VBP1', 'ATAD2', 'BUB1B', 'CCNB2', 'TMSB15A', 'EIF5B', 'MIS18A', 'C21orf58', 'PCNT', 'FDPS', 'IER2', 'RPL8', 'SRSF2', 'RACGAP1', 'SPC24', 'ASRGL1', 'MAGOH', 'RBBP4', 'NFIA', 'USP1', 'PEA15', 'KIAA1524', 'EOMES', 'SGOL2', 'GMPS', 'TOPBP1', 'KIF15', 'RFC4', 'SLBP', 'RNF168', 'H2AFZ', 'PGRMC2', 'HMGB2', 'MAD2L1', 'ANXA5', 'RHOBTB3', 'STK17A', 'PTTG1', 'CDCA7L', 'FABP5', 'RAD21', 'PSIP1', 'HNRNPK', 'MELK', 'SPTSSA', 'SKA3', 'LRR1', 'E2F7', 'PSMC3', 'CEP295', 'CKB', 'CENPN', 'MCM7', 'CENPV', 'B2M', 'FAM111A', 'KIAA0101', 'SNRPD1', 'ACAA2', 'RRM1', 'TPM4', 'CHAF1A', 'C19orf48', 'PRDX2', 'TK1', 'SRRM2', 'RPSA', 'PBK', 'RBPJ', 'GNG4', 'HIST1H1E', '44441', 'DTYMK', 'FEN1', 'STXBP6', 'HNRNPH1', 'SDC2', 'CKAP2L', 'BUB1', 'CNBP', 'HNRNPF', 'UBE2E3', 'KCNAB3', 'HNRNPA3', 'CDK1', 'UBB', 'FOS', 'EMX2', 'PA2G4', 'LSM3', 'SHCBP1', 'CHD7', 'ESCO2', 'CXXC5', 'RRM2', 'RPS7', 'ID4', 'CKS1B', 'INSM1', 'SMARCC1', 'GOLIM4', 'GNG5', 'EXO1', 'ZWILCH', 'LARP7', 'CEP135', 'RSRC1', 'UBE2C', 'CSRP2', 'CCNE2', 'BANF1', 'CCDC14', 'NR2F1', 'COX8A', 'TYMS', 'PXMP2', 'RPLP2', 'JUN', 'HNRNPA0', 'ARL6IP6', 'KDELC2', 'GEN1', 'SUZ12', 'RMI1', 'AURKB', 'RAD23A', 'SSTR2', 'NPM1', 'PENK', 'SOX2', 'ZBTB20', 'NEUROG1', 'SNRPE', 'RTKN2', 'IDH2', 'SKA2', 'HIST2H2AC', 'HIST1H1B', 'POU3F2', 'H1FX', 'NDUFA6', 'SIVA1', 'ZFP36L1', 'MYBL1', 'NKAIN3', '44449', 'NAP1L1', 'PTMA', 'HIST1H1C', 'TUBB4B', 'H2AFX', 'SUMO2', 'FAM111B', 'H1F0', 'HMGB1', 'PPIA', 'XRCC6', 'XRCC2', 'HIST1H4C', 'PCBP2', 'BLM', 'HNRNPAB', 'HES5', 'ELOVL2', 'PRIM1', 'HMGN5', 'RPL23A', 'ASPH', 'WDHD1', 'BAZ1A', 'SMOC1', 'ARHGAP11A', 'HMGN2', 'CCDC152', 'SMC5', 'PRC1', 'CCDC167', 'CENPW', 'GPANK1', 'NAP1L4', 'TMSB4X', 'HMGN1', 'HN1L', 'DNAJC9', 'MIR99AHG', 'CKLF', 'UBA52', 'FGD5-AS1', 'DHFR', 'RPL41', 'DLEU2', 'LINC01158', 'MAGI2-AS3', 'PEG10', 'SNHG6', 'TMEM158', 'PRKDC')
  

  g2m.genes <- c('CDC27', 'DBF4', 'PAX6', 'SPAG9', 'NCAPD2', 'ANLN', 'BRCA1', 'TACC3', 'DEPDC1', 'VIM', 'HMGB3', 'DEPDC1B', 'MAT2B', 'SPDL1', 'PSMA4', 'CNTLN', 'TPR', 'SLC4A8', 'POLQ', 'MPHOSPH9', 'NNAT', 'SYNE2', 'CCAR1', 'COL11A1', 'QSER1', 'SPA17', 'SUGP2', 'HMG20B', 'ASPM', 'MPPED2', 'PRR11', 'LAPTM4A', 'NUCKS1', 'SMC1A', 'HMMR', 'NDE1', 'SRI', 'GTSE1', 'ACTB', 'SPAG5', 'UBE2T', 'JADE1', 'PPP2R5C', 'PCM1', 'SLC1A3', 'KIF22', 'NDC80', 'STK17B', 'XPO1', 'REST', 'SEPHS1', 'AURKA', 'AAMDC', 'TPX2', 'DYNLL1', 'KIF4A', 'ORC6', 'G2E3', 'PHGDH', 'EZR', 'CBX5', 'SUCO', 'HNRNPH3', 'IFT74', 'HNRNPM', 'RANBP1', 'RANGAP1', 'CDKN3', 'KIAA0586', 'DHRS7', 'CEP128', 'ERH', 'VRK1', 'EMC9', 'CDC25B', 'FAM83D', 'SMARCA1', 'CMC2', 'CEP152', 'OIP5', 'MYEF2', 'SFRP1', 'EEF1D', 'HNRNPUL1', 'CARD8', 'CDK6', 'PON2', 'PTN', 'H2AFV', 'HSPB1', 'PPP1R17', 'LSM5', 'EZH2', 'RHEB', 'SMC3', 'UBE2S', 'CBX1', 'NMU', 'NEIL3', 'WHSC1', 'NCAPG', 'CCDC34', 'MDK', 'CORO1C', 'ATP5B', 'PTGES3', 'FOXM1', 'RAD51AP1', 'CDKN1B', 'TIMELESS', 'MRPL51', 'CDCA3', 'FBXO5', 'SRSF3', 'GMNN', 'QKI', 'TTK', 'BRD8', 'KIF20A', 'LMNB1', 'H2AFY', 'SMC4', 'CEP70', 'TFDP2', 'HES1', 'ECT2', 'FXR1', 'CENPA', 'GCA', 'SFPQ', 'TTF2', 'CDC20', 'PRDX1', 'STMN1', 'NEK2', 'CENPF', 'TXNDC12', 'KIF14', 'HMGN3', 'FBXL5', 'CCND2', 
    'CNTRL', 'PHF19', 'CENPL', 'ENY2', 'EXOSC8', 'EGR1', 'TMPO', 'NCAPH', 'MND1', 'PSPC1', 'KIF18A', 'DESI2', 'GPSM2', 'ZC3H7A', 'CCDC18', 'CBX3', 'HNRNPA2B1', 'NPY', 'CALD1', 'ZWINT', 'CIT', 'CDKN2C', 'DDX39A', 'CENPK', 'NEUROD4', 'TUBA1B', 'STIL', 'HJURP', 'MORF4L2', 'CKS2', 'SNRPC', 'HIST1H1D', 'HIST1H1A', 'GLO1', 'DEK', 'MT2A', 'SOX9', 'MGME1', 'HNRNPR', 'NSRP1', 'DLGAP5', 'HP1BP3', 'KNSTRN', 'PALLD', 'FAM64A', 'MIS18BP1', 'SGOL1', 'AKAP12', 'TOP2A', 'DNAJB1', 'RAN', 'PCBD2', 'NES', 'MYH10', 'CCNA1', 'CCNB1', 'PSRC1', 'LDHA', 'CDCA8', 'AKIRIN2', 'TROAP', 'HNRNPA1', 'RNASEH2B', 'CKAP2', 'BORA', 'LMO7', 'SCRN1', 'IGF2BP3', 'CALCOCO2', 'DCAF7', 'ACTL6A', 'TRA2B', 'ODF2', 'SMC2', 'CDK5RAP2', 'ANP32B', 'DCTN3', 'ARHGEF39', 'RDX', 'NUSAP1', 'KIF23', 'CASC5', 'CENPO', 'KIF11', 'CEP55', 'KIF20B', 'BARD1', 'COX17', 'CENPE', 'PHLDA1', 'NEDD1', 'ASCL1', 'GAS2L3', 'BRCA2', 'TMBIM6', 'DIAPH3', 'TMX1', 'SERF2', 'PIF1', 'TICRR', 'PLK4', 'KIF2C', 'NUF2', 'HDGF', 'ANP32E', 'RAB13', 'ILF2', 'CNIH4', 'LBR', 'HNRNPLL', 'CALM2', 'SNRPG', 'CCDC150', 'HES6', 'FANCD2', 'CENPC', 'CCNA2', 'SFRP2', 'MYO10', 'G3BP1', 'PHIP', 'CDCA5', 'NCAPG2', 'RBMX', 'PLIN2', 'ZEB1', 'ADD3', 'MKI67', 'SESN3', 'INCENP', 'HIRIP3', 'CCT5', 'SCLT1', 'CENPU', 'CENPJ', 'MZT2B', 'SPC25', 'CENPH', 'CETN3', 'SMARCA5', 'HNRNPU', 'CEP112', 'ENAH', 'BUB3', 'BTG3', 'SKA1', 'DBI', 'TMEM237', 'VBP1', 'FBXO43', 'ATAD2', 'BUB1B', 'NRG1', 'CCNB2', 'TMSB15A', 'CDC25C', 'TAGLN2', 'MIS18A', 'PTMS', 'CALM3', 'C21orf58', 'PCNT', 'SRSF2', 'RACGAP1', 'SPC24', 'CCNF', 'ASRGL1', 'PPAP2B', 'NFIA', 'USP1', 'FUBP1', 'PEA15', 'NEUROD1', 'DCAF16', 'KIAA1524', 'EOMES', 'SGOL2', 'SMIM14', 'KIF15', 'H2AFZ', 'INTU', 'HMGB2', 'SAP30', 'MAD2L1', 'ANXA5', 'CEP44', 'ITGA2', 'STK17A', 'PTTG1', 'FABP5', 'RAD21', 'PSIP1', 'HNRNPK', 'MELK', 'SKA3', 'CEP295', 'IKBIP', 'CKB', 'CENPN', 'WEE1', 'HSP90B1', 'B2M', 'FAM111A', 'KIAA0101', 'PLK1', 'TPM4', 'TUBA1C', 'CTNNB1', 'PBK', 'HIST1H1E', 'DTYMK', 'HNRNPH1', 'CKAP2L', 'BUB1', 'DCXR', 'HNRNPA3', 'CDK1', 'UBB', 'FOS', 'EMX2', 'ARL6IP1', 'NUDCD2', 'KIF5B', 'SHCBP1', 'CHD7', 'ESCO2', 'ATF7IP', 'RHNO1', 'RRM2', 'ID4', 'ZNF24', 'DCP2', 'CKS1B', 'RNF26', 'FKBP2', 'GOLIM4', 'GNG5', 'LARP7', 'CEP135', 'RSRC1', 'UBE2C', 'CKAP5', 'BANF1', 'CCDC14', 'NR2F1', 'RUVBL1', 'TUBB6', 'ACBD7', 'COX8A', 'TYMS', 'TGIF1', 'JUN', 'HNRNPA0', 'C2orf69', 'LCORL', 'GEN1', 'SUZ12', 'APOLD1', 'AURKB', 'PENK', 'SOX2', 'ZBTB20', 'RTKN2', 'FIGN', 'KPNA2', 'CEP97', 'SKA2', 'CEP57L1', 'RUVBL2', 'PTTG1IP', 'SETD8', 'HIST1H1B', 'POU3F2', 'CDCA2', 'H1FX', 'RPS27L', 'UBALD2', 'PARPBP', 'ZFP36L1', 'MYBL1', 'NKAIN3', 'SAPCD2', 'PPP1CC', '44449', 'NAP1L1', 'HIST1H1C', 'ARHGAP11B', 'TUBB4B', 'H2AFX', 'HN1', 'HMGB1', 'XRCC6', 'XRCC2', 'ZMYM1', 'HIST1H4C', 'PCBP2', 'HES5', 'HMGN5', 'HSD17B11', 'HYLS1', 'ECI2', 'SMOC1', 'ARHGAP11A', 'HMGN2', 'CCDC152', 'TOP1', 'PRC1', 'CCDC167', 'CENPW', 'HSPA1B', 'HSPA1A', 'MZT1', 'TMSB4X', 'HMGN1', 'HLA-A', 'TRIM59', 'MXD3', 'MIR99AHG', 'DLEU2', 'LINC01158', 'CRNDE')


  DefaultAssay(orgo_cirm43) <- 'GeneActivity'

  orgo_cirm43 <- CellCycleScoring(orgo_cirm43, s.features = s.genes, g2m.features = g2m.genes, set.ident = F,search=T)

  saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")


  plt<-FeaturePlot(orgo_cirm43,features=c("G2M.Score","S.Score"),cols=c("lightgrey","red"),order=T,min.cutoff="q75")
  ggsave(plt,file="cellcycle_score.pdf",width=20)
  system("slack -F cellcycle_score.pdf ryan_todo")

  summary(orgo_cirm43$G2M.Score)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.02677  0.01296  0.02688  0.02729  0.04103  0.11135

  summary(orgo_cirm43$S.Score)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.02898  0.01271  0.02555  0.02580  0.03848  0.08810

```


### Count of celltype and clusters per organoid
Generated bar plots of cell count per organoid. 
Also looking at FOXG1+ Expression per organoid to make sure they are forebrain specified.

```R
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(ComplexHeatmap)
  library(patchwork)
  library(dplyr)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

  dat<- as.data.frame(as.data.frame(orgo_cirm43@meta.data) %>% group_by(differentiation_exp,orgID,DIV,seurat_clusters) %>% summarize(count=n()))
  plt<-ggplot(dat,aes(x=as.factor(paste(differentiation_exp,orgID)),y=count,fill=paste(seurat_clusters)))+geom_bar(position="stack",stat="identity")+facet_wrap(dat$differentiation_exp ~ dat$DIV,nrow=2,scale="free_x")+theme_bw()
  ggsave(plt,file="cellcount_stacked_bar.pdf")
  system("slack -F cellcount_stacked_bar.pdf ryan_todo")

  plt<-ggplot(dat,aes(x=as.factor(paste(differentiation_exp,orgID)),y=count,fill=paste(seurat_clusters)))+geom_bar(position="fill",stat="identity")+facet_wrap(dat$differentiation_exp ~ dat$DIV,nrow=2,scale="free_x")+theme_bw()
  ggsave(plt,file="cellcount_perc_bar.pdf")
  system("slack -F cellcount_perc_bar.pdf ryan_todo")

  plt<-VlnPlot(orgo_cirm43,feat="FOXG1",group.by=paste(differentiation_exp,orgID))+facet_wrap(orgo_cirm43$differentiation_exp ~ orgo_cirm43$DIV,nrow=2,scale="free_x")+theme_bw() #feature from gene activity
  ggsave(plt,file="uniq_orgID_FOXG1.pdf")
  system("slack -F uniq_orgID_FOXG1.pdf ryan_todo")

  plt<-VlnPlot(orgo_cirm43,feat="SATB2",group.by=paste(differentiation_exp,orgID))+facet_wrap(orgo_cirm43$differentiation_exp ~ orgo_cirm43$DIV,nrow=2,scale="free_x")+theme_bw() #feature from gene activity
  ggsave(plt,file="uniq_orgID_SATB2.pdf")
  system("slack -F uniq_orgID_SATB2.pdf ryan_todo")

  plt<-VlnPlot(orgo_cirm43,feat="CUX1",group.by="cluster_ID")+facet_wrap(orgo_cirm43$differentiation_exp ~ orgo_cirm43$DIV,nrow=2,scale="free_x")+theme_bw() #feature from gene activity
  ggsave(plt,file="uniq_orgID_CUX1.pdf")
  system("slack -F uniq_orgID_CUX1.pdf ryan_todo")

  #Based on this organoid 3 looks suspect, seeing where it falls in the plot
  orgo_cirm43$suspect_organoid<-"FALSE"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$uniq_orgID==3,]$suspect_organoid<-"TRUE"
  plt<-DimPlot(orgo_cirm43,group.by="suspect_organoid",order=T)
  ggsave(plt,file="uniq_orgID3.pdf")
  system("slack -F uniq_orgID3.pdf ryan_todo")
```

### Feature values from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6544371/#SD2

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

#From https://www.nature.com/articles/s41586-020-1962-0#MOESM1 Supplementary Table 9 and markers by Supplementary Table 3
areal_signatures<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/areal_signatures.tsv",sep="\t",head=T)
organoid_cluster_markers<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/bhaduri_organoid_markers.tsv",sep="\t",head=T)
organoid_cluster_desc<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/bhaduri_organoid_marker_clusters.tsv",sep="\t",head=T)


module_membership<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/gene_modules.membership.txt",sep="\t",head=T) #Formatted from Table S4
colnames(module_membership)<-c("gene","module_name")
module_description<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/gene_module.description.txt",sep="\t",head=T) #Formatted from Table S4
module_curated<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/gene_modules_curated.txt",sep="\t",head=F)
colnames(module_curated)<-c("module_name","module_no","module_description","celltype")

module_dat<-lapply(unique(module_description$Network),function(x) {
  modulename<-strsplit(x,"[.]")[[1]][4]
  c(module_membership[module_membership$module_name==modulename,]$gene)})
names(module_dat)<-unique(module_description$Network)
lapply(module_dat, function(i) c(length(i),sum(i %in% row.names(orgo_cirm43@assays$GeneActivity@data))))

orgo_cirm43_wip<-AddModuleScore(orgo_cirm43,features=module_dat,assay="GeneActivity",name=paste0("module_",names(module_dat),"_"), search=T)
saveRDS(orgo_cirm43_wip,file="orgo_cirm43.SeuratObject.ModulesWIP.Rds")
orgo_cirm43_wip<-readRDS("orgo_cirm43.SeuratObject.ModulesWIP.Rds")

predictdat<-orgo_cirm43_wip@meta.data
predictdat<-predictdat[startsWith(colnames(predictdat),"module")| colnames(predictdat) %in% c("seurat_clusters")]
predictdat<-predictdat[!(colnames(predictdat) %in% c("module_prediction.score.max","module_predicted.id"))]

predictdat<-melt(predictdat)
predictdat<-as.data.frame(predictdat %>% group_by(seurat_clusters,variable) %>% summarize(average=median(value)))

predictdat<-dcast(predictdat,seurat_clusters~variable)
row.names(predictdat)<-predictdat$seurat_clusters
predictdat<-predictdat[,2:ncol(predictdat)]
predictdat<-as.data.frame(t(predictdat))
    clus_order<-c("5","4","0","3","2","1","6")
predictdat<-predictdat[colnames(predictdat) %in% clus_order,]
predictdat<-as.data.frame(t(scale(t(predictdat),scale=T)))


plt<-Heatmap(predictdat,
row_names_gp = gpar(fontsize = 4)
#column_order=clus_order
#row_split=unlist(lapply(strsplit(row.names(predictdat),"[.]"),"[",1))

)

pdf("predictedid.heatmap.modules.pdf")
plt
dev.off()
system("slack -F predictedid.heatmap.modules.pdf ryan_todo")

orgo_cirm43[['Bhaduri_modules']] <- CreateAssayObject(data = t(orgo_cirm43_wip@meta.data[grepl("module_",colnames(orgo_cirm43_wip@meta.data))]))
orgo_cirm43@meta.data<-orgo_cirm43_wip@meta.data[!grepl("module_",colnames(orgo_cirm43_wip@meta.data))]


saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")


```

### Transcription Factor Modules

https://www.cell.com/neuron/pdf/S0896-6273(19)30561-6.pdf Defined waves of transcription factors active is neocoritcal mid-gestation.
Supplementary table 8 details this. We looked at these in both gene activity and transcription factor motif accessibility where available.
The regulon is defined by the transcription factor (ChromVar Motif Score) and acts on the listed genes (Gene Activity Score). We can then relate these regulons within our subclusters to that of the human mid-gestational data (Supp table 9).


```R

  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(parallel) 
  library(zoo)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

  #TF Modules
  tf_modules<-read.csv("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/PUBMED31303374.SuppTable8.tsv",sep="\t",head=T)
  tf_modules<-tf_modules[,2:ncol(tf_modules)]
  modules<-lapply(1:ncol(tf_modules),function(x) unlist(tf_modules[x]))
  modules<-lapply(modules,function(x) x[x !=""])
  names(modules)<-colnames(tf_modules)
  #Unname genes in the vectors
  modules<-lapply(modules,unname)
  orgo_cirm43<-AddModuleScore(orgo_cirm43,features=modules,assay="GeneActivity",name=paste0("TF_",names(modules),"_"), search=T)
  #Save them as a separate object rather than meta data
  orgo_cirm43[['TF_modules']] <- CreateAssayObject(data = t(orgo_cirm43@meta.data[grepl("TF_",colnames(orgo_cirm43@meta.data))]))
  orgo_cirm43@meta.data<-orgo_cirm43@meta.data[!grepl("TF_",colnames(orgo_cirm43@meta.data))]
  saveRDS(orgo_cirm43,file="orgo_cirm43.QC.SeuratObject.Rds")


  #Generate Heatmap of TF ChromVar score and TF modules
  modules<-as.data.frame(t(orgo_cirm43@assays$TF_modules@data))
  module_tfs<-unlist(lapply(strsplit(colnames(modules),"-"),"[",2)) #set up TF names

  tf_chrom<-orgo_cirm43@assays$chromvar@data #set up chromvar matrix
  tfList <- getMatrixByID(JASPAR2020, ID=row.names(tf_chrom)) #correct names
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  row.names(tf_chrom)<-tfList
  tfList<-tfList[tfList %in% module_tfs]
  tf_chrom<-as.data.frame(t(tf_chrom[tfList,]))

  modules<-modules[module_tfs %in% colnames(tf_chrom)] #filter modules to chromvar matrix
  colnames(modules)<-unlist(lapply(strsplit(colnames(modules),"-"),"[",2))

  #Combine over subclusters
  tf_chrom<-split(tf_chrom,orgo_cirm43$seurat_clusters) #group by rows to seurat clusters
  tf_chrom<-lapply(tf_chrom,function(x) apply(x,2,median)) #take average across group
  tf_chrom<-do.call("rbind",tf_chrom) #condense to smaller data frame

  modules<-split(modules,paste0(orgo_cirm43$seurat_clusters)) #group by rows to seurat clusters
  modules<-lapply(modules,function(x) apply(x,2,median)) #take average across group
  modules<-do.call("rbind",modules) #condense to smaller data frame

  modules<-t(scale(modules))
  tf_chrom<-t(scale(tf_chrom))

  #This col order to be checked
    clus_order<-c("5","4","0","3","2","1","6")

  tf_chrom<-tf_chrom[,match(clus_order,colnames(tf_chrom),nomatch=0)]
  modules<-modules[,match(clus_order,colnames(modules),nomatch=0)]

  plt1<-Heatmap(tf_chrom,column_order=1:ncol(tf_chrom))
  plt2<-Heatmap(modules,column_order=1:ncol(modules))
  pdf("tf.complexheatmap.pdf",height=20)
  plt1+plt2
  dev.off()
  system("slack -F tf.complexheatmap.pdf ryan_todo")

  #Generate a tanglegram to look at gene opening before or after chromvar activity
  library(ggdendro)
  library(dendextend)
  tf_chrom_dend <- tf_chrom %>% dist("maximum") %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 5)
  modules_dend <- modules %>% dist("maximum") %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 5)
  tang<-tanglegram(tf_chrom_dend,modules_dend) %>% untangle(method = "ladderize")
  pdf("tf.tangle.pdf",width=20)
  tang %>% plot(main = paste("entanglement =", round(entanglement(tang), 2)))
  dev.off()
  system("slack -F tf.tangle.pdf ryan_todo")



```


### Marker Plotting Function

```R
library(Signac)
library(Seurat)
library(ggplot2)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")
  DefaultAssay(orgo_cirm43)<-"peaks"
  plt<-FeaturePlot(orgo_cirm43,feat=c("PAX6",'SOX2',"EOMES","PPP1R17","TBR1","NEUROD1"),col=c("white","red"),ncol=1,order=T)
  ggsave(plt,file="markers.featureplot.pdf",height=50,limitsize=F)
  system("slack -F markers.featureplot.pdf ryan_todo")

  plt<-CoveragePlot(orgo_cirm43,region=c("PAX6",'SOX2',"EOMES","PPP1R17","TBR1","NEUROD1"),extend.upstream=2000,extend.downstream=2000,show.bulk=T,group.by="seurat_clusters",ncol=1)
  ggsave(plt,file="markers.coverageplot.pdf",height=50,limitsize=F)
  system("slack -F markers.coverageplot.pdf ryan_todo")

```


## Public ATAC Comparison
Using the descartes data set of human chromatin accessibility (sciATAC in GW14-15)
https://descartes.brotmanbaty.org/bbi/human-chromatin-during-development/dataset/cerebrum

```bash
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://atlas.fredhutch.org/data/bbi/descartes/human_atac/downloads/cerebrum_filtered.seurat.for_website.RDS
```
## Integrate data 
Have to convert hg19 to hg38 genomic locations
Following https://satijalab.org/signac/articles/integration.html?q=liftover#preprocessing


```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(parallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)

  hg19_to_hg38<-rtracklayer::import.chain("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/hg19ToHg38.over.chain")
  fetal_atac<-readRDS("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/cerebrum_filtered.seurat.for_website.RDS")
  plt<-DimPlot(fetal_atac,group.by=c("cell_type","sex","day_of_pregnancy"))
  ggsave(plt,file="fetal_umap.pdf",width=30)
  system('slack -F fetal_umap.pdf ryan_todo')


  #from this filtering out cells to just RG and excitatory neurons
  fetal_atac<-subset(fetal_atac,cells=colnames(fetal_atac)[which(fetal_atac$cell_type %in% c("Excitatory neurons","Cerebrum_Unknown.3"))])

  plt<-FeaturePlot(fetal_atac,features=c("FOXG1",#forebrain
    "SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3", #RG
    "MOXD1", #vRG
    "NEUROG1","EOMES","PPP1R17","PTPRZ1","NEUROD4",
    "SLC17A7","NEUROD6","RIT1","TBR1","SLA","NHLH1",
    "DLX2","DLX1","LHX6","GAD1",
    "NEUROD1",
    "SATB2","CUX1","RELN","MEF2C",
    "SLC1A3","GLI3","MKI67","STMN2","GADD45B","MAP2","CALB2","FBXO32",
    "PGK1","GORASP2" #stress
    ),order=T,min.cutoff="q25",cols=c("lightgrey","red"))
  ggsave(plt,file="fetal_umap.markers.pdf",width=30,height=30)
  system('slack -F fetal_umap.markers.pdf ryan_todo')

  saveRDS(fetal_atac,"/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/cerebrum_subset.RDS")
  fetal_atac<-readRDS("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/cerebrum_subset.RDS")
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")


  sci_peaks_hg19 <- StringToGRanges(regions = rownames(fetal_atac@assays$peaks@counts), sep = c("-", "-"))
  sci_peaks_hg38 <- rtracklayer::liftOver(x = sci_peaks_hg19, chain = hg19_to_hg38)
  names(sci_peaks_hg38) <- rownames(fetal_atac@assays$peaks@counts)

  # discard any peaks that were mapped to >1 region in hg38
  correspondence <- elementNROWS(sci_peaks_hg38)
  sci_peaks_hg38 <- sci_peaks_hg38[correspondence == 1]
  sci_peaks_hg38 <- unlist(sci_peaks_hg38)
  fetal_atac@assays$peaks@counts <- fetal_atac@assays$peaks@counts[names(sci_peaks_hg38), ]
  # rename peaks with hg38
  rownames(fetal_atac@assays$peaks@counts) <- GRangesToString(grange = sci_peaks_hg38)
  DefaultAssay(fetal_atac)<-"peaks"
  DefaultAssay(orgo_cirm43)<-"peaks"
  fetal_atac[["peaks"]]<-CreateChromatinAssay(count=fetal_atac@assays$peaks@counts,genome="hg38")

  # find intersecting regions
  intersecting.regions <- Signac::findOverlaps(query = fetal_atac, subject = orgo_cirm43)
  intersecting.regions<-intersecting.regions[which(!duplicated(intersecting.regions@from) & !duplicated(intersecting.regions@to)),] # 70583 peaks uniquely shared
  intersections.fetal <- unique(queryHits(intersecting.regions))
  intersections.orgo <- unique(subjectHits(intersecting.regions))
  fetal_atac[["shared_peaks"]]<-CreateChromatinAssay(count=fetal_atac@assays$peaks@counts[intersections.fetal,],genome="hg38")
  orgo_cirm43[["shared_peaks"]]<-CreateChromatinAssay(count=orgo_cirm43@assays$peaks@counts[intersections.orgo,],genome="hg38")

fetal_tomerge<-CreateSeuratObject(fetal_atac[["shared_peaks"]],assay="ATAC",metadata=fetal_atac@meta.data)
fetal_tomerge$dataset<-"primary"
orgo_tomerge<-CreateSeuratObject(orgo_cirm43[["shared_peaks"]],assay="ATAC",metadata=orgo_cirm43@meta.data)
orgo_tomerge$dataset<-"orgo"

combined <- merge(
  x = fetal_tomerge,
  y = orgo_tomerge,
  add.cell.ids = c("fetal", "orgo")
)

saveRDS(combined,"orgo_primary.integration.SeuratObject.Rds")
```

### Integration: Adding chromVAR scores to fetal atlas for comparison

Limiting Chromvar analysis to shared peaks across public data set and our own.

```R
 library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  set.seed(1234)

  #lowerign cores to be used by chromvar to 10
  library(BiocParallel)
  register(MulticoreParam(5))
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  combined<-readRDS(file="orgo_primary.integration.SeuratObject.Rds")

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(x = JASPAR2020, opts = list(species =9606, all_versions = FALSE))

  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  motif.matrix <- CreateMotifMatrix(features = granges(combined), pwm = pfm, genome = 'hg38', use.counts = FALSE)

  # Create a new Mofif object to store the results
  motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)

  DefaultAssay(combined)<-"ATAC"
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(x = JASPAR2020, opts = list(species =9606, all_versions = FALSE))
  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  motif.matrix <- CreateMotifMatrix(features = granges(combined), pwm = pfm, genome = 'hg38', use.counts = FALSE)

  # Create a new Mofif object to store the results
  motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)

  # Add the Motif object to the assays and run ChromVar
  combined <- SetAssayData(object = combined, assay = 'ATAC', slot = 'motifs', new.data = motif)
  combined <- RegionStats(object = combined, genome = BSgenome.Hsapiens.UCSC.hg38)
  combined <- RunChromVAR( object = combined,genome = BSgenome.Hsapiens.UCSC.hg38)

  saveRDS(combined,file="orgo_primary.integration.SeuratObject.Rds")


```

### Integration: Differetial chromvar TF and peaks used 

```R

  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  set.seed(1234)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  combined<-readRDS(file="orgo_primary.integration.SeuratObject.Rds")
  #read in other RDS files for meta data
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  fetal_atac<-readRDS("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/cerebrum_subset.RDS")
  orgo_met<-as.data.frame(orgo_cirm43@meta.data)
  fetal_met<-as.data.frame(fetal_atac@meta.data)
  row.names(orgo_met)<-paste("orgo",row.names(orgo_met),sep="_")
  row.names(fetal_met)<-paste("fetal",row.names(fetal_met),sep="_")
  combined<-AddMetaData(combined,metadata=orgo_met)
  combined<-AddMetaData(combined,metadata=fetal_met)

  combined$assay_cluster<-paste(combined$cell_type,combined$seurat_clusters)  
  saveRDS(combined,file="orgo_primary.integration.SeuratObject.Rds")

  da_one_v_one<-function(i,obj,group,j_list,assay.="ATAC",logfc.threshold.=0.25){
      i<-as.character(i)
      da_tmp_2<-list()
      for (j in j_list){
          if ( i != j){
          if (substr(i,0,3)!=substr(j,0,3)){
          da_peaks_tmp <- FindMarkers(
              object = obj,
              ident.1 = i,
              ident.2 = j,
              group.by = group,
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              only.pos=T,
              assay=assay.,
              logfc.threshold=logfc.threshold.
              )
          da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
          if(assay.=="ATAC"){
          #closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
          #da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
          }
          da_peaks_tmp$enriched_group<-c(i)
          da_peaks_tmp$compared_group<-c(j)
          da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
          }
      }
      }
      return(da_tmp_2)
    }


  #set up an empty list for looping through
  da_peaks<-list()

  #Perform parallel application of DA test
  library(parallel)

  n.cores=8
  da_peaks<-mclapply(
      unique(combined$assay_cluster),
      FUN=da_one_v_one,
      obj=combined,
      group="assay_cluster",
      j_list=unique(combined$assay_cluster),
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_peaks<-do.call("rbind",do.call("rbind",da_peaks))
  i<-"primary_orgo.combined"
  write.table(da_peaks,file=paste0(i,".onevone.da_peaks.txt"),sep="\t",col.names=T,row.names=T,quote=F)
 
  n.cores=1
  da_tf<-mclapply(
      unique(combined$assay_cluster),
      FUN=da_one_v_one,
      obj=combined,
      group="assay_cluster",
      j_list=unique(combined$assay_cluster),
      assay.="chromvar",
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_tf<-do.call("rbind",do.call("rbind",da_tf))
  i<-"primary_orgo.combined"
  write.table(da_tf,file=paste0(i,".onevone.da_tf.txt"),sep="\t",col.names=T,row.names=T,quote=F)
 

```

### Integration: Now Clustering together with cistopic and using harmony to integrate



```R
library(harmony,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/lib_backup_210125")
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(Matrix)

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

combined<-readRDS("orgo_primary.integration.SeuratObject.Rds")
combined_cistopic_counts_frmt<-combined$ATAC@counts

#renaming row names to fit granges expectation of format
row.names(combined_cistopic_counts_frmt)<-sub("-", ":", row.names(combined_cistopic_counts_frmt))
combined_atac_cistopic<-cisTopic::createcisTopicObject(combined_cistopic_counts_frmt) #set up CisTopicObjects
combined_atac_cistopic_models<-cisTopic::runWarpLDAModels(combined_atac_cistopic,topic=c(22,24,26,28,30),nCores=5,addModels=FALSE) #Run warp LDA on objects
saveRDS(combined_atac_cistopic_models,file="orgo_primary.integration.CisTopicObject.Rds") #Saving all models for posterity


#Setting up topic count selection
pdf("orgo_primary.integration.model_selection.pdf")
par(mfrow=c(3,3))
combined_atac_cistopic_models <- selectModel(combined_atac_cistopic_models, type='derivative')
dev.off()
system("slack -F orgo_primary.integration.model_selection.pdf ryan_todo")

cisTopicObject<-cisTopic::selectModel(combined_atac_cistopic_models,type='derivative',keepModels=F)

#saving model selected RDS
saveRDS(cisTopicObject,file="orgo_primary.integration.CisTopicObject.Rds")
cisTopicObject<-readRDS(file="orgo_primary.integration.CisTopicObject.Rds")

#run UMAP on topics
topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
row.names(dims)<-colnames(topic_df)
colnames(dims)<-c("x","y")
dims$cellID<-row.names(dims)
dims<-merge(dims,combined@meta.data,by.x="cellID",by.y="row.names")

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
cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="ATAC",key="UMAP_")
combined@reductions$cistopic<-cistopic_obj
combined@reductions$umap<-cistopic_umap

n_topics<-ncol(Embeddings(combined,reduction="cistopic"))

combined[["peaks"]] <- combined[["ATAC"]] # just renaming to use default settings
DefaultAssay(combined)<-"peaks"

combined <- FindNeighbors(
  object = combined ,
  reduction = 'cistopic',
  dims = 1:n_topics)

combined <- FindClusters(
  object = combined,
  verbose = TRUE,
  resolution=0.01)

#Correcting bias with harmony
pdf("combined.convergence.pdf")
harm_mat<-HarmonyMatrix(combined@reductions$cistopic@cell.embeddings, combined@meta.data$dataset,do_pca=FALSE,nclust=4,plot_convergence=T)
dev.off();system("slack -F combined.convergence.pdf ryan_todo")
combined@reductions$harmony<-CreateDimReducObject(embeddings=as.matrix(harm_mat),assay="peaks",key="topic_")
combined<-RunUMAP(combined, reduction = "harmony",dims=2:ncol(combined@reductions$harmony))
#combined <- FindNeighbors(object = combined,reduction = 'harmony')
#combined <- FindClusters(object = combined,verbose = TRUE,resolution=0.1)


plt<-DimPlot(combined,group.by=c("dataset","seurat_clusters"))
ggsave(plt,file="combined.umap.pdf")
system("slack -F combined.umap.pdf ryan_todo")

saveRDS(combined,"orgo_primary.integration.SeuratObject.Rds")

#Doesnt really integrate but we can still look at diferences across the assays


```


## Differential Gene Activity through Clusters


```R
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(Signac)
library(ggplot2)
library(ComplexHeatmap)
library(TFBSTools)
library(ggdendro)
library(dendextend)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(circlize)

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group,assay.="GeneActivity"){
      da_ga_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T,
        assay=assay.,
        logfc.threshold=0
        )
    da_ga_tmp$da_region<-row.names(da_ga_tmp)
    da_ga_tmp$enriched_group<-c(i)
    da_ga_tmp$compared_group<-c("all_other_cells")
    return(da_ga_tmp)
  }

#peaks
da_peaks_df<-lapply(unique(orgo_cirm43$seurat_clusters),function(x) da_one_v_rest(x,obj=orgo_cirm43,group="seurat_clusters",assay.="peaks"))
da_peaks_df<-do.call("rbind",da_peaks_df)
write.table(da_peaks_df,file="orgo_cirm43.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)


#gene activity
da_ga_df<-lapply(unique(orgo_cirm43$seurat_clusters),function(x) da_one_v_rest(x,obj=orgo_cirm43,group="seurat_clusters",assay.="GeneActivity"))
da_ga_df<-do.call("rbind",da_ga_df)
write.table(da_ga_df,file="orgo_cirm43.onevrest.da_ga.txt",sep="\t",col.names=T,row.names=T,quote=F)

#now doing chromvar
da_ga_df<-lapply(unique(orgo_cirm43$seurat_clusters),function(x) da_one_v_rest(x,obj=orgo_cirm43,group="seurat_clusters",assay.="chromvar"))
da_ga_df<-do.call("rbind",da_ga_df)
write.table(da_ga_df,file="orgo_cirm43.onevrest.da_chromvar.txt",sep="\t",col.names=T,row.names=T,quote=F)

#now doing bhaduri
da_ga_df<-lapply(unique(orgo_cirm43$seurat_clusters),function(x) da_one_v_rest(x,obj=orgo_cirm43,group="seurat_clusters",assay.="Bhaduri_modules"))
da_ga_df<-do.call("rbind",da_ga_df)
write.table(da_ga_df,file="orgo_cirm43.onevrest.da_eigengenes.txt",sep="\t",col.names=T,row.names=T,quote=F)

#now doing tf modules
da_ga_df<-lapply(unique(orgo_cirm43$seurat_clusters),function(x) da_one_v_rest(x,obj=orgo_cirm43,group="seurat_clusters",assay.="TF_modules"))
da_ga_df<-do.call("rbind",da_ga_df)
write.table(da_ga_df,file="orgo_cirm43.onevrest.da_tfmodules.txt",sep="\t",col.names=T,row.names=T,quote=F)


```

## Differential Accessibility between Clusters (one v one)

<!-- To run

```R

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
#  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
#  library(parallel)
  library(ggplot2)
#  library(ggrepel)
  library(dplyr)
#  library(JASPAR2020)
#  library(TFBSTools)
#  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
#  library(parallel) 
#  library(zoo)
#  library(reshape2)
#  library(ape)
#  library(ggdendro)
#  library(dendextend)
  library(circlize)
#  library(dendsort)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)


  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

  da_one_v_one<-function(i,obj,group,j_list,assay.="peaks",logfc.threshold.=0.25){
      i<-as.character(i)
      da_tmp_2<-list()
      for (j in j_list){
          if ( i != j){
          da_peaks_tmp <- FindMarkers(
              object = obj,
              ident.1 = i,
              ident.2 = j,
              group.by = group,
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              only.pos=T,
              assay=assay.,
              logfc.threshold=logfc.threshold.
              )
          da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
          if(assay.=="peaks"){
          closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
          da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
          }
          da_peaks_tmp$enriched_group<-c(i)
          da_peaks_tmp$compared_group<-c(j)
          da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
          }
      }
      return(da_tmp_2)
    }


  #set up an empty list for looping through
  da_peaks<-list()

  #Perform parallel application of DA test
  library(parallel)

  n.cores=8
  da_peaks<-mclapply(
      unique(orgo_cirm43$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      j_list=unique(orgo_cirm43$seurat_clusters),
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_peaks<-do.call("rbind",do.call("rbind",da_peaks))
  i<-"orgo_cirm43"
  write.table(da_peaks,file=paste0(i,".onevone.da_peaks.txt"),sep="\t",col.names=T,row.names=T,quote=F)


  da_tf<-list()
  da_tf<-mclapply(
      unique(orgo_cirm43$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      assay.="chromvar",
      j_list=unique(orgo_cirm43$seurat_clusters),
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_tf<-do.call("rbind",do.call("rbind",da_tf))
  da_tf$da_tf <- unlist(lapply(unlist(lapply(da_tf$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(da_tf,file=paste0(i,".onevone.da_tf.txt"),sep="\t",col.names=T,row.names=T,quote=F)

  da_ga<-list()
  da_ga<-mclapply(
      unique(orgo_cirm43$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      assay.="GeneActivity",
      j_list=unique(orgo_cirm43$seurat_clusters),
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_ga<-do.call("rbind",do.call("rbind",da_ga))
  write.table(da_ga,file=paste0(i,".onevone.ga_tf.txt"),sep="\t",col.names=T,row.names=T,quote=F)

  da_tfmod<-list()
  da_tfmod<-mclapply(
      unique(orgo_cirm43$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      assay.="TF_modules",
      j_list=unique(orgo_cirm43$seurat_clusters),
      logfc.threshold.=0,
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_tfmod<-do.call("rbind",do.call("rbind",da_tfmod))
  write.table(da_tfmod,file=paste0(i,".onevone.tfmodule_tf.txt"),sep="\t",col.names=T,row.names=T,quote=F)


  da_module<-list()
  da_module<-mclapply(
      unique(orgo_cirm43$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      assay.="Bhaduri_modules",
      j_list=unique(orgo_cirm43$seurat_clusters),
      logfc.threshold.=0,
      mc.cores=n.cores)
  #Merge the final data frame from the list for 1vrest DA
  da_module<-do.call("rbind",do.call("rbind",da_module))
  write.table(da_module,file=paste0(i,".onevone.tfmodules.txt"),sep="\t",col.names=T,row.names=T,quote=F)


####################### PLOTTING ##########################################
#Plot out top ga for each cluster
da_ga<-read.csv(file="orgo_cirm43.onevone.ga_tf.txt",head=T,sep="\t",row.names=NULL)
da_ga$gene_name<-da_ga$da_region
da_ga<-da_ga[complete.cases(da_ga),]

da_ga$label<-""
for (x in unique(da_ga$enriched_group)){
selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x)  %>% filter(p_val < 0.05) %>% arrange(rev(desc(p_val))) %>% slice(1:10))$da_region
da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_ga<-as.data.frame(t(as.data.frame(orgo_cirm43[["GeneActivity"]]@data)))
sum_ga<-split(dat_ga,orgo_cirm43$seurat_clusters) #group by rows to seurat clusters
sum_ga<-lapply(sum_ga,function(x) apply(x,2,mean)) #take average across group
sum_ga<-do.call("rbind",sum_ga) #condense to smaller data frame

sum_ga<-t(scale(sum_ga))
sum_ga<-sum_ga[row.names(sum_ga) %in% unique(da_ga$label),]

#cluster by all marker genes
sum_da_dend <- t(sum_ga) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="orgo_cirm43.geneactivity.dend.rds") 

sum_ga_plot<-na.omit(t(sum_ga))
clus_order<-c("7","6","0","3","1","2","5","4")
sum_ga_plot<-sum_ga_plot[clus_order,]
colfun=colorRamp2(quantile(unlist(sum_ga_plot), probs=c(0.5,0.75,0.9),na.rm=T),magma(3))
plt1<-Heatmap(sum_ga_plot,
    #cluster_rows=sum_da_dend,
    #left_annotation=side_ha,
    col=colfun,
    row_order=1:nrow(sum_ga_plot),
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

pdf("orgo_cirm43.geneactivity.heatmap.pdf",height=20,width=20)
plt1
dev.off()
system("slack -F orgo_cirm43.geneactivity.heatmap.pdf ryan_todo")


########################Plot out top TF for each cluster###################
da_tf<-read.csv(file="orgo_cirm43.onevone.da_tf.txt",head=T,sep="\t",row.names=NULL)
da_tf$tf_name<-da_tf$da_tf
da_tf<-da_tf[complete.cases(da_tf),]

da_tf$label<-""
for (x in unique(da_tf$enriched_group)){
selc_genes<-as.data.frame(da_tf %>% filter(enriched_group==x)  %>% filter(p_val < 0.05) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:10))$tf_name
da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$label<- da_tf[da_tf$tf_name %in% selc_genes & da_tf$enriched_group==x,]$tf_name
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_tf<-as.data.frame(t(as.data.frame(orgo_cirm43[["chromvar"]]@data)))
sum_tf<-split(dat_tf,orgo_cirm43$seurat_clusters) #group by rows to seurat clusters
sum_tf<-lapply(sum_tf,function(x) apply(x,2,mean)) #take average across group
sum_tf<-do.call("rbind",sum_tf) #condense to smaller data frame

sum_tf<-t(scale(sum_tf))
sum_tf<-sum_tf[,!endsWith(colnames(sum_tf),"NA")] #remove NA (doublet cells)

#clustered by all marker genes ga
sum_da_dend<-readRDS(file="orgo_cirm43.geneactivity.dend.rds") 


sum_tf<-sum_tf[row.names(sum_tf) %in% unique(da_tf[da_tf$label!="",]$da_region),]
row.names(sum_tf)<-da_tf[match(row.names(sum_tf),da_tf$da_region,nomatch=0),]$tf_name
sum_tf_plot<-na.omit(t(sum_tf))
clus_order<-c("7","6","0","3","1","2","5","4")
sum_tf_plot<-sum_tf_plot[clus_order,]

colfun=colorRamp2(quantile(unlist(sum_tf_plot), probs=c(0.5,0.75,0.95)),cividis(3))
plt1<-Heatmap(sum_tf_plot,
    #cluster_rows=sum_da_dend,
    #left_annotation=side_ha,
    col=colfun,
    row_order=1:nrow(sum_tf_plot),
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

plt1<-draw(plt1)


pdf("orgo_cirm43.tf.heatmap.pdf",height=20,width=20)
draw(plt1)
dev.off()
system("slack -F orgo_cirm43.tf.heatmap.pdf ryan_todo")

#Plot motifs alongside chromvar plot
library(ggplot2)
library(patchwork)

motif_order<-names(orgo_cirm43@assays$peaks@motifs@motif.names[match(colnames(sum_tf_plot)[column_order(plt1)],unlist(orgo_cirm43@assays$peaks@motifs@motif.names),nomatch=0)])
plt<-MotifPlot(object = orgo_cirm43,motifs = motif_order,ncol=1)+theme_void()+theme(strip.text = element_blank())

ggsave(plt,file="orgo_cirm43.tf.heatmap.motif.pdf",height=100,width=2,limitsize=F)
system("slack -F orgo_cirm43.tf.heatmap.motif.pdf ryan_todo")


#Plot out top tf modules for each cluster
da_tfmod<-read.csv(file="orgo_cirm43.onevone.tfmodule_tf.txt",head=T,sep="\t",row.names=NULL)
da_tfmod$gene_name<-unlist(lapply(strsplit(da_tfmod$da_region,"-"),"[",2))
da_tfmod<-da_tfmod[complete.cases(da_tfmod),]

da_tfmod$label<-""
for (x in unique(da_tfmod$enriched_group)){
selc_genes<-as.data.frame(da_tfmod %>% filter(enriched_group==x) %>% filter(p_val < 0.05) %>% arrange(rev(desc(p_val))) %>% slice(1:10))$da_region
da_tfmod[da_tfmod$da_region %in% selc_genes & da_tfmod$enriched_group==x,]$label<- da_tfmod[da_tfmod$da_region %in% selc_genes & da_tfmod$enriched_group==x,]$da_region
}

#Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
dat_tfmod<-as.data.frame(t(as.data.frame(orgo_cirm43[["TF_modules"]]@data)))
sum_tfmod<-split(dat_tfmod,orgo_cirm43$seurat_clusters) #group by rows to seurat clusters
sum_tfmod<-lapply(sum_tfmod,function(x) apply(x,2,mean)) #take average across group
sum_tfmod<-do.call("rbind",sum_tfmod) #condense to smaller data frame

sum_tfmod<-t(scale(sum_tfmod))
sum_tfmod<-sum_tfmod[row.names(sum_tfmod) %in% unique(da_tfmod$label),]

#cluster by all marker genes
sum_da_dend <- t(sum_tfmod) %>% dist() %>% hclust %>% as.dendrogram %>% ladderize  %>% set("branches_k_color", k = 8)
saveRDS(sum_da_dend,file="orgo_cirm43.tfmodule.dend.rds") 

sum_tfmod_plot<-na.omit(t(sum_tfmod))
clus_order<-c("7","6","0","3","1","2","5","4")
sum_tfmod_plot<-sum_tfmod_plot[clus_order,]
colfun=colorRamp2(quantile(unlist(sum_tfmod_plot), probs=c(0.5,0.75,0.9),na.rm=T),c("#253494","#41b6c4","#ffffd9"))
plt1<-Heatmap(sum_tfmod_plot,
    #cluster_rows=sum_da_dend,
    #left_annotation=side_ha,
    col=colfun,
    row_order=1:nrow(sum_tfmod_plot),
    #bottom_annotation=bottom_ha,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp=gpar(fontsize=7),
    column_names_rot=90
)

pdf("orgo_cirm43.tfmodule.heatmap.pdf",height=20,width=20)
plt1
dev.off()
system("slack -F orgo_cirm43.tfmodule.heatmap.pdf ryan_todo")


```

-->

## Performing GREAT on DA peaks


```R
  #mkdir GREAT_analysis

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(GenomicRanges)

  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")

  #To perform GREAT on peaks for enrichment per cluster
  write("Performing GREAT on all enriched sites per annotation group", stderr())
  library(rGREAT)

  #format data as bed file all seurat objects have the same peak list
  write("Preparing Background Set as all called peaks.", stderr())
  orgo_bg_bed<-do.call("rbind",strsplit(unlist(orgo_cirm43@assays$peaks@counts@Dimnames[1]),"[-]"))
  orgo_bg_bed<-as.data.frame(orgo_bg_bed)
  colnames(orgo_bg_bed)<-c("chr","start","end")
  orgo_bg_bed$start<-as.numeric(as.character(orgo_bg_bed$start))
  orgo_bg_bed$end<-as.numeric(as.character(orgo_bg_bed$end))
  orgo_bg_bed<-makeGRangesFromDataFrame(orgo_bg_bed)
  cirm43_da_peaks<-read.table("orgo_cirm43.onevrest.da_peaks.txt",header=T)

  dir.create("./GREAT_analysis")

  write("Beginning loop through all annotation groups.", stderr())

  great_processing<-function(enriched_group_input,peak_dataframe,prefix,bg){
      #subset bed file to peaks enriched in input group
      orgo_bed<-as.data.frame(do.call("rbind",strsplit(peak_dataframe[peak_dataframe$enriched_group==enriched_group_input,]$da_region,"-")))
      colnames(orgo_bed)<-c("chr","start","end")
      orgo_bed$start<-as.numeric(as.character(orgo_bed$start))
      orgo_bed$end<-as.numeric(as.character(orgo_bed$end))
      nrow(orgo_bed)
      orgo_bed<-orgo_bed[!duplicated(orgo_bed),]
      row_count<-nrow(orgo_bed)
      orgo_bed$width<-orgo_bed$end-orgo_bed$start
      orgo_bed<-makeGRangesFromDataFrame(orgo_bed)

      #run GREAT using all peaks as background
      write(paste("Using",row_count, "DA peaks from",enriched_group_input), stderr())
      job = submitGreatJob(orgo_bed,bg=bg,species="hg38",request_interval=30)
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
  mclapply(unique(cirm43_da_peaks$enriched_group), FUN=function(x){great_processing(enriched_group_input=x,peak_dataframe=cirm43_da_peaks,prefix="cirm43",bg=orgo_bg_bed)},mc.cores=7)

```


## Check Accelerated Regions
```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Signac)
  library(GenomicRanges)
  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds")
  har<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/nCHAR.txt",header=T)  
  length(unique(findOverlaps(query = makeGRangesFromDataFrame(har), subject = orgo_cirm43@assays$peaks@ranges)))
  #154/2745 #all human accelerated regions

  har<-har[har$enhancer_activity=="Yes",]

  length(unique(findOverlaps(query = makeGRangesFromDataFrame(har), subject = orgo_cirm43@assays$peaks@ranges)))
  #55/773 #predicted enhancer activity

  har<-har[!(har$predicted_tissue_activity %in% c("N/A","limb","heart;limb","other")),]
  length(unique(findOverlaps(query = makeGRangesFromDataFrame(har), subject = orgo_cirm43@assays$peaks@ranges)))
  #18/288 #predicted brain/neural tube activity


```

## Monocle

Running radial glia (RG) like and excitatory neuron like (ExN) cells serparately. This is to build a more biologically meaninful trajectory for cells.

### Recluster clusters 5, 4, 0, 3
```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cisTopic)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds") #reading in QC passing cells
  dat<-subset(orgo_cirm43,seurat_clusters %in% c("5","4","0","3"))

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(20:30),nCores=11,addModels=FALSE)

      #Setting up topic count selection
      pdf(paste(prefix,"model_selection.pdf",sep="."))
      par(mfrow=c(1,3))
      cirm43_cistopic_models <- selectModel(atac_cistopic_models, type='derivative')
      dev.off()
      system(paste0("slack -F ",paste(prefix,"model_selection.pdf",sep=".")," ryan_todo"))
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=dat,prefix="orgo_cirm43.RG")
  cistopic_models<-readRDS("orgo_cirm43.RG.CisTopicObject.Rds")

  #set topics based on derivative
  cisTopicObject<-cisTopic::selectModel(cistopic_models,type="derivative",keepModels=T)

  #saving model selected RDS
  saveRDS(cisTopicObject,file="orgo_cirm43.RG.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


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
        resolution=resolution
      )

  return(object_input)}
  DefaultAssay(dat)<-"peaks"
  dat<-cistopic_wrapper(object_input=dat,cisTopicObject=cisTopicObject,resolution=0.2)

  plt<-DimPlot(dat,group.by=c("DIV","seurat_clusters"))
  ggsave(plt,file="RG.qc.umap.pdf")
  system("slack -F RG.qc.umap.pdf ryan_todo")

  saveRDS(dat,file="orgo_cirm43.RG.SeuratObject.Rds")   ###save Seurat file

```

### Recluster clusters 2 1 6
```R
 library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cisTopic)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.QC.SeuratObject.Rds") #reading in QC passing cells
  dat<-subset(orgo_cirm43,seurat_clusters %in% c("2","1","6"))

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(20:30),nCores=11,addModels=FALSE)

      #Setting up topic count selection
      pdf(paste(prefix,"model_selection.pdf",sep="."))
      par(mfrow=c(1,3))
      cirm43_cistopic_models <- selectModel(atac_cistopic_models, type='derivative')
      dev.off()
      system(paste0("slack -F ",paste(prefix,"model_selection.pdf",sep=".")," ryan_todo"))
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=dat,prefix="orgo_cirm43.ExN")
  cistopic_models<-readRDS("orgo_cirm43.ExN.CisTopicObject.Rds")

  #set topics based on derivative
  cisTopicObject<-cisTopic::selectModel(cistopic_models,type="derivative",keepModels=T)

  #saving model selected RDS
  saveRDS(cisTopicObject,file="orgo_cirm43.ExN.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


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
        resolution=resolution
      )

  return(object_input)}
  DefaultAssay(dat)<-"peaks"
  dat<-cistopic_wrapper(object_input=dat,cisTopicObject=cisTopicObject,resolution=0.2)

  plt<-DimPlot(dat,group.by=c("DIV","seurat_clusters"))
  ggsave(plt,file="ExN.qc.umap.pdf")
  system("slack -F ExN.qc.umap.pdf ryan_todo")

  saveRDS(dat,file="orgo_cirm43.ExN.SeuratObject.Rds")   ###save Seurat file
```

## Using monocle to build trajectories

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(ggplot2)
  set.seed(1234)
  library(cicero)
  library(SeuratWrappers)
  library(patchwork)
  library(monocle3)


  monocle_processing<-function(outname, object_input,min_branch=10){
      atac.cds <- as.cell_data_set(object_input)
      atac.cds <- cluster_cells(cds = atac.cds, reduction_method = "UMAP",k=10) 
      #Read in cds from cicero processing earlier and continue processing
      atac.cds<- learn_graph(atac.cds, use_partition=F,close_loop=F,learn_graph_control=list(minimal_branch_len=min_branch))
      #plot to see nodes for anchoring
      plt1<-plot_cells(cds = atac.cds, show_trajectory_graph = TRUE, color_cells_by="DIV", label_leaves=T, label_branch_points=F, label_roots=T) 
      plt2<-plot_cells(cds = atac.cds, show_trajectory_graph = TRUE, color_cells_by="seurat_clusters",label_leaves=T, label_branch_points=F, label_roots=T) 
      #Also make a plot of just node names for easier identification
      root_nodes<-as.data.frame(t(atac.cds@principal_graph_aux$UMAP$dp_mst))
      root_nodes$label<-row.names(root_nodes)
      plt3<-ggplot(root_nodes, aes(x=UMAP_1,y=UMAP_2))+ geom_text(aes(label=label),size=3)+ theme_bw() 
      plt<-(plt1+plt2)/plt3
      ggsave(plt,file=paste(outname,"DIV_trajectory.pdf",sep="_"),width=20,height=10)
      #system(paste0("slack -F ",paste(outname,"DIV_trajectory.pdf",sep="_")," ryan_todo"))
      return(atac.cds)
  }


###Excitatory Neurons####
#Rerun cistopic and recluster RG cell subset
exn<-readRDS("orgo_cirm43.ExN.SeuratObject.Rds")
DefaultAssay(exn)<-"peaks"
outname.="ExN"
exn_sub.cds<-monocle_processing(object_input=exn,outname=outname.)
#NOTE THIS SHOULD BE CHANGED AFTER FIRST PLOT
exn_sub.cds <- order_cells(exn_sub.cds, reduction_method = "UMAP", root_pr_nodes = c("Y_184")) #Chose youngest cells as root 
saveRDS(exn_sub.cds ,paste0(outname.,".pseudotime.monoclecds.Rds"))

#plot pseudotime after determining root
pdf(paste0(outname.,".trajectory.pseudotime.pdf"))
plot_cells(cds = exn_sub.cds ,show_trajectory_graph = TRUE,color_cells_by = "pseudotime", alpha=0.5, cell_stroke=0)
dev.off()
#system("slack -F ",paste0(outname.,".trajectory.pseudotime.pdf")," ryan_todo")

#Append pseudotime to meta data of seurat object
exn<-AddMetaData(object=exn,metadata=exn_sub.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name=paste0("pseudotime_",outname.))
plt1<-FeaturePlot(exn,feature=paste0("pseudotime_",outname.))
plt2<-DimPlot(exn,group.by="seurat_clusters")
ggsave(plt1+plt2,file=paste0(outname.,"_trajectory.pseudotime.pdf"))
#system("slack -F ",paste0(outname.,"_trajectory.pseudotime.pdf"))," ryan_todo")
saveRDS(exn,"orgo_cirm43.ExN.pseudotime.SeuratObject.Rds")


### Radial Glia ########
RG<-readRDS(file="orgo_cirm43.RG.SeuratObject.Rds")   ###save Seurat file
DefaultAssay(RG)<-"peaks"
outname.="RG"
RG_sub.cds<-monocle_processing(object_input=RG,outname=outname.)
#NOTE THIS SHOULD BE CHANGED AFTER FIRST PLOT
RG_sub.cds <- order_cells(RG_sub.cds, reduction_method = "UMAP", root_pr_nodes = c("Y_153")) #Chose youngest cells as root 
saveRDS(RG_sub.cds ,paste0(outname.,".pseudotime.monoclecds.Rds"))

#plot pseudotime after determining root
pdf(paste0(outname.,".trajectory.pseudotime.pdf"))
plot_cells(cds = RG_sub.cds ,show_trajectory_graph = TRUE,color_cells_by = "pseudotime", alpha=0.5, cell_stroke=0)
dev.off()
#system("slack -F ",paste0(outname.,".trajectory.pseudotime.pdf")," ryan_todo")

#Append pseudotime to meta data of seurat object
RG<-AddMetaData(object=RG,metadata=RG_sub.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name=paste0("pseudotime_",outname.))
plt1<-FeaturePlot(RG,feature=paste0("pseudotime_",outname.))
plt2<-DimPlot(RG,group.by="seurat_clusters")
ggsave(plt1+plt2,file=paste0(outname.,"_trajectory.pseudotime.pdf"))
#system("slack -F ",paste0(outname.,"_trajectory.pseudotime.pdf"))," ryan_todo")

saveRDS(RG,"orgo_cirm43.RG.pseudotime.SeuratObject.Rds")








```

### Statistical analysis across pseudotime

```R

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(cisTopic)
  library(monocle3)
  library(JASPAR2020)
  library(TFBSTools)
  library(dplyr)

  #Find changes over pseudotime
  pseudotime_statistics<-function(seurat.object=orgo_cirm43,cds=in_cds,prefix="test",ncores=1){
  cds <- cds[,is.finite(pseudotime(cds))] #exclude earliest cells not fitting into trajectory
  chromvar.cds <- as.cell_data_set(seurat.object,assay="chromvar")
  chromvar.cds@assays@data@listData$counts[which(is.na(chromvar.cds@assays@data@listData$counts),arr.ind=T)]<-0
  geneactivity.cds <- as.cell_data_set(seurat.object,assay="GeneActivity")
  tf_module.cds <- as.cell_data_set(seurat.object,assay="TF_modules")

  chromvar.cds  <- chromvar.cds[,row.names(colData(chromvar.cds)) %in% row.names(colData(cds))]
  geneactivity.cds  <- geneactivity.cds[,row.names(colData(geneactivity.cds)) %in% row.names(colData(cds))]
  tf_module.cds  <- tf_module.cds[,row.names(colData(tf_module.cds)) %in% row.names(colData(cds))]

  pData(cds)$Pseudotime <- pseudotime(cds) #add pseudotime to meta data
  pData(chromvar.cds)$Pseudotime <- pseudotime(cds) 
  pData(geneactivity.cds)$Pseudotime <- pseudotime(cds) 
  pData(tf_module.cds)$Pseudotime <- pseudotime(cds) 

  #using peak accessibility projection for all pseudotime analysis
  #colData(cds)$Size_Factor<-colData(cds)$nCount_peaks
  #cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=ncores,reduction_method="UMAP")
  #nrow(cds_pr_test_res %>% filter(q_value<0.1))
  #saveRDS(cds_pr_test_res,paste0(prefix,"_pseudotime.peaks.rds"))

  chromvar.cds@principal_graph$UMAP<-cds@principal_graph$UMAP
  chromvar.cds@principal_graph_aux$UMAP<-cds@principal_graph_aux$UMAP
  reducedDims(chromvar.cds)[["UMAP"]]<-reducedDims(cds)[["UMAP"]]
  colData(chromvar.cds)$Size_Factor<-colData(chromvar.cds)$nFeature_peaks
  chromvar_cds_pr_test_res <- graph_test(chromvar.cds, neighbor_graph="principal_graph", cores=ncores,reduction_method="UMAP")
  nrow(chromvar_cds_pr_test_res %>% filter(q_value<0.1))
  tfList <- getMatrixByID(JASPAR2020, ID=row.names(chromvar_cds_pr_test_res))
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  chromvar_cds_pr_test_res$tf_name<-tfList
  saveRDS(chromvar_cds_pr_test_res,paste0(prefix,"_pseudotime.chromvar.rds"))

  tf_module.cds@principal_graph$UMAP<-cds@principal_graph$UMAP
  tf_module.cds@principal_graph_aux$UMAP<-cds@principal_graph_aux$UMAP
  reducedDims(tf_module.cds)[["UMAP"]]<-reducedDims(cds)[["UMAP"]]
  colData(tf_module.cds)$Size_Factor<-colData(tf_module.cds)$nFeature_peaks
  tf_module_cds_pr_test_res <- graph_test(tf_module.cds, neighbor_graph="principal_graph", cores=ncores,reduction_method="UMAP")
  nrow(tf_module_cds_pr_test_res %>% filter(q_value<0.05))
  saveRDS(tf_module_cds_pr_test_res,paste0(prefix,"_pseudotime.tfmodules.rds"))

  geneactivity.cds@principal_graph$UMAP<-cds@principal_graph$UMAP
  geneactivity.cds@principal_graph_aux$UMAP<-cds@principal_graph_aux$UMAP
  reducedDims(geneactivity.cds)[["UMAP"]]<-reducedDims(cds)[["UMAP"]]
  colData(geneactivity.cds)$Size_Factor<-colData(geneactivity.cds)$nFeature_peaks
  geneactivity_cds_pr_test_res <- graph_test(geneactivity.cds, neighbor_graph="principal_graph", cores=ncores,reduction_method="UMAP")
  nrow(geneactivity_cds_pr_test_res %>% filter(q_value<0.05))
  saveRDS(geneactivity_cds_pr_test_res,paste0(prefix,"_pseudotime.geneactivity.rds"))

}



  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  DefaultAssay(orgo_cirm43)<-"peaks"
  orgo_cirm43@meta.data<-orgo_cirm43@meta.data[sort(row.names(orgo_cirm43@meta.data)),] #as.cell_data_set assumes abc sorted cell ids
  in_cds<-readRDS("cirm43.pseudotime.monoclecds.Rds")
  in_cds<-in_cds[,row.names(colData(in_cds)) %in% row.names(orgo_cirm43@meta.data)] #exclude earliest cells not fitting into trajectory
  pseudotime_statistics(seurat.object=orgo_cirm43,cds=in_cds,prefix="orgo_cirm43",ncores=5)

  # #Need to add TF modules to subset cell objects
  rg.cds<-readRDS("radial_glia.pseudotime.monoclecds.Rds")
  rg<-subset(orgo_cirm43,cells=colnames(rg.cds))
   pseudotime_statistics(seurat.object=rg,cds=rg.cds,prefix="RG",ncores=10)


  # exN<-readRDS("./subcluster/_0_SeuratObject.Rds")
  # exN.cds<-readRDS("./subcluster/excitatory_neuron.pseudotime.monoclecds.Rds")
  # exN <- subset(exN, cells=row.names(orgo_cirm43@meta.data))
  # exN[['TF_modules']] <- subset(orgo_cirm43, cells=row.names(exN@meta.data))[['TF_modules']] 
  # exN.cds<-exN.cds[,row.names(colData(exN.cds)) %in% row.names(exN@meta.data)] #exclude earliest cells not fitting into trajectory
  # pseudotime_statistics(seurat.object=exN,cds=exN.cds,prefix="./subcluster/exN",ncores=5)



  # pData(orgo_cirm43.cds)$pseudotime_bin <- cut(pseudotime(orgo_cirm43.cds), 50) #bin data into multiple bins (10 bins)
  # pData(chromvar.cds)$pseudotime_bin <- cut(pseudotime(orgo_cirm43.cds), 50) 
  # pData(geneactivity.cds)$pseudotime_bin <- cut(pseudotime(orgo_cirm43.cds), 50) 
  # pData(tf_module.cds)$pseudotime_bin <- cut(pseudotime(orgo_cirm43.cds), 50) 

  # binned_peaks.cds <- aggregate_by_cell_bin(orgo_cirm43.cds, "pseudotime_bin")
  # binned_geneactivity.cds <- aggregate_by_cell_bin(geneactivity.cds, "pseudotime_bin")

  # binned_chromvar.cds <- aggregate_mean_by_cell_bin(chromvar.cds, "pseudotime_bin")
  # binned_tfmodule.cds <- aggregate_mean_by_cell_bin(tf_module.cds, "pseudotime_bin")

  # colData(chromvar.cds)$Size_Factor<-colData(chromvar.cds)$nFeature_peaks
  # chromvar_fits <- fit_models(chromvar.cds, model_formula_str = "~Pseudotime + differentiation_exp",expression_family="quasipoisson")
  # fit_coefs <- coefficient_table(chromvar_fits)
  #   out<-fit_coefs %>% filter(term == "Pseudotime") %>% filter(p_value <0.05) %>% select(gene_id, term, p_value, estimate, std_err, normalized_effect)
  # saveRDS(chromvar_fits,"cirm43_chromvar.pseudotime.models.Rds")
  # saveRDS(fit_coefs,"cirm43_chromvar.pseudotime.coefftable.Rds")
  
  # tfmod_fits <- fit_models(binned_peaks.cds, model_formula_str = "~Pseudotime",expression_family="quasipoisson")
  # fit_coefs <- coefficient_table(tfmod_fits)
  # saveRDS(tfmod_fits,"cirm43_tfmod.pseudotime.models.Rds")
  # saveRDS(fit_coefs,"cirm43_tfmod.pseudotime.coefftable.Rds")
  
  # ga_fits <- fit_models(binned_geneactivity.cds, model_formula_str = "~ Pseudotime + nFeature_GeneActivity + differentiation_exp")
  # fit_coefs <- coefficient_table(ga_fits)
  # out<-fit_coefs %>% filter(term == "Pseudotime") %>% filter(p_value <0.05) %>% select(site_name, term, p_value, estimate, std_err, normalized_effect)

  # saveRDS(ga_fits,"cirm43_ga.pseudotime.models.Rds")
  # saveRDS(fit_coefs,"cirm43_ga.pseudotime.coefftable.Rds")

  # acc_fits <- fit_models(binned_peaks.cds, model_formula_str = "~Pseudotime + nFeature_peaks",expression_family="binomial")
  # fit_coefs <- coefficient_table(acc_fits)
  # saveRDS(acc_fits,"cirm43_acc.pseudotime.models.Rds")
  # saveRDS(fit_coefs,"cirm43_acc.pseudotime.coefftable.Rds")



#Include differentiation experiment in model as well

```


### Analysis of ChromVAR TF motifs and Gene Activity Through Pseudotime

Decided to use a binning strategy to assess pseudotime signal.


```R

  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(parallel) 
  library(zoo)
  library(viridis)
  library(dplyr)
  library(reshape2)
  library(ape)
  library(ggdendro)
  library(dendextend)
  library(circlize)
  library(dendsort)
  library(RColorBrewer)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
    atac_sub<-orgo_cirm43
    plt1<-FeaturePlot(atac_sub,feature=c("pseudotime_allcells"),order=T)
    plt2<-FeaturePlot(atac_sub,feature=c("pseudotime_rg"),order=T)
    plt3<-DimPlot(atac_sub,group.by="seurat_clusters")

    ggsave(plt1/plt2/plt3,file="test_orgo_cirm43_trajectory.pseudotime.pdf")
    system("slack -F test_orgo_cirm43_trajectory.pseudotime.pdf ryan_todo")

  cirm43_da_peaks<-read.table("orgo_cirm43.onevone.da_peaks.txt",header=T)
  gene_list<-unique(cirm43_da_peaks$)

  chromvar_motifs_per_bin<-function(x){
  bin_tmp<-row.names(atac_sub@meta.data[atac_sub@meta.data$pseudotime_bin==x,])
  chromvar<-atac_sub@assays$chromvar@data
  return(apply(chromvar[,bin_tmp],1,mean))
  }

  geneactivity_per_bin<-function(x){
  bin_tmp<-row.names(atac_sub@meta.data[atac_sub@meta.data$pseudotime_bin==x,])
  ga<-atac_sub@assays$GeneActivity@data
  return(apply(ga[,bin_tmp],1,mean))
  }

  TF_modules_per_bin<-function(x){
  bin_tmp<-row.names(atac_sub@meta.data[atac_sub@meta.data$pseudotime_bin==x,])
  tf_module<-atac_sub@assays$TF_modules@data
  return(apply(tf_module[,bin_tmp],1,mean))
  }

  #Set up different annotation stacked barplots
  annotation_bin_summary<-function(object=atac_sub,x){
  tmp_timebin<-as.data.frame(object@meta.data %>% group_by_("pseudotime_bin",x) %>% summarize(count=n()))
  colnames(tmp_timebin)<-c("bin","var","count")
  tmp_timebin<-dcast(tmp_timebin,bin~var,value.var="count",fill=0)
  row.names(tmp_timebin)<-tmp_timebin$bin
  tmp_timebin<-as.data.frame(t(tmp_timebin[2:ncol(tmp_timebin)]))
  tmp_timebin = as.data.frame(t(scale(tmp_timebin, center = FALSE, 
               scale = colSums(tmp_timebin))))
  return(tmp_timebin)}

  plot_over_pseudotime<-function(object=orgo_cirm43,pseudotime="pseudotime_cluster0",bincount=25,prefix="test",monocle_prefix){
        atac_sub<-subset(object,cells=row.names(object@meta.data)[which(!is.na(object@meta.data[,which(colnames(object@meta.data)==pseudotime)]))]) #subset to cells with pseudotime info
        atac_sub<-subset(object,cells=row.names(object@meta.data)[which(is.finite(object@meta.data[,which(colnames(object@meta.data)==pseudotime)]))]) #subset to cells with pseudotime info
        atac_sub$pseudotime_bin<-cut_number(atac_sub@meta.data[,which(colnames(atac_sub@meta.data)==pseudotime)],bincount) #bin data into equal group numbers
        atac_sub$working_pseudotime<-atac_sub@meta.data[,which(colnames(atac_sub@meta.data)==pseudotime)] #convenience for later calling
        #Plot bins by cluster ID
        clus_bin<-as.data.frame(atac_sub@meta.data %>% group_by(pseudotime_bin,seurat_clusters) %>% summarize(count=n()))
        clus_bin<-dcast(clus_bin,pseudotime_bin~seurat_clusters,fill=0)
        row.names(clus_bin)<-clus_bin[,1]
        clus_bin<-clus_bin[,2:ncol(clus_bin)]
        clus_bin<-t(clus_bin)
        clus_bin<-clus_bin/rowSums(clus_bin)
        clus_bin<-clus_bin[order(unlist(lapply(1:nrow(clus_bin),function(x) which(clus_bin[x,]==max(clus_bin[x,]))))),]
        #dend <- clus_bin %>% dist("maximum") %>% hclust %>% as.dendrogram %>% ladderize  
        colfun=colorRamp2(quantile(unlist(clus_bin), probs=c(0.01,0.99)),c("#ffffff","#000000"))

        #Set up annotation barplots for pseudotime bins
        div_timebin<-annotation_bin_summary(x="DIV")
        cluster_timebin<-annotation_bin_summary(x="seurat_clusters")
        organoid_timebin<-annotation_bin_summary(x="orgID")
        differentiation_timebin<-annotation_bin_summary(x="differentiation_exp")

        ha_barplots = HeatmapAnnotation(
        DIV = anno_barplot(div_timebin, name="DIV",bar_width=1,gp = gpar(color = 1:nrow(div_timebin),fill = 1:nrow(div_timebin)),border=F),
        cluster = anno_barplot(cluster_timebin, name="cluster",bar_width=1,gp = gpar(color = 1:nrow(cluster_timebin),fill = 1:nrow(cluster_timebin)),border=F),
        organoid = anno_barplot(organoid_timebin, name="organoid",bar_width=1,gp = gpar(color = 1:nrow(organoid_timebin),fill = 1:nrow(organoid_timebin)),border=F),
        differentiation=anno_barplot(differentiation_timebin, name="experiment",bar_width=1,gp = gpar(color = 1:nrow(differentiation_timebin),fill = 1:nrow(differentiation_timebin)),border=F))

        lgd_list = list(
        Legend(labels = colnames(div_timebin), title = "DIV", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(div_timebin))),
        Legend(labels = colnames(cluster_timebin), title = "cluster", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(cluster_timebin))),
        Legend(labels = colnames(organoid_timebin), title = "organoid", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(organoid_timebin))),
        Legend(labels = colnames(differentiation_timebin), title = "experiment", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(differentiation_timebin))))

        plt<-Heatmap(clus_bin,
        column_order=1:ncol(clus_bin),
        bottom_annotation=ha_barplots,
        row_order=1:nrow(clus_bin),
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        show_column_names=T,
        col=colfun,
        show_heatmap_legend=T)

        pdf(paste0(prefix,".seurat_clusters.pseudotime.order.heatmap.pdf"),width=30)
        plt<-draw(plt,annotation_legend_list=lgd_list)
        dev.off()
        system(paste0("slack -F ",prefix,".seurat_clusters.pseudotime.order.heatmap.pdf"," ryan_todo"))


        plt1<-ggplot(dat=as.data.frame(atac_sub@meta.data),aes(x=working_pseudotime,fill=as.factor(pseudotime_bin)))+geom_histogram(bins=bincount*2)+theme_bw()
        plt2<-ggplot(dat=as.data.frame(atac_sub@meta.data),aes(x=working_pseudotime,fill=as.factor(seurat_clusters)))+geom_histogram(bins=bincount*2)+theme_bw()

        ggsave(plt1/plt2,file=paste0(prefix,".pseudotime_resolution_membership.pdf"),width=20)
        system(paste0("slack -F ",prefix,".pseudotime_resolution_membership.pdf"," ryan_todo"))

        #Transcription Factor activity (ChromVAR TF motifs) over pseudotime bins
        motif_mean<-mclapply(unique(atac_sub$pseudotime_bin),chromvar_motifs_per_bin,mc.cores=5)
        motif_mean<-as.data.frame(do.call("rbind",motif_mean))
        row.names(motif_mean)<-unique(atac_sub$pseudotime_bin)
        motif_mean<-motif_mean[match(levels(atac_sub$pseudotime_bin),row.names(motif_mean)),]
        #Assign human readable TF motif names
        tfList <- getMatrixByID(JASPAR2020, ID=colnames(motif_mean))
        tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
        colnames(motif_mean)<-tfList
        motif_mean<-as.data.frame(t(motif_mean))
        motif_mean<-motif_mean[complete.cases(motif_mean),]
        dim(motif_mean)
        saveRDS(motif_mean,file=paste0(prefix,".pseudotime_chromVARmotifs.rds"))

        #Gene Activity over pseudotime bins
        ga_mean<-mclapply(unique(atac_sub$pseudotime_bin),geneactivity_per_bin,mc.cores=5)
        ga_mean<-as.data.frame(do.call("rbind",ga_mean))
        row.names(ga_mean)<-unique(atac_sub$pseudotime_bin)
        ga_mean<-ga_mean[match(levels(atac_sub$pseudotime_bin),row.names(ga_mean)),]
        ga_mean<-as.data.frame(t(ga_mean))
        ga_mean<-ga_mean[complete.cases(ga_mean),]
        dim(ga_mean)
        saveRDS(ga_mean,file=paste0(prefix,".pseudotime_geneactivity.rds"))

        #TF Modules over pseudotime bins
        tfmod_mean<-mclapply(unique(atac_sub$pseudotime_bin),TF_modules_per_bin,mc.cores=5)
        tfmod_mean<-as.data.frame(do.call("rbind",tfmod_mean))
        row.names(tfmod_mean)<-unique(atac_sub$pseudotime_bin)
        tfmod_mean<-tfmod_mean[match(levels(atac_sub$pseudotime_bin),row.names(tfmod_mean)),]
        tfmod_mean<-as.data.frame(t(tfmod_mean))
        tfmod_mean<-tfmod_mean[complete.cases(tfmod_mean),]
        dim(tfmod_mean)
        saveRDS(tfmod_mean,file=paste0(prefix,".pseudotime_tfmodules.rds"))

        #filter data by significant changes over pseudotime
        #chromvar_cds_pr_test_res<-readRDS(paste0(monocle_prefix,"_pseudotime.chromvar.rds"))
        #OR
        ###Filter by top 10% variances
        # calculate variance per column
        da_tf<-read.table("orgo_cirm43.onevone.da_tf.txt")
        da_tf<-da_tf[da_tf$enriched_group  %in% c("6","0","3","1","2"),]
        motif_mean<-readRDS(paste0(prefix,".pseudotime_chromVARmotifs.rds"))
        motif_mean<-motif_mean[row.names(motif_mean) %in% unique(da_tf$da_tf),]
        #motif_mean<-motif_mean[which(apply(motif_mean,1,var) < quantile(apply(motif_mean,1,var),0.25)),]
        #motif_mean<-motif_mean[row.names(motif_mean) %in% (chromvar_cds_pr_test_res %>% filter(q_value < 0.1))$tf_name,] # use that to subset the original data
        motif_mean<-na.omit(t(scale(t(motif_mean))))

        ####build up tangleogram by matching row names
        motif_mean<-motif_mean[!duplicated(unlist(lapply(strsplit(row.names(motif_mean),'[(:]'),"[",1))),] #remove similar motifs
        row.names(motif_mean)<-unlist(lapply(strsplit(row.names(motif_mean),'[(:]'),"[",1)) #simplify row names
        #subset ga and tf modules to those in motif mean
        #row.names(tfmod_mean)<-unlist(lapply(strsplit(row.names(tfmod_mean),'[-]'),"[",2))
        #tfmod_mean<-tfmod_mean[row.names(tfmod_mean) %in% row.names(motif_mean),] #ignoring tfmod mean for now
        
        #shared rows only for ga and motif mean
        #ga_mean<-ga_mean[row.names(ga_mean) %in% row.names(motif_mean),]
        #motif_mean<-motif_mean[row.names(motif_mean) %in% row.names(ga_mean),]

        #hc_motif <- motif_mean %>% dist() %>% hclust("centroid") 
        #k_waves_motif<-hc_motif %>% find_k(2:10)
        #hc_motif <- color_branches(as.dendrogram(hc_motif),k=k_waves_motif$k)
        #hc_motif<- ladderize(hc_motif,right=T)
        #dend_motif <- gc_motif %>% as.dendrogram %>% set("labels_to_char")
        #dend_ga <- ga_sub %>% match_order_by_labels(dend_ga)
        #dends_motif_ga <- dendlist(dend_motif,dend_ga)

        colfun=colorRamp2(quantile(unlist(motif_mean),c(0.5,0.75,0.99)),cividis(3))
        motif_sub<-motif_mean[order(unlist(lapply(1:nrow(motif_mean),function(x) which(motif_mean[x,]==max(motif_mean[x,]))))),]
        plt<-Heatmap(motif_sub,
        column_order=1:ncol(motif_sub),
        row_order=1:nrow(motif_sub),
        bottom_annotation=ha_barplots,
        #cluster_rows=hc_motif,
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        show_column_names=T,
        col=colfun,
        show_heatmap_legend=T)

        pdf(paste0(prefix,".pseudotime.chromvar.heatmap.pdf"),width=30)
        plt<-draw(plt,annotation_legend_list=lgd_list)
        dev.off()
        system(paste0("slack -F ",prefix,".pseudotime.chromvar.heatmap.pdf"," ryan_todo"))

        ga_mean<-readRDS(paste0(prefix,".pseudotime_geneactivity.rds"))
        ga_mean<-na.omit(t(scale(t(ga_mean))))
        da_ga<-read.table("orgo_cirm43.onevone.ga_tf.txt")
        da_ga<-da_ga[da_ga$enriched_group  %in% c("6","0","3","1","2"),]
        #ga_cds_pr_test_res<-readRDS(paste0(monocle_prefix,"_pseudotime.geneactivity.rds"))
        #ga_mean<-ga_mean[row.names(ga_mean) %in% row.names(ga_cds_pr_test_res %>% filter(q_value < 0.001)),] # use that to subset the original data
        #ga_mean<-ga_mean[which(apply(ga_mean,1,var) < quantile(apply(ga_mean,1,var),0.01)),]

        #ga_mean<-na.omit(scale(ga_mean))

        #hc_ga <- ga_mean %>% dist() %>% hclust("centroid") 
        #k_waves_ga<-hc_ga %>% find_k(2:10)
        #hc_ga<- color_branches(as.dendrogram(hc_ga),k=k_waves_ga$k)
        #hc_ga<- ladderize(hc_ga,right=T)

        ga_mean<-ga_mean[order(unlist(lapply(1:nrow(ga_mean),function(x) which(ga_mean[x,]==max(ga_mean[x,]))))),]
        colfun=colorRamp2(quantile(unlist(ga_mean),c(0,0.5,0.9)),magma(3))

        plt<-Heatmap(ga_mean,
        column_order=1:ncol(ga_mean),
        bottom_annotation=ha_barplots,
        #row_order=1:nrow(ga_mean),
        #cluster_rows=hc_ga,
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        show_column_names=T,
        col=colfun,
        show_heatmap_legend=T)

        pdf(paste0(prefix,".pseudotime.geneactivity.heatmap.pdf"),width=30)
        plt<-draw(plt,annotation_legend_list=lgd_list)
        dev.off()
        system(paste0("slack -F ",prefix,".pseudotime.geneactivity.heatmap.pdf"," ryan_todo"))

        dends_tf_ga <- dendlist(hc_motif, hc_ga)
        pdf("entanglement.test.pdf")
        tanglegram(dends_tf_ga) %>% untangle(method="step2side") %>% plot()
        dev.off()
        system("slack -F entanglement.test.pdf ryan_todo")

        tfmod_sub<-readRDS(file=paste0(prefix,".pseudotime_tfmodules.rds"))
        tfmod_sub<-tfmod_sub[which(apply(tfmod_sub,1,var) < quantile(apply(tfmod_sub,1,var),0.1)),]

        tfmod_sub<-na.omit(t(scale(t(tfmod_sub))))

        tfmod_sub<-tfmod_sub[order(unlist(lapply(1:nrow(tfmod_sub),function(x) which(tfmod_sub[x,]==max(tfmod_sub[x,]))))),]

        colfun=colorRamp2(quantile(unlist(tfmod_sub),c(0,0.5,0.9)),viridis(3))
        plt<-Heatmap(tfmod_sub,
        column_order=1:ncol(tfmod_sub),
        bottom_annotation=ha_barplots,
        row_order=1:nrow(tfmod_sub),
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        show_column_names=T,
        col=colfun,
        show_heatmap_legend=T)

        pdf(paste0(prefix,".pseudotime.tfmodule.heatmap.pdf"),width=30)
        plt<-draw(plt,annotation_legend_list=lgd_list)
        dev.off()
        system(paste0("slack -F ",prefix,".pseudotime.tfmodule.heatmap.pdf"," ryan_todo"))

  }
  
    #plot_over_pseudotime(object=orgo_cirm43,pseudotime="pseudotime_cluster0",bincount=50,
    #  prefix="./subcluster/exN",monocle_prefix="./subcluster/exN")
    #plot_over_pseudotime(object=orgo_cirm43,pseudotime="pseudotime_cluster1",bincount=25,
    #  prefix="./subcluster/rg",monocle_prefix="./subcluster/rg")
    orgo_cirm43<-subset(orgo_cirm43,cells=row.names(orgo_cirm43@meta.data)[which(!is.na(orgo_cirm43$pseudotime_allcells))])
    plot_over_pseudotime(object=orgo_cirm43,pseudotime="pseudotime_allcells",bincount=25,
      prefix="orgo_cirm43",monocle_prefix="orgo_cirm43")

    atac_sub<-subset(orgo_cirm43,cells=row.names(orgo_cirm43@meta.data)[which(!is.na(orgo_cirm43$pseudotime_rg))])
    plot_over_pseudotime(object=atac_sub,pseudotime="pseudotime_rg",bincount=25,
      prefix="rg",monocle_prefix="rg")


```

### Plotting Pseudotime Bins And Defined Regulatory Waves


```R

    setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
    library(Seurat)
    library(Signac)
    library(GenomeInfoDb)
    library(ggplot2)
    set.seed(1234)
    library(EnsDb.Hsapiens.v86)
    library(Matrix)
    library(cicero)
    library(SeuratWrappers)
    library(ComplexHeatmap)
    library(JASPAR2020)
    library(TFBSTools)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(patchwork)
    library(parallel) 
    library(dplyr)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)

    atac_sub<-readRDS("orgo_cirm43.filt.SeuratObject.Rds")

    #sig_tfs<-readRDS(file="pseudotime_cirm43_ccan_tfmotif_enrichment.rds") #excluded for now
    motif_mean<-readRDS("atac_sub.pseudotime_chromVARmotifs.rds")
    ccan_mean<-readRDS("atac_sub.pseudotime_ccan.rds")
    ccan_mean<-as.data.frame(t(scale(ccan_mean)))

    timebin<-readRDS("atac_sub.pseudotime.bins.rds")
    timebin_df<-as.data.frame(do.call("rbind",lapply(1:length(timebin),function(x) cbind(as.character(x),unlist(timebin[x])))))
    colnames(timebin_df)<-c("pseudotime_bin","cellID_idx")
    timebin_df$cellID<-atac_sub@meta.data$cellID[as.numeric(timebin_df$cellID_idx)]
    timebin_df<-merge(atac_sub@meta.data,timebin_df,by="cellID")


    #Set up different annotation stacked barplots
    annotation_bin_summary<-function(x){
    tmp_timebin<-as.data.frame(timebin_df %>% group_by_("pseudotime_bin",x) %>% summarize(count=n()))
    colnames(tmp_timebin)<-c("bin","var","count")
    tmp_timebin<-dcast(tmp_timebin,bin~var,value.var="count",fill=0)
    row.names(tmp_timebin)<-tmp_timebin$bin
    tmp_timebin<-as.data.frame(t(tmp_timebin[2:ncol(tmp_timebin)]))
    tmp_timebin = as.data.frame(t(scale(tmp_timebin, center = FALSE, 
                 scale = colSums(tmp_timebin))))
    tmp_timebin<-tmp_timebin[order(as.numeric(row.names(tmp_timebin)),decreasing=F),]
    return(tmp_timebin)}
    
    div_timebin<-annotation_bin_summary(x="DIV")
    #celltype_timebin<-annotation_bin_summary(x="celltype")
    cluster_timebin<-annotation_bin_summary(x="seurat_clusters")
    organoid_timebin<-annotation_bin_summary(x="orgID")

    ha_barplots = HeatmapAnnotation(
      DIV = anno_barplot(div_timebin, name="DIV",bar_width=1,gp = gpar(color = 1:nrow(div_timebin),fill = 1:nrow(div_timebin)),border=F),

      #celltype = anno_barplot(celltype_timebin, name="cell type",bar_width=1,gp = gpar(color = 1:nrow(celltype_timebin),fill= 1:nrow(celltype_timebin)),border=F),
      cluster = anno_barplot(cluster_timebin, name="cluster",bar_width=1,gp = gpar(color = 1:nrow(cluster_timebin),fill = 1:nrow(cluster_timebin)),border=F),
      organoid = anno_barplot(organoid_timebin, name="organoid",bar_width=1,gp = gpar(color = 1:nrow(organoid_timebin),fill = 1:nrow(organoid_timebin)),border=F)
      )

      lgd_list = list(
      Legend(labels = colnames(div_timebin), title = "DIV", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(div_timebin))),
      #Legend(labels = colnames(celltype_timebin), title = "celltype", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(celltype_timebin))),
      Legend(labels = colnames(cluster_timebin), title = "cluster", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(cluster_timebin))),
      Legend(labels = colnames(organoid_timebin), title = "organoid", type = "points", pch = 16, legend_gp = gpar(col = 1:ncol(organoid_timebin)))
  )

    ###Plotting ChromVAR Motif Accessibility Scores through Pseudotime

    tf_markers<-c("CTCF","EMX1","EMX2","LHX2","PAX6","RFX4","SOX2",
               "TBR1","EOMES","NEUROD1","NEUROD2","NEUROG1","TGIF1","TGIF2",
               "DLX1","DLX2","DLX6","GSX2","LHX6",
               "POU3F3","POU3F2","TFAP4")
    #subset to markers
    tf_markers<-tf_markers[tf_markers %in% row.names(motif_mean)]
    motif_markers<-which(row.names(motif_mean) %in% c(tf_markers))
    motif_markers<-as.data.frame(cbind(as.numeric(motif_markers),row.names(motif_mean)[motif_markers]))
    colnames(motif_markers)<-c("motif_row","motif_name")
    ha_rowmarkers = rowAnnotation(Markers = anno_mark(
      at = as.numeric(motif_markers$motif_row),
      , labels = motif_markers$motif_name))

    motif_mean_scale<-scale(motif_mean)
    col_fun = colorRamp2(quantile(motif_mean, na.rm = TRUE, prob = c(0.2,0.4,0.6,0.8,0.95)), rev(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")))

    plt3<-Heatmap(motif_mean,
      column_order=1:ncol(motif_mean),
      row_names_gp = gpar(fontsize = 3),
      clustering_distance_rows="maximum",
      #col=col_fun,
      show_column_names=F,
      bottom_annotation=ha_barplots,
      right_annotation=ha_rowmarkers,
      row_km=6,
      show_heatmap_legend=T)

    plt<-draw(plt3,annotation_legend_list=lgd_list)

    pdf(paste("orgo_cirm43.pseudotime.motifTF.heatmap.pdf",sep="."),width=30,height=10)
    plt
    dev.off()
    system("slack -F orgo_cirm43.pseudotime.motifTF.heatmap.pdf ryan_todo")
    #set up list of factors in each wave and save
    wave_list<-as.data.frame(do.call("rbind",lapply(1:length(row_order(plt)),function(x) cbind(names(row_order(plt))[x],unlist(row_order(plt)[x])))))
    colnames(wave_list)<-c("wave_group","row_index")
    wave_list$tf<-row.names(motif_mean)[as.numeric(wave_list$row_index)]
    saveRDS(wave_list,file="cirm43_pseudotime_TF_wavelist.rds")

    col_fun = colorRamp2(c(-3, -1.5,0, 1.5,3), rev(c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0")))

    ###Plotting CCANs through Pseudotime
    #Adding ccan markers
    ccan_markers$ccan<-row.names(ccan_markers)
    colnames(ccan_markers)<-c("gene","ccan")
    ccan_markers<-ccan_markers[!(ccan_markers$gene=="NA"),]
    ha_rowmarkers = rowAnnotation(Markers = anno_mark(
      at = match(as.numeric(ccan_markers$ccan),row.names(ccan_mean)),
      , labels = ccan_markers$gene))


    plt1<-Heatmap(ccan_mean,
      column_order=1:ncol(ccan_mean),
      show_row_names=F,
      clustering_distance_rows="pearson",
      row_km=6,
      bottom_annotation=ha_barplots,
      right_annotation=ha_rowmarkers,
      col=col_fun
      )
    

    #plt<-plt1+plt2
    plt<-draw(plt1,annotation_legend_list=lgd_list)

    pdf("cirm43_pseudotime.ccan.heatmap.pdf",width=30,height=10)
    plt
    dev.off()

    system("slack -F cirm43_pseudotime.ccan.heatmap.pdf ryan_todo")
    
    #set up list of ccans in each wave and save
    wave_list<-as.data.frame(do.call("rbind",lapply(1:length(row_order(plt)),function(x) cbind(names(row_order(plt))[x],unlist(row_order(plt)[x])))))
    colnames(wave_list)<-c("wave_group","row_index")
    wave_list$tf<-row.names(ccan_mean)[as.numeric(wave_list$row_index)]
    saveRDS(wave_list,file="cirm43_pseudotime_ccan_wavelist.rds")  

```
```R






    motif_mean_scale<-scale(motif_mean)
    library(circlize)
    col_fun = colorRamp2(quantile(motif_mean_scale, na.rm = TRUE, prob = c(0.2,0.4,0.6,0.8,0.95)), rev(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")))

    plt3<-Heatmap(motif_mean_scale,
      column_order=1:ncol(motif_mean),
      row_names_gp = gpar(fontsize = 3),
      column_names_gp = gpar(fontsize = 3),
      clustering_distance_rows="maximum",
      #col=col_fun,
      show_column_names=T,
      #bottom_annotation=ha_barplots,
      #right_annotation=ha_rowmarkers,
      row_km=6,
      show_heatmap_legend=T)

    pdf("test_motif.pseudotime.heatmap.pdf")
    plt<-draw(plt3)#,annotation_legend_list=lgd_list)
    dev.off()
    system("slack -F test_motif.pseudotime.heatmap.pdf ryan_todo")



    tfmod_mean_scale<-scale(tfmod_mean)
    library(circlize)
    col_fun = colorRamp2(quantile(tfmod_mean_scale, na.rm = TRUE, prob = c(0.2,0.4,0.6,0.8,0.95)), rev(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")))

    plt3<-Heatmap(tfmod_mean_scale,
      column_order=1:ncol(tfmod_mean),
      row_names_gp = gpar(fontsize = 3),
      column_names_gp = gpar(fontsize = 3),
      clustering_distance_rows="euclidean",
      #col=col_fun,
      show_column_names=T,
      #bottom_annotation=ha_barplots,
      #right_annotation=ha_rowmarkers,
      row_km=6,
      show_heatmap_legend=T)

    pdf("test_motif.pseudotime.heatmap.pdf")
    plt<-draw(plt3)#,annotation_legend_list=lgd_list)
    dev.off()
    system("slack -F test_motif.pseudotime.heatmap.pdf ryan_todo")






  #Set up gneomic annotation
  # Generate gene annotations for CCAN assignment of marker genes
  hg38_annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  pos <-as.data.frame(hg38_annotations,row.names=NULL)
  pos$chromosome<-paste0("chr",pos$seqnames)
  pos$gene<-pos$gene_id
  pos <- subset(pos, strand == "+")
  pos <- pos[order(pos$start),] 
  pos <- pos[!duplicated(pos$tx_id),] # remove all but the first exons per transcript
  pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS
  neg <-as.data.frame(hg38_annotations,row.names=NULL)
  neg$chromosome<-paste0("chr",neg$seqnames)
  neg$gene<-neg$gene_id
  neg <- subset(neg, strand == "-")
  neg <- neg[order(neg$start,decreasing=TRUE),] 
  neg <- neg[!duplicated(neg$tx_id),] # remove all but the first exons per transcript
  neg$end <- neg$end + 1 # make a 1 base pair marker of the TSS

  gene_annotation<- rbind(pos, neg)
  gene_annotation <- gene_annotation[,c("chromosome","start","end","gene_name")] # Make a subset of the TSS annotation columns containing just the coordinates and the gene name
  names(gene_annotation)[4] <- "gene" # Rename the gene symbol column to "gene"
  gene_annotation<-makeGRangesFromDataFrame(gene_annotation,keep.extra.columns=T)


  cell_accessibility_per_ccan<-function(x){
  ccan_tmp<-ccan[x]
  peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
  dat_tmp<-as.data.frame(atac_sub@assays$peaks@data[peak_subset,])
  norm_dat_tmp<-colSums(dat_tmp)/atac_sub@meta.data[match(atac_sub@meta.data$cellID,colnames(dat_tmp)),]$nCount_peaks
  norm_dat_tmp
  }

  marker_annote_per_ccan<-function(x){
  ccan_tmp<-ccan[x]
  markers<-ifelse(countOverlaps(query=ccan_tmp,subject=gene_annotation)>0,
    unique(gene_annotation[unique(findOverlaps(query=ccan_tmp,subject=gene_annotation)@to),]$gene),
    "NA")
  markers
  }

  ccan_summarize_per_bin<-function(x){
  bin_tmp<-timebin[[x]]
  return(apply(ccan_sum[,bin_tmp],1,mean))
  }


    #Assign CCAN to RNA cell type markers by promoter overlap
    markers<-c("CTCF","EMX1","EMX2","LHX2","PAX6","RFX4","SOX2",
                 "TBR1","EOMES","NEUROD1","NEUROD2","NEUROG1","TGIF1","TGIF2",
                 "DLX1","DLX2","DLX6","GSX2","LHX6",
                 "POU3F3","POU3F2","TFAP4")
    gene_annotation<-gene_annotation[gene_annotation$gene %in% markers,]
    gene_annotation<-makeGRangesFromDataFrame(do.call("rbind",
      lapply(unique(gene_annotation$gene), function(x) {
      tmp<-gene_annotation[gene_annotation$gene==x,]
      return(as.data.frame(tmp[which.max(as.data.frame(tmp@ranges)$width),]))})),keep.extra.columns=T) 
      #select longest transcript of each gene for ccan overlap

    #Set up TF markers
    tf_markers<-c("SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3",
               "NEUROG1","SSTR2","EOMES","PPP1R17","NEUROD4",
               "SLC17A7","NEUROD6","SATB2","TBR1","SLA",
               "DLX2","DLX1","LHX6","GAD1")

    #CCANs over pseudotime
      #First generating list of peaks in CCANs
        ccan<-Links(atac_sub)
        ccan<-ccan[!is.na(ccan$group),]
        ccan<-split(ccan,f=ccan$group)

      #Plotting amount of peak membership per CCAN
        ccan_membership<-unlist(lapply(ccan,length))
        plt<-ggplot()+geom_density(aes(x=ccan_membership))+theme_bw()+xlim(c(0,500))
        ggsave(plt,file="ccan_membership.pdf")
        system("slack -F ccan_membership.pdf ryan_todo")

      #Getting list of peaks
        peaks<-row.names(atac_sub@assays$peaks@data)
        seqname=unlist(lapply(strsplit(peaks,"-"),"[",1))
        start=unlist(lapply(strsplit(peaks,"-"),"[",2))
        end=unlist(lapply(strsplit(peaks,"-"),"[",3))
      peak_info<-makeGRangesFromDataFrame(data.frame(seqname=seqname,start=start,end=end,peak_id=peaks),keep.extra.columns=T)

      #find overlaps between peaks and CCAN regions (cicero merges similar peaks, so take all within overlap)
      #then from those overlapped peak per CCAN, grab peak data from the seurat object
      #then perform colSums (per cell peak accessibility summarized by the CCAN)
      #then normalize the data based on the cell specific size factor
      ccan_sum<-mclapply(1:length(ccan),cell_accessibility_per_ccan,mc.cores=10)
      ccan_sum<-as.data.frame(do.call("rbind",ccan_sum))
      row.names(ccan_sum)<-names(ccan)
      dim(ccan_sum)

      #Find overlap between CCAN peaks and gene bodies
      ccan_markers<-mclapply(1:length(ccan),marker_annote_per_ccan,mc.cores=10)
      ccan_markers<-as.data.frame(do.call("rbind",ccan_markers))
      row.names(ccan_markers)<-names(ccan)

      #summarize ccans through pseudotime cell bins
      #using mean here as the summary statistic
      ccan_mean<-mclapply(1:length(timebin),ccan_summarize_per_bin,mc.cores=10)
      ccan_mean<-as.data.frame(do.call("rbind",ccan_mean))
      row.names(ccan_mean)<-names(timebin)
      colnames(ccan_mean)<-names(ccan)
      saveRDS(ccan_mean,file="atac_sub.pseudotime_ccan.rds")


    #Question: do we see a bias in CCAN opening through pseudotime with TF motif presence?
    #Take out motif presence in overlapping peaks with CCANs
    #This is taking from the chromvar motif peak association
    #generate background set per TF motif from all peaks
    #going to use a hypergeometric score since number of peaks per ccan is variable

    #Building contingency table
    ##             TFmotif  non-TFmotif
    ## CCAN            a  b
    ## Non-CCAN        c  d
    ## Total    f (a+c) g(b+d)

   # dat_tmp_bg_f.<-colSums(orgo_cirm43@assays$peaks@motifs@data>0)
   # dat_tmp_bg_g.<-colSums(orgo_cirm43@assays$peaks@motifs@data==0)

   # fisher_tfmotifs_per_ccan<-function(ccan_idx=x,dat_tmp_bg_f=dat_tmp_bg_f.,dat_tmp_bg_g=dat_tmp_bg_g.){
   #   ccan_tmp<-ccan[ccan_idx]
   #   peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
   #   dat_tmp<-as.data.frame(orgo_cirm43@assays$peaks@motifs@data[peak_subset,])
   #   a_list<-colSums(dat_tmp>0)
   #   b_list<-colSums(dat_tmp==0)
   #   c_list<-dat_tmp_bg_f - a_list
   #   d_list<-dat_tmp_bg_g - b_list
   #   tf_dat<-unlist(
   #     lapply(colnames(dat_tmp),FUN=function(y) {
   #       a<-a_list[y]
   #       b<-b_list[y]
   #       c<-c_list[y]
   #       d<-d_list[y]
   #       contingency_table<-matrix(c(a,b,c,d),nrow=2)
   #       fisher.test(contingency_table,alternative="greater")$p.value}
   #     ))
   #   names(tf_dat)<-colnames(dat_tmp)
   #   return(tf_dat)
   # }
      
    # ccan_tf_enrich<-mclapply(1:length(ccan),function(x) fisher_tfmotifs_per_ccan(ccan_idx=x),mc.cores=10)
    # ccan_tf_enrich<-as.data.frame(do.call("rbind",ccan_tf_enrich))
    # #Assign human readable TF motif names
    # tfList <- getMatrixByID(JASPAR2020, ID=colnames(ccan_tf_enrich))
    # tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
    # colnames(ccan_tf_enrich)<-tfList
    # sig_tfs<-as.data.frame(t(apply(ccan_tf_enrich,1,p.adjust,method="bonferroni")))
    # #sig_tfs <- sig_tfs[colSums(sig_tfs<0.05)>length(ccan)/20] #limit to TF motifs that reach significance in at least 5% ccans
    # sig_tfs <- -log10(sig_tfs)
    # saveRDS(sig_tfs,file="pseudotime_cirm43_ccan_tfmotif_enrichment.rds")

    #Question: do we see bias in GWAS associated sites within CCANs through pseudotime?
    #Take out GWAS sites in overlapping peaks with CCANs
    #This is taking from https://www.ebi.ac.uk/gwas/docs/file-downloads
    #generate background set per condition from all peaks
    #going to use a hypergeometric score since number of peaks per ccan is variable

    #Building contingency table
    ##             GWAS   non-GWAS
    ## CCAN            a  b
    ## Non-CCAN        c  d
    ## Total    f (a+c) g(b+d)
    # gwas<-read.csv("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/gwas_catalog_v1.0-associations_e100_r2020-10-20.tsv",header=T,sep="\t")
    # gwas$seqname<-paste0("chr",gwas$CHR_ID)
    # gwas$start<-as.numeric(gwas$CHR_POS)
    # gwas<-gwas[colnames(gwas) %in% c("seqname","start","MAPPED_GENE","DISEASE.TRAIT")]
    # gwas$end<-gwas$start+1
    # gwas<-gwas[complete.cases(gwas),]
    # gwas<-makeGRangesFromDataFrame(gwas,keep.extra.columns=T)
    # gwas<-split(gwas,f=gwas$DISEASE.TRAIT)

  #GWAS analysis to be changed
    #dat_tmp_bg_f<-colSums(orgo_cirm43@assays$peaks@motifs@data>0)
    #dat_tmp_bg_g<-colSums(orgo_cirm43@assays$peaks@motifs@data==0)

    #tfmotifs_per_ccan<-function(x){
  #   ccan_tmp<-ccan[x]
  #   peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
  #   dat_tmp<-as.data.frame(orgo_cirm43@assays$peaks@motifs@data[peak_subset,])
  #   a_list<-colSums(dat_tmp>0)
  #   b_list<-colSums(dat_tmp==0)
  #   c_list<-dat_tmp_bg_f - a_list
  #   d_list<-dat_tmp_bg_g - b_list
  #   tf_dat<-unlist(lapply(1:ncol(dat_tmp),FUN=fisher_test_columnwise,a_list=a_list,b_list=b_list,c_list=c_list,d_list=d_list))
  #   names(tf_dat)<-colnames(dat_tmp)
  #   return(tf_dat)
  # }

  # ccan_tf_enrich<-mclapply(1:length(ccan),tfmotifs_per_ccan,mc.cores=10)
  # ccan_tf_enrich<-as.data.frame(do.call("rbind",ccan_tf_enrich))
    #Assign human readable TF motif names
  # tfList <- getMatrixByID(JASPAR2020, ID=colnames(ccan_tf_enrich))
  # tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  # colnames(ccan_tf_enrich)<-tfList
  # sig_tfs<-as.data.frame(t(apply(ccan_tf_enrich,1,p.adjust,method="bonferroni")))
  # sig_tfs <- sig_tfs[colSums(sig_tfs<0.05)>length(ccan)/20] #limit to TF motifs that reach significance in at least 5% ccans
  # sig_tfs <- -log10(sig_tfs)

  #Set up functions
  sliding_window<-function(i){
  start_location=(stepsize*i)+1
  end_location=(stepsize*i)+binwidth
  return(ord_cells[start_location:end_location])
  }


    binwidth=as.integer(length(ord_cells)/100)
    stepsize=as.integer(length(ord_cells)/100)
    timebin<-lapply(0:round(length(ord_cells)/stepsize-1),FUN=sliding_window)
    summary(unlist(lapply(timebin,FUN=length))) #membership of bins
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #324     324     324     324     324     324


```


### 3D Plotting for better trajectory visualization

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(plotly)

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  cirm43_subset<-subset(orgo_cirm43,cells=which(orgo_cirm43$celltype %in% c("radial_glia","intermediate_progenitor","excitatory_neuron")))


  #reading model selected RDS
  library(cisTopic)

  cirm43_cistopic_models<-readRDS(file="orgo_cirm43.CisTopicObject.Rds")
  #set topics based on derivative
  cirm43_selected_topic=27
  cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,select=cirm43_selected_topic,keepModels=T)

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=3))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y","z")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


      #Add cell embeddings into seurat
      cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
      colnames(cell_embeddings)<-cisTopicObject@cell.names
      cell_embeddings<-cell_embeddings[,colnames(cell_embeddings) %in% dims$cellID,]
      n_topics<-nrow(cell_embeddings)
      row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
      cell_embeddings<-as.data.frame(t(cell_embeddings))

      #Add feature loadings into seurat
      feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
      row.names(feature_loadings)<-paste0("topic_",1:n_topics)
      feature_loadings<-as.data.frame(t(feature_loadings))

      #combined cistopic results (cistopic loadings and umap with seurat object)
      cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
      umap_dims<-as.data.frame(as.matrix(dims[2:4]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2","UMAP_3")
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
        resolution=resolution
      )

  return(object_input)}

  cirm43_subset<-cistopic_wrapper(object_input=cirm43_subset,
    cisTopicObject=cirm43_cisTopicObject,
    resolution=0.5)


  monocle_processing<-function(prefix, seurat_input){
      atac.cds <- as.cell_data_set(seurat_input)
      atac.cds <- cluster_cells(cds = atac.cds, reduction_method = "UMAP") 
      #Read in cds from cicero processing earlier and continue processing
      atac.cds<- learn_graph(atac.cds, 
                             use_partition = F, 
                             learn_graph_control=list(
                                 minimal_branch_len=10,
                                 orthogonal_proj_tip=F,
                                 prune_graph=T))
      return(atac.cds)
  }

  cirm43_subset.cicero<-monocle_processing(seurat_input=cirm43_subset,prefix="cirm43")

  #Then determine root nodes via plots and assign by order cells function.
  #cirm43.cds <- order_cells(cirm43_subset.cicero, 
  #reduction_method = "UMAP", 
  #root_pr_nodes = c("Y_253")) #Chose youngest cells as root


  #Generating a 3D Plot via Plotly of the umap projection.
  #Loading in additional libraries.
    library(plotly)
    library(htmlwidgets)
    library(RColorBrewer)
    library(igraph)

    dat<-merge(cirm43_subset@reductions$umap@cell.embeddings,cirm43_subset@meta.data,by="row.names") 
    dat$DIV<-as.character(dat$DIV) 
    dat$differentiation_exp<-as.character(dat$differentiation_exp) 
    write.table(dat,file="orgo_cirm43.3dmetadata.csv",col.names=T,row.names=F,sep=",")
    system("slack -F orgo_cirm43.3dmetadata.csv ryan_todo")
    
    #igraph plotting from monocle3 function
    plotly_pseudotime<-function(p=.,x){
    #Set up node locations
    ica_space_df <- t(x@principal_graph_aux[["UMAP"]]$dp_mst) %>%
            as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = "UMAP_1",
            prin_graph_dim_2 = "UMAP_2", prin_graph_dim_3 = "UMAP_3") %>% dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))
    #Set up edge spans
        dp_mst <- x@principal_graph[["UMAP"]]
        edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from",
            target = "to") %>% dplyr::left_join(ica_space_df %>%
            dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
                source_prin_graph_dim_2 = "prin_graph_dim_2",
                source_prin_graph_dim_3 = "prin_graph_dim_3"),
            by = "source") %>% dplyr::left_join(ica_space_df %>%
            dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1",
                target_prin_graph_dim_2 = "prin_graph_dim_2",
                target_prin_graph_dim_3 = "prin_graph_dim_3"),
            by = "target")
        for (i in 1:nrow(edge_df)) {
            p <- p %>% plotly::add_trace(x = as.vector(t(edge_df[i,
                c("source_prin_graph_dim_1", "target_prin_graph_dim_1")])),
                y = as.vector(t(edge_df[i, c("source_prin_graph_dim_2",
                  "target_prin_graph_dim_2")])), z = as.vector(t(edge_df[i,
                  c("source_prin_graph_dim_3", "target_prin_graph_dim_3")])),
                color = "rgba(0,0,0,1)", line = list(color = "rgba(0,0,0,1)",
                  width = 3), mode = "lines",
                type = "scatter3d", showlegend = FALSE)
        }
        return(p)
      }


    l <- list(
      font = list(
        family = "sans-serif",
        size = 8,
        color = "#000"),
      bgcolor = "#E2E2E2",
      bordercolor = "#FFFFFF",
      borderwidth = 2,
      travedorder="group")

    #Generate 3D Plot and standalone HTML widget
    p<-plot_ly(dat,
      x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
      type="scatter3d", mode="markers",
      size=I(4), hoverinfo="none", legendgroup="Experiment",
      color=~differentiation_exp) %>%
    plotly_pseudotime(x=cirm43_subset.cicero) %>%
    add_markers(color=~DIV,legendgroup="DIV") %>%
    add_markers(color=~pseudotime,legendgroup="Pseudotime",showlegend=F) %>%
    add_markers(color=~seurat_clusters,legendgroup="Clusters") %>% layout(legend = l)

    htmlwidgets::saveWidget(as_widget(partial_bundle(p)), "cirm43_umap.3d.html",selfcontained=TRUE)

  system("slack -F cirm43_umap.3d.html ryan_todo")
```



### Plotting ChromVAR motifs through pseudotime

 The following code is exploratory but in the end wasn't included in analysis for the manuscript.
I mainly just wanted to play around with network analysis a bit.

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  #Setting up chromvar matrix
  tfList <- getMatrixByID(JASPAR2020, ID=row.names(orgo_cirm43@assays$chromvar@data)) 
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  dat_tf<-orgo_cirm43@assays$chromvar@data
  row.names(dat_tf)<-tfList
  dat_tf<-data.frame(t(dat_tf))
  dat_tf$cellID<-row.names(dat_tf)
  cirm43_tf<-dat_tf

  cirm43_pseudotime<-orgo_cirm43@meta.data[orgo_cirm43@meta.data$cell_line=="CIRM43",c("cellID","pseudotime")]
  cirm43_tf<-merge(cirm43_tf,cirm43_pseudotime,by="cellID")
  cirm43_tf<-cirm43_tf[complete.cases(cirm43_tf),]
  cirm43_cell_row_order=order(cirm43_tf$pseudotime,decreasing=T)


  cirm43_row_ha = rowAnnotation(
      seurat_clusters=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$seurat_clusters),
      organoid=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$orgID),
      DIV=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$DIV),
      experiment=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$differentiation_exp)
  )

  #what is the distribution of pseudotime?
  plt<-ggplot()+geom_density(aes(x=orgo_cirm43$pseudotime))+theme_bw()
  ggsave(plt,file="pseudotime_distribution.pdf")
  system("slack -F pseudotime_distribution.pdf ryan_todo")
  #divide pseudotime into equally sized (by cell count) bins

  cirm43_tf<-cirm43_tf[!(colnames(cirm43_tf) %in% c("cellID","pseudotime"))]

  vars <- apply(cirm43_tf, 2, var)
  cirm43_tf<-cirm43_tf[vars > quantile(vars, 0.7)]

  pdf("cirm43_tfmotif_pseudotimeordering.pdf")
  Heatmap(as.matrix(cirm43_tf),column_km=7,left_annotation = cirm43_row_ha,
          row_order=cirm43_cell_row_order,
          show_row_names=F,
          clustering_distance_columns="maximum",
          column_names_gp = gpar(fontsize = 5)
  )
  dev.off()
  system("slack -F cirm43_tfmotif_pseudotimeordering.pdf ryan_todo")

  library(igraph)
  cor_mat<-cor(cirm43_tf,method="pearson")
  cor_mat[which(cor_mat<0.3,arr.ind=T)]<-0
  graph<-graph.adjacency(cor_mat,weighted=TRUE,"directed",diag=F)
  pdf("tf_motif.igraph.pdf")
  plot.igraph(graph,edge.size=E(graph)$weight,vertex.size=3,vertex.label.cex=0.8)
  dev.off()
  system("slack -F tf_motif.igraph.pdf ryan_todo")
```



### Define TF Enrichment In CCAN Waves

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R

```
{% endcapture %} {% include details.html %} 




## Save final seurat files for UCSC cellbrowser and set up.

To view data in an interactive way, I converted the final seurat object for UCSCs cell browser.

I then hosted locally to test.
Install cellbrowser via this (guide.)[https://cellbrowser.readthedocs.io/en/master/installation.html]

```R
#Downloaded seurat objects off clusters
#Prepare seurat data as cellbrowser file

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
library(Signac)
library(Seurat)
library(Matrix)
library(R.utils)
library(dplyr)


#modified function from https://github.com/satijalab/seurat-wrappers/blob/master/R/cellbrowser.R just for export
ExportToCellbrowser <- function(
object,
reduction="umap",
assay.name="GeneActivity",
dir,
dataset.name,
marker.file,
marker.n=5,
cluster.field="seurat_clusters",
gene_list,
skip.expr.mat=FALSE) {

  reducNames = reduction
  # see https://satijalab.org/seurat/essential_commands.html
    idents <- Idents(object = object)
    cellOrder <- row.names(object@meta.data)
    counts <-  object@assays[[assay.name]]@data
    counts<-counts[,cellOrder]
    genes <- rownames(x = counts)
    dr<-object@reductions[[reduction]]
    meta.fields<-object@meta.data
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)}

  # Export expression matrix
    # we have to write the matrix to an mtx file
    matrixPath <- file.path(dir, "matrix.mtx")
    genesPath <- file.path(dir, "features.tsv")
    barcodesPath <- file.path(dir, "barcodes.tsv")
if(!skip.expr.mat){
    message("Writing expression matrix to ", matrixPath)
    message("This may take a couple of minutes...")
    writeMM(counts, matrixPath)
    # easier to load if the genes file has at least two columns. Even though seurat objects
    write.table(as.data.frame(cbind(rownames(counts), rownames(counts))), file=genesPath, sep="\t", row.names=F, col.names=F, quote=F)
    write(colnames(counts), file = barcodesPath)

    message("Gzipping expression matrix")
    gzip(matrixPath,overwrite=T,remove=T)
    gzip(genesPath,overwrite=T,remove=T)
    gzip(barcodesPath,overwrite=T,remove=T)
}
    matrixOutPath <- "matrix.mtx.gz"

    #Export embeddings
    df <- dr@cell.embeddings
    if (ncol(x = df) > 2) {
      warning('Embedding ', embedding, ' has more than 2 coordinates, taking only the first 2')
      df <- df[, 1:2]
    }
    colnames(x = df) <- c("x", "y")
    df <- data.frame(cellId = rownames(x = df), df, check.names = FALSE)
    fname <- file.path(dir,"umap.coords.tsv")
    message("Writing embeddings to ", fname)
    write.table(df[cellOrder, ], sep="\t", file=fname, quote = FALSE, row.names = FALSE)
  embeddings.conf <- sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',reduction)
  coords.string <- paste0("[",paste(embeddings.conf, collapse = ",\n"),"]")

  # Export metadata
  enum.fields=colnames(object@meta.data)[colnames(object@meta.data) %in% c("cellid","pcr_idx","tn5_idx",cluster.field,"seurat_clusters","cluster_ID","predicted.id","celltype")]
  fname <- file.path(dir, "meta.tsv")
  message("Writing meta data to ", fname)
  out.meta<-as.matrix(object@meta.data[cellOrder, enum.fields])
  write.table(out.meta, sep = "\t", file = fname, quote = FALSE, row.names = FALSE)
  enum.string <- paste0("[",paste(paste0('"', enum.fields, '"'), collapse = ", "), "]")

  # Export markers
    file <- file.path("markers.tsv")
    fname <- file.path(dir, file)
      da_ga<-read.table(file=marker.file)
        da_ga$label<-""
        for (x in unique(da_ga$enriched_group)){
        selc_genes<-as.data.frame(da_ga %>% filter(enriched_group==x) %>% dplyr::arrange(p_val_adj) %>% dplyr::slice(1:marker.n))$da_region
        da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$label<- da_ga[da_ga$da_region %in% selc_genes & da_ga$enriched_group==x,]$da_region
        }
    da_ga<-da_ga[da_ga$label!="",]
      message("Writing top ", marker.n, " cluster markers to ", fname)
    da_ga$cluster<-da_ga$enriched_group
      da_ga<-da_ga[c("cluster","label","p_val_adj","enriched_group")]
      colnames(da_ga)<-c("cluster","symbol","Adjusted p Value","Cluster Enrichment")
      write.table(x = da_ga, file = fname, quote = FALSE, sep = "\t", col.names = T,row.names=F)
    markers.string <- sprintf('markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',file)

#Export Gene list
  file <- file.path("quickGenes.csv")
  fname <- file.path(dir, file)
  write.table(gene_list,sep=",",file = fname, quote = FALSE, row.names = FALSE)

  config <- '
# This is a bare-bones cellbrowser config file auto-generated from R.
# Look at https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf
# for a full file that shows all possible options
name="%s"
shortLabel="%1$s"
exprMatrix="%s"
tags = ["sciatac"]
meta="meta.tsv"
# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"
geneIdType="auto"
# file with gene,description (one per line) with highlighted genes, called "Dataset Genes" in the user interface
quickGenesFile="quickGenes.csv"
clusterField="%s"
labelField="%s"
enumFields=%s
%s
coords=%s
geneLabel="Gene Activity"
unit="GA"'

  config <- sprintf(
    config,
    dataset.name,
    matrixOutPath,
    cluster.field,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  confPath = file.path(dir, "cellbrowser.conf")
  message("Writing cellbrowser config to ", confPath)
  cat(config, file = confPath)
  message("Prepared cellbrowser directory ", dir)
}


orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

hg38_marker_list<- c("FOXG1", "SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3", "MOXD1", "NEUROG1","EOMES","PPP1R17","PTPRZ1","NEUROD4", "SLC17A7","NEUROD6","RIT1","TBR1","SLA","NHLH1", "DLX2","DLX1","LHX6","GAD1", "NEUROD1", "SATB2","CUX1","RELN","MEF2C", "SLC1A3","GLI3","MKI67","STMN2","GADD45B","MAP2","CALB2","FBXO32", "PGK1","GORASP2" ) 

ExportToCellbrowser(
object=orgo_cirm43,
reduction="umap",
assay.name="GeneActivity",
dir="orgo_cb",
dataset.name="orgo",
marker.file="orgo_cirm43.onevrest.da_ga.txt",
marker.n=5,
cluster.field="seurat_clusters",
gene_list=hg38_marker_list,
skip.expr.mat=FALSE)
```

Locally building and hosting the cell browser.
Using windows wsl2 ubuntu environment

```bash
#install cellbrowser via conda
#conda install -c bioconda ucsc-cell-browser

cd /mnt/c/Users/mulqueen/Documents/organoid/orgo_cb
cbBuild -o orgo_cellbrowser -p 8888

#Go to http://localhost:8888/ to use the interactive viewer

```



<!--

### Recluster cells now that data has been filtered
{% capture summary %} Code {% endcapture %} {% capture details %}  
```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cisTopic)

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(15:30),nCores=15,addModels=FALSE)
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=orgo_cirm43,prefix="orgo_cirm43_filt")
  cirm43_cistopic_models<-readRDS("orgo_cirm43_filt.CisTopicObject.Rds")


  #Setting up topic count selection
  pdf("cirm43_filt_model_selection.pdf")
  par(mfrow=c(1,3))
  cirm43_cistopic_models <- selectModel(cirm43_cistopic_models, type='derivative')
  dev.off()
  system("slack -F cirm43_filt_model_selection.pdf ryan_todo")


  ###############################################
  #Loop through cistopic models
  cistopic_loop<-function(topic_number,object_input,models_input){
      models_input<-selectModel(models_input,select=topic_number)
      #perform UMAP on topics
      topic_df<-as.data.frame(models_input@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      print("Performed UMAP.")
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")
     
      #combine with seurat object    
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$umap<-cistopic_umap
      
      #finally plot
      plt<-DimPlot(object_input,group.by=c('DIV','cell_line'),size=0.1)+ggtitle(as.character(topic_number))
      return(plt)
  }

  library(patchwork)

  plt_list<-lapply(cirm43_cistopic_models@calc.params$runWarpLDAModels$topic,
                     FUN=cistopic_loop,
                     object_input=orgo_cirm43,
                     models_input=cirm43_cistopic_models)
  plt_list<-wrap_plots(plt_list)
  ggsave(plt_list,file="cirm43.filt.umap_multipleTopicModels_clustering.png",height=20,width=60,limitsize=FALSE)
  system("slack -F cirm43.filt.umap_multipleTopicModels_clustering.png ryan_todo")

  ###############################################



  #set topics based on derivative
  cirm43_selected_topic=21
  cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,select=cirm43_selected_topic,keepModels=T)

  #saving model selected RDS
  saveRDS(cirm43_cisTopicObject,file="orgo_cirm43.filt.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   

      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")

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
        resolution=resolution,
        graph.name="peaks_snn"
      )
    plt1<-DimPlot(object_input,group.by=c('seurat_clusters'),size=0.1)
    ggsave(plt1,file="orgo_cirm43.filt.umap.clus.pdf")
    system("slack -F orgo_cirm43.filt.umap.clus.pdf ryan_todo")
    return(object_input)
}

  orgo_cirm43_out<-cistopic_wrapper(object_input=orgo_cirm43,cisTopicObject=cirm43_cisTopicObject,resolution=0.1)

  saveRDS(orgo_cirm43_out,file="orgo_cirm43.filt.SeuratObject.Rds")   ###save Seurat file

```
{% endcapture %} {% include details.html %} 


### Subclustering of celltypes and Pseudotime Trajectories
{% capture summary %} Code {% endcapture %} {% capture details %}  
```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  library(cisTopic)
  library(chromVAR)
  library(dplyr)
  library(princurve)

  orgo_cirm43<-readRDS("orgo_cirm43.filt.SeuratObject.Rds")

  #Focus on "radial_glia","intermediate_progenitor","excitatory_neuron"

  #Rerun cistopic on subset of organoid cells
  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(10,20:30),nCores=11,addModels=FALSE)
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          
  lapply(unique(orgo_cirm43$seurat_clusters), function(i) {
    cirm43_subset<-subset(orgo_cirm43,cells=which(orgo_cirm43$seurat_clusters == i))
    cistopic_processing(seurat_input=cirm43_subset,prefix=i)} )

  for (i in unique(orgo_cirm43$seurat_clusters)){
    #Setting up topic count selection
    atac_cistopic_models <-readRDS(paste(i,"CisTopicObject.Rds",sep="."))
    pdf(paste(i,"CisTopicObject_model_selection.pdf",sep="."))
    par(mfrow=c(1,3))
    atac_cistopic_models <- selectModel(atac_cistopic_models, type='derivative')
    dev.off()
    system(paste0("slack -F ",paste(i,"CisTopicObject_model_selection.pdf",sep=".")," ryan_todo"))}

  #set topics based on derivative
  topic_counts<-setNames(c(27,23,28,24),c("3","1","0","2"))

  topicmodel_list<-setNames(c(27,23,28,24), paste(c("3","1","0","2"),"CisTopicObject.Rds",sep="."))

#UMAP Projection and clustering on selected cistopic model
clustering_loop<-function(topicmodel_list.=topicmodel_list,object_input=orgo_cirm43,celltype.x,topiccount_list.=topic_counts){
    #set up subset object again
    atac_sub<-subset(object_input,cells=which(orgo_cirm43$seurat_clusters == celltype.x)) 
    #select_topic
    models_input<-readRDS(paste(celltype.x,"CisTopicObject.Rds",sep="."))
    cisTopicObject<-cisTopic::selectModel(models_input,select=topiccount_list.[celltype.x],keepModels=F)
    
    #perform UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
    dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
    row.names(dims)<-colnames(topic_df)
    colnames(dims)<-c("x","y")
    dims$cellID<-row.names(dims)
    dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")

    #get cell embeddings
    cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
    colnames(cell_embeddings)<-cisTopicObject@cell.names
    n_topics<-nrow(cell_embeddings)
    row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
    cell_embeddings<-as.data.frame(t(cell_embeddings))
    
    #get feature loadings
    feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
    row.names(feature_loadings)<-paste0("topic_",1:n_topics)
    feature_loadings<-as.data.frame(t(feature_loadings))
    
    #combine with seurat object for celltype seuratobject  
    cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
    umap_dims<-as.data.frame(as.matrix(dims[2:3]))
    colnames(umap_dims)<-c("UMAP_1","UMAP_2")
    row.names(umap_dims)<-dims$cellID
    cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
    atac_sub@reductions$cistopic<-cistopic_obj
    atac_sub@reductions$umap<-cistopic_umap
    #finally recluster
    n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic"))
    #Clustering with multiple resolutions to account for different celltype complexities
    atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics)
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.1,graph.name="peaks_snn")
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.2,graph.name="peaks_snn")
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.5,graph.name="peaks_snn")
    atac_sub <- FindClusters(object = atac_sub,verbose = TRUE,resolution=0.9,graph.name="peaks_snn")
    
    saveRDS(atac_sub,paste("./subcluster/",celltype.x,"SeuratObject.Rds",sep="_"))
    plt<-DimPlot(atac_sub,group.by=c('peaks_snn_res.0.1','peaks_snn_res.0.2','peaks_snn_res.0.5','peaks_snn_res.0.9'))
    ggsave(plt,file=paste("./subcluster/",celltype.x,"clustering.pdf",sep="_"))
    system(paste0("slack -F ",paste("./subcluster/",celltype.x,"clustering.pdf",sep="_")," ryan_todo"))

}

dir.create("./subcluster")

lapply(unique(orgo_cirm43$seurat_clusters), function(i) {clustering_loop(celltype.x=i)})

for (i in unique(orgo_cirm43$seurat_clusters)){system(paste0("slack -F ",paste("./subcluster/",i,"clustering.pdf",sep="_")," ryan_todo"))}


subcluster_assignment<-function(resolution_list.=resolution_list,object_input=orgo_cirm43,celltype.x){
    #set up subset object again
    atac_sub<-readRDS(paste("./subcluster/",celltype.x,"SeuratObject.Rds",sep="_"))
    atac_sub$seurat_subcluster<-atac_sub@meta.data[,which(colnames(atac_sub@meta.data)==resolution_list[celltype.x])]
    saveRDS(atac_sub,paste("./subcluster/",celltype.x,"SeuratObject.Rds",sep="_"))

}
#Selected resolution based on subclustering plots
resolution_list<-setNames(c("peaks_snn_res.0.1","peaks_snn_res.0.2","peaks_snn_res.0.2","peaks_snn_res.0.2"),unique(orgo_cirm43$seurat_clusters))
#Seurat subcluster assignment based on resolution plots
for (i in unique(orgo_cirm43$seurat_clusters)){subcluster_assignment(celltype.x=i)}

#Add subcluster assignment and subcluster x and y coordinates to main Rds
#Assign clustering resolution based on clustering.pdf output
cell_order<-row.names(orgo_cirm43@meta.data)
celltype_list<-unique(orgo_cirm43$seurat_clusters)

#adding all subclustering info back into main RDS object
orgo_cirm43$seurat_subcluster<-"NA"
orgo_cirm43$subcluster_x<-"NA"
orgo_cirm43$subcluster_y<-"NA"

for(celltype.x in celltype_list){
    atac_sub<-readRDS(paste("./subcluster/",celltype.x,"SeuratObject.Rds",sep="_"))
    embedding_order<-row.names(atac_sub@reductions$umap@cell.embeddings)
    row_order<-match(cell_order,embedding_order,nomatch=0)
    orgo_cirm43@meta.data[match(row.names(atac_sub@meta.data),row.names(orgo_cirm43@meta.data),nomatch=0),]$seurat_subcluster<-as.character(unlist(atac_sub@meta.data[which(colnames(atac_sub@meta.data)==resolution_list[celltype.x])]))
    orgo_cirm43@meta.data[match(row.names(atac_sub@meta.data),row.names(orgo_cirm43@meta.data),nomatch=0),]$subcluster_x<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,1])
    orgo_cirm43@meta.data[match(row.names(atac_sub@meta.data),row.names(orgo_cirm43@meta.data),nomatch=0),]$subcluster_y<-as.numeric(atac_sub@reductions$umap@cell.embeddings[row_order,2])
}

as.data.frame(orgo_cirm43@meta.data %>% group_by(seurat_clusters,seurat_subcluster)%>% summarize(count=n()))
"""   seurat_clusters seurat_subcluster count
1                0                 0  3422
2                0                 1  3127
3                0                 2  2561
4                0                 3  2123
5                0                 4  1940
6                0                 5  1590
7                1                 0  4760
8                1                 1  4336
9                1                 2  2384
10               1                 3  1043
11               1                 4   521
12               1                 5   380
13               1                 6   246
14               2                 0  1322
15               2                 1  1228
16               2                 2   303
17               2                 3   103
18               3                 0   653
19               3                 1   394
"""

  #set up unique cluster IDs by cluster and subcluster
  orgo_cirm43$cluster_ID<-paste(orgo_cirm43$seurat_clusters,orgo_cirm43$seurat_subcluster,sep="_")
saveRDS(orgo_cirm43,"orgo_cirm43.filt.SeuratObject.Rds")


```
{% endcapture %} {% include details.html %} 

### Differential Accessibillity on Clusters

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(parallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)

  orgo_cirm43<-readRDS("orgo_cirm43.filt.SeuratObject.Rds")

  DefaultAssay(orgo_cirm43)<-"peaks"
  #Perform One vs. rest DA enrichment
  #set up an empty list for looping through
  cirm43_da_peaks<-list()

  #define DA functions for parallelization
  #Use LR test for atac data
  da_one_v_rest<-function(i,obj,group){
      da_peaks_tmp <- FindMarkers(
          object = obj,
          ident.1 = i,
          group.by = group,
          test.use = 'LR',
          latent.vars = 'nCount_peaks',
          only.pos=T
          )
      da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
      closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
      da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
      da_peaks_tmp$enriched_group<-c(i)
      da_peaks_tmp$compared_group<-c("all_other_cells")
      return(da_peaks_tmp)
    }

  da_one_v_one<-function(i,obj,group,j_list){
      i<-as.character(i)
      da_tmp_2<-list()
      for (j in j_list){
          if ( i != j){
          da_peaks_tmp <- FindMarkers(
              object = obj,
              ident.1 = i,
              ident.2 = j,
              group.by = group,
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              only.pos=T
              )
          da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
          closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
          da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
          da_peaks_tmp$enriched_group<-c(i)
          da_peaks_tmp$compared_group<-c(j)
          da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
          }
      }
      return(da_tmp_2)
    }

  #Perform parallel application of DA test
  library(parallel)

  n.cores=length(unique(orgo_cirm43$cluster_ID))
  cirm43_da_peaks<-mclapply(
      unique(orgo_cirm43$cluster_ID),
      FUN=da_one_v_rest,
      obj=orgo_cirm43,
      group="cluster_ID",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1vrest DA
  cirm43_da_peaks<-do.call("rbind",cirm43_da_peaks)

  write("Outputting One v Rest DA Table.", stderr())
  write.table(cirm43_da_peaks,file="cirm43.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

  #Plot out top peaks and associated gene name for each cluster
  dat<-read.table("cirm43.onevrest.da_peaks.txt",header=T,sep="\t")

  dat$label<-c("")
  for (i in unique(dat$enriched_group)){
    selc_genes<-row.names(dat %>% filter(enriched_group==i) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:3))
    dat[row.names(dat) %in% selc_genes & dat$enriched_group==i,]$label<- dat[row.names(dat) %in% selc_genes & dat$enriched_group==i,]$gene_name
  }


  plt<-ggplot(dat,aes(x=pct.1/pct.2,y=(-log(p_val_adj)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(aes(label=label),size=3,force=15)+theme_bw()+xlim(c(0,35))
  ggsave(plt,file="cirm43_da_peaks.pdf")
  system("slack -F cirm43_da_peaks.pdf ryan_todo")

  #Empty list to rerun for 1v1 comparisons
  cirm43_da_peaks<-list()

  n.cores=length(unique(orgo_cirm43@meta.data$seurat_clusters))
  cirm43_da_peaks<-mclapply(
      unique(orgo_cirm43@meta.data$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      j_list=do.call("as.character",list(unique(orgo_cirm43@meta.data$seurat_clusters))),
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1v1 DA
  cirm43_da_peaks<-do.call("rbind",do.call("rbind",cirm43_da_peaks))

  write("Outputting One v One DA Table.", stderr())
  write.table(cirm43_da_peaks,file="cirm43.onevone.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

```
{% endcapture %} {% include details.html %} 



### Differential Motif Accessibility

{% capture summary %} Code {% endcapture %} {% capture details %}  

```R
  ###Differential TF Accessibility by cluster###
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(parallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)

  orgo_cirm43<-readRDS("orgo_cirm43.filt.SeuratObject.Rds")


  #Perform One vs. rest DA enrichment

  write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

  DefaultAssay(orgo_cirm43) <- 'chromvar'

  #set up an empty list for looping through
  cirm43_tf<-list()

  #define DA functions for parallelization
  #Use LR test for atac data
  da_one_v_rest<-function(i,obj,group){
      da_peaks_tmp <- FindMarkers(
          object = obj,
          ident.1 = i,
          group.by = group,
          test.use = 'LR',
          latent.vars = 'nCount_peaks',
          only.pos=T
          )
      da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
      da_peaks_tmp$enriched_group<-c(i)
      da_peaks_tmp$compared_group<-c("all_other_cells")
      return(da_peaks_tmp)
    }

  da_one_v_one<-function(i,obj,group,j_list){
      i<-as.character(i)
      da_tmp_2<-list()
      for (j in j_list){
          if ( i != j){
          da_peaks_tmp <- FindMarkers(
              object = obj,
              ident.1 = i,
              ident.2 = j,
              group.by = group,
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              only.pos=T
              )
          da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
          da_peaks_tmp$enriched_group<-c(i)
          da_peaks_tmp$compared_group<-c(j)
          da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
          }
      }
      return(da_tmp_2)
    }

  #Perform parallel application of DA test
  library(parallel)
  n.cores=length(unique(orgo_cirm43$cluster_ID))
  cirm43_tf<-mclapply(
      unique(orgo_cirm43$cluster_ID),
      FUN=da_one_v_rest,
      obj=orgo_cirm43,
      group="cluster_ID",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1vrest DA
  cirm43_tf<-do.call("rbind",cirm43_tf)

  write("Outputting One v Rest DA Table.", stderr())
  write.table(cirm43_tf,file="cirm43.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

  dat<-read.table("cirm43.onevrest.da_tf.txt",header=T,sep="\t")
  #To convert JASPAR ID TO TF NAME
  dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(dat,file="cirm43.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)
  dat$label<-c("")
  for (i in unique(dat$enriched_group)){
    selc_genes<-row.names(dat %>% filter(enriched_group==i) %>% arrange(rev(desc(p_val_adj))) %>% slice(1:2))
    dat[row.names(dat) %in% selc_genes & dat$enriched_group==i,]$label<- dat[row.names(dat) %in% selc_genes & dat$enriched_group==i,]$da_tf
  }

  plt<-ggplot(dat,aes(x=pct.1/pct.2,y=(-log(p_val_adj)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(aes(label=label),size=3,force=15)+theme_bw()
  ggsave(plt,file="cirm43_oncevrest.da_tf.pdf")
  system("slack -F cirm43_oncevrest.da_tf.pdf ryan_todo")

  #Empty list to rerun for 1v1 comparisons
  cirm43_tf<-list()
      
  n.cores=length(unique(orgo_cirm43@meta.data$seurat_clusters))
  cirm43_tf<-mclapply(
      unique(orgo_cirm43@meta.data$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      j_list=do.call("as.character",list(unique(orgo_cirm43@meta.data$seurat_clusters))),
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1v1 DA
  cirm43_tf<-do.call("rbind",do.call("rbind",cirm43_tf))

  write("Outputting One v One DA Table.", stderr())
  write.table(cirm43_tf,file="cirm43.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

  dat<-read.table("cirm43.onevone.da_tf.txt",header=T,sep="\t")
  #To convert JASPAR ID TO TF NAME
  dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(dat,file="cirm43.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)
```
{% endcapture %} {% include details.html %} 
-->

