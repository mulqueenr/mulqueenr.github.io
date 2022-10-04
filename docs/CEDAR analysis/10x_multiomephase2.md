---
title: 10X Multiome Tumor All Together
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_phase2/
category: CEDAR
---

Multiome processing for 10X multiome data on Primary Tumors (Phase 1+2 and preliminary experiment combined)
*This code is a WIP and will be cleaned prior to manuscript finalization.*

## File Location
```bash
mkdir /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
sftp mulqueen@nix.ohsu.edu
#enter password
get -r /data/EXP220628HM
get -r /data/EXP220629HM

#and download WGS data
get -r /data/EXP220921HM
```


## Reference data
Using chipseq bed files held on cistrome for accessibility analysis.

```bash
/home/groups/CEDAR/mulqueen/ref/cistrome #stored reference files here, downloaded from site at .tar.gz file
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor_full_QC.txt #has information on each download peaks files
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor #has individual peak files
```
## Run cellranger-mkfastq

Set up indexes used for mkfastq

```bash
echo """Lane,Sample,Index
*,sample_13,SI-NA-B8
*,sample_14,SI-NA-B7
*,sample_15,SI-NA-B6
*,sample_16,SI-NA-B5
*,sample_17,SI-NA-B4
*,sample_18,SI-NA-B3
*,sample_19,SI-NA-B2
*,sample_20,SI-NA-B1""" > multiome_atac_phase2.csv

echo """Lane,Sample,Index
1,sample_13,SI-TT-E3
1,sample_14,SI-TT-F3
1,sample_15,SI-TT-G3
1,sample_16,SI-TT-H3
1,sample_17,SI-TT-E4
1,sample_18,SI-TT-F4
1,sample_19,SI-TT-G4
1,sample_20,SI-TT-H4""" > multiome_rna_phase2.csv #other lane used for nextera libraries

```

Sample Sheet for RNA demultiplexing. Index 2 is workflow b

```
[Header]
EMFileVersion,4
 
[Reads]
28
90
 
[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
1,sample_13,sample_13,ACCAGACAAC,CCTAGTTCCT,phase2,SI-TT-E3
1,sample_14,sample_14,GAGAGGATAT,CCCATTTCAA,phase2,SI-TT-F3
1,sample_15,sample_15,ATGACGTCGC,ATCCTGACCT,phase2,SI-TT-G3
1,sample_16,sample_16,CCCGTTCTCG,CCAATCCGTC,phase2,SI-TT-H3
1,sample_17,sample_17,AACCACGCAT,TAACCTGAAT,phase2,SI-TT-E4
1,sample_18,sample_18,CCCACCACAA,AAGCGGAGGT,phase2,SI-TT-F4
1,sample_19,sample_19,GCGCTTATGG,CTAGCCAGGC,phase2,SI-TT-G4
1,sample_20,sample_20,AGTTTCCTGG,CTGTGTGGCA,phase2,SI-TT-H4
```
Run bcl2fastq
```bash
run_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2"
sample_sheet="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/multiome_rna_samplesheet.csv"
out_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_rna/HCVLCDRX2"

bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
            --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --output-dir=${out_dir}\
            --sample-sheet=${sample_sheet}
```
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2
Run Cellranger-arc
```bash
#conda install -c bih-cubi bcl2fastq2

cellranger-arc mkfastq --id=phase_2_atac \
                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220629HM/220713_A01058_0247_AHFJY3DRX2 \
                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/multiome_atac_phase2.csv \
                     --localcores=20 \
                     --localmem=80

#mkfastq is failing on RNA so I used bcl2fastq above
#cellranger-arc mkfastq --id=phase_2_rna \
#                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220628HM/220713_A01058_0246_BHFMTNDRX2 \
#                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/multiome_rna_phase2.csv \
#                     --localcores=20 \
#                     --barcode-mismatches=2 \
#                     --with-failed-reads \
#                     --lanes=1 \
#                     --localmem=80
```
## Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc.

### RM Libraries

```bash
for i in 13 14 15 16 17 18 19 20; do
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_atac/outs/fastq_path/HFJY3DRX2/sample_"""${i}""",sample_"""${i}""",Chromatin Accessibility
/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase_2_rna/HCVLCDRX2/phase2,sample_"""${i}""",Gene Expression""" > sample_${i}.csv ; done
```

## Run CellRanger-ARC
```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2 #using this as analysis directory
```
### Slurm commands
Test run

```bash
 cellranger-arc testrun --id=tiny2
```
Run Cellranger per sample

```bash          
for i in sample_13.csv sample_14.csv sample_15.csv sample_16.csv sample_17.csv sample_18.csv sample_19.csv sample_20.csv ; do
  outname=${i::-4};
  cellranger-arc count --id=${outname} \
   --reference=/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
   --libraries=/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/${i} \
   --localcores=30 \
   --localmem=90 ; done &
   
#check web summaries
cd /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/
for i in 1 3 4 5 6 7 8 9 10 11 12; do
  cp ./sample_$i/outs/web_summary.html ./sample_$i/outs/$i.web_summary.html
  slack -F ./sample_$i/outs/$i.web_summary.html ryan_todo; done 
```

# Initial QC

## Perform Scrublet on Data to Ensure Single-cells

Code from tutorial here.[https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb]

Reading in h5 python matrices with code from https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb
```bash
pip install scrublet
```
Saving a python script as /home/groups/CEDAR/mulqueen/src/multiome_scrublet.py

```python
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import sys
import pandas as pd
np.random.seed(0)
x=sys.argv[1]

try:
    x = int(x)
    print(x)
except ValueError:
    # Handle the exception
    print('Passing to Preliminary Data:' + x)

if type(x) == int:
  if x < 13:
    input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_"+str(x)+"/outs"
    outname="sample_"+str(x)
  elif x >= 13:
    input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_"+str(x)+"/outs"
    outname="sample_"+str(x) 
else:
  input_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/"+x+"/outs"
  outname=x


#Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.
counts_matrix = scipy.io.mmread(input_dir + '/filtered_feature_bc_matrix/matrix.mtx.gz').T.tocsc()
cellIDs=gzip.open(input_dir + '/filtered_feature_bc_matrix/barcodes.tsv.gz',"rb").read().split()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
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

df = pd.DataFrame({'cellid':cellIDs, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv(input_dir+'/'+outname+'.scrublet.tsv', index=False, sep="\t")
print("Done with sample: "+outname)
print("Saved output to: "+input_dir+'/'+outname+'.scrublet.tsv')
```

```bash
for i in 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 RM_1 RM_2 RM_3 RM_4; 
  do python /home/groups/CEDAR/mulqueen/src/multiome_scrublet.py ${i}; 
  done
```

## Use SoupX to remove ambient RNA

```R
install.packages('SoupX')
```

```R
library(SoupX)

run_soupX_persample<-function(x,y=1){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  sc = load10X(wd)
  sc = autoEstCont(sc,tfidfMin=y) #1 is default
  out = adjustCounts(sc)
  saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
  print(paste("Finished:",outname))
}


lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,16,19,20,"RM_1","RM_2","RM_3","RM_4"),run_soupX_persample)

#14,15,17 and 18  failed to correct

#Sample 15 failed due to high homogeneity, autosupplying contamination fraction based on other samples
wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",15,"/outs")
outname<-paste0("sample_",15)
sc = load10X(wd)
sc = setContaminationFraction(sc, 0.1) #supplying contamination fraction manually based on values seen from other samples, using 10% 
out = adjustCounts(sc)
saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
print(paste("Finished:",outname))

```

# Seurat Generation and Processing

### Seurat Object Generation for Samples
Performing seurat analysis following https://satijalab.org/signac/articles/pbmc_multiomic.html

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# set up sample loop to load the RNA and ATAC data, save to seurat object
setupseurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################

  setwd(wd)
  counts <- Read10X_h5("filtered_feature_bc_matrix.h5") #count data
  fragpath <- "atac_fragments.tsv.gz" #atac fragments
  metadata_cellranger<-read.csv("per_barcode_metrics.csv") #metadata
  row.names(metadata_cellranger)<-metadata_cellranger$barcode
  soupx_output<-readRDS("soupx_corrected_counts.rds") #load SoupX contamination corrected output
  scrublet_output<-read.table(paste0(outname,".scrublet.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
  #clean up scrublet output to add to metadata columns
  #just a hold over from a python output that I'm correcting.
  if(startsWith(scrublet_output$cellid[1],"b")){ 
  scrublet_output$cellID<-unlist(lapply(scrublet_output$cellid, function(x) substr(x,2,nchar(x))))}
  row.names(scrublet_output)<-scrublet_output$cellID
  scrublet_output<-scrublet_output[,c("doublet_scores","predicted_doublets")]

  # create a Seurat object containing the RNA data
  dat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  dat[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  #Create corrected RNA data and add to object
  dat[["SoupXRNA"]]<-CreateAssayObject(
    counts=soupx_output)

  #QC cells
  DefaultAssay(dat) <- "ATAC"
  dat <- NucleosomeSignal(dat)
  dat <- TSSEnrichment(dat)
  dat<-AddMetaData(dat,metadata=metadata_cellranger)
  dat<-AddMetaData(dat,metadata=scrublet_output)

  plt<-VlnPlot(
    object = dat,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
  ggsave(plt,file=paste0(outname,".qc.pdf"))
  system(paste0("slack -F ",outname,".qc.pdf ryan_todo"))
  saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))
}

#generate all seurat objects
lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),setupseurat)

```

Initial Merged Seurat Object from all phases

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".SeuratObject.rds"))
  dat$sample<-outname #set up sample metadata
  return(dat)}

out<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),merge_seurat)


dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"), project = "all_data")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.SeuratObject.rds")
```

### Call Peaks and Dimensionality Reduction

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")

dat<-readRDS("phase2.SeuratObject.rds")
dat
table(dat$sample)
#sample_1 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15 sample_16
#     3523     20000      5575     14071       186        30      1489       576
#sample_17 sample_18 sample_19 sample_20 sample_21 sample_22 sample_23 sample_24
#       19        37      2234       722      1851      1463       876       931
# sample_3  sample_4  sample_5  sample_6  sample_7  sample_8  sample_9
#     7698     15253      1453       666      2656       724      1451

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# call peaks using MACS2
DefaultAssay(dat)<-"ATAC"
peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")
#use this set of peaks for all samples

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

DefaultAssay(dat) <- "ATAC"
saveRDS(peaks,file="combined.peakset.rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = annotation
)

saveRDS(dat,file="phase2.SeuratObject.rds")
```

## Run Dim Reduction Per Sample
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

peaks <- readRDS(file="combined.peakset.rds")

#perform initial clustering, and remove scrublet detected doublets
single_sample_qc<-function(x){
 if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".SeuratObject.rds"))
  dat$sample<-x
  outname<-x
  }
# call peaks using MACS2
DefaultAssay(dat) <- "ATAC"

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = annotation
)

#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)
p1<-DimPlot(dat,reduction="rna_umap")+ggtitle("RNA UMAP")

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)
p2<-DimPlot(dat,reduction="atac_umap")+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted_doublets")+ggtitle("Multimodal UMAP Doublets")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-FeaturePlot(dat,reduction="multimodal_umap",features="doublet_scores")+ggtitle("Multimodal UMAP Scublet Scores")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file=paste0(wd,"/",outname,".umap.pdf"))
system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
table(dat$predicted_doublets)
if((sum(dat@meta.data$predicted_doublets=="True")/ncol(dat))<0.05){
cellids<-row.names(dat@meta.data[dat@meta.data$predicted_doublets=="False",])
}else{
cellids<-row.names(dat@meta.data[dat@meta.data$doublet_scores<quantile(dat@meta.data$doublet_scores,0.95),])
}
saveRDS(cellids,paste0(wd,"/",outname,".cellids.rds"))
dat<-subset(dat,cells=cellids)
saveRDS(dat,file=paste0(wd,"/",outname,".QC.SeuratObject.rds"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_qc)

#then rerun clustering now that they are filtered.
single_sample_qc2<-function(x){
 if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds"))
  dat$sample<-x
  outname<-x
  }
# call peaks using MACS2
DefaultAssay(dat) <- "ATAC"

#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")

p1<-DimPlot(dat,reduction="rna_umap",group.by="seurat_clusters")+ggtitle("RNA UMAP")
p2<-DimPlot(dat,reduction="atac_umap",group.by="seurat_clusters")+ggtitle("ATAC UMAP")
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP")

#Finally Plot results
plt<-(p1 | p2)/(p3)
ggsave(plt,file=paste0(wd,"/",outname,".umap2.pdf"))
system(paste0("slack -F ",paste0(wd,"/",outname,".umap2.pdf")," ryan_todo"))
saveRDS(dat,file=paste0(wd,"/",outname,".QC.SeuratObject.rds"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_qc2)
```


## Run cisTopic for ATAC Dimensionality Reduction

Cistopic Per sample (Updated to include other directory folders)

```bash
nano /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/cistopic_per_sample.R
```
```R
#######CISTOPIC PROCESSING PER CELL LINE#################
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
library(RColorBrewer)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

single_sample_cistopic_generation<-function(x){
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds"))
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds"))
  }
  cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
  saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))

  sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =atac_sub@meta.data)
  pdf(paste0(wd,"/",outname,"_model_selection.pdf"))
  par(mfrow=c(3,3))
  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,"_model_selection.pdf")," ryan_todo"))
  
  saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
  print("finshed running cistopic")

  #Add cell embeddings into seurat
  cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
  colnames(cell_embeddings)<-sub_cistopic_models@cell.names
  n_topics<-nrow(cell_embeddings)
  row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
  cell_embeddings<-as.data.frame(t(cell_embeddings))

  #Add feature loadings into seurat
  feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
  row.names(feature_loadings)<-paste0("topic_",1:n_topics)
  feature_loadings<-as.data.frame(t(feature_loadings))

  #combined cistopic results (cistopic loadings and umap with seurat object)
  cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
  print("Cistopic Loading into Seurat")
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  #plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(wd,"/",outname,".cistopic.umap.pdf"),width=10)
  print(plt1)
  #print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".cistopic.umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".QC.SeuratObject.rds"))
  }
single_sample_cistopic_generation(x=args[1])

#for i in 1 3 4 5 6 7 8 9 10 11 12 15 16 19 20 "RM_1" "RM_2" "RM_3" "RM_4"; do Rscript cistopic_per_sample.R $i; done &

```

# Public Datasets for Comparison 
## Using Transfer Anchors for Cell identification.

Using Swarbrick paper labels for transfer. https://pubmed.ncbi.nlm.nih.gov/34493872/

Download data
```bash
cd /home/groups/CEDAR/mulqueen/ref/swarbrick
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
tar -xvf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```
Make Seurat Object with Metadata

```R
library(Seurat)

setwd("/home/groups/CEDAR/mulqueen/ref/swarbrick")
counts<-ReadMtx(mtx="count_matrix_sparse.mtx",cells="count_matrix_barcodes.tsv",features="count_matrix_genes.tsv",feature.column=1) #sparse matrix of counts
metadata<-read.csv("metadata.csv") #metadata
row.names(metadata)<-metadata$X
# create a Seurat object containing the RNA adata
swarbrick <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)
swarbrick<-AddMetaData(swarbrick,metadata=metadata)
saveRDS(swarbrick,"/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
```

### Using EMBO paper for transfer of signatures as well. 
https://doi.org/10.15252/embj.2020107333
Full code here: https://www.nature.com/articles/s41597-022-01236-2 data available here https://doi.org/10.6084/m9.figshare.17058077

Download data from GEO FTP server

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
wget https://figshare.com/ndownloader/articles/17058077/versions/1
unzip 1
```

# Swarbrick Paper Label Transfer
## Transfer Swarbrick cell types

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Using Label transfer to label cell types by Swarbrick paper
#seurat object made by AD
swarbrick<-readRDS("/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
swarbrick<-NormalizeData(swarbrick)
swarbrick<-FindVariableFeatures(swarbrick)
swarbrick<-ScaleData(swarbrick)

##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat@assays$peaks<-dat@assays$ATAC
  dat$sample<-paste0(x)
  }
  DefaultAssay(dat)<-"SoupXRNA"
  dat<-NormalizeData(dat)
  dat<-FindVariableFeatures(dat)
  dat<-ScaleData(dat)
  saveRDS(dat,file=file_in)

  transfer.anchors <- FindTransferAnchors(
    reference = swarbrick,
    reference.assay="RNA",
    query = dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = swarbrick$celltype_major,
  )

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)
  plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_label_transfer)
#
```

### Add sample metadata

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#from tumor sample information
meta_data_in<-as.data.frame(cbind("sample"=c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"),
  "diagnosis"= c("DCIS", "IDC", "IDC", "IDC", "IDC", "DCIS", "IDC", "IDC", "IDC", "IDC", "IDC", "NAT", "DCIS", "NAT", "IDC", "ILC", "IDC", "IDC", "NAT"), "molecular_type"=c(
"DCIS", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "DCIS", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "NA", "DCIS", "NA", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2+", "ER+/PR+/HER2-", "NA")))

sample_metadata_merged<-function(dat){
  dat_file_path=dat
  file_in=basename(dat)
  dir_in=dirname(dat)
  dat<-readRDS(dat) #read in as seurat object
  print(paste("Read in",dat_file_path))
  saveRDS(dat,paste0(dat_file_path,".backup")) #save a backup RDS file
  print("Made backup file")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$diagnosis,col.name="diagnosis")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$molecular_type,col.name="molecular_type")
  saveRDS(dat,dat_file_path)
  print("Finished Sample")
}

#run umap projections of merged samples
dat="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.QC.SeuratObject.rds"
sample_metadata_merged(dat)


sample_metadata_persample<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  print(paste("Read in",file_in))
  saveRDS(dat,paste0(file_in,".backup")) #save a backup RDS file
  print("Made backup file")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$diagnosis,col.name="diagnosis")
  dat<-AddMetaData(dat,meta_data_in[match(dat$sample,meta_data_in$sample),]$molecular_type,col.name="molecular_type")
  saveRDS(dat,file_in)
  print("Finished Sample")
}
#add metadata to each sample
lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x) sample_metadata_persample(x))


```

### Plot multimodal dimensionality reduction for cistopic embedding

```R
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
library(RColorBrewer)
library(ggplot2)



plot_reductions_persample<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }

  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
  alpha_val=0.33

  #add cistopic reduction to umap projections
  p1<-DimPlot(dat,reduction="rna_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("RNA UMAP")+theme(legend.position="none")
  p2<-DimPlot(dat,reduction="atac_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("ATAC UMAP")+theme(legend.position="none")
  p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (LSI)")+theme(legend.position="none")


  #Try multimodal with cistopic
  dat <- RunUMAP(
    object = dat,
    reduction="cistopic",
    reduction.name="cistopic_umap",
    dims=1:ncol(dat@reductions$cistopic),
    assay = "peaks",
    verbose = TRUE
  )


  # build a joint neighbor graph using both assays
  dat2 <- FindMultiModalNeighbors(
    object = dat,
    reduction.list = list("pca", "cistopic"), 
    dims.list = list(1:50, 1:ncol(dat@reductions$cistopic)), 
    modality.weight.name = "RNA.weight",
    verbose = TRUE
  )

  # build a joint UMAP visualization
  dat2 <- RunUMAP(
    object = dat2,
    nn.name = "weighted.nn",
    reduction.name="multimodal_umap",
    assay = "SoupXRNA",
    verbose = TRUE
  )

  ###########plot cistopic umap too
  p5<-DimPlot(dat2,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  p6<-DimPlot(dat,reduction="cistopic_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Cistopic UMAP")+theme(legend.position="none")

  #Finally Plot results
  plt<-(p1 | p2)/(p3)/(p6 | p5)
  ggsave(plt,file=paste0(wd,"/",outname,".multimodal.umap.pdf"),width=10,height=20)
  system(paste0("slack -F ",wd,"/",outname,".multimodal.umap.pdf ryan_todo"))

  ggsave(p1,file=paste0(wd,"/",outname,".multimodal.umap.rna.pdf"),width=10,height=10)
  system(paste0("slack -F ",wd,"/",outname,".multimodal.umap.rna.pdf ryan_todo"))

  ggsave(p2,file=paste0(wd,"/",outname,".multimodal.umap.atac.pdf"),width=10,height=10)
  system(paste0("slack -F ",wd,"/",outname,".multimodal.umap.atac.pdf ryan_todo"))

  saveRDS(dat,file=file_in)

  plt_cell_count<-dat@meta.data[,c("sample","predicted.id","diagnosis","molecular_type")]
  return(plt_cell_count)
}


#run umap projections of all dim reductions per sample
cell_count<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x) plot_reductions_persample(x))

saveRDS(cell_count,"sample_swarbrick_celltype_assignment.rds") #save nested list of cell type assignment

#plot output of celltype count per sample
out<-readRDS("sample_swarbrick_celltype_assignment.rds")
out<-do.call("rbind",out)
colnames(out)<-c("sample","celltype","diagnosis","molecular_type") #rename just for simplicity
  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
plt1<-ggplot(out,aes(x=sample,fill=celltype))+geom_bar()+theme_minimal()+scale_fill_manual(values=type_cols)+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt2<-ggplot(out,aes(x=sample,fill=celltype))+geom_bar(position="fill")+theme_minimal()+scale_fill_manual(values=type_cols)+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt<-plt1/plt2
ggsave(plt,file="sample_swarbrick_celltype_assignment.pdf")
system("slack -F sample_swarbrick_celltype_assignment.pdf ryan_todo")

library(dplyr)
out_readable<-out %>% count(sample,celltype)
write.table(out_readable,file="sample_swarbrick_celltype_assignment.tsv",col.names=T,row.names=F,sep="\t",quote=F)
system("slack -F sample_swarbrick_celltype_assignment.tsv ryan_todo") #note this was calculated per sample as well as in the merged data set,  the assumption is that they will be the same

```

### EMBO Cell Signature Transfer

Make Seurat Object with Metadata

```R
library(Seurat)

setwd("/home/groups/CEDAR/mulqueen/ref/embo")

# set up sample loop to load the RNA and ATAC data, save to seurat object
files_in<-list.files(pattern="*.rds")
in_rds<-lapply(files_in,function(x) readRDS(x))
dat <- merge(in_rds[[1]], y = as.list(unlist(in_rds[2:length(in_rds)])), project = "embo")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/ref/embo/embo.rds")

table(dat$group)
   #      B1_0023         B1_0033         B1_0090         B1_0894         ER_0001
   #        15736           20063           19347           22287           10764
   #      ER_0025      ER_0029_7C      ER_0029_9C         ER_0032      ER_0040_LN
   #        16400            3102            7261            2057            5661
   #    ER_0040_T         ER_0042      ER_0043_LN       ER_0043_T      ER_0056_LN
   #        17639            9056             745           12054            6611
   #    ER_0056_T      ER_0064_LN       ER_0064_T      ER_0114_T3         ER_0125
   #         2704            2328            3603           13428            8374
   #      ER_0151         ER_0163      ER_0167_LN       ER_0167_T      ER_0173_LN
   #         5774           11942            4392           22660            3257
   #    ER_0173_T         ER_0319         ER_0360       HER2_0031       HER2_0069
   #         7766            6965            8514            8168             910
   #    HER2_0161       HER2_0176       HER2_0308       HER2_0337     mER_0068_LN
   #         9873           21731            8179           17393            4720
   #   mER_0068_T        mER_0178    N_0019_total    N_0021_total      N_0064_epi
   #         1604            7274           17177            2346            8078
   # N_0064_total    N_0092_total      N_0093_epi    N_0093_total      N_0123_epi
   #         4011           11369           16519           16379            5648
   # N_0123_total    N_0169_total   N_0230.16_epi N_0230.17_total    N_0233_total
   #        17162           17348            3616           12873           17502
   #   N_0275_epi    N_0275_total      N_0280_epi    N_0288_total      N_0342_epi
   #        12123            2187            2178            3087           14937
   # N_0342_total      N_0372_epi    N_0372_total      N_0408_epi      N_1105_epi
   #        14736           22243            6480            5971            9797
   #   N_1469_epi         TN_0106      TN_0114_T2         TN_0126         TN_0135
   #         5505            2144            4103            6314           30424
   #   TN_B1_0131      TN_B1_0177      TN_B1_0554      TN_B1_4031
   #        17255           73030           26195           15863

  #B1: BRCA1 preneoplastic (4)
  #ER_*_T: ER+ tumor (16)
  #ER_*_LN: ER+ lymph node (7)
  #HER2: HER2+ tumor(6)
  #mER: male ER+ tumor (2)
  #N: normal pre or post menopausal (if it ends in epi it is limited to epithelial)
  #TN_B1: Triple negative BRCA1 tumor

  #Our data set most closely matches:
  #ER+ tumors (female)
  #ER_0360 ER_0151 ER_0032 ER_0125 ER_0025 ER_0001 ER_0042 ER_0319 ER_0163 ER_0114_T3 ER_0167_T ER_0040_T ER_0043_T
  #HER2+ tumors
  #HER2_0031 HER2_0069 HER2_0161 HER2_0176 HER2_0308 HER2_0337
  #Normal Total Cells (post and pre menopause, epi is limited to epithelial cells)
  #N_0019_total    N_0021_total      N_0064_epi N_0064_total    N_0092_total      N_0093_epi    N_0093_total      N_0123_epi  N_0123_total    N_0169_total   N_0230.16_epi N_0230.17_total    N_0233_total N_0275_epi    N_0275_total      N_0280_epi    N_0288_total      N_0342_epi N_0342_total      N_0372_epi    N_0372_total      N_0408_epi      N_1105_epi N_1469_epi     
```


### Additional Cell Signature
From https://github.com/yunshun/HumanBreast10X/tree/main/Signatures

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
#downloaded files from
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/Human-PosSigGenes.RData
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/ImmuneMarkers2.txt
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/PAM50.txt
```


## Use EMBO and Swarbrick Paper Cell Types to Define Signatures
```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#cell lineage
  load("/home/groups/CEDAR/mulqueen/ref/embo/Human-PosSigGenes.RData")
  ls()
  #[1] "Basal" "LP"    "ML"    "Str" #lineage types
  lineage_in=list(EMBO_Basal=Basal,EMBO_LP=LP,EMBO_ML=ML,EMBO_Str=Str)

#Immune markers
  immune_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/ImmuneMarkers2.txt",header=T,sep="\t")
  immune_in<-lapply(split(immune_in,immune_in$CellType),function(x) x$Signatures)#split up data frame to a named list of genes per cell type
  names(immune_in)<-paste0("EMBO_",names(immune_in))#rename the list just so we can track the source

#using both the given PAM50 short list and the Swarbrick supplied more extensive gene list below
  PAM50_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/PAM50.txt",header=T,sep="\t")
  PAM50_in<-lapply(split(PAM50_in,PAM50_in$Subtype),function(x) x$Gene)#split up data frame to a named list of genes per cell type
  names(PAM50_in)<-paste0("PAM50_",names(PAM50_in))
  features_in=c(immune_in,PAM50_in)   

#SCSubtype Features determined by Swarbrick manuscript (Supp Table 4)
  module_feats<-list()
  module_feats[["Basal_SC"]]=c('EMP1', 'TAGLN', 'TTYH1', 'RTN4', 'TK1', 'BUB3', 'IGLV3.25', 'FAM3C', 'TMEM123', 'KDM5B', 'KRT14', 'ALG3', 'KLK6', 'EEF2', 'NSMCE4A', 'LYST', 'DEDD', 'HLA.DRA', 'PAPOLA', 'SOX4', 'ACTR3B', 'EIF3D', 'CACYBP', 'RARRES1', 'STRA13', 'MFGE8', 'FRZB', 'SDHD', 'UCHL1', 'TMEM176A', 'CAV2', 'MARCO', 'P4HB', 'CHI3L2', 'APOE', 'ATP1B1', 'C6orf15', 'KRT6B', 'TAF1D', 'ACTA2', 'LY6D', 'SAA2', 'CYP27A1', 'DLK1', 'IGKV1.5', 'CENPW', 'RAB18', 'TNFRSF11B', 'VPS28', 'HULC', 'KRT16', 'CDKN2A', 'AHNAK2', 'SEC22B', 'CDC42EP1', 'HMGA1', 'CAV1', 'BAMBI', 'TOMM22', 'ATP6V0E2', 'MTCH2', 'PRSS21', 'HDAC2', 'ZG16B', 'GAL', 'SCGB1D2', 'S100A2', 'GSPT1', 'ARPC1B', 'NIT1', 'NEAT1', 'DSC2', 'RP1.60O19.1', 'MAL2', 'TMEM176B', 'CYP1B1', 'EIF3L', 'FKBP4', 'WFDC2', 'SAA1', 'CXCL17', 'PFDN2', 'UCP2', 'RAB11B', 'FDCSP', 'HLA.DPB1', 'PCSK1N', 'C4orf48', 'CTSC')
  module_feats[["Her2E_SC"]]=c('PSMA2', 'PPP1R1B', 'SYNGR2', 'CNPY2', 'LGALS7B', 'CYBA', 'FTH1', 'MSL1', 'IGKV3.15', 'STARD3', 'HPD', 'HMGCS2', 'ID3', 'NDUFB8', 'COTL1', 'AIM1', 'MED24', 'CEACAM6', 'FABP7', 'CRABP2', 'NR4A2', 'COX14', 'ACADM', 'PKM', 'ECH1', 'C17orf89', 'NGRN', 'ATG5', 'SNHG25', 'ETFB', 'EGLN3', 'CSNK2B', 'RHOC', 'PSENEN', 'CDK12', 'ATP5I', 'ENTHD2', 'QRSL1', 'S100A7', 'TPM1', 'ATP5C1', 'HIST1H1E', 'LGALS1', 'GRB7', 'AQP3', 'ALDH2', 'EIF3E', 'ERBB2', 'LCN2', 'SLC38A10', 'TXN', 'DBI', 'RP11.206M11.7', 'TUBB', 'CRYAB', 'CD9', 'PDSS2', 'XIST', 'MED1', 'C6orf203', 'PSMD3', 'TMC5', 'UQCRQ', 'EFHD1', 'BCAM', 'GPX1', 'EPHX1', 'AREG', 'CDK2AP2', 'SPINK8', 'PGAP3', 'NFIC', 'THRSP', 'LDHB', 'MT1X', 'HIST1H4C', 'LRRC26', 'SLC16A3', 'BACE2', 'MIEN1', 'AR', 'CRIP2', 'NME1', 'DEGS2', 'CASC3', 'FOLR1', 'SIVA1', 'SLC25A39', 'IGHG1', 'ORMDL3', 'KRT81', 'SCGB2B2', 'LINC01285', 'CXCL8', 'KRT15', 'RSU1', 'ZFP36L2', 'DKK1', 'TMED10', 'IRX3', 'S100A9', 'YWHAZ')
  module_feats[["LumA_SC"]]=c('SH3BGRL', 'HSPB1', 'PHGR1', 'SOX9', 'CEBPD', 'CITED2', 'TM4SF1', 'S100P', 'KCNK6', 'AGR3', 'MPC2', 'CXCL13', 'RNASET2', 'DDIT4', 'SCUBE2', 'KRT8', 'MZT2B', 'IFI6', 'RPS26', 'TAGLN2', 'SPTSSA', 'ZFP36L1', 'MGP', 'KDELR2', 'PPDPF', 'AZGP1', 'AP000769.1', 'MYBPC1', 'S100A1', 'TFPI2', 'JUN', 'SLC25A6', 'HSP90AB1', 'ARF5', 'PMAIP1', 'TNFRSF12A', 'FXYD3', 'RASD1', 'PYCARD', 'PYDC1', 'PHLDA2', 'BZW2', 'HOXA9', 'XBP1', 'AGR2', 'HSP90AA1') 
  module_feats[["LumB_SC"]]=c('UGCG', 'ARMT1', 'ISOC1', 'GDF15', 'ZFP36', 'PSMC5', 'DDX5', 'TMEM150C', 'NBEAL1', 'CLEC3A', 'GADD45G', 'MARCKS', 'FHL2', 'CCDC117', 'LY6E', 'GJA1', 'PSAP', 'TAF7', 'PIP', 'HSPA2', 'DSCAM.AS1', 'PSMB7', 'STARD10', 'ATF3', 'WBP11', 'MALAT1', 'C6orf48', 'HLA.DRB1', 'HIST1H2BD', 'CCND1', 'STC2', 'NR4A1', 'NPY1R', 'FOS', 'ZFAND2A', 'CFL1', 'RHOB', 'LMNA', 'SLC40A1', 'CYB5A', 'SRSF5', 'SEC61G', 'CTSD', 'DNAJC12', 'IFITM1', 'MAGED2', 'RBP1', 'TFF1', 'APLP2', 'TFF3', 'TRH', 'NUPR1', 'EMC3', 'TXNIP', 'ARPC4', 'KCNE4', 'ANPEP', 'MGST1', 'TOB1', 'ADIRF', 'TUBA1B', 'MYEOV2', 'MLLT4', 'DHRS2', 'IFITM2')
  module_feats[["proliferation_score"]]<-c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")

#Swarbrick Gene Module Classification (Supp Table 5)
gene_module<-list()
  gene_module[["gene_module_1"]]<-c('ATF3', 'JUN', 'NR4A1', 'IER2', 'DUSP1', 'ZFP36', 'JUNB', 'FOS', 'FOSB', 'PPP1R15A', 'KLF6', 'DNAJB1', 'EGR1', 'BTG2', 'HSPA1B', 'HSPA1A', 'RHOB', 'CLDN4', 'MAFF', 'GADD45B', 'IRF1', 'EFNA1', 'SERTAD1', 'TSC22D1', 'CEBPD', 'CCNL1', 'TRIB1', 'MYC', 'ELF3', 'LMNA', 'NFKBIA', 'TOB1', 'HSPB1', 'BRD2', 'MCL1', 'PNRC1', 'IER3', 'KLF4', 'ZFP36L2', 'SAT1', 'ZFP36L1', 'DNAJB4', 'PHLDA2', 'NEAT1', 'MAP3K8', 'GPRC5A', 'RASD1', 'NFKBIZ', 'CTD-3252C9.4', 'BAMBI', 'RND1', 'HES1', 'PIM3', 'SQSTM1', 'HSPH1', 'ZFAND5', 'AREG', 'CD55', 'CDKN1A', 'UBC', 'CLDN3', 'DDIT3', 'BHLHE40', 'BTG1', 'ANKRD37', 'SOCS3', 'NAMPT', 'SOX4', 'LDLR', 'TIPARP', 'TM4SF1', 'CSRNP1', 'GDF15', 'ZFAND2A', 'NR4A2', 'ERRFI1', 'RAB11FIP1', 'TRAF4', 'MYADM', 'ZC3H12A', 'HERPUD1', 'CKS2', 'BAG3', 'TGIF1', 'ID3', 'JUND', 'PMAIP1', 'TACSTD2', 'ETS2', 'DNAJA1', 'PDLIM3', 'KLF10', 'CYR61', 'MXD1', 'TNFAIP3', 'NCOA7', 'OVOL1', 'TSC22D3', 'HSP90AA1', 'HSPA6', 'C15orf48', 'RHOV', 'DUSP4', 'B4GALT1', 'SDC4', 'C8orf4', 'DNAJB6', 'ICAM1', 'DNAJA4', 'MRPL18', 'GRB7', 'HNRNPA0', 'BCL3', 'DUSP10', 'EDN1', 'FHL2', 'CXCL2', 'TNFRSF12A', 'S100P', 'HSPB8', 'INSIG1', 'PLK3', 'EZR', 'IGFBP5', 'SLC38A2', 'DNAJB9', 'H3F3B', 'TPM4', 'TNFSF10', 'RSRP1', 'ARL5B', 'ATP1B1', 'HSPA8', 'IER5', 'SCGB2A1', 'YPEL2', 'TMC5', 'FBXO32', 'MAP1LC3B', 'MIDN', 'GADD45G', 'VMP1', 'HSPA5', 'SCGB2A2', 'TUBA1A', 'WEE1', 'PDK4', 'STAT3', 'PERP', 'RBBP6', 'KCNQ1OT1', 'OSER1', 'SERP1', 'UBE2B', 'HSPE1', 'SOX9', 'MLF1', 'UBB', 'MDK', 'YPEL5', 'HMGCS1', 'PTP4A1', 'WSB1', 'CEBPB', 'EIF4A2', 'S100A10', 'ELMSAN1', 'ISG15', 'CCNI', 'CLU', 'TIMP3', 'ARL4A', 'SERPINH1', 'SCGB1D2', 'UGDH', 'FUS', 'BAG1', 'IFRD1', 'TFF1', 'SERTAD3', 'IGFBP4', 'TPM1', 'PKIB', 'MALAT1', 'XBP1', 'HEBP2', 'GEM', 'EGR2', 'ID2', 'EGR3', 'HSPD1', 'GLUL', 'DDIT4', 'CDC42EP1', 'RBM39', 'MT-ND5', 'CSNK1A1', 'SLC25A25', 'PEG10', 'DEDD2')

gene_module[["gene_module_2"]]<-c('AZGP1', 'ATP5C1', 'ATP5F1', 'NHP2', 'MGP', 'RPN2', 'C14orf2', 'NQO1', 'REEP5', 'SSR2', 'NDUFA8', 'ATP5E', 'SH3BGRL', 'PIP', 'PRDX2', 'RAB25', 'EIF3L', 'PRDX1', 'USMG5', 'DAD1', 'SEC61G', 'CCT3', 'NDUFA4', 'APOD', 'CHCHD10', 'DDIT4', 'MRPL24', 'NME1', 'DCXR', 'NDUFAB1', 'ATP5A1', 'ATP5B', 'ATOX1', 'SLC50A1', 'POLR2I', 'TIMM8B', 'VPS29', 'TIMP1', 'AHCY', 'PRDX3', 'RBM3', 'GSTM3', 'ABRACL', 'RBX1', 'PAFAH1B3', 'AP1S1', 'RPL34', 'ATPIF1', 'PGD', 'CANX', 'SELENBP1', 'ATP5J', 'PSME2', 'PSME1', 'SDHC', 'AKR1A1', 'GSTP1', 'RARRES3', 'ISCU', 'NPM1', 'SPDEF', 'BLVRB', 'NDUFB3', 'RPL36A', 'MDH1', 'MYEOV2', 'MAGED2', 'CRIP2', 'SEC11C', 'CD151', 'COPE', 'PFN2', 'ALDH2', 'SNRPD2', 'TSTD1', 'RPL13A', 'HIGD2A', 'NDUFC1', 'PYCARD', 'FIS1', 'ITM2B', 'PSMB3', 'G6PD', 'CST3', 'SH3BGRL3', 'TAGLN2', 'NDUFA1', 'TMEM183A', 'S100A10', 'NGFRAP1', 'DEGS2', 'ARPC5', 'TM7SF2', 'RPS10', 'LAMTOR5', 'TMEM256', 'UQCRB', 'TMEM141', 'KRTCAP2', 'HM13', 'NDUFS6', 'PARK7', 'PSMD4', 'NDUFB11', 'TOMM7', 'EIF6', 'UQCRHL', 'ADI1', 'VDAC1', 'C9orf16', 'ETFA', 'LSM3', 'UQCRH', 'CYB5A', 'SNRPE', 'BSG', 'SSR3', 'DPM3', 'LAMTOR4', 'RPS11', 'FAM195A', 'TMEM261', 'ATP5I', 'EIF5A', 'PIN4', 'ATXN10', 'ATP5G3', 'ARPC3', 'UBA52', 'BEX4', 'ROMO1', 'SLC25A6', 'SDCBP', 'EIF4EBP1', 'PFDN6', 'PSMA3', 'RNF7', 'SPCS2', 'CYSTM1', 'CAPG', 'CD9', 'GRHPR', 'SEPP1', 'ESF1', 'TFF3', 'ARPC1B', 'ANXA5', 'WDR83OS', 'LYPLA1', 'COMT', 'MDH2', 'DNPH1', 'RAB13', 'EIF3K', 'PTGR1', 'LGALS3', 'TPI1', 'COPZ1', 'LDHA', 'PSMD8', 'EIF2S3', 'NME3', 'EIF3E', 'MRPL13', 'ZFAND6', 'FAM162A', 'ATP6V0E1', 'TMED10', 'HNRNPA3', 'PPA1', 'SNX17', 'APOA1BP', 'TUFM', 'ECHS1', 'GLTSCR2', 'RPS27L', 'NDUFB1', 'SSBP1', 'PRDX6', 'ENO1', 'PPP4C', 'COA3', 'TCEAL4', 'MRPL54', 'LAMTOR2', 'PAIP2', 'DAP', 'RPL22L1', 'C6orf203', 'TECR', 'PEBP1', 'TMED9', 'ATP6V1F', 'ESD', 'EIF3I', 'SCO2', 'ATP5D', 'UAP1', 'TMEM258', 'COX17')

gene_module[["gene_module_3"]]<-c('HLA-B', 'HLA-A', 'VIM', 'CD74', 'SRGN', 'HLA-C', 'IFI27', 'HLA-E', 'IFITM1', 'PSMB9', 'RGCC', 'S100A4', 'HLA-DRA', 'ISG15', 'IL32', 'SPARC', 'TAGLN', 'IFITM3', 'IFITM2', 'IGFBP7', 'CALD1', 'HLA-DPB1', 'HLA-DPA1', 'B2M', 'TIMP1', 'RGS1', 'FN1', 'ACTA2', 'HLA-DRB1', 'SERPING1', 'ANXA1', 'TPM2', 'TMSB4X', 'CD69', 'CCL4', 'LAPTM5', 'GSN', 'APOE', 'STAT1', 'SPARCL1', 'IFI6', 'DUSP1', 'CXCR4', 'CCL5', 'UBE2L6', 'MYL9', 'SLC2A3', 'BST2', 'CAV1', 'CD52', 'ZFP36L2', 'HLA-DQB1', 'PDLIM1', 'TNFAIP3', 'CORO1A', 'RARRES3', 'TYMP', 'C1S', 'PTRF', 'PSME2', 'CYTIP', 'COL1A1', 'PSMB8', 'NNMT', 'HLA-DQA1', 'DUSP2', 'COL1A2', 'ARHGDIB', 'COL6A2', 'FOS', 'CCL2', 'BGN', 'ID3', 'TUBA1A', 'RAC2', 'LBH', 'HLA-DRB5', 'FCER1G', 'GBP1', 'C1QA', 'COTL1', 'LUM', 'MYL6', 'GBP2', 'BTG1', 'CD37', 'HCST', 'LIMD2', 'IFIT3', 'IL7R', 'PTPRC', 'NKG7', 'FYB', 'TAP1', 'LTB', 'S100A6', 'COL3A1', 'EMP3', 'A2M', 'JUNB', 'TPM1', 'FABP4', 'TXNIP', 'SAT1', 'FXYD5', 'CD3E', 'HLA-DMA', 'CTSC', 'TSC22D3', 'MYL12A', 'CST3', 'CNN2', 'PHLDA1', 'LYZ', 'IFI44L', 'MARCKS', 'ID1', 'DCN', 'TGFBI', 'BIRC3', 'THY1', 'LGALS1', 'GPX1', 'C1QB', 'CD2', 'CST7', 'COL6A3', 'ACAP1', 'IFI16', 'ITM2B', 'POSTN', 'LDHB', 'FLNA', 'FILIP1L', 'CDKN1A', 'IRF1', 'LGALS3', 'SERPINH1', 'EFEMP1', 'PSME1', 'SH3BGRL3', 'IL2RG', 'CD3D', 'SFRP2', 'TIMP3', 'ALOX5AP', 'GMFG', 'CYBA', 'TAGLN2', 'LAP3', 'RGS2', 'CLEC2B', 'TRBC2', 'NR4A2', 'S100A8', 'PSMB10', 'OPTN', 'CTSB', 'FTL', 'KRT17', 'AREG', 'MYH9', 'MMP7', 'COL6A1', 'GZMA', 'RNASE1', 'PCOLCE', 'PTN', 'PYCARD', 'ARPC2', 'SGK1', 'COL18A1', 'GSTP1', 'NPC2', 'SOD3', 'MFGE8', 'COL4A1', 'ADIRF', 'HLA-F', 'CD7', 'APOC1', 'TYROBP', 'C1QC', 'TAPBP', 'STK4', 'RHOH', 'RNF213', 'SOD2', 'TPM4', 'CALM1', 'CTGF', 'PNRC1', 'CD27', 'CD3G', 'PRKCDBP', 'PARP14', 'IGKC', 'IGFBP5', 'IFIT1', 'LY6E')

gene_module[["gene_module_4"]]<-c('STMN1', 'H2AFZ', 'UBE2C', 'TUBA1B', 'BIRC5', 'HMGB2', 'ZWINT', 'TUBB', 'HMGB1', 'DEK', 'CDK1', 'HMGN2', 'UBE2T', 'TK1', 'RRM2', 'RANBP1', 'TYMS', 'CENPW', 'MAD2L1', 'CKS2', 'CKS1B', 'NUSAP1', 'TUBA1C', 'PTTG1', 'KPNA2', 'PCNA', 'CENPF', 'HIST1H4C', 'CDKN3', 'UBE2S', 'CCNB1', 'HMGA1', 'DTYMK', 'SNRPB', 'CDC20', 'NASP', 'MCM7', 'PLP2', 'TUBB4B', 'PLK1', 'CCNB2', 'MKI67', 'TOP2A', 'TPX2', 'PKMYT1', 'PRC1', 'SMC4', 'CENPU', 'RAN', 'DUT', 'PA2G4', 'BUB3', 'RAD21', 'SPC25', 'HN1', 'CDCA3', 'H2AFV', 'HNRNPA2B1', 'CCNA2', 'PBK', 'LSM5', 'DNAJC9', 'RPA3', 'TMPO', 'SNRPD1', 'CENPA', 'KIF20B', 'USP1', 'H2AFX', 'PPM1G', 'NUF2', 'SNRPG', 'KIF22', 'KIAA0101', 'DEPDC1', 'RNASEH2A', 'MT2A', 'STRA13', 'ANLN', 'CACYBP', 'NCL', 'NUDT1', 'ECT2', 'LSM4', 'ASF1B', 'CENPN', 'TMEM106C', 'CCT5', 'HSPA8', 'HMMR', 'SRSF3', 'AURKB', 'GGH', 'AURKA', 'TRIP13', 'CDCA8', 'HMGB3', 'HNRNPAB', 'FAM83D', 'CDC25B', 'GGCT', 'KNSTRN', 'CCT6A', 'PTGES3', 'ANP32E', 'CENPK', 'MCM3', 'DDX21', 'HSPD1', 'SKA2', 'CALM2', 'UHRF1', 'HINT1', 'ORC6', 'MZT1', 'MIS18BP1', 'WDR34', 'NAP1L1', 'TEX30', 'SFN', 'HSPE1', 'CENPM', 'TROAP', 'CDCA5', 'RACGAP1', 'SLC25A5', 'ATAD2', 'DBF4', 'KIF23', 'CEP55', 'SIVA1', 'SAC3D1', 'PSIP1', 'CLSPN', 'CCT2', 'DLGAP5', 'PSMA4', 'SMC2', 'AP2S1', 'RAD51AP1', 'MND1', 'ILF2', 'DNMT1', 'NUCKS1', 'LMNB1', 'RFC4', 'EIF5A', 'NPM3', 'ARL6IP1', 'ASPM', 'GTSE1', 'TOMM40', 'HNRNPA1', 'GMNN', 'FEN1', 'CDCA7', 'SLBP', 'TNFRSF12A', 'TM4SF1', 'CKAP2', 'CENPE', 'SRP9', 'DDX39A', 'COMMD4', 'RBM8A', 'CALM3', 'RRM1', 'ENO1', 'ANP32B', 'SRSF7', 'FAM96A', 'TPRKB', 'FABP5', 'PPIF', 'SERPINE1', 'TACC3', 'RBBP7', 'NEK2', 'CALM1', 'GMPS', 'EMP2', 'HMG20B', 'SMC3', 'HSPA9', 'NAA20', 'NUDC', 'RPL39L', 'PRKDC', 'CDCA4', 'HIST1H1A', 'HES6', 'SUPT16H', 'PTMS', 'VDAC3', 'PSMC3', 'ATP5G1', 'PSMA3', 'PGP', 'KIF2C', 'CARHSP1')

gene_module[["gene_module_5"]]<-c('GJA1', 'SCGB2A2', 'ARMT1', 'MAGED2', 'PIP', 'SCGB1D2', 'CLTC', 'MYBPC1', 'PDZK1', 'MGP', 'SLC39A6', 'CCND1', 'SLC9A3R1', 'NAT1', 'SUB1', 'CYP4X1', 'STC2', 'CROT', 'CTSD', 'FASN', 'PBX1', 'SLC4A7', 'FOXA1', 'MCCC2', 'IDH1', 'H2AFJ', 'CYP4Z1', 'IFI27', 'TBC1D9', 'ANPEP', 'DHRS2', 'TFF3', 'LGALS3BP', 'GATA3', 'LTF', 'IFITM2', 'IFITM1', 'AHNAK', 'SEPP1', 'ACADSB', 'PDCD4', 'MUCL1', 'CERS6', 'LRRC26', 'ASS1', 'SEMA3C', 'APLP2', 'AMFR', 'CDV3', 'VTCN1', 'PREX1', 'TP53INP1', 'LRIG1', 'ANK3', 'ACLY', 'CLSTN1', 'GNB1', 'C1orf64', 'STARD10', 'CA12', 'SCGB2A1', 'MGST1', 'PSAP', 'GNAS', 'MRPS30', 'MSMB', 'DDIT4', 'TTC36', 'S100A1', 'FAM208B', 'STT3B', 'SLC38A1', 'DMKN', 'SEC14L2', 'FMO5', 'DCAF10', 'WFDC2', 'GFRA1', 'LDLRAD4', 'TXNIP', 'SCGB3A1', 'APOD', 'N4BP2L2', 'TNC', 'ADIRF', 'NPY1R', 'NBPF1', 'TMEM176A', 'GLUL', 'BMP2K', 'SLC44A1', 'GFPT1', 'PSD3', 'CCNG2', 'CGNL1', 'TMED7', 'NOVA1', 'ARCN1', 'NEK10', 'GPC6', 'SCGB1B2P', 'IGHG4', 'SYT1', 'SYNGR2', 'HSPA1A', 'ATP6AP1', 'TSPAN13', 'MT-ND2', 'NIFK', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'EVL', 'GRN', 'ERH', 'CD81', 'NUPR1', 'SELENBP1', 'C1orf56', 'LMO3', 'PLK2', 'HACD3', 'RBBP8', 'CANX', 'ENAH', 'SCD', 'CREB3L2', 'SYNCRIP', 'TBL1XR1', 'DDR1', 'ERBB3', 'CHPT1', 'BANF1', 'UGDH', 'SCUBE2', 'UQCR10', 'COX6C', 'ATP5G1', 'PRSS23', 'MYEOV2', 'PITX1', 'MT-ND4L', 'TPM1', 'HMGCS2', 'ADIPOR2', 'UGCG', 'FAM129B', 'TNIP1', 'IFI6', 'CA2', 'ESR1', 'TMBIM4', 'NFIX', 'PDCD6IP', 'CRIM1', 'ARHGEF12', 'ENTPD5', 'PATZ1', 'ZBTB41', 'UCP1', 'ANO1', 'RP11-356O9.1', 'MYB', 'ZBTB44', 'SCPEP1', 'HIPK2', 'CDK2AP1', 'CYHR1', 'SPINK8', 'FKBP10', 'ISOC1', 'CD59', 'RAMP1', 'AFF3', 'MT-CYB', 'PPP1CB', 'PKM', 'ALDH2', 'PRSS8', 'NPW', 'SPR', 'PRDX3', 'SCOC', 'TMED10', 'KIAA0196', 'NDP', 'ZSWIM7', 'AP2A1', 'PLAT', 'SUSD3', 'CRABP2', 'DNAJC12', 'DHCR24', 'PPT1', 'FAM234B', 'DDX17', 'LRP2', 'ABCD3', 'CDH1', 'NFIA') 

gene_module[["gene_module_6"]]<-c('AGR2', 'TFF3', 'SELM', 'CD63', 'CTSD', 'MDK', 'CD74', 'S100A13', 'IFITM3', 'HLA-B', 'AZGP1', 'FXYD3', 'IFITM2', 'RABAC1', 'S100A14', 'CRABP2', 'LTF', 'RARRES1', 'HLA-A', 'PPIB', 'HLA-C', 'S100A10', 'S100A9', 'TIMP1', 'DDIT4', 'S100A16', 'LGALS1', 'LAPTM4A', 'SSR4', 'S100A6', 'CD59', 'BST2', 'PDIA3', 'KRT19', 'CD9', 'FXYD5', 'SCGB2A2', 'NUCB2', 'TMED3', 'LY6E', 'CFD', 'ITM2B', 'PDZK1IP1', 'LGALS3', 'NUPR1', 'SLPI', 'CLU', 'TMED9', 'HLA-DRA', 'SPTSSB', 'TMEM59', 'KRT8', 'CALR', 'HLA-DRB1', 'IFI6', 'NNMT', 'CALML5', 'S100P', 'TFF1', 'ATP1B1', 'SPINT2', 'PDIA6', 'S100A8', 'HSP90B1', 'LMAN1', 'RARRES3', 'SELENBP1', 'CEACAM6', 'TMEM176A', 'EPCAM', 'MAGED2', 'SNCG', 'DUSP4', 'CD24', 'PERP', 'WFDC2', 'HM13', 'TMBIM6', 'C12orf57', 'DKK1', 'MAGED1', 'PYCARD', 'RAMP1', 'C11orf31', 'STOM', 'TNFSF10', 'BSG', 'TMED10', 'ASS1', 'PDLIM1', 'CST3', 'PDIA4', 'NDUFA4', 'GSTP1', 'TYMP', 'SH3BGRL3', 'PRSS23', 'P4HA1', 'MUC5B', 'S100A1', 'PSAP', 'TAGLN2', 'MGST3', 'PRDX5', 'SMIM22', 'NPC2', 'MESP1', 'MYDGF', 'ASAH1', 'APP', 'NGFRAP1', 'TMEM176B', 'C8orf4', 'KRT81', 'VIMP', 'CXCL17', 'MUC1', 'COMMD6', 'TSPAN13', 'TFPI', 'C15orf48', 'CD151', 'TACSTD2', 'PSME2', 'CLDN7', 'ATP6AP2', 'CUTA', 'MT2A', 'CYB5A', 'CD164', 'TM4SF1', 'SCGB1D2', 'GSTM3', 'EGLN3', 'LMAN2', 'IFI27', 'PPP1R1B', 'B2M', 'ANXA2', 'SARAF', 'MUCL1', 'CSRP1', 'NPW', 'SLC3A2', 'PYDC1', 'QSOX1', 'TSPAN1', 'GPX1', 'TMSB4X', 'FGG', 'GUK1', 'IL32', 'ATP6V0E1', 'BCAP31', 'CHCHD10', 'TSPO', 'TNFRSF12A', 'MT1X', 'PDE4B', 'HSPA5', 'SCD', 'SERINC2', 'PSCA', 'VAMP8', 'ELF3', 'TSC22D3', 'S100A7', 'GLUL', 'ZG16B', 'TMEM45A', 'APMAP', 'RPS26', 'CALU', 'OSTC', 'NCCRP1', 'SQLE', 'RPS28', 'SSR2', 'SOX4', 'CLEC3A', 'TMEM9', 'RPL10', 'MUC5AC', 'HLA-DPA1', 'ZNHIT1', 'AQP5', 'CAPG', 'SPINT1', 'NDFIP1', 'FKBP2', 'C1S', 'LDHA', 'NEAT1', 'RPL36A', 'S100A11', 'LCN2', 'TUBA1A', 'GSTK1', 'SEPW1', 'P4HB') 

gene_module[["gene_module_7"]]<-c('KCNQ1OT1', 'AKAP9', 'RHOB', 'SOX4', 'VEGFA', 'CCNL1', 'RSRP1', 'RRBP1', 'ELF3', 'H1FX', 'FUS', 'NEAT1', 'N4BP2L2', 'SLC38A2', 'BRD2', 'PNISR', 'CLDN4', 'MALAT1', 'SOX9', 'DDIT3', 'TAF1D', 'FOSB', 'ZNF83', 'ARGLU1', 'DSC2', 'MACF1', 'GTF2I', 'SEPP1', 'ANKRD30A', 'PRLR', 'MAFB', 'NFIA', 'ZFAS1', 'MTRNR2L12', 'RNMT', 'NUPR1', 'MT-ND6', 'RBM39', 'HSPA1A', 'HSPA1B', 'RGS16', 'SUCO', 'XIST', 'PDIA6', 'VMP1', 'SUGP2', 'LPIN1', 'NDRG1', 'PRRC2C', 'CELF1', 'HSP90B1', 'JUND', 'ACADVL', 'PTPRF', 'LMAN1', 'HEBP2', 'ATF3', 'BTG1', 'GNAS', 'TSPYL2', 'ZFP36L2', 'RHOBTB3', 'TFAP2A', 'RAB6A', 'KMT2C', 'POLR2J3', 'CTNND1', 'PRRC2B', 'RNF43', 'CAV1', 'RSPO3', 'IMPA2', 'FAM84A', 'FOS', 'IGFBP5', 'NCOA3', 'WSB1', 'MBNL2', 'MMP24-AS1', 'DDX5', 'AP000769.1', 'MIA3', 'ID2', 'HNRNPH1', 'FKBP2', 'SEL1L', 'PSAT1', 'ASNS', 'SLC3A2', 'EIF4EBP1', 'HSPH1', 'SNHG19', 'RNF19A', 'GRHL1', 'WBP1', 'SRRM2', 'RUNX1', 'ASH1L', 'HIST1H4C', 'RBM25', 'ZNF292', 'RNF213', 'PRPF38B', 'DSP', 'EPC1', 'FNBP4', 'ETV6', 'SPAG9', 'SIAH2', 'RBM33', 'CAND1', 'CEBPB', 'CD44', 'NOC2L', 'LY6E', 'ANGPTL4', 'GABPB1-AS1', 'MTSS1', 'DDX42', 'PIK3C2G', 'IAH1', 'ATL2', 'ADAM17', 'PHIP', 'MPZ', 'CYP27A1', 'IER2', 'ACTR3B', 'PDCD4', 'COLCA1', 'KIAA1324', 'TFAP2C', 'CTSC', 'MYC', 'MT1X', 'VIMP', 'SERHL2', 'YPEL3', 'MKNK2', 'ZNF552', 'CDH1', 'LUC7L3', 'DDIT4', 'HNRNPR', 'IFRD1', 'RASSF7', 'SNHG8', 'EPB41L4A-AS1', 'ZC3H11A', 'SNHG15', 'CREB3L2', 'ERBB3', 'THUMPD3-AS1', 'RBBP6', 'GPBP1', 'NARF', 'SNRNP70', 'RP11-290D2.6', 'SAT1', 'GRB7', 'H1F0', 'EDEM3', 'KIAA0907', 'ATF4', 'DNAJC3', 'DKK1', 'SF1', 'NAMPT', 'SETD5', 'DYNC1H1', 'GOLGB1', 'C4orf48', 'CLIC3', 'TECR', 'HOOK3', 'WDR60', 'TMEM101', 'SYCP2', 'C6orf62', 'METTL12', 'HIST1H2BG', 'PCMTD1', 'PWWP2A', 'HIST1H3H', 'NCK1', 'CRACR2B', 'NPW', 'RAB3GAP1', 'TMEM63A', 'MGP', 'ANKRD17', 'CALD1', 'PRKAR1A', 'PBX1', 'ATXN2L', 'FAM120A', 'SAT2', 'TAF10', 'SFRP1', 'CITED2') 

#maybe add Ecotypes as well, no given gene lists from the publication??



single_sample_cell_signature_transfer<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  dat_sub<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial"))

  #embo lineage
  dat<-AddModuleScore(dat,features=lineage_in,names=names(lineage_in),assay="SoupXRNA",seed=123,search=TRUE)
  colnames(dat@meta.data)[startsWith(prefix="Cluster",colnames(dat@meta.data))]<-c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str") #Rename them

  #Immune cell features and PAM50 canonical short list of genes
  for(i in 1:length(features_in)){
    features_in[[i]]<-features_in[[i]][features_in[[i]] %in% row.names(dat[["SoupXRNA"]])] #make sure gene names match
    dat<-MetaFeature(dat,features=c(features_in[[i]]),meta.name=names(features_in)[i],assay="SoupXRNA")}

  #SCSubype List of genes
  #run only on epithelial cells
  module_scores<-AddModuleScore(dat_sub,features=module_feats,assay="SoupXRNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
  module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-(length(module_feats)-1),ncol(module_scores@meta.data))]
  colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like
  dat<-AddMetaData(dat,metadata=module_scores)

  #Swarbrick Gene Modules
  #run only on epithelial cells
  gene_module_out<-AddModuleScore(dat_sub,features=gene_module,assay="SoupXRNA",search=TRUE,name=names(gene_module)) #use add module function to add cell scores
  gene_module_out<-gene_module_out@meta.data[seq(ncol(gene_module_out@meta.data)-(length(gene_module)-1),ncol(gene_module_out@meta.data))]#get the 7 added gene modules
  colnames(gene_module_out)<-names(gene_module) 
  dat<-AddMetaData(dat,metadata=gene_module_out)
  saveRDS(dat,file=file_in)
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_cell_signature_transfer)


single_sample_EMBO_assignment<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  met<-dat@meta.data
  met<-met[met$predicted.id %in% c("Cancer Epithelial","Normal Epithelial"),]
  embo_list<-  c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str" )
  max_embo<-lapply(1:nrow(met),function(i) embo_list[which(met[i,embo_list]==max(met[i,embo_list],na.rm=T))])
  max_embo<-unlist(lapply(1:length(max_embo),function(i) do.call("paste",as.list(max_embo[[i]]))))
  max_embo<-unlist(lapply(max_embo,function(i) gsub("EMBO_","",i)))
  names(max_embo)<-row.names(met)
  dat<-AddMetaData(dat,max_embo,col.name="EMBO_designation")
  print(paste("Saving",outname))
  saveRDS(dat,file=file_in)
}

single_sample_PAM50_assignment<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  met<-dat@meta.data
  met<-met[met$predicted.id %in% c("Cancer Epithelial","Normal Epithelial"),]
  pam_list<-  c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")
  max_pam<-lapply(1:nrow(met),function(i) pam_list[which(met[i,pam_list]==max(met[i,pam_list],na.rm=T))])
  max_pam<-unlist(lapply(1:length(max_pam),function(i) do.call("paste",as.list(max_pam[[i]]))))
  names(max_pam)<-row.names(met)
  dat<-AddMetaData(dat,max_pam,col.name="SCSubtype_designation")
  print(paste("Saving",outname))
  saveRDS(dat,file=file_in)
}



lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_PAM50_assignment)
lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_EMBO_assignment)


```
Plot Features on Cells Per Sample

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

single_sample_epcam<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".EPCAM.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  DefaultAssay(dat)<-"SoupXRNA"
  plt<-FeaturePlot(dat,features="EPCAM",reduction="multimodal_umap",order=T)
  ggsave(plt,file=out_plot,width=10,height=10)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  plt<-VlnPlot(dat,features="EPCAM",group.by="predicted.id")
  ggsave(plt,file=paste0(out_plot,"VlnPlt.pdf"),width=10,height=10)
  system(paste0("slack -F ",paste0(out_plot,"VlnPlt.pdf")," ryan_todo"))
}

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_epcam)

```

## Determine Tumor Cells via CNV Callers

### InferCNV on RNA Profiles

. Immune and endothelial cells were used to 
define the reference cell-inferred copy number profiles.
(From https://www.nature.com/articles/s41588-021-00911-1)

```R
####Run InferCNV
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

####RUNNING INFERCNV#####
infercnv_per_sample<-function(x){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
    if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }


  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  #excluding T-cells due do low level but nonexclusive cell label prediction
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL")) 

  #write out gene order list
  gene_order<-annotation[!duplicated(annotation$gene_name),]
  gene_order<-as.data.frame(gene_order[gene_order$gene_name %in% row.names(dat),])
  gene_order<-gene_order[c("gene_name","seqnames","start","end")]
  chrorder<-paste0("chr",c(1:22,"X","Y","M"))
  gene_order$seqnames<-factor(gene_order$seqnames,levels=chrorder) # set chr order
  gene_order<-with(gene_order, gene_order[order(seqnames, start),]) #order by chr and start position
  write.table(gene_order,file="inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)
  gene_order<-read.table("inferCNV.gene_order.txt")

  counts=as.matrix(dat@assays$RNA@counts[,colnames(dat)])
  write.table(counts,file=paste0(wd,"/",outname,"_inferCNV.counts.txt"),sep="\t",col.names=T,row.names=T,quote=F)
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["cnv_ref"]))
  write.table(cell_annotation,file=paste0(wd,"/",outname,"_inferCNV.annotation.txt"),sep="\t",col.names=F,row.names=F,quote=F)

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(wd,"/",outname,"_inferCNV.counts.txt"),
                                      annotations_file=paste0(wd,"/",outname,"_inferCNV.annotation.txt"),
                                      delim="\t",
                                      gene_order_file="inferCNV.gene_order.txt",
                                      ref_group_names="TRUE")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0(wd,"/",outname,"_inferCNV"), 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               resume_mode=F,
                               num_threads=20)
  saveRDS(infercnv_obj,paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
  system(paste0("slack -F ",wd,"/",outname,"_inferCNV","/","infercnv.png"," -T ","\"",outname,"\"" ," ryan_todo") )
}

lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",4,10,12),infercnv_per_sample)
#
#4 10 12 listed last because I might want to rereun them with less threads
```

## Run CaSpER on RNA 

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

###Trying CASPER for CNV Profiling

library(CaSpER) 

casper_per_sample<-function(x){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  system(paste0("mkdir ",dir_in,"/casper"))
  bam_location<-paste0(dir_in,"/gex_possorted_bam.bam")
  BAFExtract_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract"
  hg38_list_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list" #downloaded from https://github.com/akdess/BAFExtract
  hg38_folder_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/"
  baf_sample_directory<-paste0(dir_in,"/casper")

  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL")) 

  control<-names(dat$cnv_ref == "TRUE") #pulling this from the inferCNV function
  log.ge <- as.matrix(dat@assays$RNA@data)
  genes <- rownames(log.ge)
  annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
  log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
  rownames(log.ge) <- annotation$Gene
  log.ge <- log2(log.ge +1)

  system(paste0("samtools view ",bam_location," | ",BAFExtract_location," -generate_compressed_pileup_per_SAM stdin ",hg38_list_location," ",baf_sample_directory," 30 0 && wait;")) #generate BAF calls
  #example of actual call: samtools view /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs/gex_possorted_bam.bam| /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper 30 0 &
  system(paste0(BAFExtract_location," -get_SNVs_per_pileup ",hg38_list_location," ",baf_sample_directory," ",hg38_folder_location," 1 1 0.1 ",baf_sample_directory,"/test.snp")) #generage snv files from BAF
  #example of actual call: /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -get_SNVs_per_pileup /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/ 1 1 0.1 /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper/test.snp
  loh <- readBAFExtractOutput ( path=baf_sample_directory, sequencing.type="bulk") 
  names(loh) <- gsub(".snp", "", names(loh))
  load(paste0(hg38_folder_location,"/maf.rda")) ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
  loh<- list()
  loh[[1]] <- maf
  names(loh) <- sample_name
  loh.name.mapping <- data.frame (loh.name= sample_name , sample.name=colnames(log.ge))

  #analysis demonstration: https://rpubs.com/akdes/673120
  object <- CreateCasperObject(raw.data=log.ge,
    loh.name.mapping=loh.name.mapping, 
    sequencing.type="single-cell", 
    cnv.scale=3, 
    loh.scale=3, 
    expr.cutoff=0.1, 
    filter="median", 
    matrix.type="normalized",
    annotation=annotation, 
    method="iterative", 
    loh=loh, 
    control.sample.ids=control, 
    cytoband=cytoband)

  saveRDS(object,paste0(dir_in,"/casper/",sample_name,".initialobj.rds"))

  pdf(paste0(dir_in,"/casper/",sample_name,".Distrubution.pdf"))
  plot(density(as.vector(object@control.normalized[[3]])))
  plot(density(log2(object@control.normalized.noiseRemoved[[3]]+1)))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"/casper/",sample_name,".Distrubution.pdf"), " ryan_todo"))

  object<-readRDS(paste0(dir_in,"/casper/",sample_name,".initialobj.rds"))
  ## runCaSpER
  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

  ## summarize large scale events 
  finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
  final.obj <- final.objects[[9]]
  saveRDS(final.obj,paste0(dir_in,"/casper/",sample_name,".finalobj.rds"))
  saveRDS(finalChrMat,paste0(dir_in,"/casper/",sample_name,".finalchrmat.rds"))
  final.obj<-readRDS(paste0(dir_in,"/casper/",sample_name,".finalobj.rds"))

  #remove all NA columns
  #final.obj <- final.obj[!apply(is.na(final.obj[-1,]), 1, all),]

  #breakdown of function CaSpER::plotHeatmap10x
    assignInNamespace(x = "draw_matrix", value = draw_matrix2,
        ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
        ns = asNamespace("pheatmap"))
    data <- final.obj@control.normalized.noiseRemoved[[3]]
    x.center <- mean(data)
    quantiles = quantile(data[data != x.center], c(0.01, 0.99))
    delta = max(abs(c(x.center - quantiles[1], quantiles[2] -x.center)))
    low_threshold = x.center - delta
    high_threshold = x.center + delta
    x.range = c(low_threshold, high_threshold)
    data[data < low_threshold] <- low_threshold
    data[data > high_threshold] <- high_threshold
    breaks <- seq(x.range[1], x.range[2], length = 16)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:22)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
    xpos <- c(half[1], (idx[-22] + half))
    xlabel[xpos] <- 1:22

    pheatmap(t(data), cluster_cols = F, cluster_rows = T,
        gaps_col = idx, color = color, breaks = breaks, labels_col = xlabel,
        show_rownames = F, filename = paste0(dir_in,"/casper/",sample_name,".heatmap.pdf"))
  system(paste0("slack -F ",paste0(dir_in,"/casper/",sample_name,".heatmap.pdf"), " ryan_todo"))


  pdf(paste0(dir_in,"/casper/",sample_name,".LargeCNV.heatmap.pdf"))
  Heatmap(finalChrMat,
    cluster_columns=F,
    show_row_names=F,
    column_names_rot = 45,
    column_split=sort(rep(c(1:22),2)),
    column_names_gp = gpar(fontsize = 5))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"/casper/",sample_name,".LargeCNV.heatmap.pdf")," ryan_todo"))


  #### VISUALIZATION 
  chrMat <- finalChrMat
  plot.data <- melt(chrMat)
  plot.data$value2 <- "neutral"
  plot.data$value2[plot.data$value > 0] <- "amplification"
  plot.data$value2[plot.data$value < 0] <- "deletion"
  plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", 
      "deletion", "neutral"))
  plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
  p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + 
      geom_tile(colour = "white", size = 0.01) + 
      labs(x = "", y = "") +
       scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + 
      theme(legend.position = "right", legend.direction = "vertical", 
      legend.title = element_blank(), strip.text.x = element_blank(), 
      legend.text = element_text(colour = "black", size = 7, 
        face = "bold"), legend.key.height = grid::unit(0.8, 
        "cm"), legend.key.width = grid::unit(0.5, "cm"), 
      axis.text.x = element_text(size = 5, colour = "black", 
        angle = -45, hjust = 0), axis.text.y = element_text(size = 6, 
        vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
      plot.title = element_text(colour = "black", hjust = 0, 
        size = 6, face = "bold"))
  ggsave(p,file=paste0(dir_in,"/casper/",sample_name,".final_plot.pdf"))
  system(paste0("slack -F ",paste0(dir_in,"/casper/",sample_name,".final_plot.pdf"), " ryan_todo"))

}

lapply(c(16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12),function(x) casper_per_sample(x))

#1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10

```

## Run CopyKat
https://github.com/navinlabcode/copykat
```R
#library(devtools)
#install_github("navinlabcode/copykat")

library(Signac)
library(Seurat)
library(copykat)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")


copykat_per_sample<-function(x){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  system(paste0("mkdir ",dir_in,"/copykat"))
  exp.rawdata <- as.matrix(dat@assays$RNA@counts)


  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  cnv_ref<-row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",])
  copykat_out <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=sample_name, distance="euclidean", norm.cell.names=cnv_ref,output.seg="FALSE", plot.genes="TRUE", genome="hg20",n.cores=5)
  saveRDS(copykat_out,paste0(dir_in,"/copykat/",sample_name,".copykat.RDS"))
}

lapply(c(10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x) copykat_per_sample(x))

#done 1,3,4,5,6,7,8,9,
```

### Plotting of InferCNV and CaSpER Output
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(reshape2)
library(philentropy)
library(CaSpER)
library(dendextend)
library(ggalluvial)


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  pam50_colors<-c("Basal"="red","Her2"="pink","LumA"="blue","LumB"="cyan","Normal"="grey","NA"="black")
  embo_colors<-c("Basal"="green","LP"="blue","ML"="orange","Str"="red","NA"="black")
  ########################################


####RUNNING INFERCNV PLOTTING#####
infercnv_per_sample_plot<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
    #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  dat_file_path=file_in
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #this is same as initial run of inferCNV, just didn't save seurat object
  infercnv_obj<-readRDS(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
  cnv<-t(infercnv_obj@expr.data)
  cnv_ref<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",]),]
  cnv<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]
  col_fun = colorRamp2(c(min(unlist(cnv)), median(unlist(cnv)), max(unlist(cnv))), c("blue", "white", "red"))

  dist_method="manhattan"
  dist_x<-philentropy::distance(cnv,method=dist_method,as.dist.obj=T,use.row.names=T)
  dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
  k_search<-find_k(dend,krange=2:10) #search for optimal K from 2-10
  k_clus_number<-k_search$nc
  k_clus_id<-k_search$pamobject$clustering
  dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters
  saveRDS(dend,file=paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.dend.Rds")) #save dendrogram

  #set up heatmap annotation
  met<-as.data.frame(dat@meta.data)
  met_ref<-met[row.names(met) %in% row.names(cnv_ref),]
  met<-met[row.names(met) %in% row.names(cnv),]
  if(any(!(unique(met$PAM50_designation) %in% names(pam50_colors)))){
    met[met$PAM50_designation %in% unique(met$PAM50_designation)[!(unique(met$PAM50_designation) %in% names(pam50_colors))],]$PAM50_designation<-"NA"}
  if(any(!(unique(met$EMBO_designation) %in% names(embo_colors)))){
    met[met$EMBO_designation %in% unique(met$EMBO_designation)[!(unique(met$EMBO_designation) %in% names(embo_colors))],]$EMBO_designation<-"NA"}

  read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
    max(met$gex_exonic_umis+met$gex_intronic_umis)), 
    c("white","black"))

  ha = HeatmapAnnotation(which="row",
    cell_type=met$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met$gex_exonic_umis+met$gex_intronic_umis,
    pam_50=met$PAM50_designation,
    embo=met$EMBO_designation,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col,
            embo=embo_colors,
            pam_50=pam50_colors))
  plt1<-Heatmap(cnv,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      cluster_rows=dend,
      left_annotation=ha,
      column_split=infercnv_obj@gene_order$chr)
  ha_ref = HeatmapAnnotation(which="row",
    cell_type=met_ref$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met_ref$gex_exonic_umis+met_ref$gex_intronic_umis,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col))
  plt1_ref<-Heatmap(cnv_ref,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      left_annotation=ha_ref,
      column_split=infercnv_obj@gene_order$chr)
    pdf(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.heatmap.pdf"),width=20)
    print(plt1_ref)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.heatmap.pdf")," ryan_todo"))
}

####RUNNING CASPER PLOT#####
casper_per_sample_plot<-function(x){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  DefaultAssay(dat)<-"RNA"
  print(paste("Reading in:",paste0(dir_in,"_casper/",sample_name,".finalobj.rds")))
  final.obj<-readRDS(paste0(dir_in,"_casper/",sample_name,".finalobj.rds"))
  finalChrMat<-readRDS(paste0(dir_in,"/casper/",sample_name,".finalchrmat.rds"))
  cnv_ref<-finalChrMat[row.names(finalChrMat) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",]),]
  cnv<-finalChrMat[row.names(finalChrMat) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]
  col_fun = setNames(c("blue","white","red"),seq(-1, 1, by = 1))
  print("Performing distance calculation.")
  dist_method="euclidean"
  dist_x<-philentropy::distance(cnv,method=dist_method,as.dist.obj=T,use.row.names=T)
  dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
  k_search<-find_k(dend,krange=2:10) #search for optimal K from 2-10
  k_clus_number<-k_search$nc
  print(paste("Determined ",k_clus_number," of clusters."))
  k_clus_id<-k_search$pamobject$clustering
  dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters
  saveRDS(dend,file=paste0(dir_in,"_casper/",sample_name,".casper.dend.Rds")) #save dendrogram

  print("Generating heatmap of CNVs.")
  
  #set up heatmap annotation
  met<-as.data.frame(dat@meta.data)
  met_ref<-met[row.names(met) %in% row.names(cnv_ref),]
  met<-met[row.names(met) %in% row.names(cnv),]
  if(any(!(unique(met$PAM50_designation) %in% names(pam50_colors)))){
  met[met$PAM50_designation %in% unique(met$PAM50_designation)[!(unique(met$PAM50_designation) %in% names(pam50_colors))],]$PAM50_designation<-"NA"}
  if(any(!(unique(met$EMBO_designation) %in% names(embo_colors)))){
  met[met$EMBO_designation %in% unique(met$EMBO_designation)[!(unique(met$EMBO_designation) %in% names(embo_colors))],]$EMBO_designation<-"NA"}

  read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
    max(met$gex_exonic_umis+met$gex_intronic_umis)), 
    c("white","black"))

  ha = HeatmapAnnotation(which="row",
    cell_type=met$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met$gex_exonic_umis+met$gex_intronic_umis,
    pam_50=met$PAM50_designation,
    embo=met$EMBO_designation,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col,
            embo=embo_colors,
            pam_50=pam50_colors))

  plt1<-Heatmap(cnv,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      cluster_rows=dend,
      left_annotation=ha,
      column_split=factor(substr(colnames(cnv),1,nchar(colnames(cnv))-1),levels=unique(substr(colnames(cnv),1,nchar(colnames(cnv))-1)))
      )
  ha_ref = HeatmapAnnotation(which="row",
    cell_type=met_ref$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met_ref$gex_exonic_umis+met_ref$gex_intronic_umis,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col))
  plt1_ref<-Heatmap(cnv_ref,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv_ref),
      col=col_fun,
      left_annotation=ha_ref,
      column_split=factor(substr(colnames(cnv),1,nchar(colnames(cnv))-1),levels=unique(substr(colnames(cnv),1,nchar(colnames(cnv))-1)))
      )
    pdf(paste0(dir_in,"_casper/",sample_name,".casper.heatmap.pdf"),width=20)
    print(plt1_ref)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".casper.heatmap.pdf")," ryan_todo"))
}



lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) infercnv_per_sample_plot(x))
lapply(c(16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) casper_per_sample_plot(x))
#15
#Done: 1,3,5,6,7,8,9,
lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) compare_RNA_cnv_results(x))


```

## Run IntClust on Samples
Using iC10 CRAN Package https://cran.r-project.org/web/packages/iC10/iC10.pdf

```R
#install.packages("iC10")
library(iC10)
library(Seurat)
library(Signac)
#CN = ID (Sample Name) \t chromosome_name (Chr) \t loc.start (start) loc.end (end) seg.mean (log2ratio of segment)
#OR
#CN = Row (hgnc gene names) X Column (Sample)
#Exp =  Row (hgnc gene names) X Column (Sample)

#using InferCNV(gene level CNV calling) as CN matrix, and RNA data as Exp Matrix

iC10_per_sample<-function(x){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  Idents(dat)<-dat$predicted.id
  dat_ep<-subset(dat, cells = row.names(dat@meta.data[dat@meta.data$predicted.id%in% c("Normal Epithelial", "Cancer Epithelial"),]))
  infercnv_obj<-readRDS(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.Rds"))
  cnv<-log2(infercnv_obj@expr.data)
  cnv<-cnv[,colnames(cnv) %in% colnames(dat_ep)]
  exp<-dat_ep[["RNA"]]@counts

  out<-matchFeatures(CN=cnv,Exp=exp,
    CN.by.feat="gene",
    Exp.by.feat="gene",
    ref=NULL)
  out<-normalizeFeatures(out, method="scale")
  out<-iC10(out)
  saveRDS(out,paste0(wd,"/",outname,"_iC10.Rds"))
  dat<-AddMetaData(dat,out$class,col.name="ic10_class")
  table(dat$ic10_class)
  saveRDS(dat,file=file_in) #overwrite old file
  print(paste("Finished Sample:",sample_name))
}

lapply(c(1,3,5,6,7,8,9,16,19,20,"RM_1","RM_2","RM_3",11,4,10,12), function(x) iC10_per_sample(x))
#done 
#15,"RM_4" not done
```

## Epithelial Subtyping Per Sample
```R
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)

  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
  alpha_val=0.33

epithelial_class_persample<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  dat<-readRDS(file=file_in)
  atac_sub<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial"))
  if(!("ic10_class" %in% colnames(atac_sub@meta.data))){
    atac_sub$ic10_class<-"NA"
  }
  plt_cell_count<-atac_sub@meta.data[,c("sample","predicted.id","diagnosis","molecular_type","PAM50_designation","EMBO_designation","ic10_class")]
  print(outname)
  return(plt_cell_count)
}


#grab all epithelial classifications
cell_count<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),function(x)epithelial_class_persample(x))

saveRDS(cell_count,"sample_epithelial_designations.rds") #save nested list of cell type assignment

#plot output of celltype count per sample
out<-readRDS("sample_epithelial_designations.rds")
out<-do.call("rbind",out)
colnames(out)<-c("sample","predicted.id","diagnosis","molecular_type","SCSubtype_designation","EMBO_designation","ic10_class") #rename just for simplicity
#clean up for samples with equal values
out[!(out$SCSubtype_designation %in% c("Basal","Her2","LumA","LumB","Normal")),]$SCSubtype_designation<-NA
  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################
plt1<-ggplot(out,aes(x=sample,fill=SCSubtype_designation))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt2<-ggplot(out,aes(x=sample,fill=EMBO_designation))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt3<-ggplot(out,aes(x=sample,fill=ic10_class))+geom_bar(position="fill")+theme_minimal()+facet_wrap(.~diagnosis+molecular_type,scale="free_x")
plt<-plt1/plt2/plt3
ggsave(plt,file="sample_epithelial_type_assignment.pdf")
system("slack -F sample_epithelial_type_assignment.pdf ryan_todo")

library(dplyr)
write.table(out,file="sample_epithelial_type_assignment.tsv",col.names=T,row.names=F,sep="\t",quote=F)
system("slack -F sample_epithelial_type_assignment.tsv ryan_todo") #note this was calculated per sample as well as in the merged data set,  the assumption is that they will be the same


```

### ATAC CNV Calling with copyscAT
Using scATAC calling algorithm copyscAT from git repo https://github.com/spcdot/CopyscAT/
Installation...
```R
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install_github("spcdot/copyscat")
library(CopyscAT)
```
Modifies the copyscAT python script (https://github.com/spcdot/CopyscAT/blob/master/process_fragment_file.py) to filter based on a metadata table rather than read count (since I already QC cells) then posted to a subdirectory

```bash
mkdir /home/groups/CEDAR/mulqueen/ref/copyscat
```

Now Running samples

Code from https://github.com/spcdot/CopyscAT/blob/master/copyscat_tutorial.R

```R
library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Generate tile references
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = "/home/groups/CEDAR/mulqueen/ref/copyscat")

##### REGULAR WORKFLOW #####
#initialize the environment
initialiseEnvironment(genomeFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
                      cytobandFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=500,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#Set up copyscAT Loop per sample
copyscAT_per_sample<-function(x){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  system(paste0("mkdir ",dir_in,"/copyscat"))
  #do python script preprocessing (basically just count fragments per window per cell)
  system(paste0("python /home/groups/CEDAR/mulqueen/ref/copyscat/process_fragment_file.py ",
  "-i ",dir_in,"/atac_fragments.tsv.gz",
  " -o ",dir_in,"/copyscat/copyscat.1mb.tsv",
  " -b ","1000000",
  " -f ","500",
  " -g ","/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
  " -c ",dir_in,"/metadata.tsv")) #modification takes in metadata table to filter cells by name, ignores -f flag
  setOutputFile(paste0(dir_in,"/copyscat"),"copyscat_out")
  #PART 1: INITIAL DATA NORMALIZATION
  scData<-readInputTable(paste0(dir_in,"/copyscat/copyscat.1mb.tsv"))
  scData_k_norm <- normalizeMatrixN(scData,
    imputeZeros = FALSE,
    dividingFactor=1,
    blacklistProp = 0.8,
    blacklistCutoff=50,
    upperFilterQuantile = 1)
  #collapse into chromosome arm level
  summaryFunction<-cutAverage
  scData_collapse<-collapseChrom3N(scData_k_norm,
    summaryFunction=summaryFunction,
    binExpand = 1,
    minimumChromValue = 100,
    logTrans = FALSE,
    tssEnrich = 1,
    logBase=2,
    minCPG=300,
    powVal=0.73) 
  #show unscaled chromosome list
  graphCNVDistribution(scData_collapse,outputSuffix = "test_violinsn2")
  #PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
  #ALTERNATE METHOD FOR CNV CALLING (with normal cells as background)
  #Using same normal cell selection as used for CASPER and InferCNV
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  control<-names(dat$cnv_ref == "TRUE") #pulling this from the inferCNV function
  #compute central tendencies based on normal cells only
  control <- control[control %in% colnames(scData_collapse)] #filter control list to control cells that survived filter
  median_iqr <- computeCenters(scData_collapse %>% select(chrom,control),summaryFunction=summaryFunction)
  #setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes in as normalCells can improve accuracy of CNV calls
  candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,
    useDummyCells = FALSE,
    propDummy=0.25,
    minMix=0.01,
    deltaMean = 0.03,
    deltaBIC2 = 0.25,
    bicMinimum = 0.1,
    subsetSize=50,
    fakeCellSD = 0.09,
    uncertaintyCutoff = 0.65,
    summaryFunction=summaryFunction,
    maxClust = 4,
    mergeCutoff = 3,
    IQRCutoff = 0.25,
    medianQuantileCutoff = -1,
    normalCells=control) 
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)
  #to save this data you can use annotateCNV4 as per usual, using normal barcodes
  final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, expectedNormals=control, saveOutput=TRUE,
    outputSuffix = "clean_cnv_b2",sdCNV = 0.6,filterResults=FALSE,filterRange=0.4,minAlteredCellProp = 0.5)
  saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs.rds"))
  print(paste("Finished sample",sample_name))
}

lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),copyscAT_per_sample)
#Done

```


# Comparison of Clones determined by CNV Callers

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(reshape2)
library(philentropy)
library(CaSpER)
library(dendextend)
library(ggalluvial)


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  pam50_colors<-c("Basal"="red","Her2"="pink","LumA"="blue","LumB"="cyan","Normal"="grey","NA"="black")
  embo_colors<-c("Basal"="green","LP"="blue","ML"="orange","Str"="red","NA"="black")
  ########################################

plotting_copyscat_persample<-function(x){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  obj_name=basename(file_in)
  dir_in=dirname(file_in)

  copyscat_out<-readRDS(paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs.rds"))
  cnv_dat<-(copyscat_out[[3]])
  row.names(cnv_dat)<-cnv_dat[,1]
  cnv_dat<-cnv_dat[,2:ncol(cnv_dat)]
  cnv_dat<-cnv_dat[,colnames(cnv_dat) %in% c('chr1p', 'chr1q', 'chr2p', 'chr2q', 'chr3p', 'chr3q', 'chr4p', 'chr5q', 'chr6p', 'chr6q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr9p', 'chr9q', 'chr10p', 'chr10q', 'chr11p', 'chr11q', 'chr12p', 'chr12q', 'chr13p', 'chr13q', 'chr14p', 'chr14q', 'chr15p', 'chr15q', 'chr16p', 'chr16q', 'chr17p', 'chr17q', 'chr18p', 'chr18q', 'chr19p', 'chr19q', 'chr20p', 'chr20q', 'chr21p', 'chr21q', 'chr22p', 'chr22q')] 

  cnv_dat<-as.data.frame(cnv_dat)
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  DefaultAssay(dat)<-"RNA"

  cnv_ref<-cnv_dat[row.names(cnv_dat) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",]),]
  cnv<-cnv_dat[row.names(cnv_dat) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]

  col_fun = colorRamp2(c(0,2,4),c("blue","white","red"))
  print("Performing distance calculation.")

  dist_method="manhattan"
  dist_x<-philentropy::distance(cnv,method=dist_method,as.dist.obj=T,use.row.names=T)
  dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
  k_search<-find_k(dend,krange=2:10) #search for optimal K from 2-10
  k_clus_number<-k_search$nc
  print(paste("Determined ",k_clus_number," of clusters."))
  k_clus_id<-k_search$pamobject$clustering
  dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters
  saveRDS(dend,paste0(dir_in,"/copyscat/",sample_name,"copyscat_dend.rds"))

  print("Generating heatmap of CNVs.")
  
  #set up heatmap annotation
  met<-as.data.frame(dat@meta.data)
  met_ref<-met[row.names(met) %in% row.names(cnv_ref),]
  met<-met[row.names(met) %in% row.names(cnv),]
  if(any(!(unique(met$PAM50_designation) %in% names(pam50_colors)))){
  met[met$PAM50_designation %in% unique(met$PAM50_designation)[!(unique(met$PAM50_designation) %in% names(pam50_colors))],]$PAM50_designation<-"NA"}
  if(any(!(unique(met$EMBO_designation) %in% names(embo_colors)))){
  met[met$EMBO_designation %in% unique(met$EMBO_designation)[!(unique(met$EMBO_designation) %in% names(embo_colors))],]$EMBO_designation<-"NA"}

  read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
    max(met$gex_exonic_umis+met$gex_intronic_umis)), 
    c("white","black"))

  ha = HeatmapAnnotation(which="row",
    cell_type=met$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met$gex_exonic_umis+met$gex_intronic_umis,
    pam_50=met$PAM50_designation,
    embo=met$EMBO_designation,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col,
            embo=embo_colors,
            pam_50=pam50_colors))

  plt1<-Heatmap(cnv,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      cluster_rows=dend,
      left_annotation=ha#,
      #column_split=factor(substr(colnames(cnv),1,nchar(colnames(cnv))-1),levels=unique(substr(colnames(cnv),1,nchar(colnames(cnv))-1)))
      )
  ha_ref = HeatmapAnnotation(which="row",
    cell_type=met_ref$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met_ref$gex_exonic_umis+met_ref$gex_intronic_umis,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col))
  plt1_ref<-Heatmap(cnv_ref,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv_ref),
      col=col_fun,
      left_annotation=ha_ref#,
      #column_split=factor(substr(colnames(cnv),1,nchar(colnames(cnv))-1),levels=unique(substr(colnames(cnv),1,nchar(colnames(cnv))-1)))
      )
    pdf(paste0(dir_in,"/copyscat",sample_name,".copyscat.heatmap.pdf"),width=20)
    print(plt1_ref)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0(dir_in,"/copyscat",sample_name,".copyscat.heatmap.pdf")," ryan_todo"))
}

lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),plotting_copyscat_persample)
#1,3,5,8

```
### Files for Travis

Transfering to /home/groups/CEDAR/scATACcnv/Hisham_data/final_data:
- Metadata
- Counts matrix (atac)
- peaks bed file
- inferCNV folder
- Casper folder
- atac_possorted_bam.bam
- Seurat Object

```R
library(Signac)
library(Seurat)

transfer_data<-function(x){
if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
  }
  metadata<-as.data.frame(dat@meta.data)
  counts_matrix<-as.data.frame(dat[["peaks"]]@counts)
  peaks_bed_file<-as.data.frame(do.call("rbind",strsplit(row.names(dat[["peaks"]]),"-")))
  write.table(metadata,file=paste0(wd,"/","metadata.tsv"),sep="\t",col.names=T,row.names=T)
  write.table(counts_matrix,file=paste0(wd,"/","counts_matrix.tsv"),sep="\t",col.names=T,row.names=T)
  write.table(peaks_bed_file,file=paste0(wd,"/","peaks.bed"),sep="\t",col.names=F,row.names=F)
  print(paste("Finished sample:",sample_name))
}

lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12),function(x) transfer_data(x))
```

- Metadata
- Counts matrix (atac)
- peaks bed file
- inferCNV folder
- Casper folder
- atac_possorted_bam.bam
- Seurat Object

```bash
out_dir="/home/groups/CEDAR/scATACcnv/Hisham_data/final_data"
mkdir $out_dir


for i in 1 3 5 6 7 8 9 10 11 12; do
  sample="sample_"${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_"${i}"/outs"
  mkdir ${out_dir}"/"${sample};
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/peaks.bed" ${out_dir}"/sample_"${i}
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/sample_"${i}
  cp -r ${in_dir}"/casper" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/"${sample}".QC.SeuratObject.rds" ${out_dir}"/sample_"${i} & done &


for i in 15 16 19 20; do
  sample="sample_"${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_"${i}"/outs"
  mkdir ${out_dir}"/"${sample};
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/peaks.bed" ${out_dir}"/sample_"${i}
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/sample_"${i}
  cp -r ${in_dir}"/casper" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/sample_"${i}
  cp ${in_dir}"/"${sample}".QC.SeuratObject.rds" ${out_dir}"/sample_"${i} & done &

for i in "RM_1" "RM_2" "RM_3" "RM_4"; do
  sample=${i}
  in_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/"${i}"/outs"
  mkdir ${out_dir}"/"${sample};
  cp ${in_dir}"/metadata.tsv" ${out_dir}"/"${sample}
  cp ${in_dir}"/counts_matrix.tsv" ${out_dir}"/"${sample}
  cp ${in_dir}"/peaks.bed" ${out_dir}"/"${sample}
  cp -r ${in_dir}"/"${sample}"_inferCNV" ${out_dir}"/"${sample}
  cp -r ${in_dir}"/casper" ${out_dir}"/"${sample}
  cp ${in_dir}"/atac_possorted_bam.bam" ${out_dir}"/"${sample}
  cp ${in_dir}"/"${sample}".QC.SeuratObject.rds" ${out_dir}"/"${sample} & done &

```

Call peaks per cell type
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")


dat<-readRDS("phase2.QC.SeuratObject.rds")

#call peaks per predicted.id cell type
peaks <- CallPeaks(dat, 
  assay="ATAC",
  group.by="predicted.id",
  combine.peaks=FALSE,
  macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")
  #use this set of peaks for all samples

for(i in 1:length(peaks)){
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peak_name<-unique(dat$predicted.id)[i]
  peaks_out <- keepStandardChromosomes(peaks[[i]], pruning.mode = "coarse")
  peaks_out <- subsetByOverlaps(x = peaks_out, ranges = blacklist_hg38_unified, invert = TRUE)
  print(paste0("Generated peakset for ",peak_name))
  write.table(as.data.frame(peaks_out)[1:3],file=paste0(peak_name,".bed"),sep="\t",quote=F,col.names=F,row.names=F)
}

```

```bash
cp *bed /home/groups/CEDAR/scATACcnv/Hisham_data/final_data
```
<!--

## Output metadata per sample
```R

library(Seurat)
library(Signac)
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

dat<-sample_in[2]

####RUNNING INFERCNV#####
outmetadata<-function(dat){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  print(paste("Reading in:",dat))
  file_in=basename(dat)
  sample_name=substr(file_in,1,nchar(file_in)-17)
  dir_in=dirname(dat)
  dat<-readRDS(dat)
  out_metadata<-as.data.frame(dat@meta.data)
  out_metadata<-out_metadata[c("sample","predicted.id","diagnosis","molecular_type")]
  out_metadata$cellID<-row.names(out_metadata)
  write.table(out_metadata,file=paste0(dir_in,"/",sample_name,".sample_metadata.tsv"),col.names=T,sep="\t",row.names=F,quote=F)
  peaks_bed<-as.data.frame(dat@assays$peaks@ranges)[1:3]
  write.table(peaks_bed,file=paste0(dir_in,"/",sample_name,".sample_peaks.bed"),col.names=F,sep="\t",row.names=F,quote=F)
  counts_matrix<-as.data.frame(dat@assays$ATAC@counts)
  write.table(counts_matrix,file=paste0(dir_in,"/",sample_name,".counts_matrix.tsv"),col.names=T,sep="\t",row.names=T,quote=F)
}

lapply(sample_in,function(x) outmetadata(dat=x))


#/home/groups/CEDAR/scATACcnv/Hisham_data
for i in 1 3 4 5 6 7 8 9 10 11 12;
do cp /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_${i}/outs/sample_${i}.sample_peaks.bed /home/groups/CEDAR/scATACcnv/Hisham_data;
cp /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_${i}/outs/sample_${i}.sample_metadata.tsv /home/groups/CEDAR/scATACcnv/Hisham_data;
cp /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_${i}/outs/sample_${i}.counts_matrix.tsv /home/groups/CEDAR/scATACcnv/Hisham_data;
done

```
-->

# Recreate new merged Seurat object with filtered cells


Final Merged Seurat Object from all phases after QC

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".QC.SeuratObject.rds"))
  dat$sample<-outname #set up sample metadata
  return(dat)}

out<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),merge_seurat)


dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"), project = "all_data")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.QC.SeuratObject.rds")
```

### Perform Merged Object Clustering

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/")

dat<-readRDS("phase2.QC.SeuratObject.rds")


#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
#set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
  diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
  molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
  ########################################

my_cols<-colorRampPalette(brewer.pal(8, "Spectral"))(length(unique(dat$seurat_clusters)))
alpha_val=0.33
p1<-DimPlot(dat,reduction="rna_umap",group.by="predicted.id")+ggtitle("RNA UMAP")+theme_minimal()+scale_fill_manual(values=type_cols)
p2<-DimPlot(dat,reduction="atac_umap",group.by="predicted.id")+ggtitle("ATAC UMAP")+theme_minimal()+scale_fill_manual(values=type_cols)
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted.id")+ggtitle("Multimodal UMAP")+theme_minimal()+scale_fill_manual(values=type_cols)



#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file="phase2_multimodal.umap.pdf",width=10,height=10)
system("slack -F phase2_multimodal.umap.pdf ryan_todo")
saveRDS(dat,file="phase2.QC.SeuratObject.rds")

```

### Cistopic on merged samples

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

atac_sub<-readRDS("phase2.QC.SeuratObject.rds")
wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
outname<-paste0("phase2")
cistopic_counts_frmt<-atac_sub$peaks@counts
row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
print("made cistopic object")
sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))

sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =atac_sub@meta.data)
pdf(paste0(wd,"/",outname,"_model_selection.pdf"))
par(mfrow=c(3,3))
sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
dev.off()
system(paste0("slack -F ",paste0(wd,"/",outname,"_model_selection.pdf")," ryan_todo"))

saveRDS(sub_cistopic_models,file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))
print("finshed running cistopic")

#Add cell embeddings into seurat
cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
colnames(cell_embeddings)<-sub_cistopic_models@cell.names
n_topics<-nrow(cell_embeddings)
row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
cell_embeddings<-as.data.frame(t(cell_embeddings))

#Add feature loadings into seurat
feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
row.names(feature_loadings)<-paste0("topic_",1:n_topics)
feature_loadings<-as.data.frame(t(feature_loadings))

#combined cistopic results (cistopic loadings and umap with seurat object)
cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
atac_sub@reductions$cistopic<-cistopic_obj
n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("sample","seurat_clusters"))
plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
pdf(paste0(wd,"/",outname,".umap.pdf"),width=10)
print(plt1)
print(plt2)
dev.off()
system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
saveRDS(atac_sub,"phase2.QC.SeuratObject.rds")
```
### Vibe check on cell type prediction

Plot out a heatmap of cell type scores per sample and prediction. I'm trying to figure out how specific they are and if results are concordant. Also plotting top 10 genes from Swarbrick gross cell type identification as a heatmap
Genes from Supplementary Table 9 (major classification).

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(seriation)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat_in<-readRDS("phase2.QC.SeuratObject.rds")

dat<-dat_in@meta.data


swarbrick_out<-as.data.frame(dat %>% group_by(sample,predicted.id) %>% summarize(
  swarbrick_Normal.Epithelial=median(prediction.score.Normal.Epithelial,na.rm=T),
  swarbrick_Cancer.Epithelial=median(prediction.score.Cancer.Epithelial,na.rm=T),
  swarbrick_Endothelial=median(prediction.score.Endothelial,na.rm=T),
  swarbrick_CAFs=median(prediction.score.CAFs,na.rm=T),
  swarbrick_B.cells=median(prediction.score.B.cells,na.rm=T),
  swarbrick_T.cells=median(prediction.score.T.cells,na.rm=T),
  swarbrick_PVL=median(prediction.score.PVL,na.rm=T),
  swarbrick_Myeloid=median(prediction.score.Myeloid,na.rm=T),
  swarbrick_Plasmablasts=median(prediction.score.Plasmablasts,na.rm=T)))

EMBO_out<-as.data.frame(dat %>% group_by(sample,predicted.id) %>% summarize(
  EMBO_BCell=median(EMBO_BCell,na.rm=T),
  EMBO_DC=median(EMBO_DC,na.rm=T),
  EMBO_Endo=median(EMBO_Endo,na.rm=T),
  EMBO_Fibro=median(EMBO_Fibro,na.rm=T),
  EMBO_Macro=median(EMBO_Macro,na.rm=T),
  #EMBO_Mega=median(EMBO_Mega,na.rm=T), all NA values
  EMBO_NK=median(EMBO_NK,na.rm=T),
  EMBO_TCell=median(EMBO_TCell,na.rm=T),
  EMBO_TCells=median(EMBO_TCell2,na.rm=T)))

row.names(swarbrick_out)<-paste(swarbrick_out$sample,swarbrick_out$predicted.id)
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", 
#immune
"B-cells" ="#089099", "T-cells" ="#003147", 
#other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

side_ha<-rowAnnotation(df= data.frame(celltype=swarbrick_out$predicted.id, sample=swarbrick_out$sample),
                col=list(
                    celltype=setNames(type_cols,names(type_cols)),
                    cluster=setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique(swarbrick_out$sample))),unique(swarbrick_out$sample))
                    ))
#seriate order for nice lookin heatmap
o = seriate(swarbrick_out[,3:ncol(swarbrick_out)], method = "BEA_TSP")

swarbrick<-Heatmap(swarbrick_out[,3:ncol(swarbrick_out)],
  left_annotation=side_ha,
  row_order = get_order(o, 1), column_order = get_order(o, 2),
  col=colorRamp2(c(0, max(swarbrick_out[,3:ncol(swarbrick_out)])), c("white", "red")),
  show_row_names=T)

EMBO<-Heatmap(EMBO_out[,3:ncol(EMBO_out)],
  col=colorRamp2(c(0, max(EMBO_out[,3:ncol(EMBO_out)],na.rm=T)), c("white", "red")),
  row_order=get_order(o,1))

pdf("predictions.heatmap.pdf",width=30)
swarbrick+EMBO
dev.off()
system("slack -F predictions.heatmap.pdf ryan_todo")


table(dat$sample)
#     RM_1      RM_2      RM_3      RM_4  sample_1 sample_10 sample_11 sample_12
#     1845      1461       869       922      3518     18789      5575     13303
#sample_15 sample_16 sample_19 sample_20  sample_3  sample_4  sample_5  sample_6
#     1486       574      2116       718      7698     15250      1444       661
# sample_7  sample_8  sample_9
#     2656       717      1445

table(dat[dat$predicted.id%in%c("Cancer Epithelial","Normal Epithelial"),]$sample)
#     RM_1      RM_2      RM_3      RM_4  sample_1 sample_10 sample_11 sample_12
#     1367       129       313       169        53      7023      4678     11791
#sample_15 sample_16 sample_19 sample_20  sample_3  sample_4  sample_5  sample_6
#      704       356      1648       284      4466      2768       792       448
# sample_7  sample_8  sample_9
#     1860       415      1096



#########THIS PORTION TO BE FIXED#######################

gene_list<-read.table("/home/groups/CEDAR/mulqueen/ref/breast_cancer_celltype_genelist.tsv",sep="\t",head=T)
gene_list<-as.data.frame(gene_list %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC))
gene_list<-gene_list[!duplicated(gene_list$gene),]

rna_mat<-dat_in@assays$SoupXRNA@data
rna_mat<-rna_mat[row.names(rna_mat) %in% gene_list$gene,]
rna_mat<-as.data.frame(t(as.data.frame(rna_mat)))
rna_mat$samp<-paste(dat$sample,dat$predicted.id)
rna_mat<-as.data.frame(rna_mat) %>% group_by(samp) %>% summarise_at(vars(!starts_with("samp")), ~mean(as.numeric(.x), na.rm = TRUE))
#rna_mat<-AggregateExpression(dat_in, assays = "SoupXRNA", features = gene_list$gene, return.seurat = FALSE, group.by = c("sample","predicted.id"), #sample slot = "data")
rna_mat<-as.data.frame(rna_mat)
row.names(rna_mat)<-rna_mat$samp
rna_mat<-rna_mat[colnames(rna_mat) != "samp"]
rna_mat<-as.data.frame(scale(rna_mat))


gene_list<-gene_list[gene_list$gene %in% colnames(rna_mat),]
top_ha<-columnAnnotation(df= data.frame(celltype=gene_list$cluster),
                col=list(
                    celltype=setNames(type_cols,names(type_cols))
                    ))

rna<-Heatmap(rna_mat,
  column_order=1:ncol(rna_mat),
  column_split=gene_list$cluster,
  top_annotation=top_ha,
  row_order = get_order(o, 1),  
  col=colorRamp2(c(0, max(rna_mat,na.rm=T)), c("white", "black"))
)



pdf("predictions.heatmap.pdf",width=30)
swarbrick+EMBO+rna
dev.off()
system("slack -F predictions.heatmap.pdf ryan_todo")

```
For some reason T-cells are showing high expression across the board. 


## Run ChromVAR on all data

```R
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)
library(BiocParallel)
register(SerialParam()) #using single core mode

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.SeuratObject.rds")
DefaultAssay(dat)<-"ATAC"
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(dat[["ATAC"]]))) %in% main.chroms)
dat[["ATAC"]] <- subset(dat[["ATAC"]], features = rownames(dat[["ATAC"]][keep.peaks]))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
peaks<-granges(dat[["ATAC"]])

motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, 
  pwm = pfm, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  use.counts = FALSE)

motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, 
  pwm = pfm)

dat <- SetAssayData(object = dat, 
  assay = 'ATAC', 
  slot = 'motifs', 
  new.data = motif.hg38)

dat <- RegionStats(object = dat, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

dat <- RunChromVAR( object = dat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

saveRDS(dat,file="phase2.QC.SeuratObject.rds")

```
## Using Signac Gene Activity Function
This is to generate enhancer promoter linkages at genes (by proximity). Running on all data.

```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(SeuratObjects)
library(EnsDb.Hsapiens.v86)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.SeuratObject.rds")
gene_activity<-GeneActivity(dat,process_n=10000)
saveRDS(gene_activity,file="phase2.QC.GeneActivity.rds")

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file="phase2.QC.SeuratObject.rds")

```

## Pseudobulk Clustering of Stroma, Immune and Epithelial Cells

### Epithelial Clustering
```R
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cisTopic)
library(SeuratWrappers)
library(patchwork)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
library(RColorBrewer)
library(ggplot2)
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.SeuratObject.rds")


#set up colors for samples
###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")
diag_cols<-c("IDC"="red", "DCIS"="grey","ILC"="blue","NAT"="orange")
molecular_type_cols<-c("DCIS"="grey", "ER+/PR+/HER2-"="#EBC258", "ER+/PR-/HER2-"="#F7B7BB","ER+/PR-/HER2+"="#4c9173","NA"="black")
########################################
alpha_val=0.33


celltype_cistopic_generation<-function(celltype_list=c("Cancer Epithelial","Normal Epithelial"),outname="epithelial"){
  atac_sub<-subset(dat,predicted.id %in% celltype_list)

  cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
  saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))

  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  
  saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(outname,".CisTopicObject.Rds"))
  print("finshed running cistopic")

  #Add cell embeddings into seurat
  cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
  colnames(cell_embeddings)<-sub_cistopic_models@cell.names
  n_topics<-nrow(cell_embeddings)
  row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
  cell_embeddings<-as.data.frame(t(cell_embeddings))

  #Add feature loadings into seurat
  feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
  row.names(feature_loadings)<-paste0("topic_",1:n_topics)
  feature_loadings<-as.data.frame(t(feature_loadings))

  #combined cistopic results (cistopic loadings and umap with seurat object)
  cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
  print("Cistopic Loading into Seurat")
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  pdf(paste0(outname,".cistopic.umap.pdf"),width=10)
  print(plt1)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistopic.umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(outname,".SeuratObject.rds"))
  }

celltype_cistopic_generation(celltype_list=c("Cancer Epithelial","Normal Epithelial"),outname="epithelial")
celltype_cistopic_generation(celltype_list=c("B-cells","T-cells","Myeloid","Plasmablasts"),outname="immune")
celltype_cistopic_generation(celltype_list=c("CAFs","Endothelial","PVL"),outname="stromal")

#Rerun other clustering now that data is subset
celltype_clustering<-function(x,outname){
  dat<-readRDS(x)

  #set up colors for samples
  my_cols = brewer.pal(1,"Spectral")
  alpha_val=0.33

  #RNA Processing
  DefaultAssay(dat) <- "SoupXRNA"
  dat <- SCTransform(dat)
  dat <- RunPCA(dat)
  dat<- RunUMAP(object = dat, reduction.name="rna_umap", reduction="pca", assay = "SoupXRNA", verbose = TRUE, dims=1:50 ) 
  p1<-DimPlot(dat,reduction="rna_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("RNA UMAP")+theme(legend.position="none")

  #DNA Accessibility processing
  DefaultAssay(dat) <- "peaks"
  dat <- FindTopFeatures(dat, min.cutoff = 5)
  dat <- RunTFIDF(dat)
  dat <- RunSVD(dat)
  dat<- RunUMAP(object = dat, reduction.name="atac_umap", reduction="lsi", assay = "peaks", verbose = TRUE, dims=2:40 )
  p2<-DimPlot(dat,reduction="atac_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("ATAC UMAP")+theme(legend.position="none")


  # build a joint neighbor graph using both assays (ATAC LSI)
    dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:40),modality.weight.name = "RNA.weight", verbose = TRUE )
    # build a joint UMAP visualization
    dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_umap", assay = "SoupXRNA", verbose = TRUE ) 
    p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP Doublets")+theme(legend.position="none")

  #Try multimodal with cistopic
    dat <- RunUMAP(object = dat, reduction="cistopic", reduction.name="cistopic_umap", dims=1:ncol(dat@reductions$cistopic), assay = "peaks", verbose = TRUE )
    # build a joint neighbor graph using both assays
    dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("pca", "cistopic"), dims.list = list(1:50, 1:ncol(dat@reductions$cistopic)), modality.weight.name = "RNA.weight", verbose = TRUE )
    # build a joint UMAP visualization
    dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_umap", assay = "SoupXRNA", verbose = TRUE )

  #plot cistopic umap too
  p4<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  p5<-DimPlot(dat,reduction="cistopic_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Cistopic UMAP")+theme(legend.position="none")
  p6<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  #Cluster on multimodal graph
  dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")


  #Finally Plot results
  plt<-(p1 | p2)/(p3 | p4)/(p5|p6)
  ggsave(plt,file=paste0(outname,".umap.pdf"))
  system(paste0("slack -F ",paste0(outname,".umap.pdf")," ryan_todo"))
  saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))
}

celltype_clustering(x="stromal.SeuratObject.rds",outname="stromal")
celltype_clustering(x="immune.SeuratObject.rds",outname="immune")
celltype_clustering(x="epithelial.SeuratObject.rds",outname="epithelial")



```

### Integration: Now Clustering together on RNA profiles using harmony to integrate

```R
library(harmony)
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")


###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

diag_cols<-c("IDC"="red", "DCIS"="grey","NAT"="lightblue","ILC"="green")

molecular_type_cols<-c("DCIS"="grey", "ER+/PR-/HER2-"="#EBC258", "ER+/PR-/HER2+"="#F7B7BB","ER+/PR+/HER2-"="#ff6699","NA"="lightblue")
########################################


harmony_sample_integration<-function(x,outname){
  dat<-readRDS(x)
  dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_atac",assay.use="ATAC",reduction="cistopic",project.dim=F)
  dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_rna",assay.use="RNA",reduction="pca",project.dim=F)

  dat<-RunUMAP(dat,reduction.name="harmonyumap_rna",reduction = "harmony_rna",dims=1:dim(dat@reductions$harmony_rna)[2]) 
  dat<-RunUMAP(dat,reduction.name="harmonyumap_atac",reduction = "harmony_atac",dims=1:dim(dat@reductions$harmony_atac)[2]) 

  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(
    object = dat,
    reduction.list = list("harmony_rna", "harmony_atac"), 
    dims.list = list(1:dim(dat@reductions$harmony_rna)[2], 1:dim(dat@reductions$harmony_atac)[2]), 
    modality.weight.name = "multimodal.weight",
    weighted.nn.name="multimodal_harmony.nn",
    verbose = TRUE
  )
  # build a joint UMAP Harmony visualization
  dat <- RunUMAP(object = dat, nn.name = "multimodal_harmony.nn",reduction.name="multimodal_harmony_umap", assay = "SoupXRNA", verbose = TRUE ) 

  i="predicted.id"
  plt1<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i,cols=type_cols)
  ggsave(plt1,file=paste0(outname,".",i,".pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  i="diagnosis"
  plt1<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i,cols=diag_cols)
  ggsave(plt1,file=paste0(outname,".",i,".pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  i="sample"
  plt1<-DimPlot(dat,reduction="multimodal_harmony_umap",group.by=i)
  ggsave(plt1,file=paste0(outname,".",i,".pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,".",i,".pdf ryan_todo"))

  saveRDS(dat,file=x)
}

harmony_sample_integration(x="stromal.SeuratObject.rds",outname="stromal")
harmony_sample_integration(x="immune.SeuratObject.rds",outname="immune")
harmony_sample_integration(x="epithelial.SeuratObject.rds",outname="epithelial")
harmony_sample_integration(x="phase2.QC.SeuratObject.rds",outname="all_cells")

```

## Transcription Factor Expression Markers

Based on seurat tutorial https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1
Using average AUC to define markers that work across modalities (RNA, Gene Activity, and TF motifs). Doing this across cell types, and then within cell types across diagnoses.

```R
library(Signac)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(seriation)
library(viridis)
library(circlize)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")



###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

diag_cols<-c("IDC"="red", "DCIS"="grey","NAT"="lightblue","ILC"="green")

molecular_type_cols<-c("DCIS"="grey", "ER+/PR-/HER2-"="#EBC258", "ER+/PR-/HER2+"="#F7B7BB","ER+/PR+/HER2-"="#ff6699","NA"="lightblue")
########################################


#Grab top overlapping TFs
topTFs <- function(markers_list,celltype, padj.cutoff = 1e-2,rna=NA,ga=NA,motifs=NA) {
  ctmarkers_rna <- dplyr::filter(
    rna, RNA.group == celltype) %>% 
    arrange(-RNA.auc)

    if(is.data.frame(motifs)) {
    ctmarkers_motif <- dplyr::filter(
      motifs, chromvar.group == celltype) %>% 
      arrange(-chromvar.auc)
    }

    if(is.data.frame(ga)) {
    ctmarkers_ga<- dplyr::filter(
      ga, GeneActivity.group == celltype) %>% 
      arrange(-GeneActivity.auc)
    }

    if(is.data.frame(motifs) && is.data.frame(ga)){    
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs <- inner_join(
        x = top_tfs ,
        y = ctmarkers_ga [,c(2, 11, 6, 7)], by = "gene"
      )
    }else if(is.data.frame(motifs)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
    } else if (is.data.frame(ga)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_ga[,c(2, 11, 6, 7)], by = "gene"
      )
    } 
  auc_colnames<-grep(".auc$",colnames(top_tfs))
  top_tfs$avg_auc <-  rowMeans(top_tfs[auc_colnames])
  top_tfs <- arrange(top_tfs, -avg_auc)
  top_tfs$celltype<-celltype
  return(top_tfs)
}

#Identify top markers
Identify_Marker_TFs<-function(x,group_by.="predicted.id",assay.="RNA"){
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by., assay = 'data', seurat_assay = assay.)
    colnames(markers) <- paste(assay., colnames(markers),sep=".")
    if (assay. == "chromvar") {
      motif.names <- markers[,paste0(assay.,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="ATAC")
    } else {
    markers$gene <- markers[,paste0(assay.,".feature")]
    }
    return(markers) 
}

#Average markers across groups
average_features<-function(x=hg38_atac,features=da_tf_markers$motif.feature,assay="chromvar",group_by.="predicted.id"){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    dat_motif<-x[[assay]]@data[features,]
    dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
    sum_motif<-split(dat_motif,x@meta.data[,group_by.]) #group by rows to seurat clusters
    sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take average across group
    sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame

    sum_motif<-t(scale(sum_motif))
    sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

#Make a heatmap of aligned multiple modalities
plot_top_TFs<-function(x=stromal,tf_markers=da_tf_markers,prefix="stromal",group_by.="predicted.id",CHROMVAR=TRUE,GA=TRUE){
    tf_rna<-average_features(x=x,features=tf_markers$gene,assay="RNA",group_by.=group_by.)
    tf_rna<-tf_rna[row.names(tf_rna) %in% tf_markers$gene,]

  if(CHROMVAR){
    tf_motif<-average_features(x=x,features=tf_markers$chromvar.feature,assay="chromvar",group_by.=group_by.)
    tf_motif<-tf_motif[row.names(tf_motif) %in% tf_markers$chromvar.feature,]
    row.names(tf_motif)<-tf_markers[tf_markers$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
  }

  if(GA){
    tf_ga<-average_features(x=x,features=tf_markers$gene,assay="GeneActivity",group_by.=group_by.)
    tf_ga<-tf_ga[row.names(tf_ga) %in% tf_markers$gene,]
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_ga<-tf_ga[markers_list,]

  }
  if(GA&&CHROMVAR){
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
  }

    #set up heatmap seriation and order by RNA
    o = seriate(max(tf_rna) - tf_rna, method = "BEA_TSP")
    saveRDS(o,file=paste0(prefix,".geneactivity.dend.rds")) 
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o,1),]$RNA.auc)
    colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.80,0.95)),plasma(3))

  if(CHROMVAR){
    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.80,0.95)),cividis(3))
    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="ATAC",motifs = motif_list[get_order(o,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

  }
  if(GA){
    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.80,0.95)),magma(3))

  }

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))


    rna_auc<-Heatmap(side_ha_rna,
        row_order = get_order(o,1),
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))

    rna_plot<-Heatmap(tf_rna,
        row_order = get_order(o,1),
        column_order = get_order(o,2),
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90)

  if(GA){
      ga_auc<-Heatmap(side_ha_ga,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)
  }
  if(CHROMVAR){
      motif_auc<-Heatmap(side_ha_motif,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
  }

  if(all(CHROMVAR,GA)){
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot)
  } else if(CHROMVAR){
      plt1<-draw(rna_auc+rna_plot+motif_auc+motif_plot)
  } else {
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot)
  }


    pdf(paste0(prefix,".tf.heatmap.pdf"))
    print(plt1)
    dev.off()

    system(paste0("slack -F ",prefix,".tf.heatmap.pdf ryan_todo"))
    system(paste0("slack -F ",prefix,".tf.heatmap.motif.pdf ryan_todo"))
}

#Final wrapper function
run_top_TFs<-function(obj=stromal,prefix="stromal",i="predicted.id"){
  markers<-lapply(c("RNA","GeneActivity","chromvar"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay))
  names(markers)<-c("RNA","GeneActivity","chromvar")
  markers_out<-do.call("rbind",lapply(unique(obj@meta.data[,i]),function(x) head(topTFs(markers_list=markers,celltype=x,rna=markers$RNA,ga=markers$GeneActivity,motifs=markers$chromvar),n=10))) #grab top 5 TF markers per celltype
  dim(markers_out)
  markers_out<-markers_out[!duplicated(markers_out$gene),]
  dim(markers_out)
  saveRDS(markers_out,file=paste0(prefix,"_celltype_TF_markers.RDS"))
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
  plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i,CHROMVAR=TRUE,GA=TRUE)
}

#stromal celltypes
stromal<-readRDS("stromal.SeuratObject.rds")
run_top_TFs(obj=stromal,prefix="stromal",i="predicted.id")

#immune celltypes
immune<-readRDS("immune.SeuratObject.rds")
run_top_TFs(obj=immune,prefix="immune",i="predicted.id")


#all cells celltypes
dat<-readRDS("phase2.QC.SeuratObject.rds")
run_top_TFs(obj=dat,prefix="dat",i="predicted.id")

#Per cell type, run diagnosis differences
#Removing ILC cells first just because we don't have enough samples for comparison
dat<-subset(dat,diagnosis!="ILC")

for(k in unique(dat$predicted.id)){
  dat_sub<-subset(dat,predicted.id==k)
  run_top_TFs(obj=dat_sub,prefix=paste0(k,"_diagnosis"),i="diagnosis")
  print(paste("Done with",k))
}



```

## Bar plots across cells

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr) 
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
 

###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

diag_cols<-c("IDC"="red", "DCIS"="grey")

molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
########################################



dat<-readRDS("phase2.QC.SeuratObject.rds")




  # build a joint UMAP visualization
  dat2 <- RunUMAP(
    object = dat2,
    nn.name = "weighted.nn",
    reduction.name="multimodal_umap",
    assay = "SoupXRNA",
    verbose = TRUE
  )


p1<-DimPlot(dat,reduction="rna_umap",group.by="sample")+ggtitle("RNA UMAP")
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample")+ggtitle("ATAC UMAP")
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP")

plt<-(p1|p2)/p3
ggsave(plt,file="all_samples.umap.pdf")
system(paste0("slack -F ","all_samples.umap.pdf"," ryan_todo"))

#Finally Plot results
ggsave(plt,file=paste0(wd,"/",outname,".umap2.pdf"))
system(paste0("slack -F ",paste0(wd,"/",outname,".umap2.pdf")," ryan_todo"))
saveRDS(dat,file=paste0(wd,"/",outname,".QC.SeuratObject.rds"))
}


#Set up metadata and set up facet labels as factors for ordering
metadat<-as.data.frame(dat@meta.data)
metadat$diagnosis = factor(metadat$diagnosis, levels=c("NAT","DCIS","IDC","ILC"), labels=c("NAT","DCIS","IDC","ILC")) 
metadat$molecular_type = factor(metadat$molecular_type, levels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2-","ER+/PR-/HER2+"), labels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2-","ER+/PR-/HER2+")) 


#Cells PF (log10)
metadat$epi<-"Nonepi"
metadat[metadat$predicted.id %in% c("Cancer Epithelial","Normal Epithelial"),]$epi<-"Epi"
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,epi) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=epi,y=n))+geom_bar(stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="barplot_qc_cellcount.pdf")
system("slack -F barplot_qc_cellcount.pdf ryan_todo")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_fill_manual(values=type_cols)+facet_wrap(.~diagnosis+molecular_type,nrow=1,scales="free")
ggsave(plt1,file="barplot_qc_celltype.pdf")
system("slack -F barplot_qc_celltype.pdf ryan_todo")



```

## Plot of Differential Genes across Normal epithelial (NAT) DCIS and IDC

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(SeuratWrappers)
library(cisTopic)
library(patchwork)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

atac_sub<-readRDS("epithelial.SeuratObject.rds")
Idents(atac_sub)<-atac_sub$diagnosis
DefaultAssay(atac_sub)<-"SoupXRNA"
RNA_markers <- FindMarkers(atac_sub, ident.1 = "IDC", ident.2 = c("DCIS","NAT"), min.pct = 0.1)
RNA_markers$gene_name<-row.names(RNA_markers)
DefaultAssay(atac_sub)<-"ATAC"
ATAC_markers <- FindMarkers(atac_sub,ident.1 = "IDC", ident.2 = c("DCIS","NAT"), min.pct = 0.1)
ATAC_markers$da_region<-row.names(ATAC_markers)
closest_genes <- ClosestFeature(atac_sub,ATAC_markers$da_region)
ATAC_markers<-cbind(ATAC_markers,closest_genes)

rna<-RNA_markers[RNA_markers$gene_name %in% ATAC_markers$gene_name,]
atac<-ATAC_markers[ATAC_markers$gene_name %in% rna$gene_name,]



cov_plots<-function(dat=atac_sub,gene_name,idents_in){
  plt_cov <- CoveragePlot(
    object = atac_sub,
    region = gene_name,
    features = gene_name,
    assay="ATAC",
    expression.assay = "SoupXRNA",
    extend.upstream = 5000,
    extend.downstream = 5000,
    idents=idents_in)
  plt_feat <- FeaturePlot(
    object = atac_sub,
    features = gene_name,
    raster=T,
    reduction="multimodal_umap",
    order=T)
  return((plt_feat|plt_cov)+ggtitle(gene_name))
}


DefaultAssay(atac_sub)<-"SoupXRNA"
for (i in c(rna$gene_name[1:25])){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}

for (i in c("SOX10","SOX9","SOX4","SOX2","TEAD4","RUNX1")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("FOXM1","FOXA1","FOXA3","GRHL2","FOXP1","ATF3")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("HOXB13","EN1","DLX4","TBX15","SLC6A12","PAX6","FAM83A","ERICH5")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}


for (i in c("ESR1")){
  plt<-cov_plots(dat=atac_sub,gene_name=i,idents_in=c("NAT","DCIS","IDC"))
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}

#TFS
#IDC Fox Family (GRHL2) 
#ILC Sox Family TEAD RUNX EGR1 RPBJ HMGA1
#DCIS STAT3/BCL9
#DCIS more likely to be invasice (Methylation IDd) HOXB13 EN1 DLX4 TBX15 SLC6A12 PAX6 

#GENES
#FAM83A
#ERICH5
```
## Comparison of cell types across diagnoses and other factors.

```R

```
# 3D Plotting in Blender

```R
#3d umap
out_3d <- RunUMAP(object=dat, n.components=3, reduction.name="harmonyumap",reduction = "harmony", dims = 1:20)
#format
#Astrocytes    TAGGTCCGACGTACTAGGGCCTCGGTCTATGGCCTA    4.24424248742567    -1.74691044949975    -6.48374510684418    #1C7D54
#Astrocytes    ATTCAGAAGCATCGCGCAGCCAGACTCTATGGCCTA    3.60301401455387    -1.96493138894082    -6.47136162049336    #1C7D54
#Astrocytes    TCAACGAGTTCGCGATGGTCAGAGCCCGCCGATATC    5.51775913941571    -1.87741656898663    -6.76243310557264    #1C7D54
out_3d_dat<-as.data.frame(cbind(out_3d@meta.data[,c("predicted.id")],row.names(out_3d@meta.data),Embeddings(out_3d,"harmonyumap")))
colnames(out_3d_dat)[1]<-"predicted.id"
col_dat<-as.data.frame(cbind(names(type_cols),unname(type_cols)))
colnames(col_dat)<-c("predicted.id","cols")
dat_out<-merge(out_3d_dat,col_dat,by="predicted.id")
write.table(dat_out,file="multiome_tumor.tsv",sep="\t",quote=F,col.names=F,row.names=F)
system("slack -F multiome_tumor.tsv ryan_todo")

```





# Low Pass Whole Genome Sequencing Data

## Align to hg38 with bwa-mem

```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2
cat readme.txt 
#Run     Lane    Sample  I7 Index ID     Index1  I5 Index ID     Index2
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_A1-12        D712    AGCGATAG        D501    AGGCTATA
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_B1-12        D712    AGCGATAG        D502    GCCTCTAT
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_C1-12        D712    AGCGATAG        D503    AGGATAGG
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_D1-12        D712    AGCGATAG        D504    TCAGAGCC
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_E1-12        D712    AGCGATAG        D505    CTTCGCCT
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_F1-12        D712    AGCGATAG        D506    TAAGATTA
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_G1-12        D712    AGCGATAG        D507    ACGTCCTG
#220929_A01058_0265_AHNGVCDRX2   1       EXP220921HM_BC32-3_H1-12        D712    AGCGATAG        D508    GTCAGTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG01   N701    TAAGGCGA        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG03   N702    CGTACTAG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG04   N703    AGGCAGAA        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG05   N704    TCCTGAGC        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG06   N705    GGACTCCT        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG07   N706    TAGGCATG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG08   N707    CTCTCTAC        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG09   N710    CGAGGCTG        S505    CTCCTTAC
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG10   N701    TAAGGCGA        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG11   N702    CGTACTAG        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG12   N703    AGGCAGAA        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG15   N704    TCCTGAGC        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG16   N705    GGACTCCT        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG19   N706    TAGGCATG        S506    TATGCAGT
#220929_A01058_0265_AHNGVCDRX2   2       EXP220921HM_BCMM_WG20   N707    CTCTCTAC        S506    TATGCAGT

#I'm taking the BCMM samples as the low pass whole genome

cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM




```

## Batch script for Alignment
Using the bwa mem for alignment.

wgs_alignment.sbatch
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-15
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=30 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=5:00:00 ## ask for 3 hour on the node
#SBATCH --

fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
ref="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
file_in=`ls $fastq_dir/*WG*R1*.fastq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' `

#Align and output as bam file
R1=$file_in
R2=`echo $bam_in | awk '{ gsub("_R1_", "_R2_"); print $0}' `
output_name=${R1::-9}".batch.bam"
bwa mem -t 10 $ref $R1 $R2 | samtools sort -T . -@5 - | samtools view -b -@5 - > $output_name

```

```bash
sbatch wgs_alignment.sbatch
```

## Batch script for deduplication via Picard
wgs_dedup.sbatch
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-15
#SBATCH --tasks-per-node=15 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 3 hour on the node
#SBATCH --

picard_dir="/home/groups/CEDAR/tools/picard-tools-1.119/"
fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
list_files=`ls $fastq_dir/*WG*R1*.bam`
bam_in=`ls $fastq_dir/*WG*R1*.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
i=$bam_in
output_name=${i::-4}".dedup.bam"
output_metrics=${i::-4}".dedup.metrics.txt"

#picard mark duplicates
java -jar ${picard_dir}/MarkDuplicates.jar \
        I=$i \
        O=$output_name \
        M=$output_metrics
```

```bash
sbatch wgs_dedup.sbatch
```

## Install GATK4
Following https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants

```bash
cd /home/groups/CEDAR/mulqueen/src/gatk
conda install -c gatk4

```


```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-15
#SBATCH --tasks-per-node=15 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=3:00:00 ## ask for 3 hour on the node
#SBATCH --

gatk="/home/groups/CEDAR/mulqueen/src/gatk/gatk-4.2.0.0/gatk"
ref="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
fastq_dir="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/EXP220921HM/220929_A01058_0265_AHNGVCDRX2/EXP220921HM"
list_files=`ls $fastq_dir/*WG*R1*.dedup.bam`
bam_in=`ls $fastq_dir/*WG*R1*.dedup.bam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
i=$bam_in
output_name=${i::-10}



$gatk PreprocessIntervals \
    -R $ref \
    --bin-length 1000 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O preprocessed.1000.interval_list
```



#ER binding poor and good outcome from patients, overlap with ATAC data
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Poor%20outcome%20ER%20regions.bed.gz
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Good%20outcome%20ER%20regions.bed.gz

#sample specificity between cancer and normal epithelial, superenhancers are open more frequently
#plot epithelial cell subtype per sample
#use cnv to infer changes, then do lineage comparisons, plot BAF across regions
#use HMMcopy to test ATAC changes
#plot cell prediction specificity (predicted ID values density plot colored by predicted ID cell types)


-->