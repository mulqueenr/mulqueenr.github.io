---
title: 10X Multiome Tumor All Together
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_phase2/
category: CEDAR
---

Multiome processing for 10X multiome data on Primary Tumors (Phase 1)

## File Location
```bash
mkdir /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
cd /home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2
sftp mulqueen@nix.ohsu.edu
#enter password
get -r /data/EXP220628HM
get -r /data/EXP220629HM

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
"DCIS", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "DCIS", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2-", "ER+/PR-/HER2-", "ER+/PR+/HER2-", "NA", "DCIS", "NA", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR-/HER2+", "ER+/PR-NA/HER2-", "NA")))

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

## Use EMBO Paper Cell Type Signatures
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

immune_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/ImmuneMarkers2.txt",header=T,sep="\t")
immune_in<-lapply(split(immune_in,immune_in$CellType),function(x) x$Signatures)#split up data frame to a named list of genes per cell type
names(immune_in)<-paste0("EMBO_",names(immune_in))#rename the list just so we can track the source

#using both the given PAM50 and the EMBO supplied more extensive gene list
PAM50_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/PAM50.txt",header=T,sep="\t")
PAM50_in<-lapply(split(PAM50_in,PAM50_in$Subtype),function(x) x$Gene)#split up data frame to a named list of genes per cell type
names(PAM50_in)<-paste0("PAM50_",names(PAM50_in))
features_in=c(immune_in,PAM50_in)   

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
  dat<-AddModuleScore(dat,features=lineage_in,names=names(lineage_in),assay="SoupXRNA",seed=123,search=TRUE)
  colnames(dat@meta.data)[startsWith(prefix="Cluster",colnames(dat@meta.data))]<-c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str") #Rename them
  for(i in 1:length(features_in)){
    features_in[[i]]<-features_in[[i]][features_in[[i]] %in% row.names(dat[["SoupXRNA"]])] #make sure gene names match
    dat<-MetaFeature(dat,features=c(features_in[[i]]),meta.name=names(features_in)[i],assay="SoupXRNA")
  }
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
  pam_list<-  c("PAM50_Basal","PAM50_Her2","PAM50_LumA","PAM50_LumB","PAM50_Normal")
  max_pam<-lapply(1:nrow(met),function(i) pam_list[which(met[i,pam_list]==max(met[i,pam_list],na.rm=T))])
  max_pam<-unlist(lapply(1:length(max_pam),function(i) do.call("paste",as.list(max_pam[[i]]))))
  max_pam<-unlist(lapply(max_pam,function(i) gsub("PAM50_","",i)))
  names(max_pam)<-row.names(met)
  dat<-AddMetaData(dat,max_pam,col.name="PAM50_designation")
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

<!--
Plot differences between epithelial lineages

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
library(reshape2)
library(philentropy)
library(dendextend)
library(ggrepel)
  library(TFBSTools)
library(JASPAR2020)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

find_markers_epithelial_lineage<-function(x,assay_name,logfc.threshold_set=0,min.pct_set=0,latent.vars="NA",i1="epithel_1",i2="epithel_2"){
  if(!(latent.vars=="NA")){
  x_da<-FindMarkers(object=x,ident.1 = i1, ident.2=i2, test.use = 'LR', logfc.threshold=logfc.threshold_set,latent.vars = latent.vars, only.pos=F, assay=assay_name,min.pct=min.pct_set)
  } else {
  x_da<-FindMarkers(object=x,ident.1 = i1, ident.2=i2, test.use = 'LR', logfc.threshold=logfc.threshold_set, only.pos=F, assay=assay_name,min.pct=min.pct_set)
  }
  #if statements to handle formating the different modalities
  if(assay_name=="chromvar"){ #make chromvar names readable
      print("Translating Chromvar Motif Names for Plot")
      x_da$out_name <- unlist(lapply(unlist(lapply(row.names(x_da), function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  } else{
    x_da$out_name<-row.names(x_da)
  }
  #if(assay_name=="peaks"){ #assign peak output to closest features
  #  x_da_peaks<-colsplit(row.names(x_da), "\\-", names=c("chr", "start", "end"))
  #  x_da_peaks<-GRanges(seqnames = x_da_peaks$chr, ranges = IRanges(start = x_da_peaks$start, end = x_da_peaks$end))
  #  DefaultAssay(x)<-"peaks"
  #  closest_genes <- ClosestFeature(x,x_da_peaks)
  #  x_da <-cbind(x_da ,closest_genes)
  #}

  x_da$sig<-ifelse(x_da$p_val_adj<=0.05,"sig","non_sig")
  x_da$label<-""

  if(assay_name!="peaks"){ #ignore adding peaks names (too many), and add top 10 significant names to either side of comparison
    ident_1_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC>0)),] ,n=10))#significant, is ident1 specific, is top 10
    ident_2_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC<0)),] ,n=10))#significant, is ident1 specific, is top 10
    x_da[c(ident_1_labels,ident_2_labels),]$label<-x_da[c(ident_1_labels,ident_2_labels),]$out_name
        #trim cistrome cell line from name (for readability)
  }
  x_scale<-max(abs(x_da$avg_log2FC))*1.1 #find proper scaling for x axis

  plt<-ggplot(x_da,aes(x=avg_log2FC,y=-log10(p_val_adj),color=sig,label=label,alpha=0.1))+
    geom_point(size=0.5)+
    theme_bw()+
    scale_fill_manual(values=c("non_sig"="#999999", "sig"="#FF0000"))+
    xlim(c(-x_scale,x_scale))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=-log10(0.05))+
    geom_text_repel(size=2,segment.size=0.1,max.overlaps=Inf,min.segment.length = 0, nudge_y = 1, segment.angle = 20,color="black") +
    theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
    ggsave(plt,file=paste0(outname,"_",assay_name,"_",i1,"v",i2,".pdf"),width=2,units="in",height=3)
  system(paste0("slack -F ",outname,"_",assay_name,"_",i1,"v",i2,".pdf"," ryan_todo"))
  write.table(x_da,file=paste0(outname,"_",assay_name,"_",i1,"v",i2,".markers.txt"),col.names=T,sep="\t")
  system(paste0("slack -F ",paste0(outname,"_",assay_name,"_",i1,"v",i2,".markers.txt")," ryan_todo"))
}


sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))


  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

  ref_cols<-c(
  "TRUE"="grey",
  "FALSE"="red")

  lineage_cols<-c(
    "epithel_1"="#00B050",
    "epithel_2"="#00b0f0",
    "epithel_3"="#ff0000",
    "NA"="grey")
  #########################################
dat<-sample_in[7]

####RUNNING INFERCNV#####
infercnv_lineage_differences<-function(dat){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  dat_file_path=dat
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-3)
  dat<-readRDS(dat)
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  DefaultAssay(dat)<-"RNA"
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-17) #redo outname to remove "seuratobject" from filename
  dend<-readRDS(file=paste0(outname,".inferCNV.dend.Rds")) #save dendrogram
  tree_info<-data.frame(cutree(dend,3)) #i know this is 3 manually, but maybe it is saved in the object somewhere?
  colnames(tree_info)<-"epithelial_lineage"
  tree_info$epithelial_lineage<-paste0("epithel_",tree_info$epithelial_lineage)
  dat<-AddMetaData(dat,metadata=tree_info,col.name="epithelial_lineage")
  Idents(dat)<-dat$epithelial_lineage
  plt<-DimPlot(dat,group.by="epithelial_lineage",cols=lineage_cols,reduction="rna_umap")
  ggsave(plt,file=paste0(outname,".inferCNV.lineage.dim.pdf"))
  system(paste0("slack -F ",paste0(outname,".inferCNV.lineage.dim.pdf")," ryan_todo"))

  meta_dat<-as.data.frame(dat@meta.data)
  meta_dat<-meta_dat[!isNA(meta_dat$epithelial_lineage),]
  plt1<-ggplot(meta_dat,aes(x=epithelial_lineage,y=Basal_SC,color=epithelial_lineage))+geom_boxplot()+theme_minimal()+  theme(legend.position="none")
  plt2<-ggplot(meta_dat,aes(x=epithelial_lineage,y=Her2E_SC,color=epithelial_lineage))+geom_boxplot()+theme_minimal()+  theme(legend.position="none")

  plt3<-ggplot(meta_dat,aes(x=epithelial_lineage,y=LumA_SC,color=epithelial_lineage))+geom_boxplot()+theme_minimal()+  theme(legend.position="none")

  plt4<-ggplot(meta_dat,aes(x=epithelial_lineage,y=LumB_SC,color=epithelial_lineage))+geom_boxplot()+theme_minimal()+  theme(legend.position="none")
  plt5<-ggplot(meta_dat,aes(x=epithelial_lineage,y=proliferation_score,color=epithelial_lineage))+geom_boxplot()+theme_minimal()+  theme(legend.position="none")
  plt<-plt1|plt2|plt3|plt4|plt5
  ggsave(plt,file="test.pdf")
  system("slack -F test.pdf ryan_todo")


  find_markers_epithelial_lineage(dat,i1="epithel_2",i2="epithel_3",assay_name="SCT",latent.vars="nFeature_SCT")
  find_markers_epithelial_lineage(dat,i1="epithel_2",i2="epithel_3",assay_name="peaks",latent.vars="atac_peak_region_fragments",logfc.threshold_set=0.1,min.pct_set=0.1)
  find_markers_epithelial_lineage(dat,i1="epithel_1",i2="epithel_3",assay_name="chromvar")


  gene_name="SRRM2"
  Idents(dat)<-dat$seurat_clusters
  dat@assays$peaks
  frag.path=paste0(dirname(dat_file_path),"/atac_fragments.tsv.gz")
  fragments <- Fragments(dat)  # get list of fragment objects
  fragments_in<-CreateFragmentObject(path=frag.path)
  fragments <- UpdatePath(fragments, new.path = frag.path)
  Fragments(dat) <- NULL  # remove fragment information from assay

  plt_cov <- CoveragePlot(
    object = dat,
    region = gene_name,
    features = gene_name,
    expression.assay = "SCT",
    assay="peaks",
    extend.upstream=100000,
    extend_downstream=100000,
    Links=F)


  ggsave(plt_cov,file="test.pdf")
  system("slack -F test.pdf ryan_todo")



  # gather the footprinting information for sets of motifs

bone <- Footprint(
  object = dat,
  motif.name = c("MA0748.2"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="peaks"
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(bone, features = c("GATA3"))
```
--->

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

lapply(c(12),function(x) casper_per_sample(x))

#1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10

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

  dist_method="euclidean"
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


compare_RNA_cnv_results<-function(x){
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

  infercnv_dend<-readRDS(paste0(wd,"/",outname,"_inferCNV","/",outname,".inferCNV.dend.Rds"))
  casper_dend<-readRDS(paste0(dir_in,"_casper/",sample_name,".casper.dend.Rds")) 

  #use Bkplot to define cluster count to use
  dl<-dendlist(intersect_trees(infercnv_dend,casper_dend))
  pdf(paste0(dir_in,"/",sample_name,".RNA.CNVs.Bkplot.pdf"),width=20)
  Bk_plot(dl[[1]],dl[[2]],p.adjust.method="bonferroni",k=2:25,xlim=c(0,25))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"/",sample_name,".RNA.CNVs.Bkplot.pdf")," ryan_todo"))

  k_search<-find_k(infercnv_dend,krange=2:10) #search for optimal K from 2-10
  infer_cnv_k_clus_number<-k_search$nc
  print(paste("Determined ",infer_cnv_k_clus_number," of clusters for InferCNV."))
  k_search<-find_k(casper_dend,krange=2:10) #search for optimal K from 2-10
  casper_k_clus_number<-k_search$nc
  print(paste("Determined ",casper_k_clus_number," of clusters for Casper."))
  
  k_clus_number<-max(c(infer_cnv_k_clus_number,casper_k_clus_number)) #use max number of clusters
  print(paste("Using ",k_clus_number," of clusters for both."))

  dat<-AddMetaData(dat,cutree(infercnv_dend,k=k_clus_number),col.name="InferCNV_clusters")
  dat<-AddMetaData(dat,cutree(casper_dend,k=k_clus_number),col.name="CaSpER_clusters")
  saveRDS(dat,file_in)

  metadata<-dat@meta.data
  metadata<-metadata[metadata$predicted.id %in% c("Cancer Epithelial","Normal Epithelial" ),]
  metadata<-metadata[!is.na(metadata$InferCNV_clusters) & !is.na(metadata$CaSpER_clusters),]
  metadata$InferCNV_clusters<-as.character(metadata$InferCNV_clusters)
  metadata$CaSpER_clusters<-as.character(metadata$CaSpER_clusters)

  plt<-ggplot(metadata,
       aes(y =   match(row.names(metadata),labels(infercnv_dend)[order.dendrogram(infercnv_dend)]),
        axis1 = InferCNV_clusters,
        axis2 = CaSpER_clusters)) +
  geom_alluvium(aes(fill = predicted.id), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("InferCNV_clusters", "CaSpER_clusters"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Copy Number Lineages")
  ggsave(plt,file=paste0(dir_in,"/",sample_name,".RNA.CNVs.alluvial.pdf"))
  system(paste0("slack -F ",paste0(dir_in,"/",sample_name,".RNA.CNVs.alluvial.pdf")," ryan_todo"))

  if(nrow(metadata)<500){
  infercnv_dend <- color_branches(infercnv_dend, k = k_clus_number)    #split breakpoint object by clusters
  casper_dend <- color_branches(casper_dend, k = k_clus_number)    #split breakpoint object by clusters
  dl<-dendlist(intersect_trees(infercnv_dend,casper_dend))
  pdf(paste0(dir_in,"/",sample_name,".RNA.CNVs.tanglegram.pdf"),width=20)
  tanglegram(dl, 
    main_left="InferCNV",
    main_right="CaSpER",
    main="Tangled",
    sub=paste("entanglement =", round(entanglement(dl), 2)))
  dl_sorted<-untangle(dl, method = "step2side")
  tanglegram(dl_sorted,    
    main_left="InferCNV",
    main_right="CaSpER",
    main="Untangled",
    sub=paste("entanglement =", round(entanglement(dl_sorted), 2)))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"/",sample_name,".RNA.CNVs.tanglegram.pdf")," ryan_todo"))
  cor_cophenetic(dl[[1]],dl[[2]])
  cor_bakers_gamma(dl[[1]],dl[[2]])
  }
}

lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) infercnv_per_sample_plot(x))
lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) casper_per_sample_plot(x))
lapply(c(1,3,5,6,7,8,9,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) compare_RNA_cnv_results(x))


```

Compare CASPER and InferCNV Results

```R
library(Seurat)
library(Signac)
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
library(dendextend)
library(CaSpER)

cnv_comparisons<-function(x){
 
}

lapply(c("RM_1","RM_2","RM_3","RM_4",11,4,10,12), function(x) cnv_comparisons(x))
#to run 15, 16, 19, 20
#ran 1,3,5,6,7,8,9
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
                      minFrags=1e4,
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
    logNorm = FALSE,
    maxZero=2000,
    imputeZeros = FALSE,
    blacklistProp = 0.8,
    blacklistCutoff=125,
    dividingFactor=1,
    upperFilterQuantile = 0.95)
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
  #apply additional filters
  scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)
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
    useDummyCells = TRUE,
    propDummy=0.25,
    minMix=0.01,
    deltaMean = 0.03,
    deltaBIC2 = 0.25,
    bicMinimum = 0.1,
    subsetSize=800,
    fakeCellSD = 0.09,
    uncertaintyCutoff = 0.65,
    summaryFunction=summaryFunction,
    maxClust = 4,
    mergeCutoff = 3,
    IQRCutoff = 0.25,
    medianQuantileCutoff = -1,
    normalCells=control) 
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)
  #to save this data you can use annotateCNV4 as per usual
  final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4)
  saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs.rds"))
  #the other option is to use annotateCNV4B and feed in the normalBarcodes - this will set the "normal" population to 2 -- if the data is noisy it may lead to false positives so use with caution
  #you may also use this version if you have a list of normal barcodes generated elsewhere
  #final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, control, saveOutput=TRUE,outputSuffix = "clean_cnv_b2",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4,minAlteredCellProp = 0.5)
  #PART 2B: smoothing CNV calls with clusters
  #data smoothing: can provide CNV as list or as an input file (CSV)
  #smoothedCNVList<-smoothClusters(scDataSampClusters,inputCNVList = final_cnv_list[[3]],percentPositive = 0.4,removeEmpty = FALSE)
  #PART 3: identify double minutes / amplifications
  #note: this is slow, and may take ~5 minutes
  #if very large dataset, may run on subset of the data to estimate the amplifications in distinct clusters
  #option to compile this code
  #library(compiler)
  #dmRead<-cmpfun(identifyDoubleMinutes)
  #
  #minThreshold is a time-saving option that doesn't call changepoints on any cell with a maximum Z score less than 4 - you can adjust this to adjust sensitivity of double minute calls (note - lower value = slower)
  #dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 4) 
  #write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"samp_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")
  #PART 4: assess putative LOH regions
  #note: this is in beta, interpret results with caution
  #loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 3,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,lossCutoffCells = 200,quantileLimit=0.2,cpgCutoff=100,dummyQuantile=0.6,dummyPercentile=0.4,dummySd=0.1)
  #if(length(loh_regions>0)){
  #write.table(x=loh_regions[[1]],file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"samp_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")}
  #PART 5: quantify cycling cells
  #this uses the signal from a particular chromosome (we use chromosome X as typically not altered in our samples) to identify 2N versus 4N DNA content within cells
  #if there is a known alteration in X in your samples, try using a different chromosome
  #barcodeCycling<-estimateCellCycleFraction(scData,sampName="sample",cutoff=1000)
  #take max as 
  #write.table(barcodeCycling[order(names(barcodeCycling))]==max(barcodeCycling),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cycling_cells.tsv"),sep="\t",quote=FALSE,row.names=TRUE,col.names=FALSE)
  print(paste("Finished sample",sample_name))
}

lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2","RM_3","RM_4",4,10,12),copyscAT_per_sample)

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
  col=colorRamp2(c(0, max(swarbrick_out[,3:ncol(swarbrick_out)])), c("white", "blue")),
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


## Add per cell subtyping to epithelial cells
This is from SC subtyping (the method from the Swarbrick paper). Supplemental Table 4. 
```R
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

dat<-readRDS("phase2.QC.SeuratObject.rds")

#Features determined by EMBO manuscript
module_feats<-list()
module_feats[["Basal_SC"]]=c('EMP1', 'TAGLN', 'TTYH1', 'RTN4', 'TK1', 'BUB3', 'IGLV3.25', 'FAM3C', 'TMEM123', 'KDM5B', 'KRT14', 'ALG3', 'KLK6', 'EEF2', 'NSMCE4A', 'LYST', 'DEDD', 'HLA.DRA', 'PAPOLA', 'SOX4', 'ACTR3B', 'EIF3D', 'CACYBP', 'RARRES1', 'STRA13', 'MFGE8', 'FRZB', 'SDHD', 'UCHL1', 'TMEM176A', 'CAV2', 'MARCO', 'P4HB', 'CHI3L2', 'APOE', 'ATP1B1', 'C6orf15', 'KRT6B', 'TAF1D', 'ACTA2', 'LY6D', 'SAA2', 'CYP27A1', 'DLK1', 'IGKV1.5', 'CENPW', 'RAB18', 'TNFRSF11B', 'VPS28', 'HULC', 'KRT16', 'CDKN2A', 'AHNAK2', 'SEC22B', 'CDC42EP1', 'HMGA1', 'CAV1', 'BAMBI', 'TOMM22', 'ATP6V0E2', 'MTCH2', 'PRSS21', 'HDAC2', 'ZG16B', 'GAL', 'SCGB1D2', 'S100A2', 'GSPT1', 'ARPC1B', 'NIT1', 'NEAT1', 'DSC2', 'RP1.60O19.1', 'MAL2', 'TMEM176B', 'CYP1B1', 'EIF3L', 'FKBP4', 'WFDC2', 'SAA1', 'CXCL17', 'PFDN2', 'UCP2', 'RAB11B', 'FDCSP', 'HLA.DPB1', 'PCSK1N', 'C4orf48', 'CTSC')
module_feats[["Her2E_SC"]]=c('PSMA2', 'PPP1R1B', 'SYNGR2', 'CNPY2', 'LGALS7B', 'CYBA', 'FTH1', 'MSL1', 'IGKV3.15', 'STARD3', 'HPD', 'HMGCS2', 'ID3', 'NDUFB8', 'COTL1', 'AIM1', 'MED24', 'CEACAM6', 'FABP7', 'CRABP2', 'NR4A2', 'COX14', 'ACADM', 'PKM', 'ECH1', 'C17orf89', 'NGRN', 'ATG5', 'SNHG25', 'ETFB', 'EGLN3', 'CSNK2B', 'RHOC', 'PSENEN', 'CDK12', 'ATP5I', 'ENTHD2', 'QRSL1', 'S100A7', 'TPM1', 'ATP5C1', 'HIST1H1E', 'LGALS1', 'GRB7', 'AQP3', 'ALDH2', 'EIF3E', 'ERBB2', 'LCN2', 'SLC38A10', 'TXN', 'DBI', 'RP11.206M11.7', 'TUBB', 'CRYAB', 'CD9', 'PDSS2', 'XIST', 'MED1', 'C6orf203', 'PSMD3', 'TMC5', 'UQCRQ', 'EFHD1', 'BCAM', 'GPX1', 'EPHX1', 'AREG', 'CDK2AP2', 'SPINK8', 'PGAP3', 'NFIC', 'THRSP', 'LDHB', 'MT1X', 'HIST1H4C', 'LRRC26', 'SLC16A3', 'BACE2', 'MIEN1', 'AR', 'CRIP2', 'NME1', 'DEGS2', 'CASC3', 'FOLR1', 'SIVA1', 'SLC25A39', 'IGHG1', 'ORMDL3', 'KRT81', 'SCGB2B2', 'LINC01285', 'CXCL8', 'KRT15', 'RSU1', 'ZFP36L2', 'DKK1', 'TMED10', 'IRX3', 'S100A9', 'YWHAZ')
module_feats[["LumA_SC"]]=c('SH3BGRL', 'HSPB1', 'PHGR1', 'SOX9', 'CEBPD', 'CITED2', 'TM4SF1', 'S100P', 'KCNK6', 'AGR3', 'MPC2', 'CXCL13', 'RNASET2', 'DDIT4', 'SCUBE2', 'KRT8', 'MZT2B', 'IFI6', 'RPS26', 'TAGLN2', 'SPTSSA', 'ZFP36L1', 'MGP', 'KDELR2', 'PPDPF', 'AZGP1', 'AP000769.1', 'MYBPC1', 'S100A1', 'TFPI2', 'JUN', 'SLC25A6', 'HSP90AB1', 'ARF5', 'PMAIP1', 'TNFRSF12A', 'FXYD3', 'RASD1', 'PYCARD', 'PYDC1', 'PHLDA2', 'BZW2', 'HOXA9', 'XBP1', 'AGR2', 'HSP90AA1') 
module_feats[["LumB_SC"]]=c('UGCG', 'ARMT1', 'ISOC1', 'GDF15', 'ZFP36', 'PSMC5', 'DDX5', 'TMEM150C', 'NBEAL1', 'CLEC3A', 'GADD45G', 'MARCKS', 'FHL2', 'CCDC117', 'LY6E', 'GJA1', 'PSAP', 'TAF7', 'PIP', 'HSPA2', 'DSCAM.AS1', 'PSMB7', 'STARD10', 'ATF3', 'WBP11', 'MALAT1', 'C6orf48', 'HLA.DRB1', 'HIST1H2BD', 'CCND1', 'STC2', 'NR4A1', 'NPY1R', 'FOS', 'ZFAND2A', 'CFL1', 'RHOB', 'LMNA', 'SLC40A1', 'CYB5A', 'SRSF5', 'SEC61G', 'CTSD', 'DNAJC12', 'IFITM1', 'MAGED2', 'RBP1', 'TFF1', 'APLP2', 'TFF3', 'TRH', 'NUPR1', 'EMC3', 'TXNIP', 'ARPC4', 'KCNE4', 'ANPEP', 'MGST1', 'TOB1', 'ADIRF', 'TUBA1B', 'MYEOV2', 'MLLT4', 'DHRS2', 'IFITM2')
module_feats[["proliferation_score"]]<-c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")


dat_sub<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial"))
module_scores<-AddModuleScore(dat_sub,features=module_feats,assay="RNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-4,ncol(module_scores@meta.data))]
colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like

dat<-AddMetaData(dat,metadata=module_scores)
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

diag_cols<-c("IDC"="red", "DCIS"="grey")

molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
########################################






harmony_sample_integration<-function(x,outname){
dat<-readRDS(x)
dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_atac",assay.use="ATAC",reduction="cistopic",project.dim=F)
dat<-RunHarmony(dat,group.by.vars="sample",reduction.save="harmony_rna",assay.use="RNA",reduction="pca",project.dim=F)

dat<-RunUMAP(dat,reduction.name="harmonyumap_rna",reduction = "harmony_rna",dims=1:dim(dat@reductions$harmony_rna)[2]) 
dat<-RunUMAP(dat,reduction.name="harmonyumap_atac",reduction = "harmony_atac",dims=1:dim(dat@reductions$harmony_atac)[2]) 

alpha_val=0.33

plt1<-DimPlot(dat,reduction="harmonyumap_rna",group.by="sample")#+scale_fill_manual(samp_cols)
plt2<-DimPlot(dat,reduction="harmonyumap_rna",group.by="predicted.id")#,cols=alpha(type_cols,alpha_val))
plt3<-DimPlot(dat,reduction="harmonyumap_rna",group.by="diagnosis")#,cols=alpha(diag_cols,alpha_val))
plt4<-DimPlot(dat,reduction="harmonyumap_rna",group.by="molecular_type")#,cols=alpha(molecular_type_cols,alpha_val))

plt5<-DimPlot(dat,reduction="harmonyumap_atac",group.by="sample")#+scale_fill_manual(samp_cols)
plt6<-DimPlot(dat,reduction="harmonyumap_atac",group.by="predicted.id")#,cols=alpha(type_cols,alpha_val))
plt7<-DimPlot(dat,reduction="harmonyumap_atac",group.by="diagnosis")#,cols=alpha(diag_cols,alpha_val))
plt8<-DimPlot(dat,reduction="harmonyumap_atac",group.by="molecular_type")#,cols=alpha(molecular_type_cols,alpha_val))
ggsave((plt1|plt2|plt3|plt4)/(plt5|plt6|plt7|plt8),file=paste0(outname,".Harmonyumap.pdf"),width=20)
system(paste0("slack -F ",outname,".Harmonyumap.pdf ryan_todo"))

saveRDS(dat,file=x)
}

harmony_sample_integration(x="stromal.SeuratObject.rds",outname="stromal")
harmony_sample_integration(x="immune.SeuratObject.rds",outname="immune")
harmony_sample_integration(x="epithelial.SeuratObject.rds",outname="epithelial")

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

## 3D Plot in Blender

Open blender, go to python console (Shift+F4) then copy and paste the code below

```python
#1. import modules
import bpy
import math
import time
import bmesh

#set up variables
file_in="//rdsdcw.ohsu.edu/cedarX/Projects/breast_cancer_multiome/Experiment_2_Phase1/multiome_tumor.tsv" 
file_out="//rdsdcw.ohsu.edu/cedarX/Projects/breast_cancer_multiome/Experiment_2_Phase1/multiome_tumor.blend"


#Read in file and store it in memory (this doesn't take up much memory)
file_xyz=open(file_in,"r") #change path to whatever filepath you want. I got my computer refurbished and it was named Chad. I swear it wasn't me.
tabraw=file_xyz.readlines()[1:]
data_count=len(tabraw)
file_xyz.close()

#initialize an object, a sphere, for our data points.
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.05,segments=64, ring_count=32) #higher segments and ring_counts will make a smoother sphere, but I dont think its necessary
obj=bpy.context.active_object #select the sphere we just made

#set up a master shader material
mat = bpy.data.materials.new(name='mymat')
mat.use_nodes = True #use node trees, these can be seen by switching a panel to the shader editor if you want. It will look like the above shader, just not nicely placed.
mat_nodes = mat.node_tree.nodes
mat_links = mat.node_tree.links
mat = bpy.data.materials['mymat'] #Get the material you want 
node_to_delete =  mat.node_tree.nodes['Principled BSDF'] #Get the node in its node tree (replace the name below)
mat.node_tree.nodes.remove( node_to_delete ) #Remove it
#add all the nodes, using col_node as variable of each node as it is being made. then using that to modify default value fields
col_node=mat_nodes.new('ShaderNodeRGB')
col_node=mat_nodes.new('ShaderNodeFresnel')
bpy.data.materials["mymat"].node_tree.nodes['Fresnel'].inputs[0].default_value = 1.33

col_node=mat_nodes.new('ShaderNodeHueSaturation')
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[0].default_value = 1
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[1].default_value = 0.7
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[2].default_value = 2
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[3].default_value = 0

col_node=mat_nodes.new('ShaderNodeMath')
bpy.data.materials["mymat"].node_tree.nodes["Math"].operation = 'MULTIPLY'

col_node=mat_nodes.new('ShaderNodeBsdfRefraction')
bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].inputs[1].default_value = 1

col_node=mat_nodes.new('ShaderNodeBsdfGlossy')
bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].inputs[1].default_value = 1

col_node=mat_nodes.new('ShaderNodeHueSaturation')
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[0].default_value = 1
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[1].default_value = 0.4
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[2].default_value = 2

col_node=mat_nodes.new('ShaderNodeMixShader')
col_node=mat_nodes.new('ShaderNodeVolumeAbsorption')
bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].inputs[1].default_value = 0.3

col_node=mat_nodes.new('ShaderNodeBsdfTranslucent')
col_node=mat_nodes.new('ShaderNodeLightPath')
col_node=mat_nodes.new('ShaderNodeMixShader')

#build node tree links (going from left most inputs)
#sorry this is a monstrosity
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[4])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[4])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Fresnel'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Math"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Hue Saturation Value'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].inputs[0])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Hue Saturation Value'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Math"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[0])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[2])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Translucent BSDF"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Material Output"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Translucent BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[2])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Light Path"].outputs[1], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Material Output"].inputs[0])


#set up render engine and scene
bpy.context.scene.render.engine="CYCLES" #set render engine to CYCLES
bpy.data.scenes["Scene"].cycles.denoiser="NLM" #set denoiser for render
bpy.data.scenes["Scene"].cycles.samples=512 #this is a whole lotta sampling
bpy.context.scene.render.image_settings.color_depth = '16' #more color channels!
bpy.context.scene.render.resolution_x = 3840 #up the resolution
bpy.context.scene.render.resolution_y = 2160
bpy.data.objects["Sphere"].hide_render = True # hide sphere in render
bpy.data.objects["Sphere"].hide_viewport=True
bpy.data.lights["Light"].energy = 100000 # increase light wattage
bpy.data.lights["Light"].shadow_soft_size= 1
bpy.data.objects["Light"].location=(5,-5,10) #location and rotation i deteremined manually and just set up here for convenience

#set up stage by cutting up the default cube vertices and smoothing it
obj_cube=bpy.data.objects["Cube"]
obj_cube.scale=(30,30,30) #scale up the cube
#this is to cut out a vertex to make an open box
bpy.context.view_layer.objects.active = obj_cube
bpy.ops.object.mode_set(mode='EDIT')
bpy.ops.mesh.select_mode(type="VERT")  # Switch to edge select mode
bm = bmesh.from_edit_mesh(obj_cube.data)  # Create bmesh object for easy mesh evaluation
bm.verts.ensure_lookup_table()
bm.verts.remove(bm.verts[2]) # Write the mesh back
bmesh.update_edit_mesh(obj_cube.data)  # Update the mesh in edit mode
bpy.ops.object.mode_set(mode='OBJECT') #switch back to object mode when done
bpy.ops.object.modifier_add(type='SUBSURF') #make it smooth
bpy.data.objects["Cube"].modifiers["Subdivision"].render_levels=6
bpy.data.objects["Cube"].location=(-4,4.3,17.725) #change the location for more dramatic shadows

#move the camera and rotate
bpy.data.objects["Camera"].location=(34.61997604370117, -40.53969955444336, 25.66326904296875)
bpy.data.objects["Camera"].rotation_euler=(1.1093189716339111, 0.0, 0.8149281740188599)

#finally ready to start reading in our data
scene=bpy.context.scene

#set up a material per hex color, name as annotation
#this is looping through the file, grabbing the unique clusters and there color codes, then making a dictionary for look up later
start = time.time()
annot={}
for line in tabraw[1:]:
  line=line.replace('\n','')
  l=line.split('\t')
  if l[0] not in annot:
    hexcode=l[5].lstrip("#")
    rgb=[int(hexcode[i:i+2], 16) for i in (0, 2, 4)]
    r=float(rgb[0])/255 #color of spheres, blender uses 0-1 scale
    g=float(rgb[1])/255
    b=float(rgb[2])/255
    clust=str(l[0])
    annot[clust]=[r,g,b]

end = time.time()
print(end - start)

#make a custom material shader for each annotation (just changing color)
#this copies the material shader we set up earlier, and then changes the input color
master_mat=source_mat = bpy.data.materials["mymat"]
for i in annot.keys():
  copied_mat = master_mat.copy()
  copied_mat.name=i
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[0]=annot[i][0]
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[1]=annot[i][1]
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[2]=annot[i][2]


#make a custom collection for each annotation. this makes a "master sphere" to link for each cluster also
for i in annot.keys():
  collection = bpy.data.collections.new(i) #make new collection
  bpy.context.scene.collection.children.link(collection) #link new collection
  mat = bpy.data.materials.get(i) #set material properties of collection
  name=str(i)+"_master" #make name of master sphere
  new_obj = bpy.data.objects.new(name, scene.objects.get("Sphere").data) #make a new copy
  new_obj.data = scene.objects.get("Sphere").data.copy()
  bpy.data.collections[i].objects.link(new_obj) #link new object to collection
  new_obj.data.materials.append(mat) #add material
  bpy.data.objects[name].hide_render = True # hide masters
  bpy.data.objects[name].hide_viewport=True

#make a dictionary look up for copying master spheres
master_sphere={}
for i in annot.keys():
  master_sphere[i]=scene.objects.get(i+"_master").data

#define a nice function to copy data points and link them to the master spheres. also places the copies into nice cluster named collections for easier navigation.
def add_data_point(input_dat):
    line=input_dat
    line=line.replace('\n','')
    l=line.split('\t')
    #print(line)
    x=float(l[2]) #location of spheres
    y=float(l[3])
    z=float(l[4])
    name=str(l[1])
    clust=str(l[0])
    my_new_obj = bpy.data.objects.new(name,master_sphere[clust])
    my_new_obj.location = (x,y,z)       
    my_new_obj.hide_viewport=False
    my_new_obj.hide_render=False
    bpy.data.collections[clust].objects.link(my_new_obj)

n=1000
in_list = [tabraw[i * n:(i + 1) * n] for i in range((len(tabraw) + n + 1) // n )] 
for in_dat_list in in_list:
  start = time.time()
  out=[add_data_point(in_dat) for in_dat in in_dat_list] 
  end = time.time()
  print(end - start)


#some last minute tweaks, here are some convenience functions if you want to change things. I also encourage you to play around with lighting and camera positioning to get some interesting views of your data.
#to adjust size of points
for clust in annot.keys():
  for i in bpy.data.collections[clust].objects:
    i.scale=(0.8,0.8,0.8)


#to adjust alpha value and translucence of material properties
for clust in annot.keys():
  bpy.data.materials[clust].node_tree.nodes["RGB"].outputs[0].default_value[3] = 0.3
  bpy.data.materials[clust].node_tree.nodes["Volume Absorption"].inputs[1].default_value = 0.1


bpy.ops.wm.save_as_mainfile(filepath=file_out) #save blender file

bpy.context.scene.render.filepath = 'cortex.test.png'
bpy.ops.render.render(write_still=True) #render and save file

```

Load final blend file onto exacloud and render

```bash
#run on exacloud
blender -b 3dplot_multiome.blend -o //render_test/ -E CYCLES -s 20 -e 1000 -t 40 -a -x 1 -F PNG &
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