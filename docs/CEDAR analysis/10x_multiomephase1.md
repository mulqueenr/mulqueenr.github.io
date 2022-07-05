---
title: 10X Multiome Tumor Phase 1
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_phase1/
category: CEDAR
---

Multiome processing for 10X multiome data on Primary Tumors (Phase 1)

## File Location
```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1
sftp mulqueen@nix.ohsu.edu
#enter password
get -r /data/EXP220411HM
get -r /data/EXP220412HM
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
*,sample_1,SI-NA-A3
*,sample_3,SI-NA-B3
*,sample_4,SI-NA-C3
*,sample_5,SI-NA-D3
*,sample_6,SI-NA-E3
*,sample_7,SI-NA-F3
*,sample_8,SI-NA-G3
*,sample_9,SI-NA-H3
*,sample_10,SI-NA-A4
*,sample_11,SI-NA-B4
*,sample_12,SI-NA-C4""" > multiome_atac_phase1.csv

echo """Lane,Sample,Index
*,sample_1,SI-TT-A5
*,sample_3,SI-TT-B5
*,sample_4,SI-TT-C5
*,sample_5,SI-TT-D5
*,sample_6,SI-TT-E5
*,sample_7,SI-TT-F5
*,sample_8,SI-TT-G5
*,sample_9,SI-TT-H5
*,sample_10,SI-TT-A6
*,sample_11,SI-TT-B6
*,sample_12,SI-TT-C6""" > multiome_rna_phase1.csv

```

Run Cellranger-arc
```bash
#For some reason the core did not send me an RTAComplete.txt file, so just fudging it.
echo "Yep, it's done." > /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw/RTAComplete.txt
echo "Yep, it's done." > /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw/RTAComplete.txt
#add bcl2fastq to path (I'm just using an old one from the adey lab)
export PATH=/home/groups/oroaklab/src/bcl2fastq/bcl2fastq-2.19.0/:$PATH

cellranger-arc mkfastq --id=phase_1_atac \
                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw \
                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/multiome_atac_phase1.csv \
                     --localcores=20 \
                     --localmem=80

cellranger-arc mkfastq --id=phase_1_rna \
                     --run=/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw \
                     --csv=/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/multiome_rna_phase1.csv \
                      --localcores=20 \
                     --localmem=80
```
## Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc.

### RM Libraries

```bash
for i in 1 3 4 5 6 7 8 9 10 11 12; do
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase_1_atac/outs/fastq_path/HCVKLDRX2/sample_"""${i}""",sample_"""${i}""",Chromatin Accessibility
/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase_1_rna/outs/fastq_path/HCVLCDRX2,sample_"""${i}""",Gene Expression""" > sample_${i}.csv ; done
```

## Run CellRanger-ARC
```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1 #using this as analysis directory
```
### Slurm commands
Test run

```bash
 cellranger-arc testrun --id=tiny2
```
Run Cellranger per sample

```bash          
for i in sample_1.csv sample_3.csv sample_4.csv sample_5.csv sample_6.csv sample_7.csv sample_8.csv sample_9.csv ; do
  outname=${i::-4};
  cellranger-arc count --id=${outname} \
   --reference=/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
   --libraries=/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/${i} \
   --localcores=30 \
   --localmem=90 ; done &
   
#check web summaries
cd /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/
for i in 1 3 4 5 6 7 8 9 10 11 12; do
  cp ./sample_$i/outs/web_summary.html ./sample_$i/outs/$i.web_summary.html
  slack -F ./sample_$i/outs/$i.web_summary.html ryan_todo; done 
```

## Phase 1 Primary Tissue Analysis

### Seurat Object Generation for Tumors
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

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# set up sample loop to load the RNA and ATAC data, save to seurat object
setupseurat<-function(i){
  setwd(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/",i,"/outs"))
  counts <- Read10X_h5("filtered_feature_bc_matrix.h5") #count data
  fragpath <- "atac_fragments.tsv.gz" #atac fragments
  metadata_cellranger<-read.csv("per_barcode_metrics.csv") #metadata
  row.names(metadata_cellranger)<-metadata_cellranger$barcode

  # create a Seurat object containing the RNA adata
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

  #QC cells
  DefaultAssay(dat) <- "ATAC"

  dat <- NucleosomeSignal(dat)
  dat <- TSSEnrichment(dat)
  dat<-AddMetaData(dat,metadata=metadata_cellranger)

  plt<-VlnPlot(
    object = dat,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
  ggsave(plt,file=paste0(i,".qc.pdf"))
  system(paste0("slack -F ",i,".qc.pdf ryan_todo"))
  saveRDS(dat,file=paste0(i,".SeuratObject.2.rds"))
}

#generate all seurat objects
lapply(c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12))),setupseurat)

#combine seurat objects and add metadata column for sample
samp_1<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs/sample_1.SeuratObject.rds"); samp_1$sample<-"sample_1"
samp_3<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_3/outs/sample_3.SeuratObject.rds"); samp_3$sample<-"sample_3"
samp_4<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_4/outs/sample_4.SeuratObject.rds"); samp_4$sample<-"sample_4"
samp_5<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_5/outs/sample_5.SeuratObject.rds"); samp_5$sample<-"sample_5"
samp_6<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_6/outs/sample_6.SeuratObject.rds"); samp_6$sample<-"sample_6"
samp_7<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_7/outs/sample_7.SeuratObject.rds"); samp_7$sample<-"sample_7"
samp_8<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_8/outs/sample_8.SeuratObject.rds"); samp_8$sample<-"sample_8"
samp_9<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_9/outs/sample_9.SeuratObject.rds"); samp_9$sample<-"sample_9"
samp_10<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_10/outs/sample_10.SeuratObject.rds"); samp_10$sample<-"sample_10"
samp_11<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_11/outs/sample_11.SeuratObject.rds"); samp_11$sample<-"sample_11"
samp_12<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_12/outs/sample_12.SeuratObject.rds"); samp_12$sample<-"sample_12"

dat <- merge(samp_1, y = c(samp_3,samp_4,samp_5,samp_6,samp_7,samp_8,samp_9,samp_10,samp_11,samp_12), add.cell.ids = c("samp_1","samp_3","samp_4","samp_5","samp_6","samp_7","samp_8","samp_9","samp_10","samp_11","samp_12"), project = "primary_phase1")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")

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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

dat<-readRDS("phase1.SeuratObject.rds")
dat
table(dat$sample)
# sample_1 sample_10 sample_11 sample_12  sample_3  sample_4  sample_5  sample_6
#     3523     20000      5575     14071      7698     15253      1453       666
# sample_7  sample_8  sample_9
#     2656       724      1451

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# filter out low quality cells
#dat <- subset(
#  x = dat,
#  subset = nCount_ATAC < 100000 &
#    nCount_RNA < 25000 &
#    nCount_ATAC > 500 &
#    nCount_RNA > 500 &
#    nucleosome_signal < 2 &
#    TSS.enrichment > 1
#)
#dat
#table(dat$sample)

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

#set up colors for samples
my_cols = brewer.pal(11,"Spectral")
alpha_val=0.33
#RNA Processing
DefaultAssay(dat) <- "RNA"
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
p1<-DimPlot(dat,reduction="rna_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("RNA UMAP")

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
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("ATAC UMAP")


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
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("Multimodal UMAP")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file="phase1_multimodal.umap.pdf")
system("slack -F phase1_multimodal.umap.pdf ryan_todo")
saveRDS(dat,file="phase1.SeuratObject.rds")

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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

peaks <- readRDS(file="combined.peakset.rds")


single_sample_cluster<-function(x){
dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds"))
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
DefaultAssay(dat) <- "RNA"
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
p1<-DimPlot(dat,reduction="rna_umap",cols=alpha(my_cols,alpha_val))+ggtitle("RNA UMAP")

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
p2<-DimPlot(dat,reduction="atac_umap",cols=alpha(my_cols,alpha_val))+ggtitle("ATAC UMAP")


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
p3<-DimPlot(dat,reduction="multimodal_umap",cols=alpha(my_cols,alpha_val))+ggtitle("Multimodal UMAP")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,"_multimodal.umap.pdf"))
system(paste0("slack -F /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,"_multimodal.umap.pdf ryan_todo"))
saveRDS(dat,file=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds"))

}

lapply(c(1,3,4,5,6,7,8,9,10,11,12),single_sample_cluster)


```


## Run cisTopic for ATAC Dimensionality Reduction
Cistopic on all samples
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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

atac_sub<-readRDS("phase1.SeuratObject.rds")
wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")
outname<-paste0("phase1")
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
saveRDS(atac_sub,"phase1.SeuratObject.rds")
```


Cistopic Per sample
```R
#######CISTOPIC PROCESSING PER CELL LINE#################
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")
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
myargs = commandArgs(trailingOnly=TRUE)
input_sample=myargs[1]

single_sample_cistopic_generation<-function(x){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds"))
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
  print("Cistopic Loading into Seurat")
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(wd,"/",outname,".umap.pdf"),width=10)
  print(plt1)
  print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.2.rds"))
  }

single_sample_cistopic_generation(x=input_sample)

```


Helper function to load already generated models to load into Seurat Object

```R
#######CISTOPIC LOADING INTO SEURAT OBJECTS PER SAMPLE#####
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")
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


single_sample_cistopic_loading<-function(x){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds"))
  sub_cistopic_models<-readRDS(file=paste0(wd,"/",outname,".CisTopicObject.Rds"))

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
  pdf(paste0(wd,"/",outname,".umap.pdf"),width=10)
  print(plt1)
  #print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.2.rds"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12),single_sample_cistopic_loading)

```
### Batch script for per sample cistopic
Using the R script above. Saved as single_sample_cistopic.Rscript, where trailing argument is the sample number (1-12)

```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node per task
#SBATCH --array=1,3,4,5,6,7,8,9,10,11,12 #set up array of sample names
#SBATCH --tasks-per-node=20 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request 10 gigabyte per cpu (200GB total)
#SBATCH --time=3:00:00 ## ask for 3 hour on the node
#SBATCH --

srun Rscript /home/groups/CEDAR/mulqueen/src/single_sample_cistopic.Rscript $SLURM_ARRAY_TASK_ID

#batch script kept failing so i just ran it as a loop
for i in 1 3 4 5 6 7 8 9 10 11 12;
do Rscript /home/groups/CEDAR/mulqueen/src/single_sample_cistopic.Rscript $i ; done &

```

### Using Transfer Anchors for Cell identification.

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

### Using EMBO paper for label transfer as well. 
Actually waiting on metadata from authors for this one for now.
https://doi.org/10.15252/embj.2020107333
Full code here: https://www.nature.com/articles/s41597-022-01236-2 data available here https://doi.org/10.6084/m9.figshare.17058077

Download data from GEO FTP server

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
wget https://figshare.com/ndownloader/articles/17058077/versions/1
unzip 1
```
Make Seurat Object with Metadata

```R
library(Seurat)

setwd("/home/groups/CEDAR/mulqueen/ref/embo")
counts<-ReadMtx(mtx="GSM4909253_N-PM0092-Total-matrix.mtx.gz",cells="GSM4909253_N-PM0092-Total-barcodes.tsv.gz",features="GSE161529_features.tsv.gz",feature.column=1) #sparse matrix of counts
metadata<-read.csv("metadata.csv") #metadata
row.names(metadata)<-metadata$X
# create a Seurat object containing the RNA adata
x<- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)
swarbrick<-AddMetaData(swarbrick,metadata=metadata)
saveRDS(swarbrick,"/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
```

Transfer Swarbrick cell type info for now


```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

dat<-readRDS("phase1.SeuratObject.rds")
DefaultAssay(dat)<-"RNA"
dat<-NormalizeData(dat)
dat<-FindVariableFeatures(dat)
dat<-ScaleData(dat)
#Using Label transfer to label cell types by Swarbrick paper
#seurat object made by AD
swarbrick<-readRDS("/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
swarbrick<-NormalizeData(swarbrick)
swarbrick<-FindVariableFeatures(swarbrick)
swarbrick<-ScaleData(swarbrick)


transfer.anchors <- FindTransferAnchors(
  reference = swarbrick,
  reference.assay="RNA",
  query = dat,
  query.assay="RNA",
  verbose=T
)

predictions<- TransferData(
  anchorset = transfer.anchors,
  refdata = swarbrick$celltype_major,
)

dat<-AddMetaData(dat,metadata=predictions)
saveRDS(dat,file="phase1.SeuratObject.rds")


plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

plt<-(plt2|plt3)/plt1
ggsave(plt,file="phase1.predictions.umap.pdf",width=20,height=30,limitsize=F)
system("slack -F phase1.predictions.umap.pdf ryan_todo")


##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")
  dat<-readRDS(file_in)
  dat$sample<-c(outname)
  DefaultAssay(dat)<-"RNA"
  dat<-NormalizeData(dat)
  dat<-FindVariableFeatures(dat)
  dat<-ScaleData(dat)
  saveRDS(dat,file=file_in)

  transfer.anchors <- FindTransferAnchors(
    reference = swarbrick,
    reference.assay="RNA",
    query = dat,
    query.assay="RNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = swarbrick$celltype_major,
  )

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)

  plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
  plot_name<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=plot_name,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",plot_name," ryan_todo"))

  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12),single_sample_label_transfer)

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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

#from tumor sample information
meta_data_in<-as.data.frame(cbind("sample"=c(paste0("sample_",seq(1,12))),
  "diagnosis"= c("DCIS", "IDC", "IDC", "IDC", "IDC", "IDC", "DCIS", "IDC", "IDC", "IDC", "IDC", "IDC"),
  "molecular_type"=c("DCIS", "er+_pr+_her2-", "er+_pr-_her2-", "er+_pr-_her2-", "er+_pr+_her2-", "er+_pr+_her2-", "DCIS", "er+_pr+_her2-", "er+_pr+_her2-", "er+_pr-_her2-", "er+_pr-_her2-", "er+_pr+_her2-")))

sample_metadata<-function(dat){
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
dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds"
sample_metadata(dat)

#add metadata to each sample
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")))
lapply(sample_in,sample_metadata)


```

### Plot multimodal dimensionality reduction for cistopic embedding

```R
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")
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



plot_reductions<-function(dat){
  #dat is full path to seurat object
  dat_file_path=dat
  file_in=basename(dat)
  dir_in=dirname(dat)
  dat<-readRDS(dat) #read in as seurat object
  outname=paste0(dir_in,"/",file_in,".")

  #set up colors for samples
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

  diag_cols<-c("IDC"="red", "DCIS"="grey")

  molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
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
    assay = "RNA",
    verbose = TRUE
  )

  ###########plot cistopic umap too
  p5<-DimPlot(dat2,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")+theme(legend.position="none")
  p6<-DimPlot(dat,reduction="cistopic_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Cistopic UMAP")+theme(legend.position="none")

  #Finally Plot results
  plt<-(p1 | p2)/(p3)/(p6 | p5)
  ggsave(plt,file=paste0(outname,"multimodal.umap.pdf"),width=10,height=20)
  system(paste0("slack -F ",outname,"multimodal.umap.pdf ryan_todo"))

  ggsave(p1,file=paste0(outname,"multimodal.umap.rna.pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,"multimodal.umap.rna.pdf ryan_todo"))

  ggsave(p2,file=paste0(outname,"multimodal.umap.atac.pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,"multimodal.umap.atac.pdf ryan_todo"))

  saveRDS(dat,file=dat_file_path)

  plt_cell_count<-dat@meta.data[,c("sample","predicted.id","diagnosis","molecular_type")]
  return(plt_cell_count)
}


#run umap projections of merged samples
plot_reductions(dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")

#run umap projections of all dim reductions per sample
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")))

out<-lapply(sample_in,function(x) plot_reductions(dat=x))
saveRDS(out,"sample_swarbrick_celltype_assignment.rds") #save nested list of cell type assignment

#plot output of celltype count per sample
out<-readRDS("sample_swarbrick_celltype_assignment.rds")
out<-do.call("rbind",out)
colnames(out)<-c("sample","celltype","diagnosis","molecular_type") #rename just for simplicity
  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

  diag_cols<-c("IDC"="red", "DCIS"="grey")

  molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
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

### Check Cell Type Assignment By Prediction Grouping
https://ars.els-cdn.com/content/image/1-s2.0-S1097276521007954-gr5.jpg

###NOT RUN YET####

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

dat<-readRDS(file="phase1.SeuratObject.rds")

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")))

#Basing identification off of https://ars.els-cdn.com/content/image/1-s2.0-S1097276521007954-gr5.jpg
#and https://www.embopress.org/doi/full/10.15252/embj.2020107333
geneset<-list()
geneset[["luminal_Hr"]]<-c("ANKRD30A","ERBB4","SYTL2","INPP4B","AFF3")
geneset[["luminal_secretory"]]<-c("IGF2BP2","ALDH1A3","COBL","PIGR","CCL28")
geneset[["basal"]]<-c("NRG1","KRT14","ACTA2","PTPRT","MYLK","KRT5","SNAI2","NOTCH4","DKK3")
geneset[["t_cells"]]<-c("PTPRC","PARP8","FYN","STAT4","AOAH")
geneset[["b_cells"]]<-c("MZB1","SSR4","HERPUD1","DERL3") #"ENAM" excluded
geneset[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74","HLA-DRB1","HLA-DPB1")
geneset[["endothelial"]]<-c("MCTP1","ADAMTS9","ZNF385D","ADGRL4","SELE") #ELTD1 is ADGRL4
geneset[["pericyte"]]<-c("C11orf95","RGS6","PRKG1","IGFBP5") #MT1A excluded
geneset[["fibroblast"]]<-c("DCN","APOD","HPSE2","CFD","LAMA2")
geneset[["epithelial"]]<-c("EPCAM","ITGA6","KRT5")
geneset[["mature_luminal"]]<-c("ESR1","PGR","FOXA1")
geneset[["luminal_progenitor"]]<-c("TNFRSF11A","KIT","SOX10")


cov_plots<-function(dat=dat,gene_name){
    gene_name<-unname(gene_name)
    plt_cov <- CoveragePlot(
      object = dat,
      region = gene_name,
      features = gene_name,
      assay="peaks",
      expression.assay = "SCT",
      extend.upstream = 5000,
      extend.downstream = 5000)
    plt_feat <- FeaturePlot(
      object = dat,
      features = gene_name,
      raster=T,
      reduction="multimodal_umap",
      order=T)
    return((plt_feat|plt_cov)+ggtitle(gene_name))
  }


plot_genes<-function(dat){
  #dat is full path to seurat object
  dat_file_path=dat
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-3)
  dat<-readRDS(dat)

  ###LINKING PEAKS TO GENES###
  DefaultAssay(dat) <- "peaks"

  # first compute the GC content for each peak
  dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

  # link peaks to genes
  #Identify cell types in data

  #filter geneset to those in data set
  unlist(geneset)[which(!(unlist(geneset) %in% row.names(dat@assays$SCT@data)))]

  #dat <- LinkPeaks(
  #  object = dat,
  #  peak.assay = "peaks",
  #  expression.assay = "SCT",
  #  genes.use = unlist(geneset)
  #)

  Idents(dat)<-dat$predicted.id

  DefaultAssay(dat) <- "SCT"
  for (i in unique(names(geneset))){
    plt_list<-lapply(unlist(geneset[i]), function(x) cov_plots(dat=dat,gene_name=x))
    plt<-patchwork::wrap_plots(plt_list, ncol = 1)
    ggsave(plt,file=paste0(outname,i,".featureplots.pdf"),height=4*length(plt_list),width=10,limitsize=F)
    system(paste0("slack -F ",outname,i,".featureplots.pdf ryan_todo"))
  }

  saveRDS(dat,file=dat_file_path)
}

plot_genes(dat="phase1.SeuratObject.rds")
plot_genes(dat=sample_in[7])
lapply(sample_in,function(x) plot_genes(dat=x))

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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")


sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")))


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

####RUNNING INFERCNV#####
infercnv_per_sample<-function(dat){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  dat_file_path=dat
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-3)
  dat<-readRDS(dat)
  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","Cancer Epithelial","Normal Epithelial","T-cells")) #just look at epithelial

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
  write.table(counts,file=paste0(outname,"_inferCNV.counts.txt"),sep="\t",col.names=T,row.names=T,quote=F)
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["cnv_ref"]))
  write.table(cell_annotation,file=paste0(outname,"_inferCNV.annotation.txt"),sep="\t",col.names=F,row.names=F,quote=F)

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(outname,"_inferCNV.counts.txt"),
                                      annotations_file=paste0(outname,"_inferCNV.annotation.txt"),
                                      delim="\t",
                                      gene_order_file="inferCNV.gene_order.txt",
                                      ref_group_names="TRUE")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0(outname,"_inferCNV"), 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               resume_mode=F,
                               num_threads=20)
  saveRDS(infercnv_obj,paste0(outname,"_inferCNV","/",basename(outname),"inferCNV.Rds"))
  system(paste0("slack -F ",outname,"_inferCNV","/","infercnv.png"," -T ","\"",outname,"\"" ," ryan_todo") )
  #infercnv_obj<-readRDS(paste0(outname,"_inferCNV","/",basename(outname),"inferCNV.Rds"))

}
lapply(sample_in,function(x) infercnv_per_sample(dat=x))


```

### Plot inferCNV output

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

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.2.rds")))


  ###########Color Schema#################
  type_cols<-c(
  #epithelial
  "Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
  "B-cells" ="#089099", "T-cells" ="#003147", #other
  "CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

  ref_cols<-c(
  "TRUE"="grey",
  "FALSE"="red")
  #########################################

####RUNNING INFERCNV#####
infercnv_per_sample_plot<-function(dat){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  dat_file_path=dat
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-3)
  dat<-readRDS(dat)
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  DefaultAssay(dat)<-"RNA"
  infercnv_obj<-readRDS(paste0(outname,"_inferCNV","/",basename(outname),"inferCNV.Rds"))
  outname<-substr(dat_file_path,1,nchar(dat_file_path)-17) #redo outname to remove "seuratobject" from filename
  cnv<-t(infercnv_obj@expr.data)
  cnv<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]
  col_fun = colorRamp2(c(min(unlist(cnv)), median(unlist(cnv)), max(unlist(cnv))), c("blue", "white", "red"))

  dist_method="euclidean"
  dist_x<-philentropy::distance(cnv,method=dist_method,as.dist.obj=T,use.row.names=T)
  dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
  k_search<-find_k(dend,krange=2:10) #search for optimal K from 2-10
  k_clus_number<-k_search$nc
  k_clus_id<-k_search$pamobject$clustering
  dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters
  saveRDS(dend,file=paste0(outname,".inferCNV.dend.Rds")) #save dendrogram

  #set up heatmap annotation
  met<-as.data.frame(dat@meta.data)
  met<-met[row.names(met) %in% row.names(cnv),]
  read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
    max(met$gex_exonic_umis+met$gex_intronic_umis)), 
    c("white","black"))
  ha = HeatmapAnnotation(which="row",
    cell_type=met$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met$gex_exonic_umis+met$gex_intronic_umis,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col))
  plt1<-Heatmap(cnv,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      cluster_rows=dend,
      left_annotation=ha,
      column_split=infercnv_obj@gene_order$chr)
    pdf(paste0(outname,".inferCNV.heatmap.pdf"),width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0(outname,".inferCNV.heatmap.pdf")," ryan_todo"))
}

lapply(sample_in,function(x) infercnv_per_sample_plot(dat=x))

```


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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

###Trying CASPER for CNV Profiling

library(CaSpER) 

casper_per_sample<-function(dat){
  file_in=basename(dat)
  sample_name=substr(file_in,1,nchar(file_in)-17)
  dir_in=dirname(dat)
  system(paste0("mkdir ",dir_in,"_casper"))
  bam_location<-paste0(dir_in,"/gex_possorted_bam.bam")
  BAFExtract_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract"
  hg38_list_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list" #downloaded from https://github.com/akdess/BAFExtract
  hg38_folder_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/"
  baf_sample_directory<-paste0(dir_in,"_casper")
  dat<-readRDS(dat)
  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Endothelial","B-cells","Cancer Epithelial","Normal Epithelial")) #just look at epithelial

  dat<-subset(dat,predicted.id %in% c("Endothelial","B-cells","Cancer Epithelial","Normal Epithelial")) #just look at epithelial
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

  saveRDS(object,paste0(dir_in,"_casper/",sample_name,".initialobj.rds"))

  pdf(paste0(dir_in,"_casper/",sample_name,".Distrubution.pdf"))
  plot(density(as.vector(object@control.normalized[[3]])))
  plot(density(log2(object@control.normalized.noiseRemoved[[3]]+1)))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".Distrubution.pdf"), " ryan_todo"))

  object<-readRDS(paste0(dir_in,"_casper/",sample_name,".initialobj.rds"))
  ## runCaSpER
  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

  ## summarize large scale events 
  finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
  final.obj <- final.objects[[9]]
  saveRDS(final.obj,paste0(dir_in,"_casper/",sample_name,".finalobj.rds"))
  saveRDS(finalChrMat,paste0(dir_in,"_casper/",sample_name,".finalchrmat.rds"))
  final.obj<-readRDS(paste0(dir_in,"_casper/",sample_name,".finalobj.rds"))

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
        show_rownames = F, filename = paste0(dir_in,"_casper/",sample_name,".heatmap.pdf"))
  system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".heatmap.pdf"), " ryan_todo"))


  pdf(paste0(dir_in,"_casper/",sample_name,".LargeCNV.heatmap.pdf"))
  Heatmap(finalChrMat,
    cluster_columns=F,
    show_row_names=F,
    column_names_rot = 45,
    column_split=sort(rep(c(1:22),2)),
    column_names_gp = gpar(fontsize = 5))
  dev.off()
  system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".LargeCNV.heatmap.pdf")," ryan_todo"))


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
  ggsave(p,file=paste0(dir_in,"_casper/",sample_name,".final_plot.pdf"))
  system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".final_plot.pdf"), " ryan_todo"))

}

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

lapply(sample_in,function(x) casper_per_sample(dat=x))

```
### Plotting of CaSpER Output
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
library(dendextend)
library(CaSpER)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

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
  #########################################

####RUNNING INFERCNV#####
casper_per_sample_plot<-function(dat){
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  #dat is full path to seurat object
  print(paste("Reading in:",dat))
  file_in=basename(dat)
  sample_name=substr(file_in,1,nchar(file_in)-17)
  dir_in=dirname(dat)
  dat<-readRDS(dat)
  dat$cnv_ref<-"FALSE"
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  DefaultAssay(dat)<-"RNA"
  print(paste("Reading in:",paste0(dir_in,"_casper/",sample_name,".finalobj.rds")))
  final.obj<-readRDS(paste0(dir_in,"_casper/",sample_name,".finalobj.rds"))

  cnv<-t(final.obj@control.normalized.noiseRemoved[[3]])
  cnv<-cnv[row.names(cnv) %in% row.names(dat@meta.data[dat@meta.data$cnv_ref=="FALSE",]),]
  col_fun = colorRamp2(c(min(unlist(cnv)), median(unlist(cnv)), max(unlist(cnv))), c("blue", "white", "red"))

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
  met<-met[row.names(met) %in% row.names(cnv),]
  read_count_col<-colorRamp2(c(min(met$gex_exonic_umis+met$gex_intronic_umis),
    max(met$gex_exonic_umis+met$gex_intronic_umis)), 
    c("white","black"))
  ha = HeatmapAnnotation(which="row",
    cell_type=met$predicted.id,
    #cnv_ref=met$cnv_ref,
    read_count= met$gex_exonic_umis+met$gex_intronic_umis,
          col = list(cell_type = type_cols,
            #cnv_ref=ref_cols,
            read_count=read_count_col))
  plt1<-Heatmap(cnv,
      show_row_names=F,
      show_column_names=F,
      column_order=1:ncol(cnv),
      col=col_fun,
      cluster_rows=dend,
      left_annotation=ha,
      column_split=final.obj@annotation.filt$Chr)
    pdf(paste0(dir_in,"_casper/",sample_name,".casper.heatmap.pdf"),width=20)
    print(plt1)
    dev.off()
    system(paste0("slack -F ",paste0(dir_in,"_casper/",sample_name,".casper.heatmap.pdf")," ryan_todo"))
}

lapply(sample_in,function(x) casper_per_sample_plot(dat=x))




```

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
## Run HMMcopy on ATAC Data
```R

library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)
library(SCOPE)
library(HMMcopy)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(reshape2)
library(philentropy)
library(dendextend)
library(Seurat)
library(Signac)
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

dat<-sample_in[7]
#set up single cell bam files from given directory
prepare_bam_bed_obj<-function(dat=sample_in[1],resol.){
    file_in=basename(dat)
    sample_name=substr(file_in,1,nchar(file_in)-17)
    dir_in=dirname(dat)
    dat<-readRDS(dat)
    print(paste("Read in file:",file_in))
    dat$cnv_ref<-"FALSE"
    dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    DefaultAssay(dat)<-"RNA"
    sampname<-row.names(dat@meta.data)
    print("Set up sample names.")

    #Initalization
    bamFile <- list.files(dir_in, pattern = "atac_possorted_bam.bam$")
    bamdir <- file.path(dir_in, bamFile)

    #set genomic window size to 500kb
    bambedObj <- get_bam_bed(bamdir = dir_in, sampname = sampname, hgref = "hg38",resolution=resol.,sex=T)#resolution of 100 = 100kbp
    print("Setting up reference genome.")
    #Compute GC content and mappability for each reference bin.
    mapp <- get_mapp(bambedObj$ref, hgref = "hg38")
    gc <- get_gc(bambedObj$ref, hgref = "hg38")
    values(bambedObj$ref) <- cbind(values(bambedObj$ref), DataFrame(gc, mapp))
    return(bambedObj)
}


#modified scope get_coverage_scDNA
read_in_reads<-function(x=bambedObj){
    bamurl <- list.files(x$bamdir, pattern = "atac_possorted_bam.bam$")
    what <- c("rname", "pos", "mapq", "qwidth")
    flag <- scanBamFlag(isPaired = TRUE,
            isUnmappedQuery = FALSE, 
            isNotPassingQualityControls = FALSE,
            isFirstMateRead = TRUE)
    param <- ScanBamParam(what = what, 
                          flag = flag,
                          tag="CB")  #added cell barcode tag to extraction
    print("Reading in bam file, this might take a while...")
    bam <- scanBam(paste0(x$bamdir,"/atac_possorted_bam.bam"), param = param)[[1]]
    bam.ref <- GRanges(seqnames = bam$rname, ranges = IRanges(start = bam[["pos"]],width = bam[["qwidth"]]))
    values(bam.ref) <- DataFrame(cellID=bam$tag$CB)
    bam.ref<-bam.ref[bam.ref$cellID %in% x$sampname,] 

    ##################filter to just sampnames, then split, then make a count overlap for each one of the nested list#######
    bam.ref<-split(bam.ref,bam.ref$cellID)

    bam.ref<-lapply(bam.ref,function(y) suppressWarnings(y[countOverlaps(y,mask.ref) == 0]))#filter out masked regions
    Y<- lapply(bam.ref,function(y) countOverlaps(x$ref, y))
    Y_out<-do.call("cbind",Y)
    colnames(Y_out)<-names(Y)
    Y_out<-as.data.frame(Y_out)
    ref_out<-as.data.frame(x$ref)
    row.names(Y_out)<-paste0(ref_out$seqnames,":",ref_out$start,"_",ref_out$end)
    return(Y_out)
}

#standard read correction
hmm_generator<-function(i){
    cellname<-colnames(Y_out)[i]
    dat<-cbind(ref,Y_out[,i])
    colnames(dat)[ncol(dat)]<-"reads"
    dat_copy <- suppressWarnings(correctReadcount(dat))
    write.table(dat_copy,
    file=paste0("counts_",cellname,".tsv"),
    sep="\t",quote=F,row.names=T,col.names=T)
    return(dat_copy)
}


#standard segmentation
copy_seg_func<-function(i){
    dat<-copy_estimate[[i]]
    dat<-HMMsegment(dat) #run HMMsegment
    segs<-makeGRangesFromDataFrame(dat$segs,keep.extra.columns=T)   #make granges object
    out_segs_idx<-as.data.frame(findOverlaps(segs,ref_granges)) #overlap GRanges for plotting consistent matrix bins
    out_segs<-segs[out_segs_idx$queryHits,]$state #return list of states length of reference bins
    return(out_segs)
}

perform_atac_cnv_calling<-function(dat,resol=1000){
  file_in=basename(dat)
  sample_name=substr(file_in,1,nchar(file_in)-17) # setting up output names here
  dir_in=dirname(dat)

  bambedObj<-prepare_bam_bed_obj(dat=dat,resol.=resol)
  #Set up genome masking of duplicated regions and gaps
  # Get segmental duplication regions
  seg.dup <- read.table(system.file("extdata", "GRCh38GenomicSuperDup.tab", package = "WGSmapp"))
  # Get hg38 gaps
  gaps <- read.table(system.file("extdata", "hg38gaps.txt", package = "WGSmapp"))
  #Set up genome masking from SCOPE get_masked_ref function
  seg.dup <- seg.dup[!is.na(match(seg.dup[,1], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
  seg.dup <- GRanges(seqnames = seg.dup[,1], ranges = IRanges(start = seg.dup[,2], end = seg.dup[,3]))
  gaps <- gaps[!is.na(match(gaps[,2], paste('chr', c(seq_len(22), 'X', 'Y'), sep = ''))),]
  gaps <- GRanges(seqnames = gaps[,2], ranges = IRanges(start = gaps[,3], end = gaps[,4]))
  mask.ref <- sort(c(seg.dup, gaps)) # Generate mask region

  #generate raw count matrix
  Y_out<-read_in_reads(x=bambedObj)
  saveRDS(Y_out,file="test.rawcount.500kbp.rds")

  #Prepare windows with gc and mappability data
    ref <- as.data.frame(bambedObj$ref)
    ref$gc<-ref$gc/100 #fitting to HMMcopy analysis
    ref$chr<-ref$seqname
    #ref$chr<-substr(ref$seqname,4,6)
    ref<-ref[c("chr","start","end","gc","mapp")]
    colnames(ref)<-c("chr","start","end","gc","map")
    ref_granges<-makeGRangesFromDataFrame(ref,keep.extra.columns=T) #used in the next few functions for overlapping consistent windows
    Y_out<-Y_out[,colSums(Y_out)>2000] #filter out cells with less that 2000 reads

  #Normalize bins
  copy_estimate<-lapply(1:ncol(Y_out),hmm_generator) #correct bins by gc content, 
  names(copy_estimate)<-colnames(Y_out)
  saveRDS(copy_estimate,file="copyestimate.500kb.rds")
  copy_estimate<-readRDS(file="copyestimate.500kb.rds")


#segmentation for plotting
copy_segmentation<-lapply(1:length(copy_estimate),copy_seg_func) #segment genome by copy estimate (log2 of cor.map)
copy_segmentation<-lapply(copy_segmentation,as.numeric)
names(copy_segmentation)<-colnames(Y_out)
copy_segmentation<-as.data.frame(do.call("cbind",copy_segmentation))
row.names(copy_segmentation)<-paste0(ref$chr,":",ref$start,"_",ref$end)
copy_segmentation<-as.data.frame(t(copy_segmentation))
saveRDS(copy_segmentation,file="copysegmentation.500kb.rds")
copy_segmentation<-readRDS(file="copysegmentation.500kb.rds")

copy_segmentation<-copy_segmentation[,!grepl("Y",colnames(copy_segmentation))] #remove Y



#filter to qc cells and set up annotations
    Y_plot<-as.data.frame(t(Y_out))
    Y_plot<-Y_plot[,!grepl("Y",colnames(Y_plot))]
    Y_plot<-as.data.frame(Y_plot %>% mutate_all(., ~ ifelse(.!=0, log10(.+0.00000000001), log10(1))))

    col_fun_reads=colorRamp2(quantile(unlist(Y_plot),c(0.1,0.2,0.3,0.5,0.6,0.8,0.9),na.rm=T),
    c("#336699","#99CCCC","#CCCCCC","#CCCC66","#CC9966","#993333","#990000"))

    col_fun_cn=structure(c("#2166ac", "#67a9cf", "#f7f7f7","#fddbc7","#ef8a62","#b2182b","#630410"), names = c("1", "2", "3", "4", "5", "6","7"))

#optimize segmentation 
#https://github.com/shahcompbio/single_cell_pipeline/blob/master/single_cell/workflows/hmmcopy/scripts/hmmcopy.R
#https://advances.sciencemag.org/content/6/50/eabd6454
#Set up clustering of cells
    dist_method="euclidean"
    dist_x<-philentropy::distance(copy_segmentation,method=dist_method,as.dist.obj=T,use.row.names=T)
    dend <- dist_x %>%  hclust(method="ward.D2") %>% as.dendrogram(edge.root=F,h=2) 
    k_search<-find_k(dend,krange=2:10) #search for optimal K from 5-10
    k_clus_number<-k_search$nc
    k_clus_id<-k_search$pamobject$clustering
    dend <- color_branches(dend, k = k_clus_number)    #split breakpoint object by clusters

#Read count raw
plt1<-Heatmap(Y_plot,
    show_row_names=F,
    #row_split=cell_annot$bulk,
    show_column_names=F,
    column_order=1:ncol(Y_plot),
    col=col_fun_reads,
    cluster_rows=dend,
    row_title="Raw read count",
    name="Log10 Reads",#,
    column_split=factor(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(Y_plot),":"),"[",1))))
    )

#HMM segmentation
plt2<-Heatmap(copy_segmentation,
    show_row_names=F,
    show_column_names=F,
    #row_split=cell_annot$bulk,
    column_order=1:ncol(copy_segmentation),
    col=col_fun_cn,
    cluster_rows=dend,
    row_title="Segmentation",
    name="copy state",
    column_split=factor(unlist(lapply(strsplit(colnames(copy_segmentation),":"),"[",1)),levels=unique(unlist(lapply(strsplit(colnames(copy_segmentation),":"),"[",1))))
)
    #left_annotation=ha,    )


pdf("HMMcopy_test.pdf",width=10,height=3)
par(mfrow=c(2,1))
plt1
plt2
dev.off()

system("slack -F HMMcopy_test.pdf ryan_todo")

}


perform_atac_cnv_calling(dat=sample_in[7])
```

## Add per cell subtyping to epithelial cells

```R
library(Signac)
library(Seurat)
library(ggplot2)
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

#Features determined by EMBO manuscript
module_feats<-list()
module_feats[["Basal_SC"]]=c('EMP1', 'TAGLN', 'TTYH1', 'RTN4', 'TK1', 'BUB3', 'IGLV3.25', 'FAM3C', 'TMEM123', 'KDM5B', 'KRT14', 'ALG3', 'KLK6', 'EEF2', 'NSMCE4A', 'LYST', 'DEDD', 'HLA.DRA', 'PAPOLA', 'SOX4', 'ACTR3B', 'EIF3D', 'CACYBP', 'RARRES1', 'STRA13', 'MFGE8', 'FRZB', 'SDHD', 'UCHL1', 'TMEM176A', 'CAV2', 'MARCO', 'P4HB', 'CHI3L2', 'APOE', 'ATP1B1', 'C6orf15', 'KRT6B', 'TAF1D', 'ACTA2', 'LY6D', 'SAA2', 'CYP27A1', 'DLK1', 'IGKV1.5', 'CENPW', 'RAB18', 'TNFRSF11B', 'VPS28', 'HULC', 'KRT16', 'CDKN2A', 'AHNAK2', 'SEC22B', 'CDC42EP1', 'HMGA1', 'CAV1', 'BAMBI', 'TOMM22', 'ATP6V0E2', 'MTCH2', 'PRSS21', 'HDAC2', 'ZG16B', 'GAL', 'SCGB1D2', 'S100A2', 'GSPT1', 'ARPC1B', 'NIT1', 'NEAT1', 'DSC2', 'RP1.60O19.1', 'MAL2', 'TMEM176B', 'CYP1B1', 'EIF3L', 'FKBP4', 'WFDC2', 'SAA1', 'CXCL17', 'PFDN2', 'UCP2', 'RAB11B', 'FDCSP', 'HLA.DPB1', 'PCSK1N', 'C4orf48', 'CTSC')
module_feats[["Her2E_SC"]]=c('PSMA2', 'PPP1R1B', 'SYNGR2', 'CNPY2', 'LGALS7B', 'CYBA', 'FTH1', 'MSL1', 'IGKV3.15', 'STARD3', 'HPD', 'HMGCS2', 'ID3', 'NDUFB8', 'COTL1', 'AIM1', 'MED24', 'CEACAM6', 'FABP7', 'CRABP2', 'NR4A2', 'COX14', 'ACADM', 'PKM', 'ECH1', 'C17orf89', 'NGRN', 'ATG5', 'SNHG25', 'ETFB', 'EGLN3', 'CSNK2B', 'RHOC', 'PSENEN', 'CDK12', 'ATP5I', 'ENTHD2', 'QRSL1', 'S100A7', 'TPM1', 'ATP5C1', 'HIST1H1E', 'LGALS1', 'GRB7', 'AQP3', 'ALDH2', 'EIF3E', 'ERBB2', 'LCN2', 'SLC38A10', 'TXN', 'DBI', 'RP11.206M11.7', 'TUBB', 'CRYAB', 'CD9', 'PDSS2', 'XIST', 'MED1', 'C6orf203', 'PSMD3', 'TMC5', 'UQCRQ', 'EFHD1', 'BCAM', 'GPX1', 'EPHX1', 'AREG', 'CDK2AP2', 'SPINK8', 'PGAP3', 'NFIC', 'THRSP', 'LDHB', 'MT1X', 'HIST1H4C', 'LRRC26', 'SLC16A3', 'BACE2', 'MIEN1', 'AR', 'CRIP2', 'NME1', 'DEGS2', 'CASC3', 'FOLR1', 'SIVA1', 'SLC25A39', 'IGHG1', 'ORMDL3', 'KRT81', 'SCGB2B2', 'LINC01285', 'CXCL8', 'KRT15', 'RSU1', 'ZFP36L2', 'DKK1', 'TMED10', 'IRX3', 'S100A9', 'YWHAZ')
module_feats[["LumA_SC"]]=c('SH3BGRL', 'HSPB1', 'PHGR1', 'SOX9', 'CEBPD', 'CITED2', 'TM4SF1', 'S100P', 'KCNK6', 'AGR3', 'MPC2', 'CXCL13', 'RNASET2', 'DDIT4', 'SCUBE2', 'KRT8', 'MZT2B', 'IFI6', 'RPS26', 'TAGLN2', 'SPTSSA', 'ZFP36L1', 'MGP', 'KDELR2', 'PPDPF', 'AZGP1', 'AP000769.1', 'MYBPC1', 'S100A1', 'TFPI2', 'JUN', 'SLC25A6', 'HSP90AB1', 'ARF5', 'PMAIP1', 'TNFRSF12A', 'FXYD3', 'RASD1', 'PYCARD', 'PYDC1', 'PHLDA2', 'BZW2', 'HOXA9', 'XBP1', 'AGR2', 'HSP90AA1') 
module_feats[["LumB_SC"]]=c('UGCG', 'ARMT1', 'ISOC1', 'GDF15', 'ZFP36', 'PSMC5', 'DDX5', 'TMEM150C', 'NBEAL1', 'CLEC3A', 'GADD45G', 'MARCKS', 'FHL2', 'CCDC117', 'LY6E', 'GJA1', 'PSAP', 'TAF7', 'PIP', 'HSPA2', 'DSCAM.AS1', 'PSMB7', 'STARD10', 'ATF3', 'WBP11', 'MALAT1', 'C6orf48', 'HLA.DRB1', 'HIST1H2BD', 'CCND1', 'STC2', 'NR4A1', 'NPY1R', 'FOS', 'ZFAND2A', 'CFL1', 'RHOB', 'LMNA', 'SLC40A1', 'CYB5A', 'SRSF5', 'SEC61G', 'CTSD', 'DNAJC12', 'IFITM1', 'MAGED2', 'RBP1', 'TFF1', 'APLP2', 'TFF3', 'TRH', 'NUPR1', 'EMC3', 'TXNIP', 'ARPC4', 'KCNE4', 'ANPEP', 'MGST1', 'TOB1', 'ADIRF', 'TUBA1B', 'MYEOV2', 'MLLT4', 'DHRS2', 'IFITM2')
module_feats[["proliferation_score"]]<-c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")


  #Cicero processing functionq
subtyping<-function(dat){
      dat_file_path=dat
      file_in=basename(dat)
      sample_name=substr(file_in,1,nchar(file_in)-17)
      dir_in=dirname(dat)
      dat<-readRDS(dat) #read in as seurat object
      dat_sub<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial"))
      module_scores<-AddModuleScore(dat_sub,features=module_feats,assay="RNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
      module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-4,ncol(module_scores@meta.data))]
      colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like
      dat<-AddMetaData(dat,metadata=module_scores)
      saveRDS(dat,file=dat_file_path)
}

#run on all samples
lapply(sample_in,subtyping)
subtyping(dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")


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
  library(BiocParallel)
  register(MulticoreParam(5))

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

  #Read in data and modify to monocle CDS file
  #read in RDS file.

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species =9606, all_versions = FALSE))

chromvar_per_sample<-function(dat){
  dat_file_path=dat
  file_in=basename(dat)
  sample_name=substr(file_in,1,nchar(file_in)-17)
  dir_in=dirname(dat)
  dat<-readRDS(dat) #read in as seurat object

  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  # Create a new Mofif object to store the results
  #for combined
  peaks<-granges(dat[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
  motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
  motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  dat <- SetAssayData(object = dat, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  dat <- RegionStats(object = dat, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  dat <- RunChromVAR( object = dat,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(dat,file=dat_file_path)

}

lapply(sample_in,chromvar_per_sample)

#chromvar of merged samples
chromvar_per_sample(dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")


```

## Cicero

```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(cicero,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") #using old libraries because exacloud doesn't like the new ones
  library(monocle3,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") #using old libraries because exacloud doesn't like the new ones
  library(SeuratObjects)
  library(EnsDb.Hsapiens.v86)
  setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")
  
  sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

  #Cicero processing function
  cicero_processing<-function(dat){
      dat_file_path=dat
      file_in=basename(dat)
      sample_name=substr(file_in,1,nchar(file_in)-17)
      dir_in=dirname(dat)
      dat<-readRDS(dat) #read in as seurat object

      #Generate CDS format from Seurat object
      atac.cds <- as.CellDataSet(dat,assay="peaks",reduction="multimodal_umap")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = atac.cds@reducedDimS)
      saveRDS(atac.cicero,paste0(dir_in,"/",sample_name,".atac_cicero_cds.Rds"))
        

      # extract gene annotations from EnsDb
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevels(annotations)<-paste0("chr",seqlevels(annotations))

      # change to UCSC style since the data was mapped to hg38
      #seqlevelsStyle(annotations) <- 'UCSC'
      genome(annotations) <- "hg38"

      genome <- annotations@seqinfo # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = genome@seqnames, "length" = genome@seqlengths) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste0(dir_in,"/",sample_name,".atac_cicero_conns.Rds"))
      
      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste0(dir_in,"/",sample_name,".atac_cicero_ccans.Rds"))
      
      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      DefaultAssay(dat)<-"peaks"
      Links(dat) <- links

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

      atac.cds<- annotate_cds_by_site(atac.cds, gene_annotation)
      unnorm_ga <- build_gene_activity_matrix(atac.cds, conns)
      saveRDS(unnorm_ga,paste0(dir_in,"/",sample_name,".atac_cicero_unnormGA.Rds"))
      # remove any rows/columns with all zeroes
      unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                             !Matrix::colSums(unnorm_ga) == 0]
      dat[['GeneActivity']]<- CreateAssayObject(counts = unnorm_ga) 
      # normalize
      # note for CCAN comparisons, multiple files should be normalized together
      dat <- NormalizeData(
        object = dat,
        assay = 'GeneActivity',
        normalization.method = 'LogNormalize',
        scale.factor = median(dat$nFeature_GeneActivity)
      )
      saveRDS(dat,dat_file_path)
  }

lapply(sample_in,cicero_processing)

#run on all merged sample
cicero_processing(dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")

```


### Integration: Now Clustering together on RNA profiles using harmony to integrate

```R
library(harmony,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/lib_backup_210125")
library(cisTopic)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1")


###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479", "Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", "PVL" ="#F2ACCA")

diag_cols<-c("IDC"="red", "DCIS"="grey")

molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
########################################


dat<-readRDS("phase1.SeuratObject.rds")
pdf("phase1.Harmonyconvergence.pdf")
dat<-RunHarmony(dat,group.by.vars="sample",reduction="pca",plot_convergence=T)
dev.off()
system("slack -F phase1.Harmonyconvergence.pdf ryan_todo")

dat <- dat %>% 
    RunUMAP(reduction.name="harmonyumap",reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

alpha_val=0.33

plt1<-DimPlot(dat,reduction="harmonyumap",group.by="sample")#+scale_fill_manual(samp_cols)
plt2<-DimPlot(dat,reduction="harmonyumap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))
plt3<-DimPlot(dat,reduction="harmonyumap",group.by="diagnosis",cols=alpha(diag_cols,alpha_val))
plt4<-DimPlot(dat,reduction="harmonyumap",group.by="molecular_type",cols=alpha(molecular_type_cols,alpha_val))
ggsave(plt1|plt2|plt3|plt4,file="phase1.Harmonyumap.pdf",width=20)
system("slack -F phase1.Harmonyumap.pdf ryan_todo")

saveRDS(dat,"phase1.SeuratObject.rds")

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


