---
title: 10X Multiome Tumor Phase 1
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_phase1/
category: CEDAR
---

Multiome processing for 10X multiome data

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
  saveRDS(dat,file=paste0(i,".SeuratObject.rds"))
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
peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

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

single_sample_cluster<-function(x){
dat<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
# call peaks using MACS2
peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

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
saveRDS(dat,file=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))

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
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
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
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.rds"))
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
  atac_sub<-readRDS(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds"))
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
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.rds"))
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
#for i in 1 3 4 5 6 7 8 9 10 11 12;
#do Rscript /home/groups/CEDAR/mulqueen/src/single_sample_cistopic.Rscript $i ; done &

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
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")
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
  sample_name=substr(file_in,1,nchar(file_in)-17)
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

#run umap projections of all dim reductions per sample
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))
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
  sample_name=substr(file_in,1,nchar(file_in)-17)
  dir_in=dirname(dat)
  dat<-readRDS(dat) #read in as seurat object
  outname=paste0(dir_in,"/",sample_name,".")

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
  p1<-DimPlot(dat,reduction="rna_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("RNA UMAP")
  p2<-DimPlot(dat,reduction="atac_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("ATAC UMAP")
  p3<-DimPlot(dat,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (LSI)")


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
  p5<-DimPlot(dat2,reduction="multimodal_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")
  p6<-DimPlot(dat,reduction="cistopic_umap",group.by="predicted.id",cols=alpha(type_cols,alpha_val))+ggtitle("Cistopic UMAP")

  #Finally Plot results
  plt<-(p1 | p2)/(p3)/(p6 | p5)
  ggsave(plt,file=paste0(outname,"multimodal.umap.pdf"),width=10,height=20)
  system(paste0("slack -F ",outname,"multimodal.umap.pdf ryan_todo"))
  saveRDS(dat,file=dat_file_path)

  plt_cell_count<-dat@meta.data[,c("sample","predicted.id","diagnosis","molecular_type")]
  return(plt_cell_count)
}


#run umap projections of merged samples
plot_reductions(dat="/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/phase1.SeuratObject.rds")

#run umap projections of all dim reductions per sample
sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

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

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))

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

  dat <- LinkPeaks(
    object = dat,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use = unlist(geneset)
  )

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

plot_genes(dat=sample_in[1])
lapply(sample_in,function(x) plot_genes(dat=x))

```


### Determine Tumor Cells via InferCNV

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

sample_in<-unlist(lapply(c(1,3,4,5,6,7,8,9,10,11,12),function(x) paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".SeuratObject.rds")))


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
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Endothelial","B-cells","Cancer Epithelial","Normal Epithelial")) #just look at epithelial

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
                               resume_mode=T)
  saveRDS(infercnv_obj,paste0(outname,"_inferCNV","/",basename(outname),"inferCNV.Rds"))
  system(paste0("slack -F ",outname,"_inferCNV","/","infercnv.png"," -T ","\"",outname,"\"" ," ryan_todo") )
  #infercnv_obj<-readRDS(paste0(outname,"_inferCNV","/",basename(outname),"inferCNV.Rds"))

}
lapply(sample_in,function(x) infercnv_per_sample(dat=x))


```

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



# celltypes<-row_order(plt) #get row order from heatmap, this includes the slice
# metadata_cluster<-as.data.frame(cbind(cellID=rownames(mat),celltype=0))
# metadata_cluster[unlist(celltypes[[1]]),]$celltype<-"1"
# metadata_cluster[unlist(celltypes[[2]]),]$celltype<-"2"
# table(metadata_cluster$celltype)
# row.names(metadata_cluster)<-metadata_cluster$cellID
# metadata_cluster$cnv_profile<-ifelse(metadata_cluster$celltype==1,"MCF7","T47D") #this is due to the mcf7 amplification we have seen in other data

# cell_line<-setNames(metadata_cluster$cnv_profile,row.names(metadata_cluster))
# dat<-AddMetaData(dat,metadata=cell_line,col.name="cnv_profile")#assign slices to cell names

# #cluster with low resolution on peaks dataset
# dat$peaks_cluster<-ifelse(as.numeric(dat@reductions$atac_umap@cell.embeddings[,1])<0,"T47D","MCF7") #it is clear from the atac data that one cluster is T47D and one is MCF7 based on chr17 profiles
# table(dat$peaks_cluster)

# plt1<-DimPlot(dat,split.by="cnv_profile",group.by="cnv_profile",na.value=NA)+ggtitle("CNV Profile Split")
# plt2<-DimPlot(dat,group.by="sample")+ggtitle("ATAC Samples")
# plt3<-DimPlot(dat,group.by="peaks_cluster")+ggtitle("Final Cell Type Assignment")

# ggsave(plt1/(plt2|plt3),file="YW_celltype.pdf",width=10)
# system("slack -F YW_celltype.pdf ryan_todo")


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

#Run SC_subtype, proportion of cells in each, what is epigenome of each?
#Her2 Amp
#Clinical data add to metadat, plot cell type distro by diagnosis, recolor umaps by diagnosis
#Magic on T cell imputations


#ER binding poor and good outcome from patients, overlap with ATAC data
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Poor%20outcome%20ER%20regions.bed.gz
http://www.carroll-lab.org.uk/FreshFiles/Data/RossInnes_Nature_2012/Good%20outcome%20ER%20regions.bed.gz

#sample specificity between cancer and normal epithelial, superenhancers are open more frequently
#plot epithelial cell subtype per sample
#use cnv to infer changes, then do lineage comparisons, plot BAF across regions
#use HMMcopy to test ATAC changes
#plot cell prediction specificity (predicted ID values density plot colored by predicted ID cell types)



<!---
###############################################################






### Run Dimensionality Reduction on just Epithelial cells

```R
###CisTopic###
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
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

cistopic_generation<-function(x,name_out,outdir="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi"){
  wd<-outdir
  outname<-name_out
  atac_sub<-x
  cistopic_counts_frmt<-atac_sub$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(10:30),nCores=5,addModels=FALSE)
  sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =x@meta.data)
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
  pdf(paste0(wd,"/",outname,".umap.pdf"))
  print(plt1)
  print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,file=paste0(wd,"/",outname,".SeuratObject.Rds"))
  }


#Running cistopic
  dat<-readRDS("rm_merged.SeuratObject.rds")
  dat<-subset(dat,predicted.id %in% c("Normal Epithelial","Cancer Epithelial"))
  cistopic_generation(x=dat,name_out="rm_epithelial")

```



## Motif Footprinting
Based on https://satijalab.org/signac/articles/footprint.html
```R
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  set.seed(1234)
  library(motifmatchr)
  library(parallel)
  setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")

combined<-readRDS("210924_cellline.SeuratObject.Rds")
mcf7<-readRDS("210924_mcf7.SeuratObject.Rds")
t47d<-readRDS("210924_t47d.SeuratObject.Rds")

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
DefaultAssay(combined)<-"peaks"
combined <- AddMotifs(combined, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
DefaultAssay(mcf7)<-"peaks"
mcf7 <- AddMotifs(mcf7, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
DefaultAssay(t47d)<-"peaks"
t47d <- AddMotifs(t47d, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

saveRDS(combined,file="210924_cellline.SeuratObject.Rds")
saveRDS(mcf7,file="210924_mcf7.SeuratObject.Rds")
saveRDS(t47d,file="210924_t47d.SeuratObject.Rds")

#function for plotting footprints in data sets

plot_footprints<-function(x,footprints,outname){
  x <- Footprint(object = x, motif.name = footprints, genome = BSgenome.Hsapiens.UCSC.hg38, in.peaks=T)   # gather the footprinting information for sets of motifs
  p2 <- PlotFootprint(x, features = footprints,label=F)   # plot the footprint data for each group of cells, might want to change idents
  return(p2)
}

Idents(combined)<-combined$origin
out<-mclapply(c("GATA3","ESR1","FOXA1","CTCF"),FUN=function(z) plot_footprints(x=combined,footprints=z,outname="210924_cellline"),mc.cores=5) #plot and return 5 at a time

outname="210924_cellline"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out))
print(wrap_plots(out) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))

```

## Cistopic Correlation to TITAN (RNA Topic Modelling)
Titan was run by Aaron Doe on the RNA data. Located here:
```bash
/home/groups/CEDAR/doe/projects/ATACRNA

#topics run by cell line (provided by AD)
/home/groups/CEDAR/mulqueen/projects/10x_atacrna/MCF7_merged_Topics_cellTopic_table.tsv #mcf7
/home/groups/CEDAR/mulqueen/projects/10x_atacrna/T47D_merged_Topics_cellTopic_table.tsv #t47d

```

# Set up Signac objects for topic comparison. 
```R
library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(circlize)
library(viridis)
library(reshape2)

setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")
mcf7_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/MCF7_RNA_strictQC.rds")
t47d_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/T47D_RNA_strictQC.rds")

plt1<-DimPlot(mcf7_rna,reduction="umap",group.by=c("origin","seurat_clusters"))
ggsave(plt1,file="mcf7_rna_umap.pdf")
system("slack -F mcf7_rna_umap.pdf ryan_todo")

plt1<-DimPlot(t47d_rna,reduction="umap",group.by=c("origin","seurat_clusters"))
ggsave(plt1,file="t47d_rna_umap.pdf")
system("slack -F t47d_rna_umap.pdf ryan_todo")

#mcf7_bins<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/MCF7_ATAC_wBins.rds")
#t47d_bins<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/T47D_ATAC_wBins.rds")

mcf7_atac<-readRDS("210924_mcf7.SeuratObject.Rds")
t47d_atac<-readRDS("210924_t47d.SeuratObject.Rds")

t47d_cisTopicObject<-readRDS("210924_t47d.CisTopicObject.Rds")
mcf7_cisTopicObject<-readRDS("210924_mcf7.CisTopicObject.Rds")

#add cistrome data to cells
t47d_atac[['cistrome_score']] <- CreateAssayObject(data = t(t47d_cisTopicObject@cell.data[(row.names(t47d_cisTopicObject@cell.data) %in% row.names(t47d_atac@meta.data)),15:ncol(t47d_cisTopicObject@cell.data)]))
mcf7_atac[['cistrome_score']] <- CreateAssayObject(data = t(mcf7_cisTopicObject@cell.data[(row.names(mcf7_cisTopicObject@cell.data) %in% row.names(mcf7_atac@meta.data)),15:ncol(mcf7_cisTopicObject@cell.data)]))


#find DA motifs from top topic differences
#set up a named vector of topic11 and topic17 cells to add to seurat object
topic_11_score<-setNames(mcf7_rna$ImputedTopic_11,paste0(row.names(mcf7_rna@meta.data),mcf7_rna$origin))
topic_17_score<-setNames(mcf7_rna$ImputedTopic_17,paste0(row.names(mcf7_rna@meta.data),mcf7_rna$origin))
mcf7_atac<-AddMetaData(object=mcf7_atac,metadata=topic_11_score,col.name="ImputedTopic_11")
mcf7_atac<-AddMetaData(object=mcf7_atac,metadata=topic_17_score,col.name="ImputedTopic_17")
mcf7_atac$topic_bin<-"NA"
mcf7_atac@meta.data[(mcf7_atac@meta.data$ImputedTopic_11>=quantile(mcf7_atac$ImputedTopic_11,0.75)) & (mcf7_atac@meta.data$ImputedTopic_17 <= quantile(mcf7_atac$ImputedTopic_17,0.25)),]$topic_bin<-"Topic11"
mcf7_atac@meta.data[(mcf7_atac@meta.data$ImputedTopic_17>=quantile(mcf7_atac$ImputedTopic_17,0.75)) & (mcf7_atac@meta.data$ImputedTopic_11 <= quantile(mcf7_atac$ImputedTopic_11,0.25)),]$topic_bin<-"Topic17"
plt<-ggplot(mcf7_atac@meta.data,aes(x=ImputedTopic_11,y=ImputedTopic_17,color=topic_bin))+geom_point()+theme_minimal()
ggsave(plt,file="mcf7_topic_scatterplot.pdf")
system("slack -F mcf7_topic_scatterplot.pdf ryan_todo")
table(mcf7_atac$topic_bin)
mcf7_atac_subset<-subset(mcf7_atac, topic_bin %in% c("Topic11","Topic17"))
Idents(mcf7_atac_subset)<-mcf7_atac_subset$topic_bin
saveRDS(mcf7_atac_subset,file="mcf7_atac_subset.Rds")
saveRDS(mcf7_atac,file="210924_mcf7.SeuratObject.Rds")

#and for t47d
topic_11_score<-setNames(t47d_rna$ImputedTopic_11,paste0(row.names(t47d_rna@meta.data),t47d_rna$origin))
topic_17_score<-setNames(t47d_rna$ImputedTopic_17,paste0(row.names(t47d_rna@meta.data),t47d_rna$origin))
t47d_atac<-AddMetaData(object=t47d_atac,metadata=topic_11_score,col.name="ImputedTopic_11")
t47d_atac<-AddMetaData(object=t47d_atac,metadata=topic_17_score,col.name="ImputedTopic_17")
t47d_atac$topic_bin<-"NA"
t47d_atac@meta.data[(t47d_atac@meta.data$ImputedTopic_11>=quantile(t47d_atac$ImputedTopic_11,0.75)) & (t47d_atac@meta.data$ImputedTopic_17 <= quantile(t47d_atac$ImputedTopic_17,0.25)),]$topic_bin<-"Topic11"
t47d_atac@meta.data[(t47d_atac@meta.data$ImputedTopic_17>=quantile(t47d_atac$ImputedTopic_17,0.75)) & (t47d_atac@meta.data$ImputedTopic_11 <= quantile(t47d_atac$ImputedTopic_11,0.25)),]$topic_bin<-"Topic17"
plt<-ggplot(t47d_atac@meta.data,aes(x=ImputedTopic_11,y=ImputedTopic_17,color=topic_bin))+geom_point()+theme_minimal()
ggsave(plt,file="t47d_topic_scatterplot.pdf")
system("slack -F t47d_topic_scatterplot.pdf ryan_todo")
t47d_atac_subset<-subset(t47d_atac, topic_bin %in% c("Topic11","Topic17"))
Idents(t47d_atac_subset)<-t47d_atac_subset$topic_bin
saveRDS(t47d_atac_subset,file="t47d_atac_subset.Rds")
saveRDS(t47d_atac,file="210924_t47d.SeuratObject.Rds")

```

## Perform topic comparisons


```R
library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(circlize)
library(viridis)
library(reshape2)


setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")
mcf7_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/MCF7_RNA_strictQC.rds")
t47d_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/T47D_RNA_strictQC.rds")

mcf7_atac<-readRDS("210924_mcf7.SeuratObject.Rds")
t47d_atac<-readRDS("210924_t47d.SeuratObject.Rds")

#Correlate cistopic with titan topics
correlate_cistopic_titan<-function(x,y,outname){
  cistopic_dat<-x@reductions$cistopic@cell.embeddings #get cistopic from x
  row.names(cistopic_dat)<-substr(row.names(cistopic_dat),0,18)
  titan_dat<-as.data.frame(y@reductions$imputedLDA@cell.embeddings) #get titan from y
  cells<-row.names(titan_dat)[which(row.names(titan_dat) %in% row.names(cistopic_dat))] # set list of shared cell names
  sum(cells %in% row.names(titan_dat))==sum(cells %in% row.names(titan_dat)) #check to make sure same number of cells
  titan_dat<-titan_dat[cells,]#reorder data frames
  cistopic_dat<-cistopic_dat[cells,]
  cor_mat<-data.frame(matrix(ncol = ncol(cistopic_dat), nrow = ncol(titan_dat))) #make an empty data frame to populate
  colnames(cor_mat)<-colnames(cistopic_dat)
  row.names(cor_mat)<-colnames(titan_dat)
  for (yi in 1:ncol(cistopic_dat)){
    for (xj in 1:ncol(titan_dat)){
      cor_mat[xj,yi]<-cor(x=as.numeric(cistopic_dat[,yi]),y=as.numeric(titan_dat[,xj]))
    }
  }

  out_plt<-Heatmap(cor_mat)
  pdf(paste0(outname,".cistopic_titan_topiccorrelation.pdf"))
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistopic_titan_topiccorrelation.pdf")," ryan_todo"))

  #add cistopic square correlation matrix
  cor_cistopic<-cor(cistopic_dat)
  out_plt<-Heatmap(cor_cistopic)
  pdf(paste0(outname,".cistopic_squarecorrelation.pdf"))
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistopic_squarecorrelation.pdf")," ryan_todo"))

  #titan square correlation
  cor_titan<-cor(titan_dat)
  out_plt<-Heatmap(cor_titan)
  pdf(paste0(outname,".titan_squarecorrelation.pdf"))
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".titan_squarecorrelation.pdf")," ryan_todo"))

}

#x is atac library signac object with cistopic in reductions slot
#y is rna library seurat object with imputedLDA in reductions slot
correlate_cistopic_titan(x=mcf7_atac,y=mcf7_rna,outname="mcf7")
correlate_cistopic_titan(x=t47d_atac,y=t47d_rna,outname="t47d")

#now to correlate chromvar motifs with titan topics

#Correlate chromvar with titan topics
correlate_chromvar_titan<-function(x,y,outname,cellline="MCF",quantile_cutoff=0.95){
  chromvar_dat<-as.data.frame(t(x@assays$chromvar@data)) #get chromvar from x
  row.names(chromvar_dat)<-substr(row.names(chromvar_dat),0,18)
  tfList <- getMatrixByID(JASPAR2020, ID=colnames(chromvar_dat)) #set up readable chromvar names
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  colnames(chromvar_dat)<-tfList
  titan_dat<-as.data.frame(y@reductions$imputedLDA@cell.embeddings) #get titan from y  
  cells<-row.names(titan_dat)[which(row.names(titan_dat) %in% row.names(chromvar_dat))] # set list of shared cell names
  sum(cells %in% row.names(titan_dat))==sum(cells %in% row.names(titan_dat)) #check to make sure same number of cells
  titan_dat<-titan_dat[cells,]#reorder data frames
  chromvar_dat<-chromvar_dat[cells,]
  cor_mat<-data.frame(matrix(ncol = ncol(chromvar_dat), nrow = ncol(titan_dat))) #make an empty data frame to populate
  colnames(cor_mat)<-colnames(chromvar_dat)
  row.names(cor_mat)<-colnames(titan_dat)

  for (yi in 1:ncol(chromvar_dat)){
    for (xj in 1:ncol(titan_dat)){
      cor_mat[xj,yi]<-cor(x=chromvar_dat[,yi],y=titan_dat[,xj],use="complete.obs")
    }
  }


  out_plt<-Heatmap(t(cor_mat),
    row_names_gp = grid::gpar(fontsize = 4),
    row_km=5
    )
  pdf(paste0(outname,".chromvar_titan_topiccorrelation.pdf"),height=30)
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".chromvar_titan_topiccorrelation.pdf")," ryan_todo"))



  cor_mat<-as.data.frame(t(cor_mat[c("imputedLDA_11","imputedLDA_17"),]))
  boundary<-quantile(abs(cor_mat$imputedLDA_11-cor_mat$imputedLDA_17),quantile_cutoff)
  cor_mat$color_fill<-ifelse(abs(cor_mat$imputedLDA_11-cor_mat$imputedLDA_17)>boundary,"red","not_red")
  cor_mat$label<-""
  cor_mat[cor_mat$color_fill=="red",]$label<-row.names(cor_mat[cor_mat$color_fill=="red",])

  plt<-ggplot(dat=cor_mat,aes(x=imputedLDA_11,y=imputedLDA_17,label=label,color=color_fill))+
  geom_point(size=0.1)+
  scale_color_manual(values=c("red"="red","not_red"="grey50"))+
  theme_minimal()+
  geom_text_repel(force=5,max.overlaps=Inf,size=2,min.segment.length = 0,box.padding=0.3,segment.size=0.1)+
  geom_abline(intercept = 0, slope = 1, col = "black")+
  geom_abline(intercept = boundary, slope = 1, col = "black",linetype="dashed")+
  geom_abline(intercept = -boundary, slope = 1, col = "black",linetype="dashed")+
  ggtitle(paste(outname,"chromvar v. titan correlations"))+
  xlab("Topic 11")+
  ylab("Topic 17")+
  theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank(),legend.position="none")+
  scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

  ggsave(plt,file=paste0(outname,".chromvar_titan_topiccorrelation.comparison.pdf"),height=4,width=4,units="in")
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".chromvar_titan_topiccorrelation.comparison.pdf")," ryan_todo"))
}

#x is atac library signac object with chromvar in assay slot
#y is rna library seurat object with imputedLDA in reductions slot
correlate_chromvar_titan(x=mcf7_atac,y=mcf7_rna,outname="mcf7")
correlate_chromvar_titan(x=t47d_atac,y=t47d_rna,outname="t47d")

#Correlate cistrome data with titan topics
correlate_cistrome_titan<-function(x,y,outname,cellline="MCF",quantile_cutoff=0.95){
  cistrome_dat<-as.data.frame(t(x@assays$cistrome_score@data)) #get cistrome from x
  row.names(cistrome_dat)<-substr(row.names(cistrome_dat),0,18)
  titan_dat<-as.data.frame(y@reductions$imputedLDA@cell.embeddings) #get titan from y  
  cells<-row.names(titan_dat)[which(row.names(titan_dat) %in% row.names(cistrome_dat))] # set list of shared cell names
  sum(cells %in% row.names(titan_dat))==sum(cells %in% row.names(titan_dat)) #check to make sure same number of cells
  titan_dat<-titan_dat[cells,]#reorder data frames
  cistrome_dat<-cistrome_dat[cells,]
  cor_mat<-data.frame(matrix(ncol = ncol(cistrome_dat), nrow = ncol(titan_dat))) #make an empty data frame to populate
  colnames(cor_mat)<-colnames(cistrome_dat)
  row.names(cor_mat)<-colnames(titan_dat)
  for (yi in 1:ncol(cistrome_dat)){
    for (xj in 1:ncol(titan_dat)){
      cor_mat[xj,yi]<-cor(x=cistrome_dat[,yi],y=titan_dat[,xj],use="complete.obs")
    }
  }

  cor_mat<-cor_mat[startsWith(prefix=cellline,x=colnames(cor_mat))]

  out_plt<-Heatmap(t(cor_mat),
    row_names_gp = grid::gpar(fontsize = 4),
    row_km=5
    )
  pdf(paste0(outname,".cistrome_titan_topiccorrelation.pdf"),height=30)
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistrome_titan_topiccorrelation.pdf")," ryan_todo"))



  cor_mat<-as.data.frame(t(cor_mat[c("imputedLDA_11","imputedLDA_17"),]))
  boundary<-quantile(abs(cor_mat$imputedLDA_11-cor_mat$imputedLDA_17),quantile_cutoff)
  cor_mat$color_fill<-ifelse(abs(cor_mat$imputedLDA_11-cor_mat$imputedLDA_17)>boundary,"red","not_red")
  cor_mat$label<-""
  hyphen_length<-lapply(strsplit(row.names(cor_mat[cor_mat$color_fill=="red",]),"-"),length)[[1]][1] #adduming no hypens in names
  cor_mat[cor_mat$color_fill=="red",]$label<-paste(
    unlist(lapply(strsplit(row.names(cor_mat[cor_mat$color_fill=="red",]),"-"),"[",hyphen_length-1)),
    unlist(lapply(strsplit(row.names(cor_mat[cor_mat$color_fill=="red",]),"-"),"[",hyphen_length)))

  plt<-ggplot(dat=cor_mat,aes(x=imputedLDA_11,y=imputedLDA_17,label=label,color=color_fill))+
  geom_point(size=0.1)+
  scale_color_manual(values=c("red"="red","not_red"="grey50"))+
  theme_minimal()+
  geom_text_repel(force=5,max.overlaps=Inf,size=2,min.segment.length = 0,box.padding=0.3,segment.size=0.1)+
  geom_abline(intercept = 0, slope = 1, col = "black")+
  geom_abline(intercept = boundary, slope = 1, col = "black",linetype="dashed")+
  geom_abline(intercept = -boundary, slope = 1, col = "black",linetype="dashed")+
  ggtitle(paste(outname,"cistrome v. titan correlations"))+
  xlab("Topic 11")+
  ylab("Topic 17")+
  theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank(),legend.position="none")+
  scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

  ggsave(plt,file=paste0(outname,".cistrome_titan_topiccorrelation.comparison.pdf"),height=4,width=4,units="in")
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistrome_titan_topiccorrelation.comparison.pdf")," ryan_todo"))
}

#x is atac library signac object with cistrome in assay slot
#y is rna library seurat object with imputedLDA in reductions slot
correlate_cistrome_titan(x=mcf7_atac,y=mcf7_rna,outname="mcf7",cellline="MCF")
correlate_cistrome_titan(x=t47d_atac,y=t47d_rna,outname="t47d",cellline="T47D")


#find differences in chromvar usage 
  mcf7_da_tf <- FindMarkers(object = mcf7_atac_subset, ident.1 = "Topic11", ident.2="Topic17", test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="chromvar") 
  t47d_da_tf <- FindMarkers(object = t47d_atac_subset, ident.1 = "Topic11", ident.2="Topic17",  test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="chromvar") 

#set up readable motif names
  mcf7_da_tf$tf_name <- unlist(lapply(unlist(lapply(row.names(mcf7_da_tf), function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  t47d_da_tf$tf_name <- unlist(lapply(unlist(lapply(row.names(t47d_da_tf), function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))

#plot
  mcf7_da_tf$sig<-ifelse(mcf7_da_tf$p_val_adj<=0.05,"sig","non_sig")
  t47d_da_tf$sig<-ifelse(t47d_da_tf$p_val_adj<=0.05,"sig","non_sig")
  mcf7_da_tf$label<-""
  mcf7_da_tf[mcf7_da_tf$sig=="sig",]$label<-mcf7_da_tf[mcf7_da_tf$sig=="sig",]$tf_name
  t47d_da_tf$label<-""
  t47d_da_tf[t47d_da_tf$sig=="sig",]$label<-t47d_da_tf[t47d_da_tf$sig=="sig",]$tf_name

  plt<-ggplot(mcf7_da_tf,aes(x=avg_log2FC,y=-log10(p_val_adj),color=sig,label=label))+geom_point()+theme_bw()+scale_fill_manual(values=c("#999999", "#FF0000"))+geom_label_repel(label.size = 0.1,max.overlaps=20)+xlim(c(-10,10))
  ggsave(plt,file="mcf7_topic11v17_chromvar.pdf")
  system("slack -F mcf7_topic11v17_chromvar.pdf ryan_todo")
write.table(mcf7_da_tf,file="mcf7_topicbin_TF.txt",col.names=T,sep="\t")
system("slack -F mcf7_topicbin_TF.txt ryan_todo")

  plt<-ggplot(t47d_da_tf,aes(x=avg_log2FC,y=-log10(p_val_adj),color=sig,label=label))+geom_point()+theme_bw()+scale_fill_manual(values=c("#999999", "#FF0000"))+geom_label_repel(label.size = 0.1,max.overlaps=20)+xlim(c(-10,10))
  ggsave(plt,file="t47d_topic11v17_chromvar.pdf")
  system("slack -F t47d_topic11v17_chromvar.pdf ryan_todo")
write.table(t47d_da_tf,file="t47d_topicbin_TF.txt",col.names=T,sep="\t")
system("slack -F t47d_topicbin_TF.txt ryan_todo")

#Find differences in peaks
  mcf7_da_peaks <- FindMarkers(object = mcf7_atac_subset, ident.1 = "Topic11", ident.2="Topic17", test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="peaks") 
  closest_genes <- ClosestFeature(mcf7_atac_subset,row.names(mcf7_da_peaks))
  mcf7_da_peaks <-cbind(mcf7_da_peaks ,closest_genes)
write.table(mcf7_da_peaks,file="mcf7_topicbin_DApeaks.txt",col.names=T,sep="\t")
system("slack -F mcf7_topicbin_DApeaks.txt ryan_todo")

  t47d_da_peaks <- FindMarkers(object = t47d_atac_subset,  ident.1 = "Topic11", ident.2="Topic17", test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="peaks")           
  closest_genes <- ClosestFeature(t47d_atac_subset,row.names(t47d_da_peaks))
  t47d_da_peaks <-cbind(t47d_da_peaks ,closest_genes)          
write.table(t47d_da_peaks,file="t47d_topicbin_DApeaks.txt",col.names=T,sep="\t")
system("slack -F t47d_topicbin_DApeaks.txt ryan_todo")

#Find differences in cistrome
#########NOT WORKING, PCT 1 and PCT 2 are both 1 for some reason
  mcf7_da_cistrome <- FindMarkers(object = mcf7_atac_subset, ident.1 = "Topic11", ident.2="Topic17", test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="cistrome_score",logfc.threshold=0) 
  mcf7_da_cistrome$sig<-ifelse(mcf7_da_cistrome$p_val_adj<=0.01,"sig","non_sig")
  mcf7_da_cistrome$label<-""
  mcf7_da_cistrome[mcf7_da_cistrome$sig=="sig",]$label<-row.names(mcf7_da_cistrome[mcf7_da_cistrome$sig=="sig",])
write.table(mcf7_da_cistrome,file="mcf7_topicbin_cistrome.txt",col.names=T,sep="\t")
system("slack -F mcf7_topicbin_cistrome.txt ryan_todo")
  plt<-ggplot(mcf7_da_cistrome,aes(x=avg_log2FC,y=-log10(p_val_adj),color=sig,label=label))+geom_point()+theme_bw()+scale_fill_manual(values=c("#999999", "#FF0000"))+geom_label_repel(label.size = 0.1,max.overlaps=20))
  ggsave(plt,file="mcf7_topic11v17_cistrome.pdf")
  system("slack -F mcf7_topic11v17_cistrome.pdf ryan_todo")

  t47d_da_cistrome <- FindMarkers(object = t47d_atac_subset, ident.1 = "Topic11", ident.2="Topic17", test.use = 'LR', latent.vars = "nCount_peaks", only.pos=F, assay="cistrome_score",logfc.threshold=0) 
  t47d_da_cistrome$sig<-ifelse(t47d_da_cistrome$p_val_adj<=0.01,"sig","non_sig")
  t47d_da_cistrome$label<-""
  t47d_da_cistrome[t47d_da_cistrome$sig=="sig",]$label<-row.names(t47d_da_cistrome[t47d_da_cistrome$sig=="sig",])
write.table(t47d_da_cistrome,file="t47d_topicbin_cistrome.txt",col.names=T,sep="\t")
system("slack -F t47d_topicbin_cistrome.txt ryan_todo")
  plt<-ggplot(t47d_da_cistrome,aes(x=avg_log2FC,y=-log10(p_val_adj),color=sig,label=label))+geom_point()+theme_bw()+scale_fill_manual(values=c("#999999", "#FF0000"))+geom_label_repel(label.size = 0.1,max.overlaps=20))
  ggsave(plt,file="t47d_topic11v17_cistrome.pdf")
  system("slack -F t47d_topic11v17_cistrome.pdf ryan_todo")



#plot topic11 vs topic17 for cells
plt<-FeaturePlot(object = mcf7_atac, features = c("ImputedTopic_11", "ImputedTopic_17"),blend = T,col=c("white","red","blue"),reduction="umap",max.cutoff=10)
ggsave(plt,file="mcf7_topicimputation.umap.pdf",width=13)
system("slack -F mcf7_topicimputation.umap.pdf ryan_todo")

plt<-FeaturePlot(object = t47d_atac, features = c("ImputedTopic_11", "ImputedTopic_17"),blend = T,col=c("white","red","blue"),reduction="umap",max.cutoff=10)
ggsave(plt,file="t47d_topicimputation.umap.pdf",width=13)
system("slack -F t47d_topicimputation.umap.pdf ryan_todo")


#function for plotting footprints in data sets
plot_footprints<-function(x,footprints){
  x <- Footprint(object = x, motif.name = footprints, genome = BSgenome.Hsapiens.UCSC.hg38, in.peaks=T)   # gather the footprinting information for sets of motifs
  p2 <- PlotFootprint(x, features = footprints,label=F,normalization="divide")   # plot the footprint data for each group of cells, might want to change idents
  return(p2)
}

out1<-mclapply(c("ESR1","FOXM1","TP53","FOXA1","TEAD3"),FUN=function(z) plot_footprints(x=t47d_atac_subset,footprints=z),mc.cores=5) #plot and return 5 at a time
outname="210924_t47dsubset"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out1))
print(wrap_plots(out1) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))


########################Plot out top TF for each cluster###################

plot_tf_heatmap<-function(table_in="t47d_topicbin_TF.txt",obj_in=t47d_atac_subset,outname="t47d"){
  da_tf<-read.csv(file=table_in,head=T,sep="\t",row.names=NULL)
  da_tf<-da_tf[complete.cases(da_tf),]
  da_tf<-da_tf[da_tf$sig=="sig",]

  #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
  dat_tf<-as.data.frame(t(as.data.frame(obj_in[["chromvar"]]@data)))
  sum_tf<-split(dat_tf,obj_in$topic_bin) #group by rows to seurat clusters
  sum_tf<-lapply(sum_tf,function(x) apply(x,2,mean)) #take average across group
  sum_tf<-do.call("rbind",sum_tf) #condense to smaller data frame

  sum_tf<-scale(t(sum_tf))
  sum_tf<-sum_tf[,!endsWith(colnames(sum_tf),"NA")] #remove NA (doublet cells)

  sum_tf<-sum_tf[row.names(sum_tf) %in% unique(da_tf$row.names),]
  row.names(sum_tf)<-da_tf[match(row.names(sum_tf),da_tf$row.names,nomatch=0),]$tf_name
  sum_tf_plot<-na.omit(t(sum_tf))

  colfun=colorRamp2(quantile(unlist(sum_tf_plot), probs=c(0,0.5,1)),cividis(3))
  plt1<-Heatmap(sum_tf_plot,
      #cluster_rows=sum_da_dend,
      #left_annotation=side_ha,
      col=colfun,
      #bottom_annotation=bottom_ha,
      column_names_gp = gpar(fontsize = 14),
      row_names_gp=gpar(fontsize=14),
      column_names_rot=90
  )

  plt1<-draw(plt1)
  pdf(paste0(outname,"_tf.heatmap.pdf"),height=20,width=20)
  draw(plt1)
  dev.off()
  system(paste0("slack -F ",outname,"_tf.heatmap.pdf ryan_todo"))
  }

plot_tf_heatmap(table_in="t47d_topicbin_TF.txt",obj_in=t47d_atac_subset,outname="t47d")

plot_tf_heatmap(table_in="mcf7_topicbin_TF.txt",obj_in=mcf7_atac_subset,outname="mcf7")


```

## Subset DA Peaks by DE Genes
Using genes associated with topics (11 is FOXM1 and 17 is ESR1)

```R
library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(circlize)
library(viridis)

setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")

mcf7_topic_geneset<-read.table("Model_MCF7_E2_20T_CLR_5000Variable_M10_top50_genes_topics.txt",sep="\t",skip=1,header=T) #MCF7 gene topics, 2 is FOXM1 8 is ESR1
mcf7_foxm1_genes<-mcf7_topic_geneset$Topic_2
mcf7_esr1_genes<-mcf7_topic_geneset$Topic_8
t47d_topic_geneset<-read.table("Model_PEPE_T47D_20T_CLR_5000Variable_M10_top50_genes_topics.txt",sep="\t",header=T) #t47d gene topics, 11 is FOXM1 17 is ESR1
t47d_foxm1_genes<-t47d_topic_geneset$Topic_11
t47d_esr1_genes<-t47d_topic_geneset$Topic_17

mcf7_da_tf<-read.table(file="mcf7_topicbin_TF.txt",header=T,sep="\t")
t47d_da_tf<-read.table(file="t47d_topicbin_TF.txt",header=T,sep="\t")

mcf7_da<-read.table("mcf7_topicbin_DApeaks.txt",sep="\t",header=T)
t47d_da<-read.table("t47d_topicbin_DApeaks.txt",sep="\t",header=T)

mcf7_da_tf$topic_enrichment<-NA
mcf7_da_tf[which(mcf7_da_tf$tf_name %in% mcf7_foxm1_genes),]$topic_enrichment<-"FOXM1"
mcf7_da_tf[which(mcf7_da_tf$tf_name %in% mcf7_esr1_genes),]$topic_enrichment<-"ESR1"
t47d_da_tf$topic_enrichment<-NA
t47d_da_tf[which(t47d_da_tf$tf_name %in% t47d_foxm1_genes),]$topic_enrichment<-"FOXM1"
t47d_da_tf[which(t47d_da_tf$tf_name %in% t47d_esr1_genes),]$topic_enrichment<-"ESR1"

mcf7_da$topic_enrichment<-NA
mcf7_da[which(mcf7_da$gene_name %in% mcf7_foxm1_genes),]$topic_enrichment<-"FOXM1"
mcf7_da[which(mcf7_da$gene_name %in% mcf7_esr1_genes),]$topic_enrichment<-"ESR1"
t47d_da$topic_enrichment<-NA
t47d_da[which(t47d_da$gene_name %in% t47d_foxm1_genes),]$topic_enrichment<-"FOXM1"
t47d_da[which(t47d_da$gene_name %in% t47d_esr1_genes),]$topic_enrichment<-"ESR1"

write.table(t47d_da,file="t47d_topicbin_DApeaks.txt",col.names=T,sep="\t")
system("slack -F t47d_topicbin_DApeaks.txt ryan_todo")
write.table(mcf7_da,file="mcf7_topicbin_DApeaks.txt",col.names=T,sep="\t")
system("slack -F mcf7_topicbin_DApeaks.txt ryan_todo")

write.table(t47d_da_tf,file="t47d_topicbin_TF.txt",col.names=T,sep="\t")
system("slack -F t47d_topicbin_TF.txt ryan_todo")
write.table(mcf7_da_tf,file="mcf7_topicbin_TF.txt",col.names=T,sep="\t")
system("slack -F mcf7_topicbin_TF.txt ryan_todo")
```
