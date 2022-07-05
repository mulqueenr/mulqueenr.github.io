---
title: 10X Multiome Tumor Preliminary
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_primaryprelim/
category: CEDAR
---

Multiome processing for 10X multiome data on Primary Tumors (Preliminary Data)

## File Location
```bash
cd /home/groups/CEDAR/mulqueen/sequencing_runs
sftp mulqueen@nix.ohsu.edu
#enter password
get -r /data/EXP211227HM
get -r /data/EXP211228HM
```


## Reference data
Using chipseq bed files held on cistrome for accessibility analysis.

```bash
/home/groups/CEDAR/mulqueen/ref/cistrome #stored reference files here, downloaded from site at .tar.gz file
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor_full_QC.txt #has information on each download peaks files
/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor #has individual peak files
```
## Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc.

### RM Libraries

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-A2_1,EXP211227HM_RM_1,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-A2_2,EXP211227HM_RM_1,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-A2_3,EXP211227HM_RM_1,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-A2_4,EXP211227HM_RM_1,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_RM_1,Gene Expression""" > RM_1.csv
```

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-B2_2,EXP211227HM_RM_2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-B2_3,EXP211227HM_RM_2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-B2_4,EXP211227HM_RM_2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_RM_2,Gene Expression""" > RM_2.csv

```

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-C2_1,EXP211227HM_RM_3,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-C2_2,EXP211227HM_RM_3,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-C2_3,EXP211227HM_RM_3,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-C2_4,EXP211227HM_RM_3,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_RM_3,Gene Expression""" > RM_3.csv

```

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-D2_1,EXP211227HM_RM_4,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-D2_2,EXP211227HM_RM_4,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-D2_3,EXP211227HM_RM_4,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-D2_4,EXP211227HM_RM_4,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_RM_4,Gene Expression""" > RM_4.csv

```
## Run CellRanger-ARC
```bash
cd /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi #using this as analysis directory
```
### Slurm commands
Test run

```bash
 cellranger-arc testrun --id=tiny2
```
Run Cellranger per sample

```bash          
for i in `ls *csv`; do
  outname=${i::-4};
  cellranger-arc count --id=${outname} \
   --reference=/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
   --libraries=/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/${i} \
   --localcores=30 \
   --localmem=90 ; done &
```

## Ryan Primary Tissue Multiome Analysis

### Seurat Object Generation for Tumors
Performing seurat analysis following https://satijalab.org/signac/articles/pbmc_multiomic.html
Generate seuratobject for RM primary samples.

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# set up sample loop to load the RNA and ATAC data, save to seurat object
setupseurat<-function(i){
  setwd(paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",i,"/outs"))
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
lapply(c("RM_1","RM_2","RM_3","RM_4"),setupseurat)

#combine seurat objects and add metadata column for sample
rm1<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/RM_1/outs/RM_1.SeuratObject.rds"); rm1$sample<-"rm1"
rm2<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/RM_2/outs/RM_2.SeuratObject.rds"); rm2$sample<-"rm2"
rm3<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/RM_3/outs/RM_3.SeuratObject.rds"); rm3$sample<-"rm3"
rm4<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/RM_4/outs/RM_4.SeuratObject.rds"); rm4$sample<-"rm4"
dat <- merge(rm1, y = c(rm2,rm3,rm4), add.cell.ids = c("rm_1","rm_2","rm_3","rm_4"), project = "primary")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/rm_merged.SeuratObject.rds")

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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

dat<-readRDS("rm_merged.SeuratObject.rds")
dat
table(dat$sample)
# rm1  rm2  rm3  rm4
#1851 1463  876  931

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
my_cols = brewer.pal(4,"Spectral")
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
ggsave(plt,file="RM_multimodal.umap.pdf")
system("slack -F RM_multimodal.umap.pdf ryan_todo")
saveRDS(dat,file="rm_merged.SeuratObject.rds")

```

## Run cisTopic for ATAC Dimensionality Reduction

```R
#######CISTOPIC PROCESSING PER CELL LINE#################
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
library(RColorBrewer)
library(ggplot2)

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
  pdf(paste0(wd,"/",outname,".umap.pdf"),width=10)
  print(plt1)
  print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.rds"))
  }


#Running cistopic on just control experiment (FOR TITAN PAPER)
dat<-readRDS("rm_merged.SeuratObject.rds")
cistopic_generation(x=dat,name_out="rm_merged")
dat<-readRDS("rm_merged.SeuratObject.rds")


#set up colors for samples
my_cols = brewer.pal(4,"Spectral")
alpha_val=0.33

#add cistopic reduction to umap projections
p1<-DimPlot(dat,reduction="rna_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("RNA UMAP")
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("ATAC UMAP")
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("Multimodal UMAP (LSI)")

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
  dims.list = list(1:50, 1:ncol(dat@reductions$cistopic)), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
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
p5<-DimPlot(dat2,reduction="multimodal_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("Multimodal UMAP (Cistopic)")
p6<-DimPlot(dat,reduction="cistopic_umap",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("Cistopic UMAP")

#Finally Plot results
plt<-(p1 | p2)/(p3)/(p6 | p5)
ggsave(plt,file="RM_multimodal.umap.pdf")
system("slack -F RM_multimodal.umap.pdf ryan_todo")
saveRDS(dat,"rm_merged.SeuratObject.rds")

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

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

dat<-readRDS("rm_merged.SeuratObject.rds")
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
saveRDS(dat,file="rm_merged.SeuratObject.rds")


plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

plt<-(plt2|plt3)/plt1
ggsave(plt,file="rm.predictedid.umap.pdf",width=20,height=30,limitsize=F)
system("slack -F rm.predictedid.umap.pdf ryan_todo")



```

### Check Cell Type Assignment By Prediction Grouping
https://ars.els-cdn.com/content/image/1-s2.0-S1097276521007954-gr5.jpg

Continuing R Session

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

dat<-readRDS("rm_merged.SeuratObject.rds")

###LINKING PEAKS TO GENES###
DefaultAssay(dat) <- "peaks"

# first compute the GC content for each peak
dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
#Identify cell types in data
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
#filter geneset to those in data set
unlist(geneset)[which(!(unlist(geneset) %in% row.names(dat@assays$SCT@data)))]

dat <- LinkPeaks(
  object = dat,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = unlist(geneset)
)

Idents(dat)<-dat$predicted.id

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

DefaultAssay(dat) <- "SCT"
for (i in unique(names(geneset))){
  plt_list<-lapply(unlist(geneset[i]), function(x) cov_plots(dat=dat,gene_name=x))
  plt<-patchwork::wrap_plots(plt_list, ncol = 1)
  ggsave(plt,file=paste0("RM_",i,".featureplots.pdf"),height=4*length(plt_list),width=10,limitsize=F)
  system(paste0("slack -F ","RM_",i,".featureplots.pdf ryan_todo"))
}

saveRDS(dat,file="rm_merged.SeuratObject.rds")
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
library(dendsort)
library(patchwork)


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

####RUNNING INFERCNV#####
#https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
dat<-readRDS("rm_merged.SeuratObject.rds") #ensure it reads in properly
DefaultAssay(dat)<-"RNA"
dat$cnv_ref<-"FALSE"
dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells"),]$cnv_ref<-"TRUE" #set cnv ref by cell type

# get gene annotations for hg38
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
  seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead


#write out gene order list
  gene_order<-annotation[!duplicated(annotation$gene_name),]
  gene_order<-as.data.frame(gene_order[gene_order$gene_name %in% row.names(dat),])
  gene_order<-gene_order[c("gene_name","seqnames","start","end")]
  chrorder<-paste0("chr",c(1:22,"X","Y","M"))
  gene_order$seqnames<-factor(gene_order$seqnames,levels=chrorder) # set chr order
  gene_order<-with(gene_order, gene_order[order(seqnames, start),]) #order by chr and start position
  write.table(gene_order,file="RM_inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)


####RUNNING INFERCNV PER SAMPLE#####
for(i in unique(dat$sample)){
  dat_tmp<-subset(dat,sample==i)
  #https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
  counts=as.matrix(dat_tmp@assays$RNA@counts[,colnames(dat_tmp)])
  write.table(counts,file=paste0(i,"_inferCNV.counts.txt"),sep="\t",col.names=T,row.names=T,quote=F)
  cell_annotation=as.data.frame(cbind(row.names(dat_tmp@meta.data),dat_tmp@meta.data["cnv_ref"]))
  write.table(cell_annotation,file=paste0(i,"_inferCNV.annotation.txt"),sep="\t",col.names=F,row.names=F,quote=F)

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(i,"_inferCNV.counts.txt"),
                                      annotations_file=paste0(i,"_inferCNV.annotation.txt"),
                                      delim="\t",
                                      gene_order_file="RM_inferCNV.gene_order.txt",
                                      ref_group_names="TRUE")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("./",i,"_inferCNV"), 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               resume_mode=T)
  saveRDS(infercnv_obj,paste0("./",i,"_inferCNV","/",i,"_inferCNV.Rds"))

  infercnv_obj<-readRDS(paste0("./",i,"_inferCNV","/",i,"_inferCNV.Rds"))
  mat<-as.data.frame(t(infercnv_obj@expr.data))

  plt<-Heatmap(mat,
    column_order=1:ncol(mat),
    show_row_names=FALSE,
    show_column_names=FALSE)

  pdf(paste0(i,".infercnv.heatmap.pdf"))
  plt<-draw(plt)
  dev.off()
  system(paste0("slack -F ",i,".infercnv.heatmap.pdf ryan_todo"))
}


###Trying CASPER for CNV Profiling

library(CaSpER) 


casper_per_sample<-function(sample_name="rm1",bam_folder_name="RM_1"){
  system(paste0("mkdir ",sample_name,"_casper_test"))
  bam_location<-paste0("./",bam_folder_name,"/outs/gex_possorted_bam.bam")
  BAFExtract_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract"
  hg38_list_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list"
  hg38_folder_location<-"/home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/"
  baf_sample_directory<-paste0("./",sample_name,"_casper_test")
  dat_tmp<-subset(dat,sample==sample_name)
  control<-names(dat_tmp$cnv_ref == "TRUE")
  log.ge <- as.matrix(dat@assays$RNA@data)
  genes <- rownames(log.ge)
  annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
  log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
  rownames(log.ge) <- annotation$Gene
  log.ge <- log2(log.ge +1)

  system(paste0("samtools view ",bam_location," | ",BAFExtract_location," -generate_compressed_pileup_per_SAM stdin ",hg38_list_location," ",baf_sample_directory," 30 0; ",BAFExtract_location," -get_SNVs_per_pileup ",hg38_list_location," ",baf_sample_directory," ",hg38_folder_location," 1 1 0.1 ",baf_sample_directory,"/test.snp")) #this needs to be fixed, not ooutputting anything

  loh <- readBAFExtractOutput ( path=baf_sample_directory, sequencing.type="bulk") 
  names(loh) <- gsub(".snp", "", names(loh))
  load("maf.rda") ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
  loh<- list()
  loh[[1]] <- maf
  names(loh) <- i
  loh.name.mapping <- data.frame (loh.name= i , sample.name=colnames(log.ge))

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


  pdf("MM135.Distrubution.pdf")
  plot(density(as.vector(object@control.normalized[[3]])))
  plot(density(log2(object@control.normalized.noiseRemoved[[3]]+1)))
  dev.off()

  ## runCaSpER
  final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

  ## summarize large scale events 
  finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 

  obj <- final.objects[[9]]
  plotHeatmap10x(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

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
      labs(x = "", 
      y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
      deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + 
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

}







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


