---
title: 10X Multiome Tumor Test
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_primary/
category: CEDAR
---

Multiome processing for 10X multiome data

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

### Yahong Libraries
```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H10_1,EXP211227HM_si_Control_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H10_2,EXP211227HM_si_Control_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H10_3,EXP211227HM_si_Control_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H10_4,EXP211227HM_si_Control_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_si_Control_plus_E2,Gene Expression""" > YW_si_Control_plus_E2.csv
```
```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H11_1,EXP211227HM_si_KLF4_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H11_2,EXP211227HM_si_KLF4_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H11_3,EXP211227HM_si_KLF4_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H11_4,EXP211227HM_si_KLF4_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_si_KLF4_plus_Veh,Gene Expression""" > YW_si_KLF4_plus_Veh.csv
```

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H12_1,EXP211227HM_si_KLF4_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H12_2,EXP211227HM_si_KLF4_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H12_3,EXP211227HM_si_KLF4_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H12_4,EXP211227HM_si_KLF4_plus_E2,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_si_KLF4_plus_E2,Gene Expression""" > YW_si_KLF4_plus_E2.csv
```

```bash
echo """fastqs,sample,library_type
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H9_1,EXP211227HM_si_Control_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H9_2,EXP211227HM_si_Control_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H9_3,EXP211227HM_si_Control_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM/SI-NA-H9_4,EXP211227HM_si_Control_plus_Veh,Chromatin Accessibility
/home/groups/CEDAR/mulqueen/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM,EXP211228HM_si_Control_plus_Veh,Gene Expression""" > YW_si_Control_plus_Veh.csv
```
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

# Seurat Analysis for Tumors
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

## Yahong Cell Line Multiome Analysis

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
lapply(c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh","YW_si_KLF4_plus_E2","YW_si_KLF4_plus_Veh"),setupseurat)

#combine seurat objects and filter
yw1<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_E2/outs/YW_si_Control_plus_E2.SeuratObject.rds") ; yw1$sample<-"YW_si_Control_plus_E2"
yw2<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_Veh/outs/YW_si_Control_plus_Veh.SeuratObject.rds") ; yw2$sample<-"YW_si_Control_plus_Veh"
yw3<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_KLF4_plus_E2/outs/YW_si_KLF4_plus_E2.SeuratObject.rds") ; yw3$sample<-"YW_si_KLF4_plus_E2"
yw4<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_KLF4_plus_Veh/outs/YW_si_KLF4_plus_Veh.SeuratObject.rds") ; yw4$sample<-"YW_si_KLF4_plus_Veh"
dat <- merge(yw1, y = c(yw2,yw3,yw4), add.cell.ids = c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh","YW_si_KLF4_plus_E2","YW_si_KLF4_plus_Veh"), project = "YW")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/yw_merged.SeuratObject.rds")

dat<-readRDS("yw_merged.SeuratObject.rds") #ensure it reads in properly

####FILTERING AND QC####
  dat
  table(dat$sample)
  # YW_si_Control_plus_E2 YW_si_Control_plus_Veh     YW_si_KLF4_plus_E2
  #                  7357                   6435                   7727
  #   YW_si_KLF4_plus_Veh
  #                  6785
  # filter out low quality cells
  dat <- subset(
    x = dat,
    subset = nCount_ATAC < 100000 &
      nCount_RNA < 25000 &
      nCount_ATAC > 500 &
      nCount_RNA > 500 &
      nucleosome_signal < 2 &
      TSS.enrichment > 1
  )
  dat
  table(dat$sample)
  # YW_si_Control_plus_E2 YW_si_Control_plus_Veh     YW_si_KLF4_plus_E2
  #                  7030                   6023                   7552
  #   YW_si_KLF4_plus_Veh
  #                  6350
  #Not much filtered out, all passing QC

```
### Call peaks and Dimensionality Reduction
Continuing R session, now performing initial ATAC seq peak calling and RNA processing.

```R

####INITIAL PROCESSING####
  # call peaks using MACS2
  peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")

  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

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
  p1<-DimPlot(dat,reduction="rna_umap",group.by="sample")+ggtitle("RNA UMAP")

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
  p2<-DimPlot(dat,reduction="atac_umap",group.by="sample")+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(
    object = dat,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:40),
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
  p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
  ggsave(plt,file="YW_multimodal.umap.pdf",width=15)
  system("slack -F YW_multimodal.umap.pdf ryan_todo")
  saveRDS(dat,file="yw_merged.SeuratObject.rds")



plt<-VlnPlot(
    object = dat,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0,
    group.by="sample"
  )
  ggsave(plt,file="yw_merged.qc.pdf")
  system("slack -F yw_merged.qc.pdf ryan_todo")

```

### Cell line deconvolution
The data generated had a mixture of T47-D and MCF-7 Cell lines in each lane.

New R session, using CNV profiles from InferCNV to split out cell lines used in experiment.

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
library(ggdendro)
library(dendextend)
library(circlize)
library(dendsort)
library(patchwork)


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

####RUNNING INFERCNV#####
#https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
dat<-readRDS("yw_merged.SeuratObject.rds") #ensure it reads in properly

DefaultAssay(dat)<-"RNA"
#write out raw counts matrix
  counts=as.matrix(dat@assays$RNA@counts[,colnames(dat)])
  write.table(counts,file="YW_inferCNV.counts.txt",sep="\t",col.names=T,row.names=T,quote=F)
#write out cell annotation
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["sample"]))
  write.table(cell_annotation,file="YW_inferCNV.annotation.txt",sep="\t",col.names=F,row.names=F,quote=F)

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
  write.table(gene_order,file="YW_inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)

#run infercnv
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix="YW_inferCNV.counts.txt",
                                      annotations_file="YW_inferCNV.annotation.txt",
                                      delim="\t",
                                      gene_order_file="YW_inferCNV.gene_order.txt",
                                      ref_group_names=NULL)

saveRDS(infercnv_obj,file="YW_inferCNV.Rds")

#read in where inferCNV keeps hanging (step 14)
infercnv_obj<-readRDS("YW_inferCNV.Rds")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="./YW_inferCNV", 
                               denoise=TRUE,
                               HMM=FALSE,
                               num_threads=10,
                               cluster_references=FALSE,
                               resume_mode=FALSE
                               )
                              #other options (removed for now)
                                #k_obs_groups=2,
                               #cluster_by_groups=FALSE, 
                               #HMM_type="i6",
                               #
                               #output_format="pdf")


saveRDS(infercnv_obj,file="YW_inferCNV.Rds")

#Run function keeps hanging on step 14, just going to cluster it myself.

infercnv_obj<-readRDS("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_inferCNV/14_invert_log_transform.infercnv_obj")
mat<-as.data.frame(t(infercnv_obj@expr.data))
mat<-mat[,row.names(infercnv_obj@gene_order[infercnv_obj@gene_order$chr=="chr17",])] #limit to chr 17 
mat<-mat[sample(1:nrow(mat),5000),] #limit to 5000 cells, which will be enough to see cluster bias

plt<-Heatmap(mat,
  column_order=1:ncol(mat),
  show_row_names=FALSE,
  show_column_names=FALSE,
  km=2)

pdf("YW.infercnv.heatmap.pdf")
plt<-draw(plt)
dev.off()
system("slack -F YW.infercnv.heatmap.pdf ryan_todo")

celltypes<-row_order(plt) #get row order from heatmap, this includes the slice
metadata_cluster<-as.data.frame(cbind(cellID=rownames(mat),celltype=0))
metadata_cluster[unlist(celltypes[[1]]),]$celltype<-"1"
metadata_cluster[unlist(celltypes[[2]]),]$celltype<-"2"
table(metadata_cluster$celltype)
row.names(metadata_cluster)<-metadata_cluster$cellID
metadata_cluster$cnv_profile<-ifelse(metadata_cluster$celltype==1,"MCF7","T47D") #this is due to the mcf7 amplification we have seen in other data

cell_line<-setNames(metadata_cluster$cnv_profile,row.names(metadata_cluster))
dat<-AddMetaData(dat,metadata=cell_line,col.name="cnv_profile")#assign slices to cell names

#cluster with low resolution on peaks dataset
dat$peaks_cluster<-ifelse(as.numeric(dat@reductions$atac_umap@cell.embeddings[,1])<0,"T47D","MCF7") #it is clear from the atac data that one cluster is T47D and one is MCF7 based on chr17 profiles
table(dat$peaks_cluster)

plt1<-DimPlot(dat,split.by="cnv_profile",group.by="cnv_profile",na.value=NA)+ggtitle("CNV Profile Split")
plt2<-DimPlot(dat,group.by="sample")+ggtitle("ATAC Samples")
plt3<-DimPlot(dat,group.by="peaks_cluster")+ggtitle("Final Cell Type Assignment")

ggsave(plt1/(plt2|plt3),file="YW_celltype.pdf",width=10)
system("slack -F YW_celltype.pdf ryan_todo")

#finally save the results, and separate data files for MCF7 and T47D
saveRDS(dat,file="yw_merged.SeuratObject.rds")
dat_mcf7<-subset(dat,peaks_cluster=="MCF7")
dat_t47d<-subset(dat,peaks_cluster=="T47D")
saveRDS(dat_mcf7,file="yw_mcf7.SeuratObject.rds")
saveRDS(dat_t47d,file="yw_t47d.SeuratObject.rds")
```

### Cell line specific Cistopic Processing

Now that cell lines are split out, can re-run cistopic for biological interpretation.
Note: I'm also going to subset the data to just control groups, this is for the TITAN paper processing.

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
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("origin","seurat_clusters"))
  plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(wd,"/",outname,".umap.pdf"))
  print(plt1)
  print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.Rds"))
  }


#Running cistopic on just control experiment (FOR TITAN PAPER)
dat_mcf7<-readRDS(file="yw_mcf7.SeuratObject.rds")
dat_mcf7<-subset(dat_mcf7, sample %in% c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh"))
cistopic_generation(x=dat_mcf7,name_out="yw_mcf7.control")
dat_t47d<-readRDS(file="yw_t47d.SeuratObject.rds")
dat_t47d<-subset(dat_mcf7, sample %in% c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh"))
cistopic_generation(x=dat_t47d,name_out="yw_t47d.control")

```

### Cistopic Cell Line Cistrome Signal

Take cistopic data per cell line and perform bed file overlap using cistrome cell line specific ChIP-seq.

```R
#####Cistopic chip-seq peak overlaps##########
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
library(ComplexHeatmap)

#read in cistrome data for topic analysis
cistrome_db<-read.csv("/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor_full_QC.txt",sep="\t") #has information on each download peaks files
cistrome_db<-cistrome_db[cistrome_db$Cell_line %in% c("MCF-7","T47D"),] #limit to our cell types
cistrome_db<-cistrome_db[cistrome_db$PeaksFoldChangeAbove10>1000,] #set lower limit for peaks
cistrome_dir="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor" #has individual peak files
ChIP_Seq_signatures <- paste(cistrome_dir, list.files(cistrome_dir), sep='/') #grab all files in directory
signature_dcids<-unlist(lapply(strsplit(unlist(lapply(strsplit(ChIP_Seq_signatures,"human_factor/"),"[",2)),"_sort"),"[",1)) #grab dcID from file name
ChIP_Seq_signatures <- ChIP_Seq_signatures[signature_dcids %in% cistrome_db$DCid] #limit to those in cistrome_db after filter
ChIP_Seq_signatures<-ChIP_Seq_signatures[order(as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(ChIP_Seq_signatures,"human_factor/"),"[",2)),"_sort"),"[",1))))] #sort by numeric, like db
cistrome_db$label<-paste(cistrome_db$Cell_line,cistrome_db$Factor,cistrome_db$DCid,sep="_")
cistrome_db$file<-ChIP_Seq_signatures

hg38_chip<-paste(cistrome_dir, list.files(cistrome_dir,pattern=".bed"), sep='/') #grab all files in directory
hg38_chip<-hg38_chip[which(hg38_chip %in% cistrome_db$file)] #filter to our files of interest
labels<-cistrome_db[cistrome_db$file %in% hg38_chip,]$label #grab meaningful file labels


#run through further processing
#this function will perform more analysis on topics in the selected model, looking at 
#1. cell topic weights
#2. cistrome chipseq peak overlaps per topic
#3. topic annotations (if peaks fall primarily in promoter regions etc)

cistopic_processing<-function(x,y,outname){
  obj<-x
  cisTopicObject<-y
  cisTopicObject<-addCellMetadata(cisTopicObject, cell.data =obj@meta.data)
  #heatmap by topic
  cisTopicObject <- runUmap(cisTopicObject, target='cell')
  pdf(paste(outname,"celltopic_heatmap.pdf",sep="."))
  print(cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c("origin","seurat_clusters")))
  dev.off()
  system(paste("slack -F", paste(outname,"celltopic_heatmap.pdf",sep="."), "ryan_todo",sep=" "))

  #region score association
  pred.matrix <- predictiveDistribution(cisTopicObject)
  cisTopicObject <- getSignaturesRegions(cisTopicObject, hg38_chip, labels=labels, minOverlap = 0.01) #output new bed files to be read in and processed
  aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)   # Compute cell rankings
  cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.1*nrow(aucellRankings),nCores=1,plot=FALSE)   # Check signature enrichment in cells

  cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
  cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.999, plot=TRUE)


  pdf(paste(outname,"signature_heatmap.pdf",sep="."),height=40)
  signaturesHeatmap(cisTopicObject,row_names_gp = gpar(fontsize = 3))
  dev.off()
  system(paste("slack -F", paste(outname,"signature_heatmap.pdf",sep="."), "ryan_todo",sep=" "))

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')
  cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE) #cluster by regions
  saveRDS(cisTopicObject,paste0(outname,".CisTopicObject.Rds"))

  #save cistopic region data
  outdat<-lapply(cisTopicObject@binarized.cisTopics, function(x) {
    temp<-merge(x,cisTopicObject@region.data,by="row.names")
    temp$topic_assigned<-names(x)
    return(temp)})

  outdat<-data.table::rbindlist(outdat,use.names=FALSE)
  write.table(outdat,file=paste0(outname,".topic_region_assignment.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  system(paste("slack -F ",paste0(outname,".topic_region_assignment.txt")," ryan_todo" ))
  pdf(paste0(outname,".topic_annotations.pdf"))
  par(mfrow=c(1,1))
  signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
  plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
  dev.off()
  system(paste("slack -F ",paste0(outname,".topic_annotations.pdf")," ryan_todo" ))

  pdf(paste0(outname,".topic_scores.pdf"))
  par(mfrow=c(4,4))
  plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
  dev.off()
  system(paste("slack -F ",paste0(outname,".topic_scores.pdf")," ryan_todo" ))
}



#add cistopic function to all cells
dat_mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
dat_mcf7_cistopic<-readRDS(file="yw_mcf7.control.CisTopicObject.Rds")
cistopic_processing(x=dat_mcf7,y=dat_mcf7_cistopic,outname="yw_mcf7.control")

dat_t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")
dat_t47d_cistopic<-readRDS(file="yw_t47d.control.CisTopicObject.Rds")
cistopic_processing(x=dat_t47d,y=dat_t47d_cistopic,outname="yw_t47d.control")


```
### ChromVar Analysis of Cell Lines
Now run Chromvar on data for agnostic transcription factor motifs

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

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

#read in RDS file.
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")



# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
 peaks<-granges(mcf7[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  mcf7 <- SetAssayData(object = mcf7, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  mcf7 <- RegionStats(object = mcf7, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  mcf7 <- RunChromVAR( object = mcf7,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")

#now run t47d
 peaks<-granges(t47d[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  t47d <- SetAssayData(object = t47d, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  t47d <- RegionStats(object = t47d, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  t47d <- RunChromVAR( object = t47d,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(t47d,file="yw_mcf7.control.SeuratObject.rds")



```
### Motif Footprinting on cells 
Perform motif footprinting on cells.
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
  setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
DefaultAssay(mcf7)<-"peaks"
mcf7 <- AddMotifs(mcf7, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
DefaultAssay(t47d)<-"peaks"
t47d <- AddMotifs(t47d, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")
saveRDS(t47d,file="yw_t47d.control.SeuratObject.rds")

#function for plotting footprints in data sets

plot_footprints<-function(x,footprints){
  x <- Footprint(object = x, motif.name = footprints, genome = BSgenome.Hsapiens.UCSC.hg38, in.peaks=T)   # gather the footprinting information for sets of motifs
  p2 <- PlotFootprint(x, features = footprints,label=F)   # plot the footprint data for each group of cells, might want to change idents
  return(p2)
}

Idents(mcf7)<-mcf7$sample
out<-mclapply(c("GATA3","ESR1","FOXA1","CTCF"),FUN=function(z) plot_footprints(x=mcf7,footprints=z),mc.cores=5) #plot and return 5 at a time

outname="yw_mcf7.control"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out))
print(wrap_plots(out) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))

Idents(t47d)<-t47d$sample
out<-mclapply(c("GATA3","ESR1","FOXA1","CTCF"),FUN=function(z) plot_footprints(x=t47d,footprints=z),mc.cores=5) #plot and return 5 at a time

outname="yw_t47d.control"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out))
print(wrap_plots(out) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))

```


































<!--

## Initial Seurat Processing of Data

Following https://satijalab.org/signac/articles/pbmc_multiomic.html

RM data
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
p1<-DimPlot(dat,reduction="rna_umap",group.by="sample")+ggtitle("RNA UMAP")

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
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample")+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
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
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-DimPlot(dat,reduction="multimodal_umap",group.by="seurat_clusters")+ggtitle("Multimodal UMAP Clusters")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file="RM_multimodal.umap.pdf")
system("slack -F RM_multimodal.umap.pdf ryan_todo")
saveRDS(dat,file="rm_merged.SeuratObject.rds")

###LINKING PEAKS TO GENES###
DefaultAssay(dat) <- "peaks"

# first compute the GC content for each peak
dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
#Identify cell types in data
#Basing identification off of https://ars.els-cdn.com/content/image/1-s2.0-S1097276521007954-gr5.jpg
geneset<-list()
geneset[["luminal_Hr"]]<-c("ANKRD30A","ERBB4","SYTL2","INPP4B","AFF3")
geneset[["luminal_secretory"]]<-c("IGF2BP2","ALDH1A3","COBL","PIGR","CCL28")
geneset[["basal"]]<-c("NRG1","KRT14","ACTA2","PTPRT","MYLK")
geneset[["t_cells"]]<-c("PTPRC","PARP8","FYN","STAT4","AOAH")
geneset[["b_cells"]]<-c("MZB1","SSR4","HERPUD1","DERL3") #"ENAM" excluded
geneset[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74","HLA-DRB1","HLA-DPB1")
geneset[["endothelial"]]<-c("MCTP1","ADAMTS9","ZNF385D","ADGRL4","SELE") #ELTD1 is ADGRL4
geneset[["pericyte"]]<-c("C11orf95","RGS6","PRKG1","IGFBP5") #MT1A excluded
geneset[["fibroblast"]]<-c("DCN","APOD","HPSE2","CFD","LAMA2")

#filter geneset to those in data set
unlist(geneset)[which(!(unlist(geneset) %in% row.names(dat@assays$SCT@data)))]

dat <- LinkPeaks(
  object = dat,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = unlist(geneset)
)

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
    reduction="multimodal_umap")
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


####RUNNING INFERCNV#####
#https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
counts=as.matrix(dat@assays$RNA@counts[,colnames(dat)])
write.table(counts,file="RM_inferCNV.counts.txt",sep="\t",col.names=T,row.names=T,quote=F)
cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["sample"]))
write.table(cell_annotation,file="RM_inferCNV.annotation.txt",sep="\t",col.names=F,row.names=F,quote=F)
gene_order<-annotation[!duplicated(annotation$gene_name),]
gene_order<-as.data.frame(gene_order[gene_order$gene_name %in% row.names(dat),])
gene_order<-gene_order[c("gene_name","seqnames","start","end")]
write.table(gene_order,file="RM_inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="RM_inferCNV.counts.txt",
                                    annotations_file="RM_inferCNV.annotation.txt",
                                    delim="\t",
                                    gene_order_file="RM_inferCNV.gene_order.txt",
                                    ref_group_names=NULL)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=".", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

saveRDS(infercnv_obj,file="RM_inferCNV.Rds")
#####################Try CaSpER also##############################
library(CaSpER)

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
  cistopic_generation(x=dat,name_out="rm_merged")


###ChromVar###
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  set.seed(1234)
  library(BiocParallel)
  register(MulticoreParam(5))

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

dat<-readRDS("rm_merged.SeuratObject.rds")

  #Read in data and modify to monocle CDS file
  #read in RDS file.

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species =9606, all_versions = FALSE))

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
  saveRDS(dat,file="rm_merged.SeuratObject.rds")



```

## Cistopic on ATAC data

```R
setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")
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

#read in RNA data to subset same cells
mcf7_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/MCF7_RNA_strictQC.rds") #updated with new file loc
t47d_rna<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/T47D_RNA_strictQC.rds")
mcf7_rna<-RenameCells(mcf7_rna,new.names=paste0(row.names(mcf7_rna@meta.data),mcf7_rna$origin))
t47d_rna<-RenameCells(t47d_rna,new.names=paste0(row.names(t47d_rna@meta.data),t47d_rna$origin))

#read in pipeline generated atac data
mcf7<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/MCF7_justATAC_SO.rds") #generated by AD
t47d<-readRDS("/home/groups/CEDAR/doe/projects/ATACRNA/hg38/T47D_justATAC_SO.rds") #generated by AD

#renaming cells to fit atac data, have to add origin so cell names are unique
mcf7<-RenameCells(mcf7,new.names=paste0(unlist(lapply(strsplit(row.names(mcf7@meta.data),"_"),"[",1)),mcf7$origin))
t47d<-RenameCells(t47d,new.names=paste0(unlist(lapply(strsplit(row.names(t47d@meta.data),"_"),"[",1)),t47d$origin))

  peaks<-granges(combined[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
#subset cells
mcf7<-subset(x=mcf7,
  cells=row.names(mcf7_rna@meta.data),
  features=row.names(mcf7[["peaks"]])[which(seqnames(granges(mcf7[["peaks"]])) %in% c(paste0("chr",c(1:22,"X","Y"))))]
  )

t47d<-subset(x=t47d,
  cells=row.names(t47d_rna@meta.data),
  features=row.names(t47d[["peaks"]])[which(seqnames(granges(t47d[["peaks"]])) %in% c(paste0("chr",c(1:22,"X","Y"))))]
)

combined<- merge(mcf7, y = t47d, add.cell.ids = c("mcf7", "t47d"), project = "atac_rnasubsetted")

#a priori filter, not using to stay consistent with rna data
#combined<-subset(
#  x = combined,
#  subset = nCount_peaks > 3000 &
#  nFeature_peaks > 10000 &
#  TSS.enrichment > 4
#)

saveRDS(combined,"210924_cellline.SeuratObject.Rds")
saveRDS(mcf7,"210924_mcf7.SeuratObject.Rds")
saveRDS(t47d,"210924_t47d.SeuratObject.Rds")

cistopic_generation<-function(x,name_out,outdir="/home/groups/CEDAR/mulqueen/projects/10x_atacrna"){
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
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("origin","seurat_clusters"))
  plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(wd,"/",outname,".umap.pdf"))
  print(plt1)
  print(plt2)
  dev.off()
  system(paste0("slack -F ",paste0(wd,"/",outname,".umap.pdf")," ryan_todo"))
  saveRDS(atac_sub,paste0(wd,"/",outname,".SeuratObject.Rds"))
  }


#Running cistopic
  combined<-readRDS("210924_cellline.SeuratObject.Rds")
  cistopic_generation(x=combined,name_out="210924_cellline")

  #Run just for t47D
  t47d<-readRDS("210924_t47d.SeuratObject.Rds")
  cistopic_generation(x=t47d,name_out="210924_t47d")

  #Run just for mcf7
  mcf7<-readRDS("210924_mcf7.SeuratObject.Rds")
  cistopic_generation(x=mcf7,name_out="210924_mcf7")

```

Take cistopic data and perform bed file overlap using cistrome cell line specific ChIP-seq.


```R
setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")
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
library(ComplexHeatmap)

#read in cistrome data for topic analysis
cistrome_db<-read.csv("/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor_full_QC.txt",sep="\t") #has information on each download peaks files
cistrome_db<-cistrome_db[cistrome_db$Cell_line %in% c("MCF-7","T47D"),] #limit to our cell types
cistrome_db<-cistrome_db[cistrome_db$PeaksFoldChangeAbove10>1000,] #set lower limit for peaks
cistrome_dir="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor" #has individual peak files
ChIP_Seq_signatures <- paste(cistrome_dir, list.files(cistrome_dir), sep='/') #grab all files in directory
signature_dcids<-unlist(lapply(strsplit(unlist(lapply(strsplit(ChIP_Seq_signatures,"human_factor/"),"[",2)),"_sort"),"[",1)) #grab dcID from file name
ChIP_Seq_signatures <- ChIP_Seq_signatures[signature_dcids %in% cistrome_db$DCid] #limit to those in cistrome_db after filter
ChIP_Seq_signatures<-ChIP_Seq_signatures[order(as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(ChIP_Seq_signatures,"human_factor/"),"[",2)),"_sort"),"[",1))))] #sort by numeric, like db
cistrome_db$label<-paste(cistrome_db$Cell_line,cistrome_db$Factor,cistrome_db$DCid,sep="_")
cistrome_db$file<-ChIP_Seq_signatures

###No longer need liftover since files were aligned to hg38
  # #need to liftover files since cistrome is hg38 and this data is hg19, following https://www.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
  # path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  # ch = import.chain(path)

  # liftover_bed<-function(x){
  #   label<-cistrome_db[cistrome_db$file==x,]$label
  #   print(paste(x,label))
  #   x_hg38<-read.table(x,sep="\t")
  #   colnames(x_hg38)[1:3]<-c("chr","start","end")  
  #   x_hg38<-x_hg38[x_hg38$chr %in% paste0("chr",c(1:23,"X","Y")),]
  #   x_hg38<-makeGRangesFromDataFrame(x_hg38) #read in bed file and make granges
  #   seqlevelsStyle(x_hg38) = "UCSC"  # necessary
  #   x_hg19<-liftOver(x_hg38, ch)
  #   x_hg19<-unlist(x_hg19)
  #   genome(x_hg19) = "hg19"
  #   export.bed(object=x_hg19,con=paste(cistrome_dir,paste0(label,".hg19liftOver.bed"),sep="/"),format="bed")
  # }

  # mclapply(as.character(cistrome_db$file),liftover_bed,mc.cores=10) #parallelize chipseq to 20 cores, saving bed files after liftover with new name


hg38_chip<-paste(cistrome_dir, list.files(cistrome_dir,pattern=".bed"), sep='/') #grab all files in directory
hg38_chip<-hg38_chip[which(hg38_chip %in% cistrome_db$file)] #filter to our files of interest
labels<-cistrome_db[cistrome_db$file %in% hg38_chip,]$label #grab meaningful file labels


#run through further processing
#this function will perform more analysis on topics in the selected model, looking at 
#1. cell topic weights
#2. cistrome chipseq peak overlaps per topic
#3. topic annotations (if peaks fall primarily in promoter regions etc)

cistopic_processing<-function(x,y,outname){
  obj<-x
  cisTopicObject<-y
  cisTopicObject<-addCellMetadata(cisTopicObject, cell.data =obj@meta.data)
  #heatmap by topic
  cisTopicObject <- runUmap(cisTopicObject, target='cell')
  pdf(paste(outname,"celltopic_heatmap.pdf",sep="."))
  print(cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c("origin","seurat_clusters")))
  dev.off()
  system(paste("slack -F", paste(outname,"celltopic_heatmap.pdf",sep="."), "ryan_todo",sep=" "))

  #region score association
  pred.matrix <- predictiveDistribution(cisTopicObject)
  cisTopicObject <- getSignaturesRegions(cisTopicObject, hg38_chip, labels=labels, minOverlap = 0.01) #output new bed files to be read in and processed
  aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)   # Compute cell rankings
  cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.1*nrow(aucellRankings),nCores=1,plot=FALSE)   # Check signature enrichment in cells

  cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
  cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.999, plot=TRUE)


  pdf(paste(outname,"signature_heatmap.pdf",sep="."),height=40)
  signaturesHeatmap(cisTopicObject,row_names_gp = gpar(fontsize = 3))
  dev.off()
  system(paste("slack -F", paste(outname,"signature_heatmap.pdf",sep="."), "ryan_todo",sep=" "))

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')
  cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE) #cluster by regions
  saveRDS(cisTopicObject,paste0(outname,".CisTopicObject.Rds"))

  #save cistopic region data
  outdat<-lapply(cisTopicObject@binarized.cisTopics, function(x) {
    temp<-merge(x,cisTopicObject@region.data,by="row.names")
    temp$topic_assigned<-names(x)
    return(temp)})

  outdat<-data.table::rbindlist(outdat,use.names=FALSE)
  write.table(outdat,file=paste0(outname,".topic_region_assignment.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  system(paste("slack -F ",paste0(outname,".topic_region_assignment.txt")," ryan_todo" ))
  pdf(paste0(outname,".topic_annotations.pdf"))
  par(mfrow=c(1,1))
  signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
  plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
  dev.off()
  system(paste("slack -F ",paste0(outname,".topic_annotations.pdf")," ryan_todo" ))

  pdf(paste0(outname,".topic_scores.pdf"))
  par(mfrow=c(4,4))
  plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
  dev.off()
  system(paste("slack -F ",paste0(outname,".topic_scores.pdf")," ryan_todo" ))
}



#add cistopic function to all cells
  combined<-readRDS("210924_cellline.SeuratObject.Rds")
  combined_cisTopicObject<-readRDS("210924_cellline.CisTopicObject.Rds")
  cistopic_processing(x=combined,y=combined_cisTopicObject,outname="210924_cellline")

  t47d<-readRDS("210924_t47d.SeuratObject.Rds")
  t47d_cisTopicObject<-readRDS("210924_t47d.CisTopicObject.Rds")
  cistopic_processing(x=t47d,y=t47d_cisTopicObject,outname="210924_t47d")

  mcf7<-readRDS("210924_mcf7.SeuratObject.Rds")
  mcf7_cisTopicObject<-readRDS("210924_mcf7.CisTopicObject.Rds")
  cistopic_processing(x=mcf7,y=mcf7_cisTopicObject,outname="210924_mcf7")

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

setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")

combined<-readRDS("210924_cellline.SeuratObject.Rds")
mcf7<-readRDS("210924_mcf7.SeuratObject.Rds")
t47d<-readRDS("210924_t47d.SeuratObject.Rds")

  #Read in data and modify to monocle CDS file
  #read in RDS file.

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species =9606, all_versions = FALSE))

  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  # Create a new Mofif object to store the results
  #for combined
  peaks<-granges(combined[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  combined <- SetAssayData(object = combined, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  combined <- RegionStats(object = combined, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  combined <- RunChromVAR( object = combined,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(combined,file="210924_cellline.SeuratObject.Rds")

 peaks<-granges(mcf7[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  mcf7 <- SetAssayData(object = mcf7, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  mcf7 <- RegionStats(object = mcf7, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  mcf7 <- RunChromVAR( object = mcf7,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(mcf7,file="210924_mcf7.SeuratObject.Rds")

 peaks<-granges(t47d[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  t47d <- SetAssayData(object = t47d, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  t47d <- RegionStats(object = t47d, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  t47d <- RunChromVAR( object = t47d,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(t47d,file="210924_t47d.SeuratObject.Rds")



```

## Cicero

```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(cicero)
  library(SeuratObjects)
  library(EnsDb.Hsapiens.v86)
  setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")
  #Cicero processing function
  cicero_processing<-function(object_input,prefix){

      #Generate CDS format from Seurat object
      atac.cds <- as.CellDataSet(object_input,assay="peaks",reduction="umap")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = atac.cds@reducedDimS)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      
      # extract gene annotations from EnsDb
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevels(annotations)<-paste0("chr",seqlevels(annotations))
      #seqnames(annotations)<-paste0("chr",seqnames(annotations))

      # change to UCSC style since the data was mapped to hg38
      #seqlevelsStyle(annotations) <- 'UCSC'
      genome(annotations) <- "hg38"

      genome <- annotations@seqinfo # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = genome@seqnames, "length" = genome@seqlengths) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      
      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      
      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      DefaultAssay(object_input)<-"peaks"
      Links(object_input) <- links
      return(object_input)
  }

t47d<-readRDS("210924_t47d.SeuratObject.Rds")
t47d<-cicero_processing(object_input=t47d,prefix="210924_t47d")
saveRDS(t47d,"210924_t47d.SeuratObject.unnormGA.Rds")

mcf7<-readRDS("210924_mcf7.SeuratObject.Rds")
mcf7<-cicero_processing(object_input=mcf7,prefix="210924_mcf7")
saveRDS(mcf7,"210924_mcf7.SeuratObject.unnormGA.Rds")

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

#combined
  conns<-as.data.frame(readRDS("210924_cellline_atac_cicero_conns.Rds"))
  obj.cicero<-readRDS("210924_cellline_atac_cicero_cds.Rds")
  geneactivity_processing(cds_input=as.CellDataSet(combined,assay="peaks",reduction="umap"),conns_input=conns,prefix="210924_cellline")

  #Read in unnormalized GA
  cicero_gene_activities<-readRDS("210924_cellline.unnorm_GA.Rds")
  combined[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 

  # normalize
  combined <- NormalizeData(
    object = combined,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(combined$nCount_GeneActivity)
  )
  saveRDS(combined,"210924_cellline.SeuratObject.Rds")


#t47d
  conns<-as.data.frame(readRDS("210924_t47d_atac_cicero_conns.Rds"))
  obj.cicero<-readRDS("210924_t47d_atac_cicero_cds.Rds")
  geneactivity_processing(cds_input=as.CellDataSet(t47d,assay="peaks",reduction="umap"),conns_input=conns,prefix="210924_t47d")

  #Read in unnormalized GA
  cicero_gene_activities<-readRDS("210924_t47d.unnorm_GA.Rds")
  t47d[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 

  # normalize
  t47d <- NormalizeData(
    object = t47d,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(t47d$nCount_GeneActivity)
  )
  saveRDS(t47d,"210924_t47d.SeuratObject.Rds")


#mcf7
  conns<-as.data.frame(readRDS("210924_mcf7_atac_cicero_conns.Rds"))
  obj.cicero<-readRDS("210924_mcf7_atac_cicero_cds.Rds")
  geneactivity_processing(cds_input=as.CellDataSet(mcf7,assay="peaks",reduction="umap"),conns_input=conns,prefix="210924_mcf7")

  #Read in unnormalized GA
  cicero_gene_activities<-readRDS("210924_mcf7.unnorm_GA.Rds")
  mcf7[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 

  # normalize
  mcf7 <- NormalizeData(
    object = mcf7,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(mcf7$nCount_GeneActivity)
  )
  saveRDS(mcf7,"210924_mcf7.SeuratObject.Rds")

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

```R

## Get aggregate window of cistrome peak by atac data

Change this to get cellID as well to assign topic values
```bash
ref="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor/46096_sort_peaks.narrowPeak.bed"
in_1="/home/groups/CEDAR/doe/projects/cellRanger/MCF7_ATACRNA/MCF7_C_hg38/outs/atac_fragments.tsv.gz"
in_2="/home/groups/CEDAR/doe/projects/cellRanger/MCF7_ATACRNA/MCF7_E_hg38/outs/atac_fragments.tsv.gz"
in_3="/home/groups/CEDAR/doe/projects/cellRanger/T47D_ATACRNA/T47D_C_hg38/outs/atac_fragments.tsv.gz"
in_4="/home/groups/CEDAR/doe/projects/cellRanger/T47D_ATACRNA/T47D_E_hg38/outs/atac_fragments.tsv.gz"
#chr start end cellID count

#awk 'OFS="\t" {if ($1 ~ /^chr/) print $1,$2; else print "chr"$1,$2}' genome.fa.fai > hg38.genome
#also added chrMT to hg3.genome as a dummy line

genome="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/fasta/hg38.genome"

#get central point of each reference bed record and set to uniform size (the assumption here is that it is a TF that binds in the middle of the reported narrowPeaks)

#1. sort cistrome data
#2. get midpoint of bed file (supposedly the TF binding site)
#3. slop it to 100bp bin width (50bp on either side)
#4. intersect the new windows with cellranger tabix fragment files
#5. take the minimal distance (start or end of read, assuming it was aligned as PE) to midpoint of new 100bp windows
#6. Use sort to count instances of this. 


#FOXM1 39443, 39495, 46311, 57117

#ESR1   2294  2295  2297  2298  2301  2303  2304  2305  2725  6549  6550  6551 6552  6553  6555  6556  6557  6558  6560  6561  6562  6563  6564  6565 6577  6578  6579  6580  6581 33100 33101 33127 33135 33139 33145 33146 33156 33216 33221 33223 33224 33491 33504 33508 33513 33520 33525 33526 35044 35046 35051 36819 36825 36834 38474 44437 46369 46374 46375 48168 48169 48170 48171 48455 49611 49612 49654 52608 52609 52612 52613 52632 52633 52636 52637 53955 53956 53957 53958 53959 53960 53961 53962 53963 53986 53987 53988 53989 53990 53991 53992 53993 53994 54023 54025 54027 54029 54031 54033 54078 54082 54084 54086 54087 54088 54089 54621 54645 54646 54648 56904 56905 59373 59374 59375 59376 59377 59378 68056 68057 68844 68845 68846 68847 68872 68873 68875 68974 71065 71855 71856 71857 71858 71859 71863 71952 72992 72993 72998 72999 74120 74121 74140 74141 74144 74145 74157 74158 74159 74381 74383 74384 74394 74395 74396 74397 74398 74399 74400 74401 76100 76101 76102 76103 76104 76105 76106 76107 76108 81409 82124 82326 82327 82328 82330 82425 82427 83180 83182 83183 83186 83427 83428 83429 83430 83431 83432 83433 83434 83816 83817 83818 85954 86212 86213 86282 86283 87114 87115 87116 87200 87688 87948 88335 88336 88345 88349 88372

#GREB1 36802 36810 36811
tf_dcis[39443]=FOXM1
tf_dcis[39495]=FOXM1
tf_dcis[46311]=FOXM1
tf_dcis[57117]=FOXM1
tf_dcis[34995]=FOXM1
tf_dcis[2294]=ESR1
tf_dcis[2295]=ESR1
tf_dcis[2297]=ESR1
tf_dcis[2298]=ESR1
tf_dcis[2301]=ESR1
tf_dcis[36802]=GREB1
tf_dcis[36810]=GREB1
tf_dcis[36811]=GREB1

for key in "${!tf_dcis[@]}"; do
    echo "$key ${tf_dcis[$key]}"
done

for key in "${!tf_dcis[@]}"; do
for i in $in_1 $in_2 $in_3 $in_4; do
dcis=$key;
tf_name=${tf_dcis[$key]};
ref="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor/${dcis}_sort_peaks.narrowPeak.bed";
outname=${i:56:-27};
sort -T . --parallel 5 -k1,1 -k2,2n -k3,3n $ref | 
awk 'OFS="\t"{mid=int(($3+$4)/2); print $1,mid,mid+1}' | 
bedtools slop -i stdin -l 1000 -r 1000 -g $genome | 
bedtools intersect -a - -b $i -wb -wa | 
awk 'OFS="\t" {midpoint=int(($3+$2)/2);ins_1=($5-midpoint);ins_2=($6-midpoint);if (sqrt(ins_1^2) < sqrt(ins_2^2)) print ins_1,$7; else print ins_2,$7}' | sort -T . --parallel=5 -k1,1n | uniq -c > ${outname}.${tf_name}.${dcis}.cistrome.txt ; done ; done &

```

```R
library(ggplot2)
library(patchwork)
library(Signac)
library(dplyr)
setwd("/home/groups/CEDAR/mulqueen/projects/10x_atacrna")

mcf7_atac<-readRDS("210924_mcf7.SeuratObject.Rds")
t47d_atac<-readRDS("210924_t47d.SeuratObject.Rds")

#plot genomic regions # focus on DE genes between topics
Idents(mcf7_atac)<-mcf7_atac$topic_bin
plt<-CoveragePlot(mcf7_atac, region = c("GREB1","ESR1","FOXA1","FOXM1","CENPF","HMGB1"),links=T,idents=c("Topic11","Topic17"),show.bulk=T)
ggsave(plt,file="mcf7_covplots.pdf",width=10,height=10,limitsize=F)
system("slack -F mcf7_covplots.pdf ryan_todo")

Idents(t47d_atac)<-t47d_atac$topic_bin
plt<-CoveragePlot(t47d_atac, region = c("GREB1","ESR1","FOXA1","FOXM1","CENPF","HMGB1"),links=T,idents=c("Topic11","Topic17"),show.bulk=T)
ggsave(plt,file="t47d_covplots.pdf",width=10,height=10,limitsize=F)
system("slack -F t47d_covplots.pdf ryan_todo")

plot_cistrome<-function(obj_in=t47d_atac_subset,cellline="T47D"){
sum_factor<-as.data.frame(obj_in@meta.data %>% group_by(topic_bin) %>% summarize(sum_total=sum(nCount_peaks)))
sum_factor$sample<-c("Topic11","Topic17")
f_in<-list.files(path="/home/groups/CEDAR/mulqueen/projects/10x_atacrna",pattern="cistrome.txt$")
obj_in$cellnames<-substr(names(obj_in$origin),1,nchar(names(obj_in$origin))-nchar(obj_in$origin))

dat<-lapply(f_in,function(x) {
  tmp<-read.table(x)
  colnames(tmp)<-c("count","dist","cellID")
  tmp<-tmp[tmp$cellID %in% obj_in$cellnames,]
  tmp$dcis<-strsplit(x,"[.]")[[1]][3]
  tmp$tf<-strsplit(x,"[.]")[[1]][2]
  return(tmp)})

dat<-do.call("rbind",dat)
dat$sample<-obj_in@meta.data[match(dat$cellID,obj_in$cellnames),]$topic_bin
dat<-dat %>% group_by(dist,dcis,tf,sample) %>% summarize(sum=sum(count))
dat<-merge(dat,sum_factor,by="sample",all.x=T)
dat$norm_count<-dat$sum/dat$sum_total
plt<-ggplot(dat,aes(x=dist,y=norm_count,color=sample))+geom_line()+theme_minimal()+xlim(c(-1000,1000))+facet_grid(dcis+tf ~ . ,scales="free")
ggsave(plt,file=paste0(cellline,".cistrome.enrichment.pdf"),height=50,limitsize=F)
system(paste0("slack -F ",cellline,".cistrome.enrichment.pdf ryan_todo"))
}

plot_cistrome(obj_in=t47d_atac_subset,cellline="T47D")
plot_cistrome(obj_in=mcf7_atac_subset,cellline="MCF7")

```


-->