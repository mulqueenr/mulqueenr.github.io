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
Make full genome plot for TM

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
  #limit to just control samples
  counts<-counts[,startsWith(colnames(counts),"YW_si_Control")] 
  write.table(counts,file="YW_inferCNV.control.counts.txt",sep="\t",col.names=T,row.names=T,quote=F)
#write out cell annotation
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["sample"]))
  cell_annotation<-cell_annotation[startsWith(row.names(cell_annotation),"YW_si_Control"),] 
  write.table(cell_annotation,file="YW_inferCNV.annotation.control.txt",sep="\t",col.names=F,row.names=F,quote=F)

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
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix="YW_inferCNV.control.counts.txt",
                                      annotations_file="YW_inferCNV.annotation.control.txt",
                                      delim="\t",
                                      gene_order_file="YW_inferCNV.gene_order.txt",
                                      ref_group_names=NULL)

saveRDS(infercnv_obj,file="YW_inferCNV.control.Rds")

#read in where inferCNV keeps hanging (step 14)
infercnv_obj<-readRDS("YW_inferCNV.control.Rds")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="./YW_inferCNV.control", 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads=10,
                               cluster_references=FALSE,
                               resume_mode=FALSE,
                               k_obs_groups=2,
                               output_format="pdf"
                               )



saveRDS(infercnv_obj,file="YW_inferCNV.control.Rds")
infercnv_obj<-readRDS("YW_inferCNV.control.Rds")

dat_control<-subset(dat,sample %in% c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh"))
cluster_assignment<-read.table("./YW_inferCNV.control/infercnv.observation_groupings.txt")

cluster_assignment$cnv_profile<-ifelse(cluster_assignment$Dendrogram.Group==1,"MCF7","T47D") #this is due to the mcf7 amplification on chr17 and 20 we have seen in other data

cell_line<-setNames(cluster_assignment$cnv_profile,row.names(cluster_assignment))
dat_control<-AddMetaData(dat_control,metadata=cell_line,col.name="cnv_profile")#assign slices to cell names

#cluster with low resolution on peaks dataset
dat_control$peaks_cluster<-ifelse(as.numeric(dat_control@reductions$atac_umap@cell.embeddings[,1])<0,"T47D","MCF7") #it is clear from the atac data that one cluster is T47D and one is MCF7 based on chr17 profiles
table(dat_control$peaks_cluster)
table(paste(dat_control$peaks_cluster, dat_control$cnv_profile))

#MCF7 MCF7 MCF7 T47D T47D MCF7 T47D T47D
#     6262       176        26      6589
plt1<-DimPlot(dat_control,split.by="cnv_profile",group.by="cnv_profile",na.value=NA)+ggtitle("CNV Profile Split")
plt2<-DimPlot(dat_control,group.by="sample")+ggtitle("ATAC Samples")
plt3<-DimPlot(dat_control,group.by="peaks_cluster")+ggtitle("Final Cell Type Assignment")

ggsave(plt1/(plt2|plt3),file="YW_control.celltype.pdf",width=10)
system("slack -F YW_control.celltype.pdf ryan_todo")

```
## Yahong Data TITAN Analysis
### Cell line specific Cistopic Processing

Now that cell lines are split out, can re-run cistopic for biological interpretation.
Note: I'm also going to subset the data to just control groups, this is for the TITAN paper processing.
After cistopic models are made, the best is selected. We then associate signature scores (from chip-seq bed files) to each cistopic topic.

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
library(ComplexHeatmap)

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
dat_mcf7_full<-readRDS(file="yw_mcf7.SeuratObject.rds")
dat_mcf7<-subset(dat_mcf7_full, sample %in% c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh"))
cistopic_generation(x=dat_mcf7,name_out="yw_mcf7.control")

dat_t47d_full<-readRDS(file="yw_t47d.SeuratObject.rds")
dat_t47d<-subset(dat_t47d_full, sample %in% c("YW_si_Control_plus_E2","YW_si_Control_plus_Veh"))
cistopic_generation(x=dat_t47d,name_out="yw_t47d.control")

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



#Take cistopic data per cell line and perform bed file overlap using cistrome cell line specific ChIP-seq.

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
  print(cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c("sample","seurat_clusters")))
  dev.off()
  system(paste("slack -F", paste(outname,"celltopic_heatmap.pdf",sep="."), "ryan_todo",sep=" "))

  #region score association
  pred.matrix <- predictiveDistribution(cisTopicObject)
  cisTopicObject <- getSignaturesRegions(cisTopicObject, hg38_chip, labels=labels, minOverlap = 0.01) #output new bed files to be read in and processed
  aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)   # Compute cell rankings
  cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.1*nrow(aucellRankings),plot=FALSE)   # Check signature enrichment in cells

  cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
  cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.999, plot=TRUE)

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

  obj[['cistrome_score']] <- CreateAssayObject(
    data = t(cisTopicObject@cell.data[(row.names(cisTopicObject@cell.data) %in% row.names(obj@meta.data)),
      (startsWith(colnames(cisTopicObject@cell.data),"MCF") | startsWith(colnames(cisTopicObject@cell.data),"T47D"))]
    ))
  saveRDS(obj,file=paste0(outname,".SeuratObject.rds"))

}



#add cistopic function to all cells
dat_mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
dat_mcf7_cistopic<-readRDS(file="yw_mcf7.control.CisTopicObject.Rds")
cistopic_processing(x=dat_mcf7,y=dat_mcf7_cistopic,outname="yw_mcf7.control")
#x=dat_mcf7;y=dat_mcf7_cistopic;outname="yw_mcf7.control"

dat_t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds") 
dat_t47d_cistopic<-readRDS(file="yw_t47d.control.CisTopicObject.Rds") 
cistopic_processing(x=dat_t47d,y=dat_t47d_cistopic,outname="yw_t47d.control")
#x=dat_t47d;y=dat_t47d_cistopic;outname="yw_t47d.control"


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
register(SerialParam()) #using single core mode

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
 peaks<-granges(mcf7[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  mcf7 <- SetAssayData(object = mcf7, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  mcf7 <- RegionStats(object = mcf7, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  mcf7 <- RunChromVAR( object = mcf7,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")


#now run t47d
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")
 peaks<-granges(t47d[["peaks"]])
  peaks<-peaks[seqnames(peaks) %in% c(paste0("chr",c(1:22,"X","Y"))),]
    motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, pwm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
    motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, pwm = pfm)
  t47d <- SetAssayData(object = t47d, assay = 'peaks', slot = 'motifs', new.data = motif.hg38)
  t47d <- RegionStats(object = t47d, genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  t47d <- RunChromVAR( object = t47d,genome = BSgenome.Hsapiens.UCSC.hg38,assay="peaks")
  saveRDS(t47d,file="yw_t47d.control.SeuratObject.rds")



```

## Cicero Enhancer-Promoter Linkage
This is to generate enhancer promoter linkages at genes. 

```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  #library(monocle3)
  library(cicero)
  library(SeuratObjects)
  library(EnsDb.Hsapiens.v86)
  setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
  mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
  t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

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

t47d<-cicero_processing(object_input=t47d,prefix="yw_t47d.control")
saveRDS(t47d,"yw_t47d.control.SeuratObject.unnormGA.Rds")

mcf7<-cicero_processing(object_input=mcf7,prefix="yw_mcf7.control")
saveRDS(mcf7,"yw_mcf7.control.SeuratObject.unnormGA.Rds")

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

#t47d
   conns<-as.data.frame(readRDS("yw_t47d.control_atac_cicero_conns.Rds"))
   obj.cicero<-readRDS("yw_t47d.control_atac_cicero_cds.Rds")
   geneactivity_processing(cds_input=as.CellDataSet(t47d,assay="peaks",reduction="umap"),conns_input=conns,prefix="yw_t47d.control")

#mcf7
   conns<-as.data.frame(readRDS("yw_mcf7.control_atac_cicero_conns.Rds"))
   obj.cicero<-readRDS("yw_mcf7.control_atac_cicero_cds.Rds")
   geneactivity_processing(cds_input=as.CellDataSet(mcf7,assay="peaks",reduction="umap"),conns_input=conns,prefix="yw_mcf7.control")

#Read in unnormalized GA
cicero_gene_activities<-readRDS("yw_t47d.control.unnorm_GA.Rds")
t47d[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
# normalize
   t47d <- NormalizeData(
     object = t47d,
     assay = 'GeneActivity',
     normalization.method = 'LogNormalize',
     scale.factor = median(t47d$nCount_GeneActivity)
   )
   saveRDS(t47d,"yw_t47d.control.SeuratObject.rds")


#Read in unnormalized GA
cicero_gene_activities<-readRDS("yw_mcf7.control.unnorm_GA.Rds")
mcf7[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 
# normalize
   mcf7 <- NormalizeData(
     object = mcf7,
     assay = 'GeneActivity',
     normalization.method = 'LogNormalize',
     scale.factor = median(mcf7$nCount_GeneActivity)
   )
   saveRDS(mcf7,"yw_mcf7.control.SeuratObject.rds")



```

### Add TITAN Topics to Data

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


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
mcf7_titan<-readRDS("/home/groups/CEDAR/doe/projects/Lydia_ATACRNA/justMCF7_wTopics.rds")
t47d_titan<-readRDS("/home/groups/CEDAR/doe/projects/Lydia_ATACRNA/justT47D_wTopics.rds")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

#Add titan reduction to full seurat object
mcf7[["titan_lda"]]<-mcf7_titan[["lda"]]
t47d[["titan_lda"]]<-t47d_titan[["lda"]]

saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")
saveRDS(t47d,file="yw_t47d.control.SeuratObject.rds")

```
### Rerun Projections on Cell Lines

```R
library(Signac)
library(Seurat)
set.seed(1234)
library(stringr)
library(ggplot2)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

cluster_everything<-function(dat,outname){
  #RNA Processing
  DefaultAssay(dat) <- "RNA"
  dat <- SCTransform(dat)
  dat <- RunPCA(dat)
  dat<- RunUMAP(object = dat, reduction.name="rna_umap", reduction="pca", assay = "SCT", verbose = TRUE, dims=1:50 )
  #DNA Accessibility processing
  DefaultAssay(dat) <- "peaks"
  dat <- FindTopFeatures(dat, min.cutoff = 5)
  dat <- RunTFIDF(dat)
  dat <- RunSVD(dat)
  dat<- RunUMAP(object = dat, reduction.name="atac_umap", reduction="lsi", assay = "peaks", verbose = TRUE, dims=2:40 )
  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:40), modality.weight.name = c("SCT.weight","peaks.weight"), verbose = TRUE ) # build a joint neighbor graph using both assays
  dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_umap", assay = "RNA", verbose = TRUE ) # build a joint UMAP visualization
  
  #run the topics UMAPs too
  dat<- RunUMAP(object = dat, reduction.name="titan_umap", reduction="titan_lda", assay = "peaks", verbose = TRUE, dims=1:ncol(dat@reductions$titan_lda@cell.embeddings) )
  dat<- RunUMAP(object = dat, reduction.name="cistopic_umap", reduction="cistopic", assay = "peaks", verbose = TRUE, dims=1:ncol(dat@reductions$cistopic@cell.embeddings) ) 
  # build a joint neighbor graph using both topic reductions
  dat <- FindMultiModalNeighbors(object = dat, reduction.list = list("titan_lda", "cistopic"), dims.list = list(1:ncol(dat@reductions$titan_lda@cell.embeddings), 1:ncol(dat@reductions$cistopic@cell.embeddings)), modality.weight.name = c("SCT.weight","peaks.weight"), verbose = TRUE )
  dat <- RunUMAP(object = dat, nn.name = "weighted.nn", reduction.name="multimodal_topic_umap", assay = "peaks", verbose = TRUE ) # build a joint UMAP visualization

#Finally Plot results
p1<-DimPlot(dat,reduction="rna_umap",group.by="sample")+ggtitle("RNA UMAP")
p2<-DimPlot(dat,reduction="atac_umap",group.by="sample")+ggtitle("ATAC UMAP")
p3<-DimPlot(dat,reduction="multimodal_umap",group.by="sample")+ggtitle("Multimodal UMAP")

p4<-DimPlot(dat,reduction="titan_umap",group.by="sample")+ggtitle("TITAN UMAP")
p5<-DimPlot(dat,reduction="cistopic_umap",group.by="sample")+ggtitle("CISTOPIC UMAP")
p6<-DimPlot(dat,reduction="multimodal_topic_umap",group.by="sample")+ggtitle("Multimodal Topic UMAP")
plt<-(p1 | p2)/(p3 | p4)/(p5 | p6)

  ggsave(plt,file=paste0(outname,".umap.pdf"),width=10,height=10)
  system(paste0("slack -F ",outname,".umap.pdf ryan_todo"))
return(dat)
}


mcf7<-cluster_everything(dat=mcf7,outname="yw_mcf7.control")
t47d<-cluster_everything(dat=t47d,outname="yw_t47d.control")
saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")
saveRDS(t47d,file="yw_t47d.control.SeuratObject.rds")
```

## Perform topic comparisons
Using TITAN topics generated by Aaron Doe to compare to epigenetic topics (cistopic)

- Perform Correlation Analysis
- Plot linear model of differences with significant changes highlighted
- Bin data for pairwise DA/DE

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
library(cisTopic)
library(ComplexHeatmap)

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

#ESR1 and FOXM1 Topics identified by aysegul
#MCF7: FOXM1=topic4 , ESR1 = topic37
#T47D: FOXM1=topic11 , ESR1 = topic23
mcf7_foxm1_topic="lda_4";mcf7_esr1_topic="lda_37"
t47d_foxm1_topic="lda_11";t47d_esr1_topic="lda_23"

#Function to fit a linear regression to the relevant topics, and then bin data by quantiles
#maybe change binning to be based on ratio?

#Plot cistopic by cell matrix
dat_mcf7_cistopic<-readRDS(file="yw_mcf7.control.CisTopicObject.Rds")
dat_t47d_cistopic<-readRDS(file="yw_t47d.control.CisTopicObject.Rds")


make_heatmap<-function(seurat_object=mcf7,
  object=dat_mcf7_cistopic,
  method='Probability',
  colorBy=c('sample'),
  topic_order=mcf7_foxm1_topic,
  other_topic=mcf7_esr1_topic,
  outname){

    colorPal <- grDevices::colorRampPalette(c("white","red"))
    col_sample<-c("YW_si_Control_plus_Veh"="red","YW_si_Control_plus_E2"="blue")
    foxm1_list<-seurat_object@reductions$titan_lda@cell.embeddings[,topic_order]
    esr1_list<-seurat_object@reductions$titan_lda@cell.embeddings[,other_topic]
    col_foxm1<-colorRamp2(breaks=c(min(foxm1_list),max(foxm1_list)),colors=c("white","blue"))
    col_esr1<-colorRamp2(breaks=c(min(esr1_list),max(esr1_list)),colors=c("white","red"))
    #col_foxm1<-grDevices::colorRampPalette(c("white","black"))
    topic.mat <- modelMatSelection(object, "cell", method)
    rownames(topic.mat) <- paste("Topic", seq(1, nrow(topic.mat)))
    colnames(topic.mat) <- object@cell.names
    object.cell.data <- object@cell.data
    #set up ordering
    object.cell.data$sample_order<-ifelse(object.cell.data$sample=="YW_si_Control_plus_Veh",0,1)
    object.cell.data$foxm1_topic<-foxm1_list
    object.cell.data$esr1_topic<-esr1_list
    col_order<-with(object.cell.data,order(sample_order,foxm1_topic))
    object.cell.data<-object.cell.data[col_order,]
    topic.mat<-topic.mat[,col_order]


        heatmap <- ComplexHeatmap::Heatmap(topic.mat,
            col = colorPal(20), 
            top_annotation = HeatmapAnnotation(treatment = anno_simple(object.cell.data$sample,col=col_sample),
              foxm1_topic=anno_simple(object.cell.data$foxm1_topic,col=col_foxm1),
              esr1_topic=anno_simple(object.cell.data$esr1_topic,col=col_esr1)),
            column_split = object.cell.data[,"sample"],
            column_order=1:ncol(topic.mat),
            name = method,
            show_column_names = FALSE, 
            show_row_names = TRUE,
            column_title = "Topic contribution per cell", 
            column_title_gp = gpar(fontface = "bold"))
        pdf(outname)
        print(heatmap)
        dev.off()
        system(paste0("slack -F ",outname, " ryan_todo"))
}

make_heatmap(object=dat_mcf7_cistopic,
  topic_order=mcf7_foxm1_topic,
  other_topic=mcf7_esr1_topic,
  seurat_object=mcf7,
  outname="yw_mcf7_cellbyTopic_heatmap.pdf")

make_heatmap(object=dat_t47d_cistopic,
  topic_order=t47d_foxm1_topic,
  other_topic=t47d_esr1_topic,
  seurat_object=t47d,
  outname="yw_t47d_cellbyTopic_heatmap.pdf")




scatterplot_linearfit<-function(x,esr1_topic,foxm1_topic,outname){
  titan_dat<-as.data.frame(x@reductions$titan_lda@cell.embeddings) #get titan from y
  cor_mat<-titan_dat[,c(esr1_topic,foxm1_topic)]
  cor_mat$sample<-x$sample
  object<-x

  #generate a linear model
  x=as.numeric(cor_mat[,esr1_topic])
  y=as.numeric(cor_mat[,foxm1_topic])
  model_fit<-summary(lm(y~x,data=cor_mat))
  intercept<-model_fit$coefficients[1,1]
  slope<-model_fit$coefficients[2,1]
  std_err<-model_fit$coefficients[2,2]
  residuals<-model_fit$residuals
  print(model_fit)

  plt<-ggplot()+
    geom_point(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$sample),size=0.1)+
    theme_minimal()+
    geom_abline(intercept = intercept, slope = slope, col = "black")+
    ggtitle(paste(outname,"ESR1 and FOXM1 Topic Correlations"))+
    xlab(paste("ESR1",esr1_topic))+
    ylab(paste("FOXM1",foxm1_topic))+
    theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank())+
    scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

    ggsave(plt,file=paste0(outname,".ESR1_FOXM1.scatter.sample.pdf"),height=4,width=7,units="in")
    print(plt)
    dev.off()
    system(paste0("slack -F ",paste0(outname,".ESR1_FOXM1.scatter.sample.pdf")," ryan_todo"))

  #bins defined as top of topic distribution only in treated sample, to be used in later analysis
  cor_mat$topic_bin<-""
  #cor_mat[cor_mat$sample=="YW_si_Control_plus_E2" & cor_mat[,foxm1_topic]>quantile(cor_mat[,foxm1_topic],0.95) & cor_mat[,esr1_topic]<quantile(cor_mat[,esr1_topic],0.75),]$topic_bin<-"FOXM1" #was >0.95 and <0.75 for other, new is 0.75 and 0.5
  #cor_mat[cor_mat$sample=="YW_si_Control_plus_E2" &cor_mat[,esr1_topic]>quantile(cor_mat[,esr1_topic],0.95) & cor_mat[,foxm1_topic]<quantile(cor_mat[,foxm1_topic],0.75),]$topic_bin<-"ESR1" 
  cor_mat[cor_mat[,foxm1_topic]>quantile(cor_mat[,foxm1_topic],0.8) & cor_mat[,esr1_topic]<quantile(cor_mat[,esr1_topic],0.75),]$topic_bin<-"FOXM1" #was >0.95 and <0.75 for other, new is 0.75 and 0.5
  cor_mat[cor_mat[,esr1_topic]>quantile(cor_mat[,esr1_topic],0.8) & cor_mat[,foxm1_topic]<quantile(cor_mat[,foxm1_topic],0.75),]$topic_bin<-"ESR1" 
  print(table(cor_mat$topic_bin))


  plt<-ggplot()+
    geom_point(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$topic_bin),size=0.1)+
    theme_minimal()+
    scale_color_manual(values=c("NA"=NA,"FOXM1"="red","ESR1"="blue"))+
    geom_abline(intercept = intercept, slope = slope, col = "black")+
    ggtitle(paste(outname,"ESR1 and FOXM1 Topic Correlations"))+
    xlab(paste("ESR1",esr1_topic))+
    ylab(paste("FOXM1",foxm1_topic))+
    theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank())+
    scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

    ggsave(plt,file=paste0(outname,".ESR1_FOXM1.scatter.binned.pdf"),height=4,width=7,units="in")
    print(plt)
    dev.off()
    system(paste0("slack -F ",paste0(outname,".ESR1_FOXM1.scatter.binned.pdf")," ryan_todo"))

  bin_metadata<-setNames(as.character(cor_mat$topic_bin),rownames(cor_mat))
  object<-AddMetaData(object,metadata=bin_metadata,col.name="topic_bin")
  return(object)
}

mcf7<-scatterplot_linearfit(x=mcf7,outname="yw_mcf7.control.newbin",foxm1_topic=mcf7_foxm1_topic,esr1_topic=mcf7_esr1_topic) #p val 0.2161 #Adjusted R-squared:  8.235e-05
t47d<-scatterplot_linearfit(x=t47d,outname="yw_t47d.control.newbin",foxm1_topic=t47d_foxm1_topic,esr1_topic=t47d_esr1_topic) # pval p-value: < 2.2e-16 #Adjusted R-squared:  0.02529

#save now that there is titan data and its binned
saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.rds")
saveRDS(t47d,file="yw_t47d.control.SeuratObject.rds")

table(mcf7$topic_bin)
table(t47d$topic_bin)
Idents(mcf7)<-mcf7$topic_bin
#      ESR1 FOXM1
# 4159  1152  1127
Idents(t47d)<-t47d$topic_bin
#       ESR1 FOXM1
# 4797   914   904

find_markers_topic_bin<-function(x,assay_name,outname,logfc.threshold_set=0,min.pct_set=0,cistrome_cellline="MCF",latent.vars="NA"){
  if(!(latent.vars=="NA")){
  x_da<-FindMarkers(object=x,ident.1 = "ESR1", ident.2="FOXM1", test.use = 'LR', logfc.threshold=logfc.threshold_set,latent.vars = latent.vars, only.pos=F, assay=assay_name,min.pct=min.pct_set)
  } else {
  x_da<-FindMarkers(object=x,ident.1 = "ESR1", ident.2="FOXM1", test.use = 'LR', logfc.threshold=logfc.threshold_set, only.pos=F, assay=assay_name,min.pct=min.pct_set)
  }
  #if statements to handle formating the different modalities
  if(assay_name=="cistrome_score"){ #if assay is cistrome limit to correct cell line
    x_da<-x_da[startsWith(row.names(x_da),cistrome_cellline),] #subset to same cell line cistrome data
    label_length<-length(strsplit(row.names(x_da),"-")[[1]])
    row.names(x_da)<-paste(unlist(lapply(strsplit(row.names(x_da),"-"),"[",label_length-1)),
                          unlist(lapply(strsplit(row.names(x_da),"-"),"[",label_length)),sep="-")
    }
  if(assay_name=="chromvar"){ #make chromvar names readable
      print("Translating Chromvar Motif Names for Plot")
      x_da$out_name <- unlist(lapply(unlist(lapply(row.names(x_da), function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  } else{
    x_da$out_name<-row.names(x_da)
  }
  if(assay_name=="peaks"){ #assign peak output to closest features
    closest_genes <- ClosestFeature(x,row.names(x_da))
    x_da <-cbind(x_da ,closest_genes)
  }

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
    ggsave(plt,file=paste0(outname,"_",assay_name,"_ESR1vFOXM1bins.pdf"),width=2,units="in",height=3)
  system(paste0("slack -F ",outname,"_",assay_name,"_ESR1vFOXM1bins.pdf"," ryan_todo"))
  write.table(x_da,file=paste0(outname,"_",assay_name,"_ESR1vFOXM1bins.markers.txt"),col.names=T,sep="\t")
  system(paste0("slack -F ",paste0(outname,"_",assay_name,"_ESR1vFOXM1bins.markers.txt")," ryan_todo"))
}

#Run DA/DE based on titan topic bins
find_markers_topic_bin(x=mcf7,outname="yw_mcf7.control",assay_name="chromvar")
find_markers_topic_bin(x=mcf7,outname="yw_mcf7.control",assay_name="cistrome_score",cistrome_cellline="MCF")
find_markers_topic_bin(x=mcf7,outname="yw_mcf7.control",assay_name="SCT")

find_markers_topic_bin(x=t47d,outname="yw_t47d.control",assay_name="chromvar")
find_markers_topic_bin(x=t47d,outname="yw_t47d.control",assay_name="cistrome_score",cistrome_cellline="T47")
find_markers_topic_bin(x=t47d,outname="yw_t47d.control",assay_name="SCT")

find_markers_topic_bin(x=mcf7,outname="yw_mcf7.control",assay_name="peaks",latent.vars="atac_peak_region_fragments")
find_markers_topic_bin(x=t47d,outname="yw_t47d.control",assay_name="peaks",latent.vars="atac_peak_region_fragments")


#Plot blend of foxm1 and esr1 topics over the umap
plt<-FeaturePlot(object = mcf7, features = c(mcf7_foxm1_topic, mcf7_esr1_topic),blend = T,col=c("white","red","blue"),reduction="multimodal_topic_umap",order=T)
ggsave(plt,file="mcf7_titantopic.umap.pdf",width=13)
system("slack -F mcf7_titantopic.umap.pdf ryan_todo")

plt<-FeaturePlot(object = t47d, features = c(t47d_foxm1_topic, t47d_esr1_topic),blend = T,col=c("white","red","blue"),reduction="multimodal_topic_umap",order=T)
ggsave(plt,file="t47d_titantopic.umap.pdf",width=13)
system("slack -F t47d_titantopic.umap.pdf ryan_todo")

#Correlate cistopic with titan topics
correlate_cistopic_titan<-function(x,outname,esr1_topic,foxm1_topic){
  cistopic_dat<-x@reductions$cistopic@cell.embeddings #get cistopic from x
  titan_dat<-as.data.frame(x@reductions$titan_lda@cell.embeddings) #get titan from y
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

  #perform a weighted linear regression, TITAN topic is weight and cistopic topics are comparisons
  list_plots<-c()
  for (xi in colnames(cistopic_dat)){ 
    for (yi in colnames(cistopic_dat)){ 
      for (zj in c(foxm1_topic)){

      #generate a linear model
      x=as.numeric(cistopic_dat[,xi])
      y=as.numeric(cistopic_dat[,yi])
      x2_titan=as.numeric(titan_dat[,zj])
      model_fit<-summary(lm(y~x+x2_titan))
      intercept<-model_fit$coefficients[1,1]
      slope<-model_fit$coefficients[2,1]
      std_err<-model_fit$coefficients[2,2]
      residuals<-model_fit$residuals

    plt_dat<-data.frame(x=x,y=y,x2_titan=x2_titan)
    plt<-ggplot(plt_dat,aes(x=x,y=y,color=x2_titan))+
    geom_abline(intercept = intercept, slope = slope, col = "red")+
    geom_point(size=0.5,alpha=0.5)+
    scale_color_gradient2()+
    theme_void()+
    xlab(xi)+
    ylab(yi)
    list_plots[[paste(xi,yi)]]<-plt
    }
  }
}

  plt_list<-wrap_plots(list_plots,nrow=length(colnames(cistopic_dat)),guides="collect")+theme(legend.position = "none")
  ggsave(plt_list,file=paste0(outname,"_cistopic_scatterplot.png"),height=50,width=50,limitsize=FALSE)
  system(paste0("slack -F ",outname,"_cistopic_scatterplot.png ryan_todo"))

  out_plt<-Heatmap(cor_mat)
  columnorder<-column_order(out_plt)
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

  cor_mat_topics<-cor_mat[c(esr1_topic,foxm1_topic),]
  #cor_mat_topics<-rbind(cor_mat_topics,abs(cor_mat_topics[esr1_topic,]-cor_mat_topics[foxm1_topic,]))
  #row.names(cor_mat_topics)<-c(esr1_topic,foxm1_topic,"difference")
  cor_mat_topics<-rbind(abs(cor_mat_topics[esr1_topic,]-cor_mat_topics[foxm1_topic,]))
  row.names(cor_mat_topics)<-c("difference")
  col_fun = colorRamp2(c(0, 0.5), c("white", "red"))

  out_plt<-Heatmap(cor_mat_topics,
    col=col_fun,
    column_order=columnorder)
  pdf(paste0(outname,".cistopic_titan_topiccorrelation_specific.pdf"))
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistopic_titan_topiccorrelation_specific.pdf")," ryan_todo"))

}

#x is seurat object with redcutions for cistopic and titan_lda
correlate_cistopic_titan(x=mcf7,outname="yw_mcf7.control",foxm1_topic=mcf7_foxm1_topic,esr1_topic=mcf7_esr1_topic)
correlate_cistopic_titan(x=t47d,outname="yw_t47d.control",foxm1_topic=t47d_foxm1_topic,esr1_topic=t47d_esr1_topic)

#now to correlate chromvar motifs with titan topics

#Correlate chromvar with titan topics
correlate_chromvar_titan<-function(x,outname,cellline="test",quantile_cutoff=0.95,esr1_topic,foxm1_topic){
  chromvar_dat<-as.data.frame(t(x@assays$chromvar@data)) #get chromvar from x
  tfList <- getMatrixByID(JASPAR2020, ID=colnames(chromvar_dat)) #set up readable chromvar names
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  colnames(chromvar_dat)<-tfList
  titan_dat<-as.data.frame(x@reductions$titan_lda@cell.embeddings) #get titan from y  
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

  cor_mat<-as.data.frame(t(cor_mat[c(esr1_topic,foxm1_topic),]))
 
  #generate a linear model
  x=as.numeric(cor_mat[,esr1_topic])
  y=as.numeric(cor_mat[,foxm1_topic])
  model_fit<-summary(lm(y~x,data=cor_mat))
  intercept<-model_fit$coefficients[1,1]
  slope<-model_fit$coefficients[2,1]
  std_err<-model_fit$coefficients[2,2]
  residuals<-model_fit$residuals

  #highlight points which were determined significant
  assay_name="chromvar"
  x_da<-read.table(paste0(outname,"_",assay_name,"_ESR1vFOXM1bins.markers.txt"))
  cor_mat$color_fill<-ifelse(abs(residuals)>std_err,"red","not_red")

  ident_1_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC>0)),] ,n=20))#significant, is ident1 specific, is top 20
  ident_2_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC<0)),] ,n=20))#significant, is ident1 specific, is top 20
  x_da[c(ident_1_labels,ident_2_labels),]$label<-x_da[c(ident_1_labels,ident_2_labels),]$out_name

  cor_mat$label<-unlist(lapply(1:nrow(cor_mat),function(y) {ifelse(row.names(cor_mat)[y] %in% x_da$label,row.names(cor_mat)[y],"")})) #add label to significant features

  plt<-ggplot()+
  geom_point(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$color_fill),size=0.1)+
  scale_color_manual(values=c("red"="red","not_red"="grey50"))+
  theme_minimal()+
  geom_text_repel(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$color_fill,label=cor_mat$label),force=5,max.overlaps=Inf,size=2,min.segment.length = 0,box.padding=0.3,segment.size=0.1)+
  geom_abline(intercept = intercept, slope = slope, col = "black")+
  ggtitle(paste(outname,"chromvar v. titan correlations"))+
  xlab(paste("ESR1",esr1_topic))+
  ylab(paste("FOXM1",foxm1_topic))+
  theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank(),legend.position="none")+
  scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

  ggsave(plt,file=paste0(outname,".chromvar_titan_topiccorrelation.comparison.pdf"),height=4,width=4,units="in")
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".chromvar_titan_topiccorrelation.comparison.pdf")," ryan_todo"))
}

correlate_chromvar_titan(x=mcf7,outname="yw_mcf7.control",foxm1_topic=mcf7_foxm1_topic,esr1_topic=mcf7_esr1_topic)
correlate_chromvar_titan(x=t47d,outname="yw_t47d.control",foxm1_topic=t47d_foxm1_topic,esr1_topic=t47d_esr1_topic)

#Correlate cistrome data with titan topics
correlate_cistrome_titan<-function(x,outname,cellline="MCF",quantile_cutoff=0.95,esr1_topic,foxm1_topic){
  cistrome_dat<-as.data.frame(t(x@assays$cistrome_score@data)) #get cistrome from x
  titan_dat<-as.data.frame(x@reductions$titan_lda@cell.embeddings) #get titan from y  
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

  cor_mat<-cor_mat[,startsWith(colnames(cor_mat),cellline)] #subset to same cell line cistrome data
  out_plt<-Heatmap(t(cor_mat),
    row_names_gp = grid::gpar(fontsize = 4),
    row_km=5
    )
  pdf(paste0(outname,".cistrome_titan_topiccorrelation.pdf"),height=30)
  print(out_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistrome_titan_topiccorrelation.pdf")," ryan_todo"))


  cor_mat<-as.data.frame(t(cor_mat[c(esr1_topic,foxm1_topic),]))
  
  #generate a linear model
  x=as.numeric(cor_mat[,esr1_topic])
  y=as.numeric(cor_mat[,foxm1_topic])
  model_fit<-summary(lm(y~x,data=cor_mat))
  intercept<-model_fit$coefficients[1,1]
  slope<-model_fit$coefficients[2,1]
  std_err<-model_fit$coefficients[2,2]
  residuals<-model_fit$residuals

  #highlight points which were determined significant
  assay_name="cistrome_score"
  x_da<-read.table(paste0(outname,"_",assay_name,"_ESR1vFOXM1bins.markers.txt"))
  ident_1_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC>0)),] ,n=10))#significant, is ident1 specific, is top 20
  ident_2_labels<- row.names(head(x_da[which((x_da$sig=="sig") & (x_da$avg_log2FC<0)),] ,n=10))#significant, is ident1 specific, is top 20
  sig_labels<-row.names(x_da[which(x_da$sig=="sig"),])#significant, is ident1 specific, is top 20
  x_da[c(ident_1_labels,ident_2_labels),]$label<-x_da[c(ident_1_labels,ident_2_labels),]$out_name
  cor_mat$label<-unlist(lapply(1:nrow(cor_mat),function(y) {ifelse(row.names(cor_mat)[y] %in% x_da$label,row.names(cor_mat)[y],"")})) #add label to significant features
  cor_mat$color_fill<-unlist(lapply(1:nrow(cor_mat),function(y) {ifelse(row.names(cor_mat)[y] %in% sig_labels,"red","not_red")})) #add label to significant features

  plt<-ggplot()+
  geom_point(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$color_fill),size=0.1)+
  scale_color_manual(values=c("red"="red","not_red"="grey50"))+
  theme_minimal()+
  geom_text_repel(aes(x=cor_mat[,esr1_topic],y=cor_mat[,foxm1_topic],color=cor_mat$color_fill,label=cor_mat$label),force=5,max.overlaps=Inf,size=2,min.segment.length = 0,box.padding=0.3,segment.size=0.1)+
  geom_abline(intercept = intercept, slope = slope, col = "black")+
  ggtitle(paste(outname,"cistrome v. titan correlations"))+
  xlab(paste("ESR1",esr1_topic))+
  ylab(paste("FOXM1",foxm1_topic))+
  theme(axis.text = element_text(size = 6),axis.title=element_text(size=8),axis.ticks=element_blank(),legend.position="none")+
  scale_x_continuous(expand = expansion(mult = 0.5))+  scale_y_continuous(expand = expansion(mult = 0.5))

  ggsave(plt,file=paste0(outname,".cistrome_titan_topiccorrelation.comparison.pdf"),height=4,width=4,units="in")
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,".cistrome_titan_topiccorrelation.comparison.pdf")," ryan_todo"))
}

correlate_cistrome_titan(x=mcf7,outname="yw_mcf7.control",foxm1_topic=mcf7_foxm1_topic,esr1_topic=mcf7_esr1_topic,cellline="MCF")
correlate_cistrome_titan(x=t47d,outname="yw_t47d.control",foxm1_topic=t47d_foxm1_topic,esr1_topic=t47d_esr1_topic,cellline="T47")



plt1<-FeaturePlot(mcf7,features=c("topic_20"),reduction="multimodal_topic_umap",cols=c("white","blue"))
plt2<-FeaturePlot(mcf7,features=c("topic_13"),reduction="multimodal_topic_umap",cols=c("white","red"))
ggsave(plt1|plt2,file="yw_mcf7.control.cistopic_differences.pdf",width=10)
system("slack -F yw_mcf7.control.cistopic_differences.pdf ryan_todo")


plt1<-FeaturePlot(t47d,features=c("topic_3"),reduction="multimodal_topic_umap",cols=c("white","blue"))
plt2<-FeaturePlot(t47d,features=c("topic_2"),reduction="multimodal_topic_umap",cols=c("white","red"))
ggsave(plt1|plt2,file="yw_t47d.control.cistopic_differences.pdf",width=10)
system("slack -F yw_t47d.control.cistopic_differences.pdf ryan_todo")

```

### De Novo Motif Analysis on Topics
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
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(cisTopic)
library(org.Hs.eg.db)
library(plyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
mcf7_cistopic<-readRDS(file="yw_mcf7.control.CisTopicObject.Rds")

t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")
t47d_cistopic<-readRDS(file="yw_t47d.control.CisTopicObject.Rds") 

# Helper functions

.doGREAT <- function(coord, genome="hg38", fold_enrichment=1, geneHits=1, sign=0.1, request_interval=10 ) {
  coord <- coord[!duplicated(names(coord)),]
  coord <- sortSeqlevels(coord)
  coord <- sort(coord)
  job <- submitGreatJob(coord, species=genome, request_interval = request_interval)
  tb <- getEnrichmentTables(job, ontology=availableOntologies(job), request_interval = request_interval)
  for (i in 1:length(tb)){
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Binom_Fold_Enrichment'] < fold_enrichment),]
    if (length(which(tb[[i]][,'Hyper_Observed_Gene_Hits'] < geneHits)) > 1){
      tb[[i]] <- tb[[i]][-which(tb[[i]][,'Hyper_Observed_Gene_Hits'] < geneHits),]
    }
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Binom_Adjp_BH'] > sign),]
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Hyper_Adjp_BH'] > sign),]
  }
  tb <- tb[sapply(tb, function(x) dim(x)[1]) > 0]
  return(tb)
}


motif_generation<-function(x,cistopic_x,topic_name){
  topic<-rownames(cistopic_x@binarized.cisTopics[[topic_name]])
  topic_feats<-paste(
    unlist(lapply(strsplit(topic,"[:-]"),"[",1)),
    unlist(lapply(strsplit(topic,"[:-]"),"[",2)),
    unlist(lapply(strsplit(topic,"[:-]"),"[",3)),
    sep="-")

  enriched.motifs <- FindMotifs(object = x,features = topic_feats)
  plt<-MotifPlot(object = x,motifs = head(rownames(enriched.motifs)))
  return(plt)
}


cistopic_continued_processing<-function(x,outname,esr1_topic,foxm1_topic,seurat_object){
  #write bed files of binarized cisTopic regions
  getBigwigFiles(x, path=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",outname,"_cistopic"), seqlengths=seqlengths(txdb))
  getBedFiles(x, path=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",outname,"_cistopic"))
  #add annotations and plot
  x <- annotateRegions(x, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')
  plt<-signaturesHeatmap(x, selected.signatures = 'annotation')
  pdf(paste0(outname,"_cistopic.annotation.heatmap.pdf"))
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,"_cistopic.annotation.heatmap.pdf")," ryan_todo"))
  #use GREAT for gene ontology
  #x <- GREAT(x, genome='hg38', fold_enrichment=1, geneHits=1, sign=0.1, request_interval=5,species="hg38")
  #have to use modified cistopic function to use hg38
  object.binarized.rGREAT <- llply(1:length(x@binarized.cisTopics),
            function(i) .doGREAT(x@region.ranges[rownames(x@binarized.cisTopics[[i]])]))

  names(object.binarized.rGREAT) <- names(x@binarized.cisTopics)
  x@binarized.rGREAT <- object.binarized.rGREAT

  plt<-ontologyDotPlot(x, top=5, topics=1:length(x@binarized.cisTopics), var.y='name', order.by='Binom_Adjp_BH')
  pdf(paste0(outname,"_cistopic.Gontology.plot.pdf"))
  print(plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,"_cistopic.Gontology.plot.pdf")," ryan_todo"))

  esr1_plt<-motif_generation(x=seurat_object,cistopic_x=x,topic_name=esr1_topic)
  pdf(paste0(outname,"_cistopic.motifs.",esr1_topic,".esr1.pdf"))
  print(esr1_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,"_cistopic.motifs.",esr1_topic,".esr1.pdf")," ryan_todo"))


  foxm1_plt<-motif_generation(x=seurat_object,cistopic_x=x,topic_name=foxm1_topic)
  pdf(paste0(outname,"_cistopic.motifs.",foxm1_topic,".foxm1.pdf"))
  print(foxm1_plt)
  dev.off()
  system(paste0("slack -F ",paste0(outname,"_cistopic.motifs.",foxm1_topic,".foxm1.pdf")," ryan_todo"))

  return(x)
}

mcf7_cistopic<-cistopic_continued_processing(seurat_object=mcf7,x=mcf7_cistopic,outname="yw_mcf7.control",esr1_topic="Topic20",foxm1_topic="Topic13")
t47d_cistopic<-cistopic_continued_processing(seurat_object=t47d,x=t47d_cistopic,outname="yw_t47d.control",esr1_topic="Topic3",foxm1_topic="Topic2")


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
  library(ggplot2)
  library(motifmatchr)
  library(chromVAR)

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

Idents(mcf7)<-mcf7$topic_bin
out<-lapply(c("GATA3","ESR1","FOXA1","CTCF"),FUN=function(z) plot_footprints(x=mcf7,footprints=z)) #plot and return 5 at a time
outname="yw_mcf7.control"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out))
print(wrap_plots(out) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))

Idents(t47d)<-t47d$topic_bin
out<-lapply(c("GATA3","ESR1","FOXA1","KLF4","CTCF"),FUN=function(z) plot_footprints(x=t47d,footprints=z)) #plot and return 5 at a time
outname="yw_t47d.control"
pdf(paste0(outname,".motif_footprints.pdf"),height=5*length(out))
print(wrap_plots(out) + patchwork::plot_layout(ncol = 1))
dev.off()
system(paste("slack -F ",paste0(outname,".motif_footprints.pdf")," ryan_todo" ))

```
### Marker Plots
Set up FOXM1 PWM based on Homer data
```bash
wget http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif120.motif
```

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
  library(ggplot2)
  library(motifmatchr)
  library(chromVAR)
  library(universalmotif)
  library(EnsDb.Hsapiens.v86)
  library(cicero,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") 
  library(monocle3,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") 
  library(SeuratWrappers)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")

#Run cicero for binned cells
  #Cicero processing function
  cicero_processing<-function(object_input,prefix){

      #Generate CDS format from Seurat object
      atac.cds <- as.cell_data_set(object_input,assay="peaks")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = atac.cds@int_colData$reducedDims$MULTIMODAL_UMAP)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      
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

t47d_foxm1<-subset(t47d,topic_bin=="FOXM1")
t47d_foxm1<-cicero_processing(object_input=t47d_foxm1,prefix="yw_t47d.control.FOXM1")
saveRDS(t47d_foxm1,"yw_t47d.control.SeuratObject.unnormGA.FOXM1.Rds")
t47d_foxm1<-readRDS("yw_t47d.control.SeuratObject.unnormGA.FOXM1.Rds")

t47d_esr1<-subset(t47d,topic_bin=="ESR1")
t47d_esr1<-cicero_processing(object_input=t47d_esr1,prefix="yw_t47d.control.ESR1")
saveRDS(t47d_esr1,"yw_t47d.control.SeuratObject.unnormGA.ESR1.Rds")
t47d_esr1<-readRDS("yw_t47d.control.SeuratObject.unnormGA.ESR1.Rds")

mcf7_foxm1<-subset(mcf7,topic_bin=="FOXM1")
mcf7_foxm1<-cicero_processing(object_input=mcf7_foxm1,prefix="yw_mcf7.control.FOXM1")
saveRDS(mcf7_foxm1,"yw_mcf7.control.SeuratObject.unnormGA.FOXM1.Rds")
mcf7_foxm1<-readRDS("yw_mcf7.control.SeuratObject.unnormGA.FOXM1.Rds")

mcf7_esr1<-subset(mcf7,topic_bin=="ESR1")
mcf7_esr1<-cicero_processing(object_input=mcf7_esr1,prefix="yw_mcf7.control.ESR1")
saveRDS(mcf7_esr1,"yw_mcf7.control.SeuratObject.unnormGA.ESR1.Rds")
mcf7_esr1<-readRDS("yw_mcf7.control.SeuratObject.unnormGA.ESR1.Rds")

```
Plot genome tracks

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
  library(ggplot2)
  library(motifmatchr)
  library(chromVAR)
  library(universalmotif)
  library(EnsDb.Hsapiens.v86)
  library(cicero,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") 
  library(monocle3,lib.loc="/home/groups/oroaklab/nishida/R_4.0.0_arsn") 
  library(SeuratWrappers)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")
mcf7$topic_bin<-factor(mcf7$topic_bin,levels=c("ESR1","FOXM1"))
t47d$topic_bin<-factor(t47d$topic_bin,levels=c("ESR1","FOXM1"))
t47d_foxm1<-readRDS("yw_t47d.control.SeuratObject.unnormGA.FOXM1.Rds")
t47d_esr1<-readRDS("yw_t47d.control.SeuratObject.unnormGA.ESR1.Rds")
mcf7_foxm1<-readRDS("yw_mcf7.control.SeuratObject.unnormGA.FOXM1.Rds")
mcf7_esr1<-readRDS("yw_mcf7.control.SeuratObject.unnormGA.ESR1.Rds")

#use link plot to generate cicero links

#Generate Coverage Plots Across Genes
#Pick DE genes from the topic_binned data
#Pull gene list from DE genes, and RIME paper https://ars.els-cdn.com/content/image/1-s2.0-S221112471300017X-figs3.jpg
geneset<-c()
#MA0112.3 is ESR1 can add more
#
#geneset<-c("GREB1","CUX2","CENPE","PGR","NEAT1","ERBB4","TOP2A","KIF14","FMN1","ABCC12","ABCC11","TMEM164","CENPF","AURKA","FMN1","PIK3R3","GREB1","STMN1")
geneset<-c("PGR","CENPF")
#download foxm1 motif from homer data base
system("wget http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif120.motif")
foxm1<-read_matrix(file="motif120.motif",header=">",positions="rows")
foxm1<-convert_motifs(foxm1,class="TFBSTools-PWMatrix")
foxm1@ID<-"foxm1"


mcf7_peaks<-read.table("yw_mcf7.control_peaks_ESR1vFOXM1bins.markers.txt")
mcf7_rna<-read.table("yw_mcf7.control_SCT_ESR1vFOXM1bins.markers.txt")
t47d_peaks<-read.table("yw_t47d.control_peaks_ESR1vFOXM1bins.markers.txt")
t47d_rna<-read.table("yw_t47d.control_SCT_ESR1vFOXM1bins.markers.txt")

mcf7_peaks[mcf7_peaks$gene_name=="PGR",]

t47d_peaks[t47d_peaks$gene_name=="PGR",]


#####Using a custom Plotting function###############
#PGR chr11-101020000-101140000
#CENPE chr1-214600000-214660000


#####################Old style of plotting######################
cov_plots<-function(dat=mcf7,gene_range="chr11-101020000-101140000",gene_name="PGR",motif.name=c("MA0112.3"),add.foxm1=TRUE,outname="test",foxm1_link=mcf7_foxm1,esr1_link=mcf7_esr1,ymax=50){
  
  peak_annot<-granges(dat@assays$peaks) #set up peak motif overlap as annotation track
  peak_annot$motif_overlap<-""
  for (i in length(motif.name)){ #jarspar formatted names
    overlap_motif_idx<-findOverlaps(granges(peak_annot),dat@assays$peaks@motifs@positions[motif.name[i]])@from
    peak_annot[overlap_motif_idx,]$motif_overlap<-paste(peak_annot[overlap_motif_idx,]$motif_overlap,motif.name[i])
  }
  if(add.foxm1){ #add foxm1 (not included in Jaspar but used in homer)
    motif.positions <- motifmatchr::matchMotifs(pwms = foxm1, subject = granges(dat@assays$peaks), out = 'positions', genome = BSgenome.Hsapiens.UCSC.hg38 )
    overlap_motif_idx<-findOverlaps(granges(peak_annot),motif.positions[[1]])@from
    peak_annot[overlap_motif_idx,]$motif_overlap<-paste(peak_annot[overlap_motif_idx,]$motif_overlap,"foxm1")
  }

  dat@assays$peaks@meta.features$motif_overlap<-peak_annot$motif_overlap
  
  annot_plot<-AnnotationPlot(object=dat, region=gene_range)
  peak_plot<-PeakPlot(object=dat,region=gene_range,group.by="motif_overlap")

  esr1_cov <- CoveragePlot(object = dat, region = gene_range, assay="peaks", ident=c("ESR1"),annotation=FALSE,peaks=FALSE,links=FALSE ,ymax=ymax) 
  esr1_links<-LinkPlot(object=esr1_link, region=gene_range,min.cutoff=0.1)+scale_colour_gradient(low="white",high="red",limits=c(0,0.5))
  foxm1_cov <- CoveragePlot(object = dat, region = gene_range, assay="peaks", ident=c("FOXM1"),annotation=FALSE,peaks=FALSE,links=FALSE ,ymax=ymax ) 
  foxm1_links<-LinkPlot(object=foxm1_link, region=gene_range,min.cutoff=0.1)+scale_colour_gradient(low="white",high="red",limits=c(0,0.5))

  esr1_expr_plot<-ExpressionPlot(object=dat, features=gene_name, assay="SCT", ident=c("ESR1"))+xlim(c(0,5))
  foxm1_expr_plot<-ExpressionPlot(object=dat, features=gene_name, assay="SCT", ident=c("FOXM1"))+xlim(c(0,5)) 
  layout<-"
  HHH#
  AAAE
  FFF#
  BBB#
  CCCI
  GGG#
  DDD#
  "

  plt<-wrap_plots(A=esr1_cov,B=esr1_links,C=foxm1_cov,D=foxm1_links,E=esr1_expr_plot,F=peak_plot,G=peak_plot,H=annot_plot,I=foxm1_expr_plot,design=layout,heights = c(1,3,1,2,3,1,2))+ggtitle(outname)
  return(plt)
}

#PGR
mcf7_plot<-cov_plots(dat=mcf7,gene_range="chr11-100920000-101200000",gene_name="PGR",motif.name=c("MA0112.3"),add.foxm1=TRUE,outname="mcf7",foxm1_link=mcf7_foxm1,esr1_link=mcf7_esr1,ymax=30)
t47d_plot<-cov_plots(dat=t47d,gene_range="chr11-100920000-101200000",gene_name="PGR",motif.name=c("MA0112.3"),add.foxm1=TRUE,outname="t47d",foxm1_link=t47d_foxm1,esr1_link=t47d_esr1,ymax=30)
ggsave(mcf7_plot/t47d_plot,file=paste0("yw_control.PGR.featureplots.pdf"),height=40,width=20,limitsize=F)
system(paste0("slack -F ","yw_control.PGR.featureplots.pdf"," ryan_todo"))
#CENPF
mcf7_plot<-cov_plots(dat=mcf7,gene_range="chr1-214550000-214700000",gene_name="CENPF",motif.name=c("MA0112.3"),add.foxm1=TRUE,outname="mcf7",foxm1_link=mcf7_foxm1,esr1_link=mcf7_esr1,ymax=30)
t47d_plot<-cov_plots(dat=t47d,gene_range="chr1-214550000-214700000",gene_name="CENPF",motif.name=c("MA0112.3"),add.foxm1=TRUE,outname="t47d",foxm1_link=t47d_foxm1,esr1_link=t47d_esr1,ymax=30)
ggsave(mcf7_plot/t47d_plot,file=paste0("yw_control.CENPF.featureplots.pdf"),height=40,width=20,limitsize=F)
system(paste0("slack -F ","yw_control.CENPF.featureplots.pdf"," ryan_todo"))
```

## Fragment overlap with chip-seq peaks

Get aggregate window of cistrome peak by atac data
```bash

cd /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi

ref="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor/46096_sort_peaks.narrowPeak.bed"
in_1="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_E2/outs/atac_fragments.tsv.gz"
in_2="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_Veh/outs/atac_fragments.tsv.gz"

#chr start end cellID count

#set up genome
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
for i in $in_1 $in_2; do
dcis=$key;
tf_name=${tf_dcis[$key]};
ref="/home/groups/CEDAR/mulqueen/ref/cistrome/human_factor/${dcis}_sort_peaks.narrowPeak.bed";
outname=${i:59:-27}; #cut to cellranger output directory name
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
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")


plot_cistrome<-function(obj_in=t47d_atac_subset,cellline="T47D"){
sum_factor<-as.data.frame(obj_in@meta.data %>% group_by(topic_bin) %>% summarize(sum_total=sum(nCount_peaks)))
sum_factor$sample<-c("","FOXM1","ESR1")
f_in<-list.files(path="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi",pattern="cistrome.txt$")

dat<-lapply(f_in,function(x) {
  tmp<-read.table(x)
  colnames(tmp)<-c("count","dist","cellID")
  tmp<-tmp[tmp$cellID %in% obj_in$barcode,]
  tmp$dcis<-strsplit(x,"[.]")[[1]][3]
  tmp$tf<-strsplit(x,"[.]")[[1]][2]
  return(tmp)})

dat<-do.call("rbind",dat)
dat$sample<-obj_in@meta.data[match(dat$cellID,obj_in$barcode),]$topic_bin
dat<-dat %>% group_by(dist,dcis,tf,sample) %>% summarize(sum=sum(count))
dat<-merge(dat,sum_factor,by="sample",all.x=T)
dat$norm_count<-dat$sum/dat$sum_total
plt<-ggplot(dat,aes(x=dist,y=norm_count,color=sample))+geom_line(alpha=0.3)+geom_smooth()+theme_minimal()+xlim(c(-1000,1000))+facet_grid(dcis+tf ~ . ,scales="free")
ggsave(plt,file=paste0(cellline,".cistrome.enrichment.pdf"),height=50,limitsize=F)
system(paste0("slack -F ",cellline,".cistrome.enrichment.pdf ryan_todo"))
}

plot_cistrome(obj_in=t47d,cellline="T47D")
plot_cistrome(obj_in=mcf7,cellline="MCF7")

```


<!--

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

-->

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


<!--

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
