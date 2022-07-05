---
title: 10X Multiome Cell Lines
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /10xmultiome_cellline/
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

### Multiome Analysis with scREG
Per reviewer suggestion going to try multiome analysis with alternative tools to Signac.
First going to use scREG https://github.com/Durenlab/RegNMF

Installing package...

```bash
#trying not to obliterate my R libraries
conda create -n regnmf
conda activate regnmf
conda uninstall r-base #default is R 3.6
conda install r-base #reinstall R v4
conda install r-stringi
conda install r-rcppparallel
conda install r-rgeos
conda install r-hdf5r
conda install macs2
```
Installing in R...
```R
#install.packages(c("withr","devtools"),lib=.libPaths()[1])
#install.packages("Signac",lib=.libPaths()[1])
#install.packages("Seurat",lib=.libPaths()[1])
#install.packages("hdf5r",lib=.libPaths()[1])
#install.packages("tidyverse",lib=.libPaths()[1])

library(withr)
setRepositories(ind=1:3)
devtools::install_github("Durenlab/RegNMF",ref="main") #retry this with bioconda installed rccp parallel
```

Variable descriptions from the github page:

in_foldername (Character) : Path of unziped Filtered feature barcode matrix MEX("filtered_feature_bc_matrix")
out_foldername (Character) : Path of folder contain result.
fragment (Character) : Path of unziped ATAC Per fragment information file("XXXatac_fragments.tsv")
macs2path (Character) : Path of macs2
awkpath (Character) : Path of bedtools
chr (Character) : Which chromatin you want to see in the result(ex. "chr16").
from (int), to (int) : Which region of the chromasome in the result.
core (int) : How many core you want to use. You can use detectCores() function to check how many core you can use in R.
width (int), height (int) : Figure size of result.

First set up fragment files
```bash
#screg wants uncompressed fragment files for a later step. merging the two runs
zcat /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_E2/outs/atac_fragments.tsv.gz | awk 'OFS="\t" {print $1,$2,$3,"YW_si_Control_plus_E2_"$4,$5}' > /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/atac_fragments.tsv
zcat /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_Veh/outs/atac_fragments.tsv.gz | awk 'OFS="\t" {print $1,$2,$3,"YW_si_Control_plus_Veh_"$4,$5}'>> /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/atac_fragments.tsv

#set up features
zcat /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_E2/outs/filtered_feature_bc_matrix/features.tsv.gz /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/YW_si_Control_plus_Veh/outs/filtered_feature_bc_matrix/features.tsv.gz > /home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/features.tsv
```
```R
#conda activate regnmf
library(RegNMF)
library(Seurat)
library(Signac)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

my_cols = brewer.pal(4,"Spectral")
alpha_val=0.33

#read in Seurat Objects
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.rds")


#read in cellranger RNA reference file to matching genes to chr
#modified from https://www.biostars.org/p/272889/ user acvill
read_gtf <- function(file) {
  cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
  # read in raw gtf as tsv and remove comment rows
  messy <- read.table(gzfile(file), header=F, sep="\t") %>% `colnames<-`(cnames)
  
  # get the unique attribute types
  # this assumes there are no spaces in the attribute names
  att_names <- messy %>%
    select(attribute) %>%
    apply(., MARGIN = 1, FUN = str_split, pattern = '; ') %>%
    unlist() %>% trimws() %>% trimws(whitespace = ";") %>%
    sub(" .*$", "", .) %>% unique()
  att_names <- att_names[att_names != ""]
  # for each attribute type, create column
  # apply over gtf to fill in rows where attribute type is found
  for (att in c("gene_name")) { #use att_names for all attributes
    colatt <- apply(messy, MARGIN = 1, function(x) {
      var_loc<-unlist(lapply(strsplit(x[9],";"),function(x) startsWith(trimws(x),att)))
      var<-unlist(strsplit(x[9],";"))[which(var_loc)] %>%
      gsub(pattern=att,replacement="") %>%
      trimws(whitespace = '[; ]')
    })
    messy <- messy %>% add_column(colatt)
    colnames(messy)[ncol(messy)]<-att 
  }
  # remove original attribute column
  messy<- messy %>% select(-c(attribute))  
}

#note this function can take ~10-20 min if it is a large file, so I'm just saving the output to skip if I rerun this script
#ref<-read_gtf(file="/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
#ref<-ref[!duplicated(ref$gene_name),] #remove gene duplicates (this is taking the first instance)
#saveRDS(ref,file="/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.filtered.gtf.rds")#save for posterity
ref<-readRDS(file="/home/groups/CEDAR/mulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.filtered.gtf.rds")

# set up sample loop to load the RNA and ATAC data, pulling from unnormalized Seurat Object so it is consistent with other dimensionality approaches

run_regNMF_onSeurat<-function(object){
  obj <- object #count data

  #split the data
  rna_counts <- obj$RNA@counts
  atac_counts <- obj$ATAC@counts
  #filter by the barcode (commented out, using all barcodes in final seurat objects)
  #barcode_use<-read.table("~/pbmc_10k/barcode_use.txt")
  E=rna_counts #dgcmatrix
  O=atac_counts #dgcmatrix

  #log normalization
  #rna
  E=log2(1+rna_counts)
  #atac
  O=log10(1+atac_counts)

  Symbol<-row.names(E) #list
  PeakName=row.names(O) #list
  barcode=data.frame(V1=colnames(E)) #1 column data frame
  chr=unique(unlist(lapply(strsplit(PeakName,"-"),"[",1)))
  Peak_location=cbind(as.numeric(match(unlist(lapply(strsplit(PeakName,"-"),"[",1)),chr)),
  as.numeric(unlist(lapply(strsplit(PeakName,"-"),"[",2)))) 
  #matrix with index of chromosome in col 1, and start location in col 2
  ref_tmp<-ref[match(Symbol,ref$gene_name),]
  Symbol_location=cbind(as.numeric(match(ref_tmp$seqname,chr)),as.numeric(ref_tmp$start))
  #matrix with index of chromosome in col 1, and start location in col 2

  out=RegNMF(E=E, 
               O=O, 
               Symbol=Symbol, 
               PeakName=PeakName, 
               Symbol_location=Symbol_location, 
               Peak_location=Peak_location,
               core=core)
  return(out)
}

t47d_out<-run_regNMF_onSeurat(t47d)
mcf7_out<-run_regNMF_onSeurat(mcf7)
saveRDS(mcf7_out,file="yw_mcf7.control.scREGout.rds")
saveRDS(t47d_out,file="yw_t47d.control.scREGout.rds")

#Add dimensionality reduction and plot cells

integrate_regnmf_with_seurat<-function(object,regnmf_in){
cell.embeddings<-t(regnmf_in$H)
row.names(cell.embeddings)<-colnames(object$RNA)
colnames(cell.embeddings)<-paste0("regnmf_",1:ncol(cell.embeddings))
dim_reduc<-CreateDimReducObject(embeddings = cell.embeddings, key = "Dim",assay = "RNA")
object[["scREG"]]<-dim_reduc
object <- RunUMAP(object, reduction = "scREG", dims = 1:100, reduction.name = "umap.RegNMF")
object <- FindNeighbors(object, reduction = "scREG", dims = 1:100)###
object <- FindClusters(object,graph.name = "RNA_snn", resolution = 0.5)
return(object)
}


t47d<-integrate_regnmf_with_seurat(object=t47d,regnmf_in=t47d_out)
p1<-DimPlot(t47d,reduction="umap.RegNMF",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("scREG")
ggsave(p1,file=paste0("t47d.screg.umap.pdf"),width=10,height=10)
system(paste0("slack -F t47d.screg.umap.pdf ryan_todo"))


mcf7<-integrate_regnmf_with_seurat(object=mcf7,regnmf_in=mcf7_out)
p1<-DimPlot(mcf7,reduction="umap.RegNMF",group.by="sample",cols=alpha(my_cols,alpha_val))+ggtitle("scREG")
ggsave(p1,file=paste0("mcf7.screg.umap.pdf"),width=10,height=10)
system(paste0("slack -F mcf7.screg.umap.pdf ryan_todo"))


#save Seurat Objects 
saveRDS(mcf7,file="yw_mcf7.control.SeuratObject.screg.rds")
saveRDS(t47d,file="yw_t47d.control.SeuratObject.screg.rds")
```

Now going to run the cluster specific analysis for both sample (control v E2+) and binned cells (ESR1 v FOXM1 cells)

```R
#conda activate regnmf
library(RegNMF)
library(Seurat)
library(Signac)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(rtracklayer)
library(Gviz)
#modified splitgroup function
SplitGroup.default<-function(foldername,barcord,W3,H,Reg_symbol_name,Reg_peak_name,cluster){
  if(!dir.exists(foldername)){
    dir.create(foldername,recursive = TRUE )
  }
  clustern=length(unique(cluster))
  cluster_names=unique(cluster)
  barcord_cluster=data.frame(barcord=barcord,cluster=cluster)
  bfilename=paste0(foldername,"barcord_cluster.bed") #should specify that folder name should have a hanging "/"
  write.table(barcord_cluster,bfilename,sep="\t",col.names = F,row.names = F,quote = FALSE)
  chr=c()
  peaks=c()
  peake=c()
  for (i in 1:length(Reg_peak_name)) { #my peak names are split by "-" only per seurat style
    a=strsplit(Reg_peak_name[i],'-')
    chr[i]=a[[1]][1]
    peaks[i]=a[[1]][2]
    peake[i]=a[[1]][3]
  }
  df=data.frame(chr=chr,peaks=peaks,peake=peake,symbolName=Reg_symbol_name)

  H_norm=H/sqrt(rowSums(H*H))
  W3_norm=t(t(W3)*sqrt(rowSums(H*H)))
  H_w=matrix(nrow = nrow(H_norm),ncol = clustern)

  for (i in 1:clustern) {
    H_w[,i]=rowMeans(H_norm[,cluster==cluster_names[i]]) #fixed here to call cluster names rather than number
  }

  W3_cluster=W3_norm%*%H_w
  RegFolderName=paste0(foldername,"old_Reg_cluster/")
  if(!dir.exists(RegFolderName)){
    dir.create(RegFolderName,recursive = TRUE )
  } else if (length(dir(RegFolderName))) {
    file.remove(paste0(RegFolderName,dir(RegFolderName)))
  }

  for (i in 1:clustern) {
    cluster_name=cluster_names[i] #fixed here to call cluster names rather than number
    topk=order(W3_cluster[,i])[1:10000]
    outdf=df[topk,]
    outdf$Reg=W3_cluster[topk,i]
    filename=paste0(RegFolderName,"Reg_cluster",cluster_name,".bed") #fixed this to output cluster name
    write.table(outdf,filename,sep="\t",col.names = F,row.names = F,quote = FALSE)
  }


  return(a=c(barcordFileName=bfilename,RegFolderName=RegFolderName))
}


my_cols = brewer.pal(4,"Spectral")
alpha_val=0.33

#use regNMF to cluster and find differences amongst cells
mcf7_regnmf<-readRDS(file="yw_mcf7.control.scREGout.rds")
t47d_regnmf<-readRDS(file="yw_t47d.control.scREGout.rds")

#read in Seurat Objects
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.screg.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.screg.rds")

#set up topic bin annotation (replace NA values)
mcf7@meta.data$topic_bin[mcf7@meta.data$topic_bin==""]<-"rest"
t47d@meta.data$topic_bin[t47d@meta.data$topic_bin==""]<-"rest"

screg_furtherprocessing<-function(regnmf_output=mcf7_regnmf,seurat_object=mcf7,outname="mcf7",group.by="sample"){
  #regnmf_output=mcf7_regnmf;seurat_object=mcf7;outname="mcf7"
  ans<-seurat_object@meta.data[,group.by]
  out_foldername=paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/",
    outname,"_screg/")
  fragments<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/",
    outname,"_atac_fragments.tsv")

  groupName=SplitGroup.default(foldername=out_foldername,
                     barcord=colnames(seurat_object$RNA),
                     W3=regnmf_output$W3,
                     H=regnmf_output$H,
                     Reg_symbol_name=regnmf_output$Reg_gene_name,
                     Reg_peak_name=regnmf_output$Reg_peak_name,
                     cluster=ans)

  visual_need=callpeak(outfolder=out_foldername,
                     fragment=fragments,
                     barcord_cluster_whole=groupName["barcordFileName"],
                     oldRegFolder=groupName["RegFolderName"],
                     macs2path="macs2",
                     awkpath="awk",
                     cluster=unique(ans),
                     clusterL=length(ans))

plt<-Visualization(
              wholef="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/screg_fragments/",#generated earlier
              peakf=visual_need["peak_clusterF"],
              regf=visual_need["RE_clusterF"],
              chr="chr11",
              from=101010000,
              to=101150000,
              clusterlist=unique(ans),
              width=20,
              height=10)
system("slack -F result.pdf ryan_todo")
}


screg_furtherprocessing(regnmf_output=mcf7_regnmf,seurat_object=mcf7,outname="mcf7")
screg_furtherprocessing(regnmf_output=t47d_regnmf,seurat_object=t47d,outname="t47d")
screg_furtherprocessing(regnmf_output=mcf7_regnmf,seurat_object=mcf7,outname="mcf7",group.by="topic_bin")
screg_furtherprocessing(regnmf_output=t47d_regnmf,seurat_object=t47d,outname="t47d",group.by="topic_bin")
```

Alternative multiome joint embedding strategies
https://github.com/openproblems-bio/neurips2021_multimodal_topmethods/tree/main/src/joint_embedding/methods/jae


### JAE Model (Amateur)

```bash
pip install numpy scipy anndata sklearn pickle-mixin umap-learn tensorflow scanpy
cd /home/groups/CEDAR/mulqueen/ref
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz
```
Prepare seurat objects to be read in
```R
#conda activate regnmf
#remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(Signac)
library(SeuratDisk)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

#read in Seurat Objects
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.screg.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.screg.rds")
DefaultAssay(mcf7)<-"SCT"
DefaultAssay(t47d)<-"SCT"

#add cell cycle scoring to metadata
s.genes <- cc.genes$s.genes #loaded with seurat
g2m.genes <- cc.genes$g2m.genes
mcf7 <- CellCycleScoring(mcf7, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
t47d <- CellCycleScoring(t47d, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

#Generate simple seurat object for amateur
#Create a Seurat object containing the RNA adata

foldername="/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/amateur"
if(!dir.exists(foldername)){
  dir.create(foldername,recursive = TRUE )
}
setwd(foldername)

writeout_data_for_amateur<-function(x,prefix){
  dat_rna <- CreateSeuratObject(
    counts = x[["RNA"]]@counts,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  dat_atac<- CreateSeuratObject(
    counts = x[["peaks"]]@counts,
    assay = "ATAC"
  )

  dat_rna<-AddMetaData(dat_rna,metadata=x@meta.data)
  dat_atac<-AddMetaData(dat_atac,metadata=x@meta.data)
  #Save as h5ad to be read in through scanpy
  SaveH5Seurat(dat_rna,filename=paste0("yw_",prefix,".control.SeuratObject.screg.rna.h5Seurat"),overwrite=T)
  Convert(paste0("yw_",prefix,".control.SeuratObject.screg.rna.h5Seurat"), dest = "h5ad",overwrite=T)
  SaveH5Seurat(dat_atac,filename=paste0("yw_",prefix,".control.SeuratObject.screg.atac.h5Seurat"),overwrite=T)
  Convert(paste0("yw_",prefix,".control.SeuratObject.screg.atac.h5Seurat"), dest = "h5ad",overwrite=T)
}

setwd(foldername)
writeout_data_for_amateur(x=mcf7,prefix="mcf7")
writeout_data_for_amateur(x=t47d,prefix="t47d")

```

SVD preprocessing of data
Script taken from https://github.com/kimmo1019/JAE/blob/main/1_pca_pretrain.py

```python
import numpy as np
import scipy.sparse as sp
import anndata as ad
from sklearn.decomposition import TruncatedSVD
import pickle as pk
from sklearn.preprocessing import normalize
import pandas as pd
'''SVD applied to all available data (aggregate all) for multiome and cite
NOTE:
    For cite-seq, phase1 is a strict subset of phase2.
    For multiome, phase1 contains 1511 cells that is not in phase2. We intent to merge them.
    The raw phase1 and phase2 data is under the DATA_DIR folder, respectively.
    Please change DATA_DIR to the right folder.
Output:
    Two SVDs for multiome mod1 and mod2, and SVD for cite mod1 will be saved to disk under the current path (./).
'''

# TODO: change this to the directory of your own
DATA_DIR = "/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/amateur"

####Doing T47D####
phase1_multiome_mod1 = ad.read_h5ad("%s/yw_t47d.control.SeuratObject.screg.rna.h5ad"%DATA_DIR)
phase1_multiome_mod2 = ad.read_h5ad("%s/yw_t47d.control.SeuratObject.screg.atac.h5ad"%DATA_DIR)
print(phase1_multiome_mod1.shape, phase1_multiome_mod2.shape)
phase1_multiome_cells = list(phase1_multiome_mod1.obs.index)

#no need to aggregrate but copying variable
agg_multiome_mod1_count=phase1_multiome_mod1.X
agg_multiome_mod2_count=phase1_multiome_mod2.X

#scale and log transform
random_seed = 123
scale = 1e4
mod1_data = scale * normalize(agg_multiome_mod1_count,norm='l1', axis=1)
mod1_data = sp.csr_matrix.log1p(mod1_data) / np.log(10)

mod2_data = scale * normalize(agg_multiome_mod2_count,norm='l1', axis=1)
mod2_data = sp.csr_matrix.log1p(mod2_data) / np.log(10)

#apply SVD to both modalities of multiome
for n_components_mod1, n_components_mod2 in [(100,100)]:
  mod1_reducer = TruncatedSVD(n_components=n_components_mod1, random_state=random_seed)
  mod1_reducer.fit(mod1_data)
  pca_data_mod1 = mod1_reducer.transform(mod1_data)
  print('multiome 1 done',pca_data_mod1.shape)
  pk.dump(mod1_reducer, open("t47d.multiome_svd1.pkl","wb"))
  mod2_reducer = TruncatedSVD(n_components=n_components_mod2, random_state=random_seed)
  mod2_reducer.fit(mod2_data)
  pca_data_mod2 = mod2_reducer.transform(mod2_data)
  print('multiome 2 done',pca_data_mod2.shape)
  pk.dump(mod2_reducer, open("t47d.multiome_svd2.pkl","wb"))



####Doing MCF7####
phase1_multiome_mod1 = ad.read_h5ad("%s/yw_mcf7.control.SeuratObject.screg.rna.h5ad"%DATA_DIR)
phase1_multiome_mod2 = ad.read_h5ad("%s/yw_mcf7.control.SeuratObject.screg.atac.h5ad"%DATA_DIR)
print(phase1_multiome_mod1.shape, phase1_multiome_mod2.shape)
phase1_multiome_cells = list(phase1_multiome_mod1.obs.index)

#no need to aggregrate but copying variable
agg_multiome_mod1_count=phase1_multiome_mod1.X
agg_multiome_mod2_count=phase1_multiome_mod2.X

#scale and log transform
random_seed = 123
scale = 1e4
mod1_data = scale * normalize(agg_multiome_mod1_count,norm='l1', axis=1)
mod1_data = sp.csr_matrix.log1p(mod1_data) / np.log(10)

mod2_data = scale * normalize(agg_multiome_mod2_count,norm='l1', axis=1)
mod2_data = sp.csr_matrix.log1p(mod2_data) / np.log(10)

#apply SVD to both modalities of multiome
for n_components_mod1, n_components_mod2 in [(100,100)]:
  mod1_reducer = TruncatedSVD(n_components=n_components_mod1, random_state=random_seed)
  mod1_reducer.fit(mod1_data)
  pca_data_mod1 = mod1_reducer.transform(mod1_data)
  print('multiome 1 done',pca_data_mod1.shape)
  pk.dump(mod1_reducer, open("mcf7.multiome_svd1.pkl","wb"))
  mod2_reducer = TruncatedSVD(n_components=n_components_mod2, random_state=random_seed)
  mod2_reducer.fit(mod2_data)
  pca_data_mod2 = mod2_reducer.transform(mod2_data)
  print('multiome 2 done',pca_data_mod2.shape)
  pk.dump(mod2_reducer, open("mcf7.multiome_svd2.pkl","wb"))

```

Running script https://github.com/kimmo1019/JAE/blob/main/2_ae_pretrain_multiome.py
Adding some modifications and removing if then flow control where it isn't relevant to our data

```python

import argparse
parser = argparse.ArgumentParser(description='multi-modal')
parser.add_argument('-tf_seed', dest='tf_seed', type=int, default=46, help='tf random seed')
parser.add_argument('-np_seed', dest='np_seed', type=int, default=56, help='np random seed')
args = parser.parse_args()
tf_seed = args.tf_seed
np_seed = args.np_seed
suffix = str(tf_seed)+'_'+str(np_seed)
import logging
import numpy as np
import anndata as ad
from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split
import pickle as pk
import scipy
import scanpy as sc
import sys
np.random.seed(np_seed)
import tensorflow as tf
tf.random.set_seed(tf_seed)
use_label = True
logging.basicConfig(level=logging.INFO)

'''Pretrain with only exploration data (with cell type label, cell cycle scores, etc)
NOTE:
    The loss function for each epoch will be printed.
    Please change the par path to the right location of explration data
Output:
    The best pretrained model (multiome.h5 or cite.h5) will be saved to disk under the current path (./).
    The par['output'] recorded the joint embedding of exploration data.
'''

def preprocess(mod1_data, mod2_data, scale=1e4):
  nb_feats_multiome_mod1, nb_feats_multiome_mod2 = 36601, 219129 #features per matrix
  n_components_mod1, n_components_mod2 = 100, 100
  mod1_reducer = pk.load(open(meta['resources_dir'] + '/t47d.multiome_svd1.pkl','rb'))
  mod2_reducer = pk.load(open(meta['resources_dir'] + '/t47d.multiome_svd2.pkl','rb'))
  if mod1_data.shape[1] != nb_feats_multiome_mod1 and mod1_data.shape[1] != nb_feats_cite_mod1:
      print('Fake data to pass sample data test')
      mod1_data = np.zeros((mod1_data.shape[0], nb_feats_multiome_mod1))
      mod2_data = np.zeros((mod2_data.shape[0], nb_feats_multiome_mod2))
      mod1_data = scipy.sparse.csc_matrix(mod1_data)
      mod2_data = scipy.sparse.csc_matrix(mod2_data)
  mod1_data = scale * normalize(mod1_data,norm='l1', axis=1)
  mod2_data = scale * normalize(mod2_data,norm='l1', axis=1)
  mod1_data = scipy.sparse.csr_matrix.log1p(mod1_data) / np.log(10)
  mod2_data = scipy.sparse.csr_matrix.log1p(mod2_data) / np.log(10)
  pca_data_mod1 = mod1_reducer.transform(mod1_data)
  pca_data_mod2 = mod2_reducer.transform(mod2_data)
  return pca_data_mod1, pca_data_mod2


class EarlyStoppingAtMinLoss(tf.keras.callbacks.Callback):
    def __init__(self, patience=0):
        super(EarlyStoppingAtMinLoss, self).__init__()
        self.patience = patience
        self.best_weights = None
    def on_train_begin(self, logs=None):
        self.wait = 0
        self.stopped_epoch = 0
        self.best = np.Inf
    def on_epoch_end(self, epoch, logs=None):
        current = logs.get("val_loss")
        if np.less(current, self.best):
            self.best = current
            self.wait = 0
            self.best_weights = self.model.get_weights()
        else:
            self.wait += 1
            if self.wait >= self.patience:
                self.stopped_epoch = epoch
                self.model.stop_training = True
                print("Restoring model weights from the end of the best epoch.")
                self.model.set_weights(self.best_weights)
    def on_train_end(self, logs=None):
        if self.stopped_epoch > 0:
            print("Epoch %05d: early stopping" % (self.stopped_epoch + 1))


class Autoencoder(tf.keras.Model):
    def __init__(self, params, name=None):
        super(Autoencoder, self).__init__(name=name)
        self.params = params
        self.encoder = self.create_encoder()
        self.decoder = self.create_decoder()
        self.classifier = self.create_classifier()
    def get_config(self):
        return {
                "params": self.params,
        }
    def call(self, inputs, training):
        encoded = self.encoder(inputs)
        decoded = self.decoder(encoded)
        digits_cell_type, digits_batch, digits_phase = self.classifier(encoded)
        if self.params['use_batch']:
            return decoded, digits_cell_type, digits_batch, digits_phase
        else:
            return decoded, digits_cell_type
    def create_encoder(self, use_resnet=True):
        if use_resnet:
            inputs = tf.keras.layers.Input(shape=(self.params['dim'],))
            for i, n_unit in enumerate(self.params['hidden_units'][:-1]):
                if i==0:
                    x_init = tf.keras.layers.Dense(n_unit, activation='relu')(inputs)
                else:
                    x_init = tf.keras.layers.Dense(n_unit, activation='relu')(x)
                x = tf.keras.layers.Dropout(0.1)(x_init)
                x = tf.keras.layers.BatchNormalization()(x)
                x = tf.keras.layers.Dense(n_unit)(x)
                x = tf.keras.layers.Add()([x,x_init])
                x = tf.keras.layers.Activation(activation='relu')(x)
            encoded = tf.keras.layers.Dense(self.params['hidden_units'][-1], activation='relu')(x)
        else:
            inputs = tf.keras.layers.Input(shape=(self.params['dim'],))
            for i, n_unit in enumerate(self.params['hidden_units'][:-1]):
                if i==0:
                    x = tf.keras.layers.Dense(n_unit, activation='relu')(inputs)
                else:
                    x = tf.keras.layers.Dense(n_unit, activation='relu')(x)
                x = tf.keras.layers.Dropout(0.1)(x)
                x = tf.keras.layers.BatchNormalization()(x)
            encoded = tf.keras.layers.Dense(self.params['hidden_units'][-1], activation='relu')(x)
        return tf.keras.Model(inputs=inputs, outputs=encoded, name='encoder')
    def create_decoder(self):
        inputs = tf.keras.layers.Input(shape=(self.params['hidden_units'][-1],))
        for i, n_unit in enumerate(self.params['hidden_units'][:-1][::-1]):
            if i==0:
                x = tf.keras.layers.Dense(n_unit, activation='relu')(inputs)
            else:
                x = tf.keras.layers.Dense(n_unit, activation='relu')(x)
        decoded = tf.keras.layers.Dense(self.params['dim'], activation='relu')(x)
        return tf.keras.Model(inputs=inputs, outputs=decoded, name='decoder')
    def create_classifier(self):
        inputs = tf.keras.layers.Input(shape=(self.params['hidden_units'][-1],))
        digits_cell_type = inputs[:,:self.params['nb_cell_types']]
        digits_batch = inputs[:,self.params['nb_cell_types']:(self.params['nb_cell_types']+self.params['nb_batches'])]
        digits_phase = inputs[:,(self.params['nb_cell_types']+self.params['nb_batches']):(self.params['nb_cell_types']+self.params['nb_batches']+self.params['nb_phases'])]
        return tf.keras.Model(inputs=inputs, outputs=[digits_cell_type, digits_batch, digits_phase], name='classifier')


def random_classification_loss(y_true, y_pred):
    return tf.keras.metrics.categorical_crossentropy(tf.ones_like(y_pred)/nb_batches, y_pred, from_logits=True)


meta = { 'resources_dir': '.' }

DATA_DIR = "/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/amateur"
n_dim=200
cell_cycle_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2',\
                    'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', \
                    'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',\
                    'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', \
                    'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', \
                    'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', \
                    'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', \
                    'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8', \
                    'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', \
                    'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', \
                    'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', \
                    'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', \
                    'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', \
                    'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', \
                    'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', \
                    'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', \
                    'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']


def run_jae(x):
  output="logs/"+x+"output_multiome.h5ad"
  method_id = "python_starter_kit"
  logging.info('Reading `h5ad` files...')
  ad_mod1 = ad.read_h5ad("%s/yw_"%DATA_DIR+x+".control.SeuratObject.screg.rna.h5ad")
  ad_mod2 = ad.read_h5ad("%s/yw_"%DATA_DIR+x+".control.SeuratObject.screg.atac.h5ad")
  mod1_obs = ad_mod1.obs
  mod1_uns = ad_mod1.uns
  ad_mod2_var = ad_mod2.var
  mod1_mat = ad_mod1.X
  mod2_mat = ad_mod2.X
  mod1_mat = scipy.sparse.csr_matrix.log1p(mod1_mat)
  print(np.max(mod1_mat))
  print(mod1_mat.shape, mod2_mat.shape)
  meta['resources_dir']=DATA_DIR
  mod1_pca, mod2_pca = preprocess(mod1_mat, mod2_mat)
  del mod1_mat, mod2_mat
  print('load data and pca done', mod1_pca.shape, mod2_pca.shape)
  pca_combined = np.concatenate([mod1_pca, mod2_pca],axis=1)
  print(pca_combined.shape)
  del mod1_pca, mod2_pca
  cell_type_labels = mod1_obs['topic_bin']
  batch_ids = mod1_obs['sample']
  phase_labels = mod1_obs['Phase'] #from seurat cell cycle scoring function
  nb_cell_types = len(np.unique(cell_type_labels))
  nb_batches = len(np.unique(batch_ids))
  nb_phases = len(np.unique(phase_labels))-1 # 2
  c_labels = np.array([list(np.unique(cell_type_labels)).index(item) for item in cell_type_labels])
  b_labels = np.array([list(np.unique(batch_ids)).index(item) for item in batch_ids])
  p_labels = np.array([list(np.unique(phase_labels)).index(item) for item in phase_labels]) #0:G1, 1:G2M, 2: S, only consider the last two
  s_genes = cell_cycle_genes[:43]
  g2m_genes = cell_cycle_genes[43:]
  sc.pp.log1p(ad_mod1)
  sc.pp.scale(ad_mod1)   #sc.tl.score_genes_cell_cycle(ad_mod1, s_genes=s_genes, g2m_genes=g2m_genes)
  S_scores = ad_mod1.obs['S.Score'].values #use seurat values
  G2M_scores = ad_mod1.obs['G2M.Score'].values
  phase_scores = np.stack([S_scores,G2M_scores]).T #(nb_cells, 2)
  pca_train, pca_test, c_train_labels, c_test_labels, b_train_labels, b_test_labels, p_train_labels, p_test_labels, phase_train_scores, phase_test_scores = train_test_split(pca_combined, c_labels, b_labels, p_labels,phase_scores, test_size=0.1, random_state=42)
  print(pca_train.shape, c_train_labels.shape, b_train_labels.shape, p_train_labels.shape, phase_train_scores.shape)
  print(pca_test.shape, c_test_labels.shape, b_test_labels.shape, p_test_labels.shape, phase_test_scores.shape)
  X_train = pca_train
  Y_train = [pca_train, c_train_labels, b_train_labels, phase_train_scores]
  X_test = pca_test
  Y_test = [pca_test, c_test_labels, b_test_labels, phase_test_scores]
  print(nb_cell_types, nb_batches, nb_phases)
  hidden_units = [150, 120, 100, nb_cell_types+nb_batches+nb_phases+5]
  params = {
    'dim' : pca_combined.shape[1],
    'lr': 1e-4,
    'hidden_units' : hidden_units,
    'nb_layers': len(hidden_units),
    'nb_cell_types': nb_cell_types,
    'nb_batches': nb_batches,
    'nb_phases': nb_phases,
    'use_batch': True
  }
  print('Model hyper parameters:', params)
  model = Autoencoder(params)
  model.compile(tf.keras.optimizers.Adam(learning_rate = params["lr"]), 
            loss = [tf.keras.losses.MeanSquaredError(), 
                    tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                    random_classification_loss,
                    tf.keras.losses.MeanSquaredError()
                    ],
            loss_weights=[0.7, 0.2, 0.05, 0.05], run_eagerly=True)
  filepath = x+'.multiome.h5'
  callbacks = [EarlyStoppingAtMinLoss(patience=5), tf.keras.callbacks.ModelCheckpoint(filepath = filepath, monitor='val_loss', save_weights_only=True)] 
  model.fit(x=X_train, y=Y_train,
                epochs = 500,
                batch_size = 32,
                shuffle=True,
                callbacks = callbacks,
                validation_data=(X_test, Y_test),
                max_queue_size = 100, workers = 28, use_multiprocessing = True)
  print('Start evaluation')
  eval_results = model.evaluate(X_test, Y_test, batch_size=128)
  print('Total loss, loss1, loss2, loss3, loss4:',eval_results)
  f_out = open(x+'.multiome.log','a+')
  f_out.write('%s\t%.4f\t%.4f\t%.4f\t%.4f\n'%(suffix, eval_results[1], eval_results[2], eval_results[3], eval_results[4]))
  f_out.close()
  joint_embeds = model.encoder.predict(pca_combined)
  dataset_id="yw."+x
  method_id="jae"
  output="yw."+x+".jae.h5ad"
  adata = ad.AnnData(X=joint_embeds,obs=mod1_obs,uns={'dataset_id': dataset_id,'method_id': method_id,},)
  adata.write_h5ad(output, compression="gzip")
  np.savetxt(fname="yw."+x+"jae.table.tsv",X=joint_embeds,delimiter="\t")


run_jae(x="t47d")
run_jae(x="mcf7")

```
Now add back into seurat object for plotting

```R
#conda activate regnmf
#remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(Signac)
library(SeuratDisk)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/")

#read in Seurat Objects
mcf7<-readRDS(file="yw_mcf7.control.SeuratObject.screg.rds")
t47d<-readRDS(file="yw_t47d.control.SeuratObject.screg.rds")

mcf7_jae<-read.table("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/amateur/yw.mcf7jae.table.tsv",header=F)
t47d_jae<-read.table("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/amateur/yw.t47djae.table.tsv",header=F)


plt<-ggplot(data=t47d_jae,aes(x=V1,y=V2))+geom_point()
ggsave(plt,file="test.pdf")
system("slack -F test.pdf ryan_todo")
```
