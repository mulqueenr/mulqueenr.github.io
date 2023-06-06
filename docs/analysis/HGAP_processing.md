---
title: HGAP Processing
layout: cedar_analysis
author: Ryan Mulqueen
category: Adey Lab
---

# Analysis for HGAP 
## Integration across public data sets.
1. Find scATAC and RNA to integrate with: small, >=5K cells (HGAP dataset: frontal cortex)
2. Integrate 1 ATAC, 1 RNA > confirm cell type assignment
3. With the ATAC or RNA (cooler) integration > define COPs vs. OPCs (module enrichment)
4. scIB - quanitfy Harmony integration

##First finding Datasets for integration
| Notes | Source | Assay | Cell Count | Data source |
|------ | ------ | ------ | -------- | ----------- | 
|MS Data set| Portmortem Brain | snRNA | 17799 | https://cells.ucsc.edu/?ds=oligodendrocyte-ms |
| GWAS ATAC comparison | Postmortem 39 healthy| scATAC | 70631 | https://cells.ucsc.edu/?ds=neuro-degen-atac+peaks |
| Allen Brain Map | Postmortem brain (MTG, ACC, V1C, M1C, S1C) |  snRNA (SMART-seq) | 49495 | https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq |
| Human Hippo Axis | Hippocampus | snRNA | 129,908 | https://cells.ucsc.edu/?ds=human-hippo-axis |
| Human Hippocampus Lifespan | Hippocampus | snRNA| 224,464 | https://cells.ucsc.edu/?bp=brain&org=Human+(H.+sapiens)&ds=hippo-lifespan |

# Plan for public data set integrations: 
* Cortex: 
** Allen Brain Map: focus on M1C (Primary Motor Cortex) (RNA) #Done https://knowledge.brain-map.org/data/HPAW0I2JNX5P35OPOPL/summary
** GWAS ATAC Comparison: Focus on Middle Frontal Gyrus (ATAC) #Done

* Hippocampus:
** Human Hippo Axis (RNA)
** Human Hippocampus Lifespan (RNA) 
** GWAS ATAC Comparison: Focus on Hippocampus (ATAC) #Done

```bash
mkdir /home/groups/CEDAR/mulqueen/human_brain_ref
```

### Set up SRA-toolkit
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

```bash
cd /home/groups/CEDAR/mulqueen/src
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.0-mac64/bin
vdb-config --prefetch-to-cwd
```

# RNA
For RNA we can match gene names for integration. This means expression matrices and metadata is a sufficient starting point.
## Cortex RNA

### Allen Brain Span M1 Cortex (RNA)
https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
<!-- Done -->
```bash
mkdir /home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex
#Human download
cd /home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip
wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds
```

### Process Data into Seurat Object
Following https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

```R
library(Seurat)
library(ggplot2)
library(data.table)


allenbrainspan_to_seurat<-function(exprMat,meta,outname){
	expr<-as.data.frame(fread(exprMat,sep=",",header=T))
	row.names(expr)<-expr[,1]
	expr<-expr[,2:ncol(expr)]
	meta<-as.data.frame(fread(meta,header=T))
	row.names(meta)<-meta$Cell
	dat <- CreateSeuratObject(counts = expr, project = outname, min.cells = 3, min.features = 500)
	dat<-AddMetaData(dat,meta)
	saveRDS(dat, file = paste0(outname,".rds"))
	dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
	dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(dat)
	dat <- ScaleData(dat, features = all.genes)
	dat <- RunPCA(dat, features = VariableFeatures(object = dat))
	plt<-ElbowPlot(dat)
	ggsave(plt,file=paste0(outname,".elbowplot.pdf"))
	system(paste0("slack -F ",outname,".elbowplot.pdf"," ryan_todo"))
	dat <- FindNeighbors(dat, dims = 1:14)
	dat <- FindClusters(dat, resolution = 0.5)
	dat <- RunUMAP(dat, dims = 1:14)
	saveRDS(dat, file = paste0(outname,".rds"))
	return(dat)
}

setwd("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex")

dat<-allenbrainspan_to_seurat(
	exprMat="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv",
	meta="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv",
	outname="allen_brainspan_humancortex"
	)
plt<-DimPlot(dat, reduction = "umap",group.by=c("class_label","subclass_label"))
ggsave(plt,file=paste0(outname,".dimplot.pdf"),width=30)
system(paste0("slack -F ",outname,".dimplot.pdf"," ryan_todo"))

```

## Hippocampus RNA Data Sets

```R
library(Seurat)
library(ggplot2)
library(data.table)

ucsc_cellbrowser_to_seurat<-function(exprMat,meta,outname){
	expr<-as.data.frame(fread(exprMat))
	row.names(expr)<-expr[,1]
	expr<-expr[,2:ncol(expr)]
	meta<-as.data.frame(fread(meta))
	row.names(meta)<-meta$Cell
	dat <- CreateSeuratObject(counts = expr, project = outname, min.cells = 3, min.features = 500)
	dat<-AddMetaData(dat,meta)
	saveRDS(dat, file = paste0(outname,".rds"))
	dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
	dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(dat)
	dat <- ScaleData(dat, features = all.genes)
	dat <- RunPCA(dat, features = VariableFeatures(object = dat))
	plt<-ElbowPlot(dat)
	ggsave(plt,file=paste0(outname,".elbowplot.pdf"))
	system(paste0("slack -F ",outname,".elbowplot.pdf"," ryan_todo"))
	dat <- FindNeighbors(dat, dims = 1:14)
	dat <- FindClusters(dat, resolution = 0.5)
	dat <- RunUMAP(dat, dims = 1:14)
	saveRDS(dat, file = paste0(outname,".rds"))
	return(dat)
}

#Human Hippo Axis
system("mkdir /home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis")
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis")
dat<-ucsc_cellbrowser_to_seurat(
	exprMat="https://cells.ucsc.edu/human-hippo-axis/exprMatrix.tsv.gz",
	meta="https://cells.ucsc.edu/human-hippo-axis/meta.tsv",
	outname="hippo_axis")
outname="hippo_axis"
plt<-DimPlot(dat, reduction = "umap",group.by=c("Cluster"))
ggsave(plt,file=paste0(outname,".dimplot.pdf"),width=10)
system(paste0("slack -F ",outname,".dimplot.pdf"," ryan_todo"))

#Human Hippo Lifespan
system("mkdir /home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan")
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan")
dat<-ucsc_cellbrowser_to_seurat(
	exprMat="https://cells.ucsc.edu/hippo-lifespan/all/exprMatrix.tsv.gz",
	meta="https://cells.ucsc.edu/hippo-lifespan/all/meta.tsv",
	outname="hippo_lifespan")
outname="hippo_lifespan"
plt<-DimPlot(dat, reduction = "umap",group.by=c("MajorCellTypes"))
ggsave(plt,file=paste0(outname,".dimplot.pdf"),width=10)
system(paste0("slack -F ",outname,".dimplot.pdf"," ryan_todo"))
```

### Integrating Public RNA data with our Gene Activity
```bash
mkdir /home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses
```

Make full data set seurat objects
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
input_dir="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/Regions/"
fragments="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/coverage.plots.BAM_DataFreeze_Cortex_and_Hippocampus/BAM_DataFreeze_Cortex_and_Hippocampus.bbrd.q10.filt.bam.fragments.tsv.gz"
#(the maga directory part was about microglia, not political affiliation)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

#function to read in sparse matrix format from atac-count
read_in_sparse<-function(x){ #x is character file prefix followed by .bbrd.q10.500.counts.sparseMatrix.values.gz
    IN<-as.matrix(read.table(paste0(x,".filt.sparseMatrix.values.gz")))
    IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
    COLS<-read.table(paste0(x,".filt.sparseMatrix.cols.gz"))
    colnames(IN)<-COLS$V1
    ROWS<-read.table(paste0(x,".filt.sparseMatrix.rows.gz"))
    row.names(IN)<-ROWS$V1
    return(IN)
}


make_seurat_object<-function(input_sparse,outname,in_meta=NULL){
	counts<-read_in_sparse(paste0(input_dir,input_sparse))
	counts<-counts[,order(colnames(counts))]

	#Generate ChromatinAssay Objects
	chromatinassay <- CreateChromatinAssay(
	  counts = counts,
	  #genome="hg38",
	  min.cells = 1,
	  annotation=annotation,
	  sep=c("_","_"),
	  fragments=fragments)

	#Create Seurat Objects
	dat <- CreateSeuratObject(
	  counts = chromatinassay,
	  assay = "peaks")

	if (in_meta!=NULL){
	dat<-AddMetaData(dat,metadata=in_meta)
	}
	#saving unprocessed SeuratObjects
	saveRDS(dat,file=paste0(outname,".SeuratObject.Rds"))
}


celltype_umap_plot<-function(in_dat,color_by="",outname){
	in_dat<-RunTFIDF(in_dat)
	in_dat<-FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)
	in_dat <- RunUMAP(in_dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

	plt1<-DimPlot(in_dat,group.by=color_by)#,reduction="umap.atac")
	ggsave(plt1,file=paste0(outname,".umap.png"),limitsize=F)
	system(paste0("slack -F ",outname,".umap.png"," ryan_todo"))
	saveRDS(in_dat,file=paste0(outname,".SeuratObject.Rds"))
	return(in_dat)
}

#cortex
meta.data<-read.table(file="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/Regions/cortex/cortex.metadata",row.names=1,header=T)
meta.data<-meta.data[order(row.names(meta.data)),]
dat<-make_seurat_object(input_sparse="cortex/cortex", outname="cortex",in_meta=meta.data)
#Forgot to add celltypes to metadata!
celltype_annot<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/Regions/cortex/cortex.putative.celltype.annot",header=F)
colnames(celltype_annot)<-c("cellID","Putative_Celltype")
row.names(celltype_annot)<-celltype_annot$cellID
celltype_annot<-celltype_annot[order(row.names(celltype_annot)),]
dat<-AddMetaData(dat,setNames(celltype_annot$Putative_Celltype,nm=row.names(celltype_annot)),col.name="Putative_Celltype")
#Finally plotting to make sure metadata is in correct order (cell types should cluster separately)
dat<-celltype_umap_plot(in_dat=dat,color_by="Putative_Celltype",outname="cortex")

#hippocampus
meta.data<-read.table(file="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/Regions/hippocampus/hippocampus.putative.celltype.annot",header=F)
colnames(meta.data)<-c("cellID","Putative_Celltype")
row.names(meta.data)<-meta.data$cellID
meta.data<-meta.data[order(row.names(meta.data)),]
dat<-make_seurat_object(input_sparse="hippocampus/hippocampus", outname="hippocampus",in_meta=meta.data)
dat<-AddMetaData(dat,setNames(meta.data$Putative_Celltype,nm=row.names(meta.data)),col.name="Putative_Celltype")
dat<-celltype_umap_plot(in_dat=dat,color_by="Putative_Celltype",outname="hippocampus")

```

<!--
Subset the data to 5% of cells per cell type. This is going to greatly increase the speed for gene activity calculation.

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#Ours
hgap_cortex<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/cortex.SeuratObject.Rds")
hgap_hippo<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/hippocampus.SeuratObject.Rds")

subset_5perc<-function(in_dat,prefix){
	Idents(in_dat)<-in_dat$Putative_Celltype
	cells_to_keep=unlist(lapply(unique(Idents(in_dat)),function(x) {
			ident_count=sum(Idents(in_dat)==x)
			cells_to_return=row.names(in_dat@meta.data[Idents(in_dat)==x,])
			cells_to_return=sample(cells_to_return,size=as.integer((ident_count/100)*5))
			return(cells_to_return)#sample to 5% per cell type
			}))
	print("Writing out cells.")
	write.table(cells_to_keep,file=paste0(prefix,".cellIDs.tsv"),sep="\t",col.names=F,row.names=F,quote=F)
}

subset_5perc(in_dat=hgap_cortex,prefix="cortex.5perc")
subset_5perc(in_dat=hgap_hippo,prefix="hippocampus.5perc")

```
Use python to subset the fragments file to just those in the 5 % per cell subtype subset.

```python
import gzip
import pandas as pd

chunk_size=500000000
read_lines=range(0,9000000000,chunk_size)

#run for cortex
cellid_in=pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/cortex.5perc.cellIDs.tsv",header=None,sep=" ") #read in cellIDs
for i in range(0,len(read_lines)-1):
	df = pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/coverage.plots.BAM_DataFreeze_Cortex_and_Hippocampus/BAM_DataFreeze_Cortex_and_Hippocampus.bbrd.q10.filt.bam.fragments.tsv.gz",sep="\t",header=None,nrows=chunk_size,skiprows=i) #read in fragments file
	df_filt=df[df[3].isin(cellid_in[0])] #filter to cell IDs in input cell id in
	print(df_filt)
	df_filt.to_csv('/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/cortex.5perc.q10.filt.bam.fragments.'+str(read_lines[i])+'tsv.gz',compression="gzip",sep="\t",index=False) #write out compressed opc fragment file

#run again for hippo, trying just without chunking
cellid_in=pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/hippocampus.5perc.cellIDs.tsv",header=None,sep=" ") #read in cellIDs
for i in range(0,len(read_lines)-1):
	df = pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/coverage.plots.BAM_DataFreeze_Cortex_and_Hippocampus/BAM_DataFreeze_Cortex_and_Hippocampus.bbrd.q10.filt.bam.fragments.tsv.gz",sep="\t",header=None,nrows=chunk_size,skiprows=i) #read in fragments file
	df_filt=df[df[3].isin(cellid_in[0])] #filter to cell IDs in input cell id in
	print(df_filt)
	df_filt.to_csv('/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/hippocampus.5perc.q10.filt.bam.fragments.'+str(read_lines[i])+'tsv.gz',compression="gzip",sep="\t",index=False) #write out compressed opc fragment file

```

Merge all the output files and index

```bash
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
cat cortex.5perc.q10.*tsv.gz > cortex.5perc.q10.filt.bam.fragments.tsv.gz &
zcat cortex.5perc.q10.filt.bam.fragments.tsv.gz | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > cortex.5perc.q10.filt.bam.fragments.sorted.tsv.gz; wait ;
zcat cortex.5perc.q10.filt.bam.fragments.sorted.tsv.gz | tail -n +18 | $bgzip > cortex.5perc.q10.filt.bam.fragments.sorted2.tsv.gz; wait ; #
$tabix -p bed cortex.5perc.q10.filt.bam.fragments.sorted2.tsv.gz &
#this is to remove header from concatenation

cat hippocampus.5perc.q10.*tsv.gz > hippocampus.5perc.q10.filt.bam.fragments.tsv.gz ; wait;
zcat hippocampus.5perc.q10.filt.bam.fragments.tsv.gz | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > hippocampus.5perc.q10.filt.bam.fragments.sorted.tsv.gz; wait ;
zcat hippocampus.5perc.q10.filt.bam.fragments.sorted.tsv.gz | tail -n +18 | $bgzip > hippocampus.5perc.q10.filt.bam.fragments.sorted2.tsv.gz; wait ; 
$tabix -p bed hippocampus.5perc.q10.filt.bam.fragments.sorted2.tsv.gz &
#this is to remove header from concatenation
```

-->


### ATAC-RNA Integration of Hippocampus and Cortex
```R
library(Seurat)
library(Signac) 
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(Matrix)
library(monocle3,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/") #using old install of monocle, just need for as.cell_data_set conversion
library(cicero,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/")
library(SeuratWrappers)
library(rtracklayer)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")


split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

make_sparse_matrix <- function(data,
                               i.name = "Peak1",
                               j.name = "Peak2",
                               x.name = "value") {
  if(!i.name %in% names(data) |
     !j.name %in% names(data) |
     !x.name %in% names(data)) {
    stop('i.name, j.name, and x.name must be columns in data')
  }
  
  data$i <- as.character(data[,i.name])
  data$j <- as.character(data[,j.name])
  data$x <- data[,x.name]
  
  if(!class(data$x) %in%  c("numeric", "integer"))
    stop('x.name column must be numeric')
  
  peaks <- data.frame(Peak = unique(c(data$i, data$j)),
                      index = seq_len(length(unique(c(data$i, data$j)))))
  
  data <- data[,c("i", "j", "x")]
  
  data <- rbind(data, data.frame(i=peaks$Peak, j = peaks$Peak, x = 0))
  data <- data[!duplicated(data[,c("i", "j", "x")]),]
  data <- data.table::as.data.table(data)
  peaks <- data.table::as.data.table(peaks)
  data.table::setkey(data, "i")
  data.table::setkey(peaks, "Peak")
  data <- data[peaks]
  data.table::setkey(data, "j")
  data <- data[peaks]
  data <- as.data.frame(data)
  
  data <- data[,c("index", "i.index", "x")]
  data2 <- data
  names(data2) <- c("i.index", "index", "x")
  
  data <- rbind(data, data2)
  
  data <- data[!duplicated(data[,c("index", "i.index")]),]
  data <- data[data$index >= data$i.index,]
  
  sp_mat <- Matrix::sparseMatrix(i=as.numeric(data$index),
                                 j=as.numeric(data$i.index),
                                 x=data$x,
                                 symmetric = TRUE)
  
  colnames(sp_mat) <- peaks[order(peaks$index),]$Peak
  row.names(sp_mat) <- peaks[order(peaks$index),]$Peak
  return(sp_mat)
}

build_composite_gene_activity_matrix <- function(input_cds,
                                                 site_weights,
                                                 cicero_cons_info,
                                                 dist_thresh=250000,
                                                 coaccess_cutoff=0.25) {
    accessibility_mat <- exprs(input_cds)
    promoter_peak_table <- fData(input_cds)
    promoter_peak_table$peak <- as.character(row.names(promoter_peak_table))
    promoter_peak_table <-
        promoter_peak_table[!is.na(promoter_peak_table$gene),]
    promoter_peak_table <- promoter_peak_table[,c("peak", "gene")]
    promoter_peak_table$gene <- as.character(promoter_peak_table$gene)

    # Make site_weight matrix
    site_names <- names(site_weights)
    site_weights <- as(Matrix::Diagonal(x=as.numeric(site_weights)),
                      "sparseMatrix")
    row.names(site_weights) <- site_names
    colnames(site_weights) <- site_names

    # Find distance between cicero peaks. If distance already calculated, skip
    if ("dist" %in% colnames(cicero_cons_info) == FALSE) {
        Peak1_cols <- split_peak_names(cicero_cons_info$Peak1)
        Peak2_cols <- split_peak_names(cicero_cons_info$Peak2)
        Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                          as.integer(Peak1_cols[,2])) / 2)
        Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                          as.integer(Peak2_cols[,2])) / 2)
        cicero_cons_info$dist <- abs(Peak2_bp - Peak1_bp)
    }

    # Get connections between promoters and distal sites above coaccess
    # threshold
    nonneg_cons <-
        cicero_cons_info[(cicero_cons_info$Peak1 %in%
                          promoter_peak_table$peak |
                          cicero_cons_info$Peak2 %in%
                          promoter_peak_table$peak) &
                          cicero_cons_info$coaccess >= coaccess_cutoff &
                          cicero_cons_info$dist < dist_thresh,]
    nonneg_cons <- nonneg_cons[,c("Peak1", "Peak2", "coaccess")]
    nonneg_cons <- nonneg_cons[!duplicated(nonneg_cons),]

    nonneg_cons$Peak1 <- as.character(nonneg_cons$Peak1)
    nonneg_cons$Peak2 <- as.character(nonneg_cons$Peak2)

    nonneg_cons <- rbind(nonneg_cons,
                        data.frame(Peak1=unique(promoter_peak_table$peak),
                                   Peak2=unique(promoter_peak_table$peak),
                                   coaccess=0))

    # Make square matrix of connections from distal to proximal
    distal_connectivity_matrix <- make_sparse_matrix(nonneg_cons,
                                                    x.name="coaccess")

    # Make connectivity matrix of promoters versus all
    promoter_conn_matrix <-
        distal_connectivity_matrix[unique(promoter_peak_table$peak),]

    # Get list of promoter and distal sites in accessibility mat
    promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

    # Get accessibility info for promoters
    promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

    # Get accessibility for distal sites
    distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

    # Scale connectivity matrix by site_weights
    scaled_site_weights <- site_weights[distal_safe_sites,distal_safe_sites, drop=FALSE]
    total_linked_site_weights <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    total_linked_site_weights <- 1/Matrix::rowSums(total_linked_site_weights,
                                                na.rm=TRUE)
    total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] <- 0
    total_linked_site_weights[is.na(total_linked_site_weights)] <- 0
    total_linked_site_weights[is.nan(total_linked_site_weights)] <- 0
    site_names <- names(total_linked_site_weights)
    total_linked_site_weights <- Matrix::Diagonal(x=total_linked_site_weights)
    row.names(total_linked_site_weights) <- site_names
    colnames(total_linked_site_weights) <- site_names
    scaled_site_weights <- total_linked_site_weights %*%
        promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

    # Multiply distal accessibility by site weights
    distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

    distal_activity_scores <-
        distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

    # Sum distal and promoter scores
    promoter_activity_scores <- distal_activity_scores +
        promoter_access_mat_in_cicero_map

    # Make and populate final matrix
    promoter_gene_mat <-
        Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                             i=as.numeric(factor(promoter_peak_table$gene)),
                             x=1)
    colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
    row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
    promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
    gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores

    return(gene_activity_scores)
}

# gene annotation sample for gene activity
get_annot<-function(){
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
	return(gene_annotation)
	}


#Cicero processing function  
cicero_processing<-function(object_input=hgap_hippo,prefix="hgap_hippo_5perc",k=500){
      #Generate CDS format from Seurat object
      #atac.cds <- as.CellDataSet(object_input,assay="peaks",reduction="umap.atac")
      object_input<-subset(object_input,features=row.names(object_input@assays$peaks@counts)[which(rowSums(object_input@assays$peaks@counts)>=1000)])
      atac.cds <- as.cell_data_set(object_input,assay="peaks",reduction="umap.atac")
      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = object_input@reductions$umap.atac@cell.embeddings,k=k)

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      atac.cicero<-readRDS(paste(prefix,"atac_cicero_cds.Rds",sep="_"))

      genome<-SeqinfoForUCSCGenome("hg38")
      genome.df<-data.frame("chr"=seqnames(genome),"length"=seqlengths(genome))
      genome.df<-genome.df[genome.df$chr %in% paste0("chr",seq(1:22)),]

      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      conns<-readRDS(paste(prefix,"atac_cicero_conns.Rds",sep="_"))

      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      ccans<-readRDS(paste(prefix,"atac_cicero_ccans.Rds",sep="_"))

      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      saveRDS(object_input,file=paste0(prefix,".","gene_activity.Rds"))
      return(object_input)
}

geneactivity_processing<-function(object_input=hgap_hippo,conns_input=conns,prefix="hgap_hippo_5perc"){
			anno<-get_annot()
			atac.cds <- as.cell_data_set(object_input,assay="peaks",reduction="umap.atac")
			#atac.cds<-newCellDataSet(cellData=object_input@assays$peaks@counts,featureData=atac.cds@featureData,phenoData=atac.cds@phenoData))
			#cds <- new("CellDataSet", assayData = assayDataNew("environment", 
      #  exprs = object_input@assays$peaks@counts), phenoData = atac.cds@phenoData, featureData = atac.cds@featureData, 
      #  lowerDetectionLimit =  0.1, expressionFamily = VGAM::negbinomial.size(), 
      #  dispFitInfo = new.env(hash = TRUE))
      atac.cds<- annotate_cds_by_site(atac.cds, anno)
      fData(atac.cds)<-cbind(fData(atac.cds),
      		site_name=row.names(fData(atac.cds)),
      		chr=gsub(pattern="chr",replace="",unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",1))),
      		bp1=unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",2)),
      		bp2=unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",3)))

			input_cds=atac.cds
			cicero_cons_info=conns
			site_weights=NULL
			dist_thresh=250000
			coaccess_cutoff=0.25

			accessibility_mat <- exprs(input_cds)
      site_weights <- Matrix::rowMeans(accessibility_mat) / Matrix::rowMeans(accessibility_mat)
      site_weights[names(site_weights)] <- 1
			gene_promoter_activity <- build_composite_gene_activity_matrix(input_cds,
                                             site_weights,
                                             cicero_cons_info,
                                             dist_thresh=dist_thresh,
                                             coaccess_cutoff=coaccess_cutoff)
			unnorm_ga<-gene_promoter_activity
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
      object_input[['GeneActivity']]<- CreateAssayObject(counts = unnorm_ga) 

		  # normalize
		  object_input <- NormalizeData(
		    object = object_input,
		    assay = 'GeneActivity',
		    normalization.method = 'LogNormalize',
		    scale.factor = median(object_input$nCount_GeneActivity)
		  )
		  saveRDS(object_input,file=paste0(prefix,".","gene_activity.Rds"))
      return(object_input)
  		
}

generate_transfer_anchor<-function(in_dat,ref_dat,feat,prefix){
	#generate LSI matrix for downstream normalization
	#downsample cells to 5% per identity

	#in_dat <- RunTFIDF(in_dat)
	#in_dat <- FindTopFeatures(in_dat, min.cutoff = 'q0')
	#in_dat <- RunSVD(in_dat)

	#generate gene activity for just variable genes, add to object and normalize
	#ga_mat<-GeneActivity(in_dat, features = feat)
	#in_dat[['RNA']] <- CreateAssayObject(counts = ga_mat)
	#in_dat<- NormalizeData(
	#  object = in_dat,
	#  assay = 'RNA',
	#  normalization.method = 'LogNormalize',
	#  scale.factor = median(in_dat$nCount_RNA)
	#)

	DefaultAssay(in_dat)<-"GeneActivity"
	transfer.anchors <- FindTransferAnchors(
	  reference = ref_dat,
	  query = in_dat,
	  reduction = 'cca'
	)
	saveRDS(transfer.anchors,paste0(prefix,".transferanchors.rds"))
}

sample_label_transfer<-function(in_dat,ref_dat,transfer.anchors.,prefix="Tcell_",transfer_label="celltype"){
  predictions<- TransferData(
    anchorset = transfer.anchors.,
    refdata = ref_dat@meta.data[,transfer_label],
    weight.reduction = in_dat[["lsi"]],
    dims=2:30
  )
  colnames(predictions)<-paste0(prefix,colnames(predictions))
  in_dat<-AddMetaData(in_dat,metadata=predictions)
  return(in_dat)
}

coembed_data<-function(in_dat,ref_dat,transfer.anchors.,feat,prefix,assay_name){
	imputation<- TransferData(
	  anchorset = transfer.anchors.,
	  refdata = GetAssayData(ref_dat, assay = "RNA", slot = "data")[feat,],
	  weight.reduction = in_dat[["lsi"]],
	  dims=1:30)

	in_dat[[assay_name]]<-imputation
	coembed<-merge(x=ref_dat,y=in_dat)
	coembed <- ScaleData(coembed, features = feat, do.scale = FALSE)
	coembed <- RunPCA(coembed, features = feat, verbose = FALSE)
	coembed <- RunUMAP(coembed, dims = 1:30)
	saveRDS(in_dat,file=paste0(prefix,".SeuratObject.Rds"))
	return(coembed)
}

generate_confusion_matrix<-function(in_dat,metadat_prefix,filter_val=0){
	#Generate Confusion Matrix
	#Add row scaling
	metdat<-in_dat@meta.data

	metdat<-metdat[metdat[paste0(metadat_prefix,"prediction.score.max")]>filter_val,]#########FILTER TO HIGHER VALUES BEFORE PLOTTING
	conf_dat<-as.data.frame.matrix(table(metdat$Putative_Celltype,metdat[,paste0(metadat_prefix,"predicted.id")]))
	col_names<-colnames(conf_dat)
	row_names<-row.names(conf_dat)
	conf_dat<-do.call("cbind",lapply(1:ncol(conf_dat),function(x) conf_dat[,x]/sum(conf_dat[,x])))
	conf_dat<-do.call("rbind",lapply(1:nrow(conf_dat),function(x) conf_dat[x,]/sum(conf_dat[x,])))
	colnames(conf_dat)<-col_names
	row.names(conf_dat)<-row_names

	colfun<-colorRamp2(c(0, 0.5, 1), c("white","grey","black"))
	pdf(paste0("hgap_",metadat_prefix,".confmat.pdf"))
	plt<-Heatmap(conf_dat,col=colfun)
	print(plt)
	dev.off()
	system(paste0("slack -F ",paste0("hgap_",metadat_prefix,".confmat.pdf")," ryan_todo"))
}

#Full set processing on long request node 
#Ours
	hgap_hippo<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/hippocampus.SeuratObject.Rds")
#Gene activity processing of HGAP Hippo
	hgap_hippo<-cicero_processing(object_input=hgap_hippo,prefix="hgap_hippo_100perc",k=500)
	conns<-readRDS("hgap_hippo_100perc_atac_cicero_conns.Rds")
	hgap_hippo<-geneactivity_processing(object_input=hgap_hippo,conns_input=conns,prefix="hgap_hippo_100perc")


#Process HGAP hippocampus and Hippo Axis integration
	hgap_hippo<-readRDS("hgap_hippo_100perc.gene_activity.Rds")
	hippo_axis<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.rds")
	Idents(hippo_axis)<-hippo_axis$Cluster
	hippo_axis_subset<-subset(hippo_axis,downsample=min(table(Idents(hippo_axis))))
	#markers<-FindAllMarkers(hippo_axis,only.pos=TRUE) 
	#saveRDS(markers,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.markers.rds")
	#markers<-readRDS(file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.markers.rds")
	#features=row.names(markers)
	features<-VariableFeatures(hippo_axis)
	generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_axis,prefix="hippo.100perc.axis",feat=features)
	generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_axis_subset,prefix="hippo.100perc.axis_subset",feat=features)

	#transfer_anchors=readRDS("hippo.axis.transferanchors.rds")
	hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,ref_dat=hippo_axis,
		prefix="hippoaxis_cluster_",
		transfer.anchors.=readRDS("hippo.100perc.axis.transferanchors.rds"),
		transfer_label="Cluster")
	hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,ref_dat=hippo_axis_subset,
		prefix="hippoaxis_subset_cluster_",
		transfer.anchors.=readRDS("hippo.100perc.axis_subset.transferanchors.rds"),
		transfer_label="Cluster")

	generate_confusion_matrix(in_dat=hgap_hippo,metadat_prefix="hippoaxis_cluster_")
	generate_confusion_matrix(in_dat=hgap_hippo,metadat_prefix="hippoaxis_subset_cluster_")

	saveRDS(hgap_hippo,file="hippo.axis.100perc.SeuratObject.Rds")



#Process HGAP hippocampus and Hippo Lifespan integration
	hgap_hippo<-readRDS("hgap_hippo_100perc.gene_activity.Rds")
	hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.rds")
	Idents(hippo_lifespan)<-hippo_lifespan$MajorCellTypes
	hippo_lifespan_subset<-subset(hippo_lifespan,downsample=min(table(Idents(hippo_lifespan))))
	#markers<-FindAllMarkers(hippo_lifespan,only.pos=TRUE) 
	#saveRDS(markers,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.markers.rds")	
	#markers<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.markers.rds")	
	#features<-row.names(markers)
	features<-VariableFeatures(hippo_lifespan)

	hgap_hippo<-readRDS("hgap_hippo_100perc.gene_activity.Rds")
	generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_lifespan,prefix="hippo.100perc.lifespan",feat=features)
	#transfer_anchors=readRDS("hippo.100perc.lifespan.transferanchors.rds")
	generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_lifespan_subset,prefix="hippo.100perc.lifespan_subset",feat=features)
	#transfer_anchors=readRDS("hippo.100perc.lifespan_subset.transferanchors.rds")

	# Downsample the number of cells per identity class
	hgap_hippo<-readRDS("hippo.axis.100perc.SeuratObject.Rds")
	hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,ref_dat=hippo_lifespan,prefix="hippo.100perc.lifespan",transfer.anchors.=transfer_anchors,transfer_label="MajorCellTypes")

		hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,ref_dat=hippo_lifespan,prefix="hippo.100perc.lifespan_subset",transfer.anchors.=transfer_anchors,transfer_label="MajorCellTypes")

generate_confusion_matrix(in_dat=hgap_hippo,metadat_prefix="hippo.100perc.lifespan")
generate_confusion_matrix(in_dat=hgap_hippo,metadat_prefix="hippo.100perc.lifespan_subset")

saveRDS(hgap_hippo,file="hippo.lifespan.100perc.SeuratObject.Rds")


#Process HGAP cortex and Allen Brainspan integration
hgap_cortex<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/cortex.SeuratObject.Rds")
#Gene activity processing of HGAP Cortex
	hgap_cortex<-cicero_processing(object_input=hgap_cortex,prefix="hgap_cortex_100perc")
	conns<-readRDS("hgap_cortex_100perc_atac_cicero_conns.Rds")
	hgap_cortex<-geneactivity_processing(object_input=hgap_cortex,conns_input=conns,prefix="hgap_cortex_100perc")
	hgap_cortex<-readRDS("hgap_cortex_100perc.gene_activity.Rds")


#Process HGAP Cortex and Allen Brainmap cortex data
	hgap_cortex<-readRDS("hgap_cortex_100perc.gene_activity.Rds")
	cortex_brainspan<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
	Idents(cortex_brainspan)<-cortex_brainspan$subclass_label
	features<-VariableFeatures(hgap_cortex)

	generate_transfer_anchor(in_dat=hgap_cortex,ref_dat=cortex_brainspan,prefix="cortex.100perc.brainspan",feat=features)
	#transfer_anchors=readRDS("cortex.100perc.brainspan.transferanchors.rds")

	hgap_cortex<-sample_label_transfer(in_dat=hgap_cortex,ref_dat=cortex_brainspan,
		prefix="cortex_brainspan_cluster_",
		transfer.anchors.=readRDS("cortex.100perc.brainspan.transferanchors.rds"),
		transfer_label="subclass_label")

generate_confusion_matrix(in_dat=hgap_cortex,metadat_prefix="cortex_brainspan_cluster_")

saveRDS(hgap_cortex,file="cortex.100perc.brainspan.SeuratObject.Rds")




```

### Identify OPC cell states in RNA data

```R
library(Seurat)
library(Signac) 
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(Matrix)
library(monocle3,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/") #using old install of monocle, just need for as.cell_data_set conversion
library(cicero,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/")
library(SeuratWrappers)
library(rtracklayer)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")


generate_transfer_anchor_atac_to_rna<-function(in_dat,ref_dat,feat,prefix){
	#generate LSI matrix for downstream normalization

	#in_dat <- RunTFIDF(in_dat)
	#in_dat <- FindTopFeatures(in_dat, min.cutoff = 'q0')
	#in_dat <- RunSVD(in_dat)

	#generate gene activity for just variable genes, add to object and normalize
	#ga_mat<-GeneActivity(in_dat, features = feat)
	#in_dat[['RNA']] <- CreateAssayObject(counts = ga_mat)
	#in_dat<- NormalizeData(
	#  object = in_dat,
	#  assay = 'RNA',
	#  normalization.method = 'LogNormalize',
	#  scale.factor = median(in_dat$nCount_RNA)
	#)

	DefaultAssay(ref_dat)<-"GeneActivity"
	transfer.anchors <- FindTransferAnchors(
	  reference = ref_dat,
	  query = in_dat,
	  reduction = 'cca'
	)
	saveRDS(transfer.anchors,paste0(prefix,".transferanchors.rds"))
}

sample_label_transfer<-function(in_dat,ref_dat,transfer.anchors.,prefix="Tcell_",transfer_label="celltype"){
  predictions<- TransferData(
    anchorset = transfer.anchors.,
    refdata = ref_dat@meta.data[,transfer_label],
    weight.reduction = in_dat[["pca"]],
    dims=2:30
  )
  colnames(predictions)<-paste0(prefix,colnames(predictions))
  in_dat<-AddMetaData(in_dat,metadata=predictions)
  return(in_dat)
}

#seurat clusters 0: OPC_BCL11B; 1:OPC_MAG
opc_cellstate<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.OPCs.cistopic.seurat_clusters.annot",header=F)
opc_cellstate<-setNames(opc_cellstate$V2,nm=opc_cellstate$V1)
opc_cellstate[opc_cellstate==0]<-"OPC_BCL11B"
opc_cellstate[opc_cellstate==1]<-"OPC_MAG"

hgap_hippo<-readRDS("hgap_hippo_100perc.gene_activity.Rds")
hgap_cortex<-readRDS("hgap_cortex_100perc.gene_activity.Rds")

hgap_hippo<-AddMetaData(hgap_hippo,opc_cellstate,col.name="opc_subtype")
hgap_cortex<-AddMetaData(hgap_cortex,opc_cellstate,col.name="opc_subtype")

hgap_cortex<-subset(hgap_cortex,opc_subtype %in% c("OPC_BCL11B","OPC_MAG"))
hgap_hippo<-subset(hgap_hippo,opc_subtype %in% c("OPC_BCL11B","OPC_MAG"))
DefaultAssay(hgap_cortex)<-"GeneActivity"
DefaultAssay(hgap_hippo)<-"GeneActivity"

Idents(hgap_cortex)<-hgap_cortex$opc_subtype
Idents(hgap_hippo)<-hgap_hippo$opc_subtype

hgap_cortex<-FindVariableFeatures(hgap_cortex) 
hgap_hippo<-FindVariableFeatures(hgap_hippo) 

#OPC subtype labelling of brainspan cortex
	cortex_brainspan<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
	cortex_brainspan<-subset(cortex_brainspan,subclass_label=="OPC")
	generate_transfer_anchor_atac_to_rna(ref_dat=hgap_cortex,in_dat=cortex_brainspan,prefix="cortex_brainspan.opc",feat=VariableFeatures(hgap_cortex))
	cortex_brainspan<-sample_label_transfer(ref_dat=hgap_cortex,in_dat=cortex_brainspan,
			prefix="opc_cellstate_",
			transfer.anchors.=readRDS("cortex_brainspan.opc.transferanchors.rds"),
			transfer_label="opc_subtype")
	plt<-FeaturePlot(cortex_brainspan,feature=c("BCL11B","MAG","opc_cellstate_prediction.score.OPC_BCL11B","opc_cellstate_prediction.score.OPC_MAG"))
	ggsave(plt,file="cortex_brainspan.opc.pdf",width=10)
	system("slack -F cortex_brainspan.opc.pdf ryan_todo")
	plt<-ggplot(cortex_brainspan@meta.data,aes(x=opc_cellstate_prediction.score.OPC_BCL11B,y=opc_cellstate_prediction.score.OPC_MAG,color=subclass_label))+geom_point()+theme_bw()
	ggsave(plt,file="cortex_brainspan.opc.pdf",width=10)
	system("slack -F cortex_brainspan.opc.pdf ryan_todo")
	#table(cortex_brainspan$opc_cellstate_predicted.id)
	#OPC_BCL11B 
	#       283 
	#looks like no MAG in this, but its also a really low cell count, so it is consistent with the percentages we see

#OPC subtype labelling of hippocampus lifespan data
	hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.rds")
	hippo_lifespan<-subset(hippo_lifespan,MajorCellTypes=="OPC") #11614 cells
	generate_transfer_anchor_atac_to_rna(ref_dat=hgap_hippo,in_dat=hippo_lifespan,prefix="hippo_lifespan.opc",feat=VariableFeatures(hgap_hippo))
	hippo_lifespan<-sample_label_transfer(ref_dat=hgap_hippo,in_dat=hippo_lifespan,
			prefix="opc_cellstate_",
			transfer.anchors.=readRDS("hippo_lifespan.opc.transferanchors.rds"),
			transfer_label="opc_subtype")
	plt<-FeaturePlot(hippo_lifespan,feature=c("BCL11B","MAG","opc_cellstate_prediction.score.OPC_BCL11B","opc_cellstate_prediction.score.OPC_MAG"))
	ggsave(plt,file="hippo_lifespan.opc.pdf",width=10)
	system("slack -F hippo_lifespan.opc.pdf ryan_todo")
	plt<-ggplot(hippo_lifespan@meta.data,aes(x=opc_cellstate_prediction.score.OPC_BCL11B,y=opc_cellstate_prediction.score.OPC_MAG,color=MajorCellTypes))+geom_point()+theme_bw()
	ggsave(plt,file="hippo_lifespan.opc.pdf",width=10)
	system("slack -F hippo_lifespan.opc.pdf ryan_todo")
#All assigned to BCL11B as well...

#OPC subtype labelling of hippocampus atlas data
	hippo_axis<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.rds")
	hippo_axis<-subset(hippo_axis,cells=row.names(hippo_axis@meta.data[which(startsWith(hippo_axis$Cluster,prefix="OPC")),])) #cells
	#generate_transfer_anchor_atac_to_rna(ref_dat=hgap_hippo,in_dat=hippo_axis,prefix="hippo_axis.opc",feat=VariableFeatures(hgap_hippo))
	hippo_axis<-sample_label_transfer(ref_dat=hgap_hippo,in_dat=hippo_axis,
			prefix="opc_cellstate_",
			transfer.anchors.=readRDS("hippo_axis.opc.transferanchors.rds"),
			transfer_label="opc_subtype")
	plt<-FeaturePlot(hippo_axis,feature=c("BCL11B","MAG","opc_cellstate_prediction.score.OPC_BCL11B","opc_cellstate_prediction.score.OPC_MAG"))
	ggsave(plt,file="hippo_axis.opc.pdf",width=10)
	system("slack -F hippo_axis.opc.pdf ryan_todo")
	plt<-ggplot(hippo_axis@meta.data,aes(x=opc_cellstate_prediction.score.OPC_BCL11B,y=opc_cellstate_prediction.score.OPC_MAG,color=Cluster))+geom_point()+theme_bw()
	ggsave(plt,file="hippo_axis.opc.pdf",width=10)
	system("slack -F hippo_axis.opc.pdf ryan_todo")


```


# ATAC
For ATAC data, we will make a new counts matrix using our peak set for integration. This means we will go back to FASTQ format and make sure data sets are aligned on the same reference genome. Also we will limit our download to just the regions we are interested in. 

Cell specific ID and cluster assignment. 
```bash 
cd /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7606627/bin/NIHMS1630836-supplement-1630836_Supp_DataSet2.xlsx

#Cluster breakdown (saving as cluster_id.tsv)
"""
Cluster	Cell_Type_Group	Cluster_Description
Cluster1	ExcitatoryNeurons	Isocortical Excitatory
Cluster2	InhibitoryNeurons	Striatal Inhibitory (Major)
Cluster3	ExcitatoryNeurons	Hippocampal Excitatory
Cluster4	ExcitatoryNeurons	Hippocampal Excitatory
Cluster5	NigralNeurons	Nigral Neurons
Cluster6	NigralNeurons	Nigral Neurons
Cluster7	UnknownNeurons	Neurons (Unclassified)
Cluster8	OPCs	OPCs
Cluster9	OPCs	OPCs
Cluster10	OPCs	Nigral OPCs
Cluster11	InhibitoryNeurons	Isocortical Inhibitory
Cluster12	InhibitoryNeurons	Striatal Inhibitory (Minor)
Cluster13	Astrocytes	Astrocytes (Unclassified)
Cluster14	Astrocytes	Nigral Astrocytes
Cluster15	Astrocytes	Isocortical Astrocytes
Cluster16	Astrocytes	Striatal Astrocytes
Cluster17	Astrocytes	Astrocytes (Unclassified)
Cluster18	Doublets	Doublets
Cluster19	Oligodendrocytes	Oligodendrocytes
Cluster20	Oligodendrocytes	Oligodendrocytes
Cluster21	Oligodendrocytes	Oligodendrocytes
Cluster22	Oligodendrocytes	Oligodendrocytes
Cluster23	Oligodendrocytes	Oligodendrocytes
Cluster24	Microglia	Microglia"""


```
Getting the metadata ready to go!
```R
library(readxl)
dat<-read_excel("NIHMS1630836-supplement-1630836_Supp_DataSet2.xlsx",sheet=1,skip=21)
dat<-as.data.frame(dat)
cluster_id<-read.table("cluster_id.tsv",sep="\t",header=T)
dat$Cell_Type_Group<-cluster_id[match(dat$Cluster,cluster_id$Cluster) ,]$Cell_Type_Group
dat$Cluster_Description<-cluster_id[match(dat$Cluster, cluster_id$Cluster) ,]$Cluster_Description
write.table(dat,file="metadata.tsv",sep="\t",col.names=T,row.names=F,quote=F)

table(paste(dat$Cluster,dat$Cluster_Description))
```

```bash
awk 'OFS="\t" {if($2=="MDFG") print $3,$2}' metadata.tsv > mfg.annot
awk 'OFS="\t" {if($2=="HIPP") print $3,$2}' metadata.tsv > hippo.annot
```
Testing if we get the expected cell count from barcode 4 (i.e. all barcode 4 is unique in 10x run and corrected). We don't looks like the read 4 fastq files are not corrected. Just using exact matches.

```bash
zcat SRR11442503.sra_4.fastq.gz | grep -v "[@+;F:,]" | sort -T . | uniq -c > barc2_test.txt
#Looks like the data is not read corrected. but I'm just going to go with exact matches for now. I think it will be sufficient for integration.
```

Subset to just Middle Frontal Gyrus and Hippocampus instead.
Downloaded metadata from SRA Run Selector
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147672
```bash
ref_dir="/home/groups/CEDAR/mulqueen/human_brain_ref"
mkdir ${ref_dir}/corces_gwas
cd ${ref_dir}/corces_gwas
mkdir hippo
mkdir mfg

#Metadata from SRA Run Selector: Metadata_SRA_AccList.txt (copy and pasted from internet)
#Filter to just Hippocampus and Middle Frontal Gyrus:
#read in list of SRR from acc list and prefetch
#barc3 is PCR index for 10x run (not that important), doesn't need to be corrected. 
#barc4 is cell specific index

#hippo
prefetch --type fastq -X 100G -O ${ref_dir}/corces_gwas/hippo SRR11442501
prefetch --type fastq -X 100G -O ${ref_dir}/corces_gwas/hippo SRR11442502
#mfg
prefetch --type fastq -X 100G -O ${ref_dir}/corces_gwas/mfg SRR11442503

#make fastq
cd ${ref_dir}/corces_gwas/hippo
for i in *sra; do fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 20 -m 100G -O ${ref_dir}/corces_gwas/hippo $i; done
cd ${ref_dir}/corces_gwas/mfg
for i in *sra; do fasterq-dump --include-technical -t . -p -c 1G -b 1G -S -e 30 -m 100G -O ${ref_dir}/corces_gwas/mfg $i; done

#gzip fastq output
gzip ${ref_dir}/corces_gwas/hippo/*fastq &
gzip ${ref_dir}/corces_gwas/mfg/*fastq &

```

Adjusting fastq readname format to jive with scitools.

Note the script changes relative to sci-platform. Indexes are already in the fastq files 3 and 4, just need to rearrange format. 
barcode_to_scitools.py
```python
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

fq1=sys.argv[1] #read argument fq1="SRR11442503.sra_1.fastq.gz"
fq2=sys.argv[2] #read argument fq2="SRR11442503.sra_2.fastq.gz"
idx3=sys.argv[3] #SRR11442503.sra_3.fastq.gz
idx4=sys.argv[4] #SRR11442503.sra_4.fastq.gz
i=0
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1:
	with gzip.open(fq2, "rt") as handle2:
		with gzip.open(idx3, "rt") as handle3:
			with gzip.open(idx4, "rt") as handle4:
				with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
					with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
						for (title1, seq1, qual1), (title2, seq2, qual2), (title3,seq3,qual3), (title3,seq4,qual4) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2),FastqGeneralIterator(handle3),FastqGeneralIterator(handle4)):
							i+=1
							readname=seq3+seq4
							outfile_fq1.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq1, qual1))
							outfile_fq2.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq2, qual2))
```

Now running this python script for both fastq 1 and fastq 2.

```bash
python ./barcode_to_scitools.py ./mfg/SRR11442503.sra_1.fastq.gz ./mfg/SRR11442503.sra_2.fastq.gz ./mfg/SRR11442503.sra_3.fastq.gz ./mfg/SRR11442503.sra_4.fastq.gz &

python ./barcode_to_scitools.py ./hippo/SRR11442501.sra_1.fastq.gz ./hippo/SRR11442501.sra_2.fastq.gz ./hippo/SRR11442501.sra_3.fastq.gz ./hippo/SRR11442501.sra_4.fastq.gz &

python ./barcode_to_scitools.py ./hippo/SRR11442502.sra_1.fastq.gz ./hippo/SRR11442502.sra_2.fastq.gz ./hippo/SRR11442502.sra_3.fastq.gz ./hippo/SRR11442502.sra_4.fastq.gz &
```

Alignment
```bash
scitools="/home/groups/oroaklab/src/scitools/scitools-dev/scitools"
mfg_dir="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg"
hippo_dir="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo"
ref_outdir="/home/groups/CEDAR/mulqueen/human_brain_ref"

bwa mem -t 10 /home/groups/oroaklab/refs/hg38/hg38.fa \
$mfg_dir/SRR11442503.sra_1.barc.fastq.gz \
$mfg_dir/SRR11442503.sra_2.barc.fastq.gz 2>> $mfg_dir/mfg.align.log | samtools view -bSu - > mfg.nsrt.bam 2>> $mfg_dir/mfg.align.log

bwa mem -t 40 /home/groups/oroaklab/refs/hg38/hg38.fa \
$hippo_dir/SRR11442501.sra_1.barc.fastq.gz \
$hippo_dir/SRR11442501.sra_2.barc.fastq.gz 2>> $hippo_dir/hippo1.align.log | samtools view -bSu - > hippo1.nsrt.bam 2>> $hippo_dir/hippo1.align.log &

bwa mem -t 20 /home/groups/oroaklab/refs/hg38/hg38.fa \
$hippo_dir/SRR11442502.sra_1.barc.fastq.gz \
$hippo_dir/SRR11442502.sra_2.barc.fastq.gz 2>> $hippo_dir/hippo2.align.log | samtools view -bSu - > hippo2.nsrt.bam 2>> $hippo_dir/hippo2.align.log &
```
Sort by readname
```bash
samtools sort -@ 10 -m 5G -T . -n hippo1.nsrt.bam > hippo1.nsrt.sort.bam &
samtools sort -@ 10 -m 5G -T . -n hippo2.nsrt.bam > hippo2.nsrt.sort.bam &
samtools sort -@ 10 -m 5G -T . -n mfg.nsrt.bam > mfg.nsrt.sort.bam &
```

Dedup
```bash
scitools="/home/groups/oroaklab/src/scitools/scitools-dev/scitools"
mfg_dir="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg"
hippo_dir="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo"

$scitools bam-rmdup -n -r $mfg_dir/mfg.nsrt.sort.bam &
$scitools bam-rmdup -n -r $hippo_dir/hippo1.nsrt.sort.bam &
$scitools bam-rmdup -n -r $hippo_dir/hippo2.nsrt.sort.bam &

```

Reorder read name for easier processing (throwing PCR index after colon)
```bash
bam_in="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo1.sort.bbrd.q10.nsrt.bam"
((samtools view -H $bam_in) & (samtools view $bam_in | awk 'OFS="\t" {split($1,a,":");pcr_idx=substr(a[1],1,8);gem_idx=substr(a[1],9);$1=gem_idx":"pcr_idx" "NR; print $0}')) | samtools view -b - > ${bam_in::-4}.readname.bam &

bam_in="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo2.sort.bbrd.q10.nsrt.bam"
((samtools view -H $bam_in) & (samtools view $bam_in | awk 'OFS="\t" {split($1,a,":");pcr_idx=substr(a[1],1,8);gem_idx=substr(a[1],9);$1=gem_idx":"pcr_idx" "NR; print $0}')) | samtools view -b - > ${bam_in::-4}.readname.bam &

bam_in="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg/mfg.sort.bbrd.q10.nsrt.bam"
((samtools view -H $bam_in) & (samtools view $bam_in | awk 'OFS="\t" {split($1,a,":");pcr_idx=substr(a[1],1,8);gem_idx=substr(a[1],9);$1=gem_idx":"pcr_idx" "NR; print $0}')) | samtools view -b - > ${bam_in::-4}.readname.bam &
```
Filter bam files to those passing qc via the metadata
```bash
$scitools bam-filter -A /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg.annot -a MDFG $mfg_dir/mfg.sort.bbrd.q10.nsrt.readname.bam &
$scitools bam-filter -A /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo.annot -a HIPP $hippo_dir/hippo1.sort.bbrd.q10.nsrt.readname.bam &
$scitools bam-filter -A /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo.annot -a HIPP $hippo_dir/hippo2.sort.bbrd.q10.nsrt.readname.bam &
```

Make counts
```bash
input_bed="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/PEAKS_DataFreeze_Cortex_and_Hippocampus.500.bed"
scitools="/home/groups/oroaklab/src/scitools/scitools-dev/scitools"

#Make counts matrix on our peaks
input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo1.sort.bbrd.q10.nsrt.readname.filt.bam"
$scitools atac-counts \
-O ${input_bam::-4} \
$input_bam \
$input_bed &

input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo2.sort.bbrd.q10.nsrt.readname.filt.bam"
$scitools atac-counts \
-O ${input_bam::-4} \
$input_bam \
$input_bed &

input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg/mfg.sort.bbrd.q10.nsrt.readname.filt.bam"
$scitools atac-counts \
-O ${input_bam::-4} \
$input_bam \
$input_bed &
```
Make Seurat Objects
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/")
input_dir="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/"
meta.data<-read.table(paste0(input_dir,"metadata.tsv"),header=T,sep="\t")

mfg_met<-meta.data[meta.data$Region=="MDFG",]
mfg_met<-mfg_met[!isNA(mfg_met$X10x_SingleCell_Barcode),]
row.names(mfg_met)<-mfg_met$X10x_SingleCell_Barcode

hipp_met<-meta.data[meta.data$Region=="HIPP",]
hipp_met<-hipp_met[!isNA(hipp_met$X10x_SingleCell_Barcode),]
hipp_met_1<-hipp_met[hipp_met$Donor_ID=="11_0393",]
row.names(hipp_met_1)<-hipp_met_1$X10x_SingleCell_Barcode
hipp_met_2<-hipp_met[hipp_met$Donor_ID=="14_0586",]
row.names(hipp_met_2)<-hipp_met_2$X10x_SingleCell_Barcode


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead


#function to read in sparse matrix format from atac-count
read_in_sparse<-function(x){ #x is character file prefix followed by .bbrd.q10.500.counts.sparseMatrix.values.gz
    IN<-as.matrix(read.table(paste0(x,".counts.sparseMatrix.values.gz")))
    IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
    COLS<-read.table(paste0(x,".counts.sparseMatrix.cols.gz"))
    colnames(IN)<-COLS$V1
    ROWS<-read.table(paste0(x,".counts.sparseMatrix.rows.gz"))
    row.names(IN)<-ROWS$V1
    return(IN)
}

#function to make seurat object
make_seurat_object<-function(input_sparse="/mfg/mfg.sort.bbrd.q10.nsrt.readname.filt",outname="mfg",meta_in){
	counts<-read_in_sparse(paste0(input_dir,input_sparse))
	counts<-counts[,order(colnames(counts))]

	#Generate ChromatinAssay Objects
	chromatinassay <- CreateChromatinAssay(
	  counts = counts,
	  #genome="hg38",
	  min.cells = 1,
	  annotation=annotation,
	  sep=c("_","_")
	  )

	#Create Seurat Objects
	dat <- CreateSeuratObject(
	  counts = chromatinassay,
	  assay = "peaks"
	)

	meta_in<-meta_in[row.names(meta_in) %in% colnames(counts),]
	dat<-AddMetaData(dat,metadata=meta_in)

	#saving unprocessed SeuratObjects
	saveRDS(dat,file=paste0(outname,".SeuratObject.Rds"))
	print(paste("Finished",outname))
}

make_seurat_object(input_sparse="/mfg/mfg.sort.bbrd.q10.nsrt.readname.filt",outname="mfg",meta_in=mfg_met)
make_seurat_object(input_sparse="/hippo/hippo1.sort.bbrd.q10.nsrt.readname.filt",outname="hippo1",meta_in=hipp_met_1)
make_seurat_object(input_sparse="/hippo/hippo2.sort.bbrd.q10.nsrt.readname.filt",outname="hippo2",meta_in=hipp_met_2)

celltype_umap_plot<-function(in_dat,color_by="",outname){
	in_dat<-RunTFIDF(in_dat)
	in_dat<-FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)
	in_dat <- RunUMAP(in_dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

	plt1<-DimPlot(in_dat,group.by=color_by)#,reduction="umap.atac")
	ggsave(plt1,file=paste0(outname,".umap.png"),limitsize=F)
	system(paste0("slack -F ",outname,".umap.png"," ryan_todo"))
}

celltype_umap_plot(in_dat=readRDS("mfg.SeuratObject.Rds"),color_by="Cluster_Description",outname="mfg")
celltype_umap_plot(in_dat=readRDS("hippo1.SeuratObject.Rds"),color_by="Cluster_Description",outname="hippo1")
celltype_umap_plot(in_dat=readRDS("hippo2.SeuratObject.Rds"),color_by="Cluster_Description",outname="hippo2")
```

Integrate ATAC profiles with Public Data Sets
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
library(harmony)

#Function to normalize merged data
normalize_and_umap<-function(in_dat){
	in_dat<-RunTFIDF(in_dat)
	in_dat<-FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)
	in_dat <- RunUMAP(in_dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
	return(in_dat)
}

#Read in Corces et al Data
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/")
mfg_atac<-readRDS("mfg.SeuratObject.Rds");mfg_atac$orig.ident<-"Corces_MFG"
#Read in HGAP data
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
hgap_cortex<-readRDS("cortex.SeuratObject.Rds")
#Merge Seurat Objects
cortex<-merge(hgap_cortex,y=mfg_atac,add.cell.ids=c("HGAP","Corces_etal"))
#Perform Harmony Integration
cortex<-normalize_and_umap(cortex)
cortex<-RunHarmony(cortex,group.by.vars="orig.ident",reduction.save="harmony_atac",assay.use="peaks",reduction="lsi",project.dim=F)
cortex<-RunUMAP(cortex,reduction.name="harmonyumap_rna",reduction = "harmony_atac",dims=2:dim(cortex@reductions$harmony_atac)[2]) 
cortex$assigned_celltype<-cortex$Putative_Celltype
cortex@meta.data[cortex$orig.ident!="SeuratProject",]$assigned_celltype<-cortex@meta.data[cortex$orig.ident!="SeuratProject",]$Cluster_Description
#Plot
plt1<-DimPlot(cortex,group.by="orig.ident",reduction="umap.atac",raster=F)
plt2<-DimPlot(cortex,group.by="orig.ident",reduction="harmonyumap_rna",raster=F)
plt3<-DimPlot(cortex,group.by="Putative_Celltype",reduction="harmonyumap_rna",raster=F)
plt4<-DimPlot(cortex,group.by="Cluster_Description",reduction="harmonyumap_rna",raster=F)
plt5<-DimPlot(cortex,group.by="assigned_celltype",reduction="harmonyumap_rna",split.by="orig.ident",raster=FALSE)
ggsave((plt1|plt2)/(plt3|plt4),file="cortex.atac.integrated.umap.png",width=15,,height=15,limitsize=F)
system(paste0("slack -F ","cortex.atac.integrated.umap.png"," ryan_todo"))

ggsave(plt1,file="cortex.atac.preintegrated.umap.pdf",width=10,,height=10,limitsize=F)
system(paste0("slack -F ","cortex.atac.preintegrated.umap.pdf"," ryan_todo"))

ggsave(plt2,file="cortex.atac.integrated.umap.pdf",width=10,,height=10,limitsize=F)
system(paste0("slack -F ","cortex.atac.integrated.umap.pdf"," ryan_todo"))

ggsave(plt5,file="cortex.atac.integrated_assignedcelltypes.umap.pdf",width=20,,height=10,limitsize=F)
system(paste0("slack -F ","cortex.atac.integrated_assignedcelltypes.umap.pdf"," ryan_todo"))


saveRDS(cortex,"cortex.atac_integrated.SeuratObject.Rds")


#Hippocampus
#Read in Corces et al Data
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/")
hippo1_atac<-readRDS("hippo1.SeuratObject.Rds");hippo1_atac$orig.ident<-"Corces_Hippo1"
hippo2_atac<-readRDS("hippo2.SeuratObject.Rds");hippo2_atac$orig.ident<-"Corces_Hippo2"
#Read in HGAP data
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
hgap_hippo<-readRDS("hippocampus.SeuratObject.Rds")
#Merge Seurat Objects
hippo<-merge(hgap_hippo,y=c(hippo1_atac,hippo2_atac),add.cell.ids=c("HGAP","Corces_hip1","Corces_hip2"))
#Perform Harmony Integration
hippo<-normalize_and_umap(hippo)
hippo<-RunHarmony(hippo,group.by.vars="orig.ident",reduction.save="harmony_atac",assay.use="peaks",reduction="lsi",project.dim=F)
hippo<-RunUMAP(hippo,reduction.name="harmonyumap_rna",reduction = "harmony_atac",dims=2:dim(hippo@reductions$harmony_atac)[2]) 
hippo$assigned_celltype<-hippo$Putative_Celltype
hippo@meta.data[hippo$orig.ident!="SeuratProject",]$assigned_celltype<-hippo@meta.data[hippo$orig.ident!="SeuratProject",]$Cluster_Description
hippo@meta.data[hippo$orig.ident!="SeuratProject",]$orig.ident<-"Corces_Hippo"

#Plot
plt1<-DimPlot(hippo,group.by="orig.ident",reduction="umap.atac",raster=F)
plt2<-DimPlot(hippo,group.by="orig.ident",reduction="harmonyumap_rna",raster=F)
plt3<-DimPlot(hippo,group.by="Putative_Celltype",reduction="harmonyumap_rna",raster=F)
plt4<-DimPlot(hippo,group.by="Cluster_Description",reduction="harmonyumap_rna",raster=F)
plt5<-DimPlot(hippo,group.by="assigned_celltype",reduction="harmonyumap_rna",split.by="orig.ident",raster=FALSE)

ggsave((plt1|plt2)/(plt3|plt4),file="hippo.atac.integrated.umap.png",width=15,,height=15,limitsize=F)
system(paste0("slack -F ","hippo.atac.integrated.umap.png"," ryan_todo"))

ggsave(plt1,file="hippo.atac.preintegrated.umap.pdf",width=10,,height=10,limitsize=F)
system(paste0("slack -F ","hippo.atac.preintegrated.umap.pdf"," ryan_todo"))

ggsave(plt2,file="hippo.atac.integrated.umap.pdf",width=10,,height=10,limitsize=F)
system(paste0("slack -F ","hippo.atac.integrated.umap.pdf"," ryan_todo"))

ggsave(plt5,file="hippo.atac.integrated_assignedcelltypes.umap.pdf",width=20,,height=10,limitsize=F)
system(paste0("slack -F ","hippo.atac.integrated_assignedcelltypes.umap.pdf"," ryan_todo"))


saveRDS(hippo,"hippo.atac_integrated.SeuratObject.Rds")
```

Integration of ATAC-ATAC data and uncovering OPC and astrocyte biology

Prepare Files
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
library(harmony)
library(cisTopic)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs")

#Using same cistopic set up of dimensions as original OPC distinctions
cistopic_generation<-function(in_dat,outname="opc.integrated"){
	atac_sub<-in_dat
   cistopic_counts_frmt<-atac_sub@assays$peaks@counts
   row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
   sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
   print("made cistopic object")
   sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(25,28,30,33,35,40),nCores=6,addModels=FALSE)
   saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  
   saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(outname,".CisTopicObject.Rds"))
   sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
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
  return(atac_sub)
  }

#Function to normalize merged data
normalize_and_umap<-function(in_dat){
	in_dat<-RunTFIDF(in_dat)
	in_dat<-FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)
	in_dat <- RunUMAP(in_dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
	return(in_dat)
}

#seurat clusters 0: OPC_BCL11B; 1:OPC_MAG
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
hgap_hippo<-readRDS("hippocampus.SeuratObject.Rds");hgap_hippo$orig.ident<-"HGAP_hippo"
hgap_cortex<-readRDS("cortex.SeuratObject.Rds");hgap_cortex$orig.ident<-"HGAP_cortex"

opc_cellstate<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.OPCs.cistopic.seurat_clusters.annot",header=F)
opc_cellstate<-setNames(opc_cellstate$V2,nm=opc_cellstate$V1)
opc_cellstate[opc_cellstate==0]<-"OPC_BCL11B"
opc_cellstate[opc_cellstate==1]<-"OPC_MAG"

hgap_hippo<-AddMetaData(hgap_hippo,opc_cellstate,col.name="opc_subtype")
hgap_cortex<-AddMetaData(hgap_cortex,opc_cellstate,col.name="opc_subtype")
hgap_hippo<-subset(hgap_hippo,opc_subtype!="NA")
hgap_cortex<-subset(hgap_cortex,opc_subtype!="NA")

#Read in Corces et al Data
mfg_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg.SeuratObject.Rds");mfg_atac$orig.ident<-"Corces_MFG"
hippo1_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo1.SeuratObject.Rds");hippo1_atac$orig.ident<-"Corces_Hippo1"
hippo2_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo2.SeuratObject.Rds");hippo2_atac$orig.ident<-"Corces_Hippo2"
mfg_atac<-subset(mfg_atac,Cluster_Description=="OPCs")
hippo1_atac<-subset(hippo1_atac,Cluster_Description=="OPCs")
hippo2_atac<-subset(hippo2_atac,Cluster_Description=="OPCs")

merged<-merge(hgap_hippo,y=c(hgap_cortex,mfg_atac,hippo1_atac,hippo2_atac),add.cell.ids=c("HGAP_hip","HGAP_cortex","Corces_mfg","Corces_hip1","Corces_hip2"))
merged$Region<-"Cortex"
merged@meta.data[grepl("hip",row.names(merged@meta.data)),]$Region<-"Hippocampus"
merged@meta.data[grepl("HGAP",row.names(merged@meta.data)),]$orig.ident<-"HGAP"
merged$assay<-"HGAP"
merged@meta.data[grepl("Corces",row.names(merged@meta.data)),]$assay<-"Corces"
saveRDS(merged,file="opc.atac.preintegration.Rds")

#astrocytes
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
hgap_hippo<-readRDS("hippocampus.SeuratObject.Rds");hgap_hippo$orig.ident<-"HGAP_hippo"
hgap_cortex<-readRDS("cortex.SeuratObject.Rds");hgap_cortex$orig.ident<-"HGAP_cortex"

astro_cellstate_label<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/astro.seurat_clusters_marker.rename.annot",header=F)
astro_cellstate<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/combined.astroependymocytes.harmony.seurat_clusters.annot",header=F)
matches<-data.frame(cellID=astro_cellstate$V1,astro_subtype=astro_cellstate_label[match(astro_cellstate$V2, astro_cellstate_label$V1),]$V2)
astro_cellstate<-setNames(matches$astro_subtype,nm=matches$cellID)

hgap_hippo<-AddMetaData(hgap_hippo,astro_cellstate,col.name="astro_subtype")
hgap_cortex<-AddMetaData(hgap_cortex,astro_cellstate,col.name="astro_subtype")

hgap_hippo<-subset(hgap_hippo,astro_subtype!="NA")
hgap_cortex<-subset(hgap_cortex,astro_subtype!="NA")

#Read in Corces et al Data
mfg_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg.SeuratObject.Rds");mfg_atac$orig.ident<-"Corces_MFG"
hippo1_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo1.SeuratObject.Rds");hippo1_atac$orig.ident<-"Corces_Hippo1"
hippo2_atac<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo2.SeuratObject.Rds");hippo2_atac$orig.ident<-"Corces_Hippo2"
mfg_atac<-subset(mfg_atac,cells=row.names(mfg_atac@meta.data[grepl("Astro",mfg_atac$Cluster_Description),]))
hippo1_atac<-subset(hippo1_atac,cells=row.names(hippo1_atac@meta.data[grepl("Astro",hippo1_atac$Cluster_Description),]))
hippo2_atac<-subset(hippo2_atac,cells=row.names(hippo2_atac@meta.data[grepl("Astro",hippo2_atac$Cluster_Description),]))

merged<-merge(hgap_hippo,y=c(hgap_cortex,mfg_atac,hippo1_atac,hippo2_atac),add.cell.ids=c("HGAP_hip","HGAP_cortex","Corces_mfg","Corces_hip1","Corces_hip2"))
merged$Region<-"Cortex"
merged@meta.data[grepl("hip",row.names(merged@meta.data)),]$Region<-"Hippocampus"
merged@meta.data[grepl("HGAP",row.names(merged@meta.data)),]$orig.ident<-"HGAP"
merged$assay<-"HGAP"
merged@meta.data[grepl("Corces",row.names(merged@meta.data)),]$assay<-"Corces"
saveRDS(merged,file="astro.atac.preintegration.Rds")



```


```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(Matrix)
library(harmony)
library(cisTopic)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#Using same cistopic set up of dimensions as original OPC distinctions
cistopic_generation<-function(in_dat,outname="opc.integrated"){
	atac_sub<-in_dat
   cistopic_counts_frmt<-atac_sub@assays$peaks@counts
   row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
   sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
   print("made cistopic object")
   sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(25,28,30,33,35,40),nCores=6,addModels=FALSE)
   saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  
   saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  sub_cistopic_models<-readRDS(file=paste0(outname,".CisTopicObject.Rds"))
   sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
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
  return(atac_sub)
  }

#Function to normalize merged data
normalize_and_umap<-function(in_dat){
	in_dat<-RunTFIDF(in_dat)
	in_dat<-FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)
	in_dat <- RunUMAP(in_dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
	return(in_dat)
}

opc<-readRDS("opc.atac.preintegration.Rds")
opc<-normalize_and_umap(opc) #standard TFIDF and LSI normalization pipeline
opc<-cistopic_generation(opc) #cistopic dimensionality reduction
#Run projection before integration
opc<-readRDS("opc.integrated.SeuratObject.rds")


integrate_dat<-function(dat,prefix,group_by="opc_subtype",dims=c(2,3),harmony_vars=c("assay")){
	dat<-RunHarmony(dat,group.by.vars=harmony_vars,reduction.save="harmony_atac",assay.use="peaks",reduction="cistopic",dims.use=1:dim(dat@reductions$cistopic)[2],project.dim=F,max.iter.harmony=20)
	dat<-RunUMAP(dat,reduction.name="harmonyumap_cistopic",reduction = "harmony_atac",dims=1:dim(dat@reductions$harmony_atac)[2],n.components=3) 
	dat$assigned_celltype<-dat$Putative_Celltype
	dat@meta.data[dat$orig.ident!="SeuratProject",]$assigned_celltype<-dat@meta.data[dat$orig.ident!="SeuratProject",]$Cluster_Description

	#Plot
	plt1<-DimPlot(dat,dims=dims,group.by=group_by,reduction="harmonyumap_cistopic",raster=F)
	plt2<-DimPlot(dat,dims=dims,group.by="orig.ident",reduction="harmonyumap_cistopic",raster=F)
	plt3<-DimPlot(dat,dims=dims,group.by="Putative_Celltype",reduction="harmonyumap_cistopic",raster=F)
	plt4<-DimPlot(dat,dims=dims,group.by="Cluster",reduction="harmonyumap_cistopic",raster=F)
	plt5<-DimPlot(dat,dims=dims,group.by="Region",reduction="harmonyumap_cistopic",split.by="orig.ident",raster=FALSE)
	plt6<-DimPlot(dat,dims=dims,group.by=group_by,reduction="harmonyumap_cistopic",split.by="orig.ident",raster=FALSE)
	ggsave((plt1|plt2)/(plt3|plt4)/(plt5),file=paste0(prefix,".atac.integrated.umap.png"),width=15,height=20,limitsize=F)
	system(paste0("slack -F ",prefix,".atac.integrated.umap.png"," ryan_todo"))

	#Plot PDFs
	plot_name="HGAP_defined_subtype"
	ggsave(plt1,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"))
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	plot_name="dataset"
	ggsave(plt2,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"))
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	plot_name="HGAP_defined_celltype"
	ggsave(plt3,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"))
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	plot_name="Corces_cluster"
	ggsave(plt4,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"))
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	plot_name="dataset_by_region"
	ggsave(plt5,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"),width=20)
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	plot_name="dataset_by_subtype"
	ggsave(plt6,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"),width=20)
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	saveRDS(dat,paste0(prefix,".integrated.SeuratObject.rds"))
	dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))
	DefaultAssay(dat)<-"peaks"
	hga_dat<-subset(dat,assay=="HGAP")
	corces_dat<-subset(dat,assay!="HGAP")
	corces_dat<- RunTFIDF(corces_dat)
	corces_dat<- FindTopFeatures(corces_dat, min.cutoff = 'q0')
	corces_dat<- RunSVD(corces_dat)
	Idents(hga_dat)<-hga_dat@meta.data[,group_by]
	hga_dat<-FindVariableFeatures(hga_dat)

	transfer.anchors <- FindTransferAnchors(
	  reference = hga_dat,
	  query = corces_dat,
	  features=VariableFeatures(hga_dat),
	  reference.assay='peaks',
	  query.assay='peaks',
	  reduction = 'cca'
	)

	predicted.labels <- TransferData(
	  anchorset = transfer.anchors,
	  refdata = hga_dat@meta.data[,group_by],
	  weight.reduction = corces_dat[['lsi']],
	  dims = 2:30
	)

	corces_dat <- AddMetaData(object = corces_dat, metadata = predicted.labels)

	plot_name="CorcesCells_HGAP_defined_celltype"
	plt3<-DimPlot(corces_dat,dims=dims,group.by="predicted.id",reduction="harmonyumap_cistopic",raster=F)
	ggsave(plt3,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"))
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

	saveRDS(corces_dat,paste0(prefix,".corces.integrated.SeuratObject.rds"))
	table(corces_dat$predicted.id)

	dat<-subset(dat,assay=="HGAP")
	names_out<-colnames(dat)
	names_out<-names_out[grepl(names_out,pattern="HGAP")]
	names_out_cellid<-unlist(lapply(strsplit(names_out,"_"),"[[",3))
	name_conversion_table<-as.data.frame(cbind(names_out,names_out_cellid))
	write.table(name_conversion_table,file=paste0(prefix,".cellIDs.tsv"),col.names=F,quote=F,row.names=F,sep="\t")
}

integrate_dat(dat=opc,prefix="opc")
#integrate_dat(dat=astro,prefix="astro",group_by="astro_subtype",dims=c(1,2),harmony_vars=c("assay","Samples"))

prefix="opc"

dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))
	DefaultAssay(dat)<-"peaks"
	hga_dat<-subset(dat,assay=="HGAP")
	corces_dat<-subset(dat,assay!="HGAP")
	hga_dat<- RunUMAP(hga_dat, reduction = "harmony_atac", dims = 1:30, return.model = TRUE)

	Idents(hga_dat)<-hga_dat$opc_subtype
	# find transfer anchors
	transfer.anchors <- FindTransferAnchors(
	  reference = hga_dat,
	  query = corces_dat,
	  reference.reduction = "lsi",
	  reduction = "lsiproject",
	  dims = 1:30
	)

	# map query onto the reference dataset
	out <- MapQuery(
	  anchorset = transfer.anchors,
	  reference = hga_dat,
	  query = corces_dat,
	  refdata = hga_dat$opc_subtype)

	corces_dat<-AddMetaData(corces_dat,out$predicted.id,col.name="predicted.opc_subtype")

	plot_name="CorcesCells_HGAP_defined_celltype"
	plt3<-DimPlot(corces_dat,dims=dims,group.by="predicted.opc_subtype",reduction="harmonyumap_cistopic",raster=F,split.by="orig.ident")
	ggsave(plt3,file=paste0(prefix,".atac.integrated.",plot_name,".umap.pdf"),width=15)
	system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))


library(dplyr)
dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))
dat<-AddMetaData(dat,out$predicted.id,col.name="predicted.opc_subtype")
saveRDS(dat,paste0(prefix,".integrated.SeuratObject.rds"))
dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))

dat@meta.data %>% group_by(orig.ident,opc_subtype,predicted.opc_subtype) %>% summarize(count=n())

dat@meta.data[dat@meta.data$orig.ident!="HGAP",]$opc_subtype<-dat@meta.data[dat@meta.data$orig.ident!="HGAP",]$predicted.opc_subtype
plot_name="CorcesCells_HGAP_defined_celltype"

plt<-ggplot(dat@meta.data,aes(x=orig.ident,fill=opc_subtype))+geom_bar(position="fill")
ggsave(plt,file=paste0(prefix,".atac.integrated.",plot_name,".barplot.pdf"),width=15)
system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".barplot.pdf") ," ryan_todo"))

```


### Tabix fragment file generation

Making tabix fragment files for final coverage plots and gene activity calculations.
Tabix file format is a tab separated multicolumn data structure.

| Column Number | Name | Description |
|:--------|:-------:|:--------|
|1 |chrom |  Reference genome chromosome of fragment |
|2 |chromStart | Adjusted start position of fragment on chromosome. |
|3 |chromEnd   | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
|4 |barcode | The 10x (or sci) cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment. |
|5 |duplicateCount |The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeq instrument are excluded from this count. |


```bash
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"

#Make counts matrix on our peaks
cd /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/
input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo1.sort.bbrd.q10.nsrt.readname.filt.bam"
output_name="hippo1"
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $4,$5,$9,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
$tabix -p bed $output_name.fragments.tsv.gz &

input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo2.sort.bbrd.q10.nsrt.readname.filt.bam"
output_name="hippo2"
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $4,$5,$9,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
$tabix -p bed $output_name.fragments.tsv.gz &


cd /home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg/
input_bam="/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg/mfg.sort.bbrd.q10.nsrt.readname.filt.bam"
output_name="mfg"
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $4,$5,$9,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
$tabix -p bed $output_name.fragments.tsv.gz &

```

Calculate Gene Activity for Final Violin Plots
```R
library(Seurat)
library(Signac) 
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
set.seed(1234)
library(Matrix)
library(monocle3,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/") #using old install of monocle, just need for as.cell_data_set conversion
library(cicero,lib.loc="/home/groups/oroaklab/src/R/R-4.0.0/library/")
library(SeuratWrappers)
library(rtracklayer)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

prefix="opc"
dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))


split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

make_sparse_matrix <- function(data,
                               i.name = "Peak1",
                               j.name = "Peak2",
                               x.name = "value") {
  if(!i.name %in% names(data) |
     !j.name %in% names(data) |
     !x.name %in% names(data)) {
    stop('i.name, j.name, and x.name must be columns in data')
  }
  
  data$i <- as.character(data[,i.name])
  data$j <- as.character(data[,j.name])
  data$x <- data[,x.name]
  
  if(!class(data$x) %in%  c("numeric", "integer"))
    stop('x.name column must be numeric')
  
  peaks <- data.frame(Peak = unique(c(data$i, data$j)),
                      index = seq_len(length(unique(c(data$i, data$j)))))
  
  data <- data[,c("i", "j", "x")]
  
  data <- rbind(data, data.frame(i=peaks$Peak, j = peaks$Peak, x = 0))
  data <- data[!duplicated(data[,c("i", "j", "x")]),]
  data <- data.table::as.data.table(data)
  peaks <- data.table::as.data.table(peaks)
  data.table::setkey(data, "i")
  data.table::setkey(peaks, "Peak")
  data <- data[peaks]
  data.table::setkey(data, "j")
  data <- data[peaks]
  data <- as.data.frame(data)
  
  data <- data[,c("index", "i.index", "x")]
  data2 <- data
  names(data2) <- c("i.index", "index", "x")
  
  data <- rbind(data, data2)
  
  data <- data[!duplicated(data[,c("index", "i.index")]),]
  data <- data[data$index >= data$i.index,]
  
  sp_mat <- Matrix::sparseMatrix(i=as.numeric(data$index),
                                 j=as.numeric(data$i.index),
                                 x=data$x,
                                 symmetric = TRUE)
  
  colnames(sp_mat) <- peaks[order(peaks$index),]$Peak
  row.names(sp_mat) <- peaks[order(peaks$index),]$Peak
  return(sp_mat)
}

build_composite_gene_activity_matrix <- function(input_cds,
                                                 site_weights,
                                                 cicero_cons_info,
                                                 dist_thresh=250000,
                                                 coaccess_cutoff=0.25) {
    accessibility_mat <- exprs(input_cds)
    promoter_peak_table <- fData(input_cds)
    promoter_peak_table$peak <- as.character(row.names(promoter_peak_table))
    promoter_peak_table <-
        promoter_peak_table[!is.na(promoter_peak_table$gene),]
    promoter_peak_table <- promoter_peak_table[,c("peak", "gene")]
    promoter_peak_table$gene <- as.character(promoter_peak_table$gene)

    # Make site_weight matrix
    site_names <- names(site_weights)
    site_weights <- as(Matrix::Diagonal(x=as.numeric(site_weights)),
                      "sparseMatrix")
    row.names(site_weights) <- site_names
    colnames(site_weights) <- site_names

    # Find distance between cicero peaks. If distance already calculated, skip
    if ("dist" %in% colnames(cicero_cons_info) == FALSE) {
        Peak1_cols <- split_peak_names(cicero_cons_info$Peak1)
        Peak2_cols <- split_peak_names(cicero_cons_info$Peak2)
        Peak1_bp <- round((as.integer(Peak1_cols[,3]) +
                          as.integer(Peak1_cols[,2])) / 2)
        Peak2_bp <- round((as.integer(Peak2_cols[,3]) +
                          as.integer(Peak2_cols[,2])) / 2)
        cicero_cons_info$dist <- abs(Peak2_bp - Peak1_bp)
    }

    # Get connections between promoters and distal sites above coaccess
    # threshold
    nonneg_cons <-
        cicero_cons_info[(cicero_cons_info$Peak1 %in%
                          promoter_peak_table$peak |
                          cicero_cons_info$Peak2 %in%
                          promoter_peak_table$peak) &
                          cicero_cons_info$coaccess >= coaccess_cutoff &
                          cicero_cons_info$dist < dist_thresh,]
    nonneg_cons <- nonneg_cons[,c("Peak1", "Peak2", "coaccess")]
    nonneg_cons <- nonneg_cons[!duplicated(nonneg_cons),]

    nonneg_cons$Peak1 <- as.character(nonneg_cons$Peak1)
    nonneg_cons$Peak2 <- as.character(nonneg_cons$Peak2)

    nonneg_cons <- rbind(nonneg_cons,
                        data.frame(Peak1=unique(promoter_peak_table$peak),
                                   Peak2=unique(promoter_peak_table$peak),
                                   coaccess=0))

    # Make square matrix of connections from distal to proximal
    distal_connectivity_matrix <- make_sparse_matrix(nonneg_cons,
                                                    x.name="coaccess")

    # Make connectivity matrix of promoters versus all
    promoter_conn_matrix <-
        distal_connectivity_matrix[unique(promoter_peak_table$peak),]

    # Get list of promoter and distal sites in accessibility mat
    promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- intersect(colnames(promoter_conn_matrix),
                                     row.names(accessibility_mat))
    distal_safe_sites <- setdiff(distal_safe_sites, promoter_safe_sites)

    # Get accessibility info for promoters
    promoter_access_mat_in_cicero_map <- accessibility_mat[promoter_safe_sites,, drop=FALSE]

    # Get accessibility for distal sites
    distal_activity_scores <- accessibility_mat[distal_safe_sites,, drop=FALSE]

    # Scale connectivity matrix by site_weights
    scaled_site_weights <- site_weights[distal_safe_sites,distal_safe_sites, drop=FALSE]
    total_linked_site_weights <- promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    total_linked_site_weights <- 1/Matrix::rowSums(total_linked_site_weights,
                                                na.rm=TRUE)
    total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] <- 0
    total_linked_site_weights[is.na(total_linked_site_weights)] <- 0
    total_linked_site_weights[is.nan(total_linked_site_weights)] <- 0
    site_names <- names(total_linked_site_weights)
    total_linked_site_weights <- Matrix::Diagonal(x=total_linked_site_weights)
    row.names(total_linked_site_weights) <- site_names
    colnames(total_linked_site_weights) <- site_names
    scaled_site_weights <- total_linked_site_weights %*%
        promoter_conn_matrix[,distal_safe_sites, drop=FALSE] %*%
        scaled_site_weights
    scaled_site_weights@x[scaled_site_weights@x > 1] <- 1

    # Multiply distal accessibility by site weights
    distal_activity_scores <- scaled_site_weights %*% distal_activity_scores

    distal_activity_scores <-
        distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),, drop=FALSE]

    # Sum distal and promoter scores
    promoter_activity_scores <- distal_activity_scores +
        promoter_access_mat_in_cicero_map

    # Make and populate final matrix
    promoter_gene_mat <-
        Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                             i=as.numeric(factor(promoter_peak_table$gene)),
                             x=1)
    colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
    row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
    promoter_gene_mat <- promoter_gene_mat[,row.names(promoter_activity_scores)]
    gene_activity_scores <- promoter_gene_mat %*% promoter_activity_scores

    return(gene_activity_scores)
}

# gene annotation sample for gene activity
get_annot<-function(){
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
	return(gene_annotation)
	}


#Cicero processing function  
cicero_processing<-function(object_input=hgap_hippo,prefix="hgap_hippo_5perc",k=500){
      #Generate CDS format from Seurat object
      #atac.cds <- as.CellDataSet(object_input,assay="peaks",reduction="umap.atac")
      atac.cds <- as.cell_data_set(object_input,assay="peaks",reduction="umap.atac")
      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = object_input@reductions$umap.atac@cell.embeddings,k=k)

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      atac.cicero<-readRDS(paste(prefix,"atac_cicero_cds.Rds",sep="_"))

      genome<-SeqinfoForUCSCGenome("hg38")
      genome.df<-data.frame("chr"=seqnames(genome),"length"=seqlengths(genome))
      genome.df<-genome.df[genome.df$chr %in% paste0("chr",seq(1:22)),]

      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      conns<-readRDS(paste(prefix,"atac_cicero_conns.Rds",sep="_"))

      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      ccans<-readRDS(paste(prefix,"atac_cicero_ccans.Rds",sep="_"))

      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      saveRDS(object_input,file=paste0(prefix,".","gene_activity.Rds"))
      return(object_input)
}

geneactivity_processing<-function(object_input=hgap_hippo,conns_input=conns,prefix="hgap_hippo_5perc"){
			anno<-get_annot()
			atac.cds <- as.cell_data_set(object_input,assay="peaks",reduction="umap.atac")
			#atac.cds<-newCellDataSet(cellData=object_input@assays$peaks@counts,featureData=atac.cds@featureData,phenoData=atac.cds@phenoData))
			#cds <- new("CellDataSet", assayData = assayDataNew("environment", 
      #  exprs = object_input@assays$peaks@counts), phenoData = atac.cds@phenoData, featureData = atac.cds@featureData, 
      #  lowerDetectionLimit =  0.1, expressionFamily = VGAM::negbinomial.size(), 
      #  dispFitInfo = new.env(hash = TRUE))
      atac.cds<- annotate_cds_by_site(atac.cds, anno)
      fData(atac.cds)<-cbind(fData(atac.cds),
      		site_name=row.names(fData(atac.cds)),
      		chr=gsub(pattern="chr",replace="",unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",1))),
      		bp1=unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",2)),
      		bp2=unlist(lapply(strsplit(row.names(fData(atac.cds)),"-"),"[",3)))

			input_cds=atac.cds
			cicero_cons_info=conns
			site_weights=NULL
			dist_thresh=250000
			coaccess_cutoff=0.25

			accessibility_mat <- exprs(input_cds)
      site_weights <- Matrix::rowMeans(accessibility_mat) / Matrix::rowMeans(accessibility_mat)
      site_weights[names(site_weights)] <- 1
			gene_promoter_activity <- build_composite_gene_activity_matrix(input_cds,
                                             site_weights,
                                             cicero_cons_info,
                                             dist_thresh=dist_thresh,
                                             coaccess_cutoff=coaccess_cutoff)
			unnorm_ga<-gene_promoter_activity
      saveRDS(unnorm_ga,paste(prefix,"unnorm_GA.Rds",sep="."))
      object_input[['GeneActivity']]<- CreateAssayObject(counts = unnorm_ga) 

		  # normalize
		  object_input <- NormalizeData(
		    object = object_input,
		    assay = 'GeneActivity',
		    normalization.method = 'LogNormalize',
		    scale.factor = median(object_input$nCount_GeneActivity)
		  )
		  saveRDS(object_input,file=paste0(prefix,".","gene_activity.Rds"))
      return(object_input)
  		
}

#Full set processing on long request node 
#Gene activity processing
dat<-cicero_processing(object_input=dat,prefix="opc.integrated",k=50)
conns<-readRDS("opc.integrated_atac_cicero_conns.Rds")
dat<-geneactivity_processing(object_input=dat,conns_input=conns,prefix="opc.integrated")
saveRDS(dat,paste0(prefix,".integrated.SeuratObject.rds"))



#Add in fragments files

#Plot Coverage plot of MAG gene
prefix="opc"
dat<-readRDS(paste0(prefix,".integrated.SeuratObject.rds"))
dat<-subset(dat,orig.ident!="HGAP") #subset to just corces data
dat<-RenameCells(dat,new.names=unlist(lapply(strsplit(row.names(dat@meta.data),"_"),"[[",3))) #rename cells to exclude orig.ident
fpath <- "/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/mfg/mfg.fragments.tsv.gz"
cells <-row.names(dat@meta.data[dat@meta.data$orig.ident=="Corces_MFG",])
cells<-unlist(lapply(strsplit(cells,"_"),"[[",3))
frags_mfg <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)

fpath <- "/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo1.fragments.tsv.gz"
cells <-row.names(dat@meta.data[dat@meta.data$orig.ident=="Corces_Hippo1",])
cells<-unlist(lapply(strsplit(cells,"_"),"[[",3))
frags_hippo1 <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)

fpath <- "/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/hippo/hippo2.fragments.tsv.gz"
cells <-row.names(dat@meta.data[dat@meta.data$orig.ident=="Corces_Hippo2",])
cells<-unlist(lapply(strsplit(cells,"_"),"[[",3))
frags_hippo2 <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)

Fragments(dat@assays$peaks) <- c(frags_mfg,frags_hippo1,frags_hippo2)

Idents(dat)<-dat$predicted.opc_subtype
DefaultAssay(dat)<-"GeneActivity"
markers<-FindAllMarkers(dat,only.pos=TRUE) 
plot_name="CorcesCells_HGAP_defined_celltype"


DefaultAssay(dat)<-"peaks"
plt<-CoveragePlot(dat,region=c("MAG"),extend.upstream=5000,extend.downstream=5000,group.by="predicted.opc_subtype")
ggsave(plt,file=paste0(prefix,".atac.integrated.",plot_name,".vlnplot.pdf"),width=30)
system(paste0("slack -F ",paste0(prefix,".atac.integrated.",plot_name,".vlnplot.pdf") ," ryan_todo"))
```

## Analysis of clustering resolution effects.

Breakdown of clustering, showing resolution relationship to clusters
Local seurat object locations 

```
Seurat Objects:
1. Oligodendrocytes: /home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.oligodendrocytes/combined.olig.updated.SeuratObject.Rds
2. OPCs: /home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.opc.updated.SeuratObject.Rds
3. Astrocytes: /home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/combined.astro.updated.SeuratObject.Rds
4. Microglia: /home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.microglia/combined.microglia.updated.SeuratObject.Rds
```

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
set.seed(1234)
library(dplyr)
library(ggrepel)
library(clustree)

setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#UMAP Projection and clustering on selected cistopic model
clustering_loop<-function(in_path,output_name,res_count){
		dat<-readRDS(in_path)
		res_list<-colnames(dat@meta.data)[startsWith(prefix="peaks_snn_res",colnames(dat@meta.data))]
    plt1<-DimPlot(dat,group.by=res_list,combine=FALSE,label=T,raster=T)
    #system(paste0("slack -F ",paste(output_name,"clustering.dimplot.pdf",sep=".")," ryan_todo"))
    plt2<-clustree(dat, prefix = "peaks_snn_res.",node_size = 10)

    if(res_count==5){
		layout <- "
		AZZ
		BZZ
		CZZ
		DZZ
		EZZ"
    plt<-wrap_plots(
    	A=plt1[[1]]+theme_void()+ theme(legend.position = "none"),
    	B=plt1[[2]]+theme_void()+ theme(legend.position = "none"),
    	C=plt1[[3]]+theme_void()+ theme(legend.position = "none"),
    	D=plt1[[4]]+theme_void()+ theme(legend.position = "none"),
    	E=plt1[[5]]+theme_void()+ theme(legend.position = "none"),
    	Z=plt2+ theme(legend.position = "none"),
    	design=layout) + plot_layout(guides="collect")
    legend <- cowplot::get_legend(plt2)
		} else if (res_count==6) {
		layout <- "
		AZZ
		BZZ
		CZZ
		DZZ
		EZZ
		FZZ"
    plt<-wrap_plots(
    	A=plt1[[1]]+theme_void()+ theme(legend.position = "none"),
    	B=plt1[[2]]+theme_void()+ theme(legend.position = "none"),
    	C=plt1[[3]]+theme_void()+ theme(legend.position = "none"),
    	D=plt1[[4]]+theme_void()+ theme(legend.position = "none"),
    	E=plt1[[5]]+theme_void()+ theme(legend.position = "none"),
    	F=plt1[[6]]+theme_void()+ theme(legend.position = "none"),
    	Z=plt2+ theme(legend.position = "none"),
    	design=layout) + plot_layout(guides="collect")
    legend <- cowplot::get_legend(plt2)
		} else if (res_count==7) {
		layout <- "
		AZZ
		BZZ
		CZZ
		DZZ
		EZZ
		FZZ
		GZZ"
    plt<-wrap_plots(
    	A=plt1[[1]]+theme_void()+ theme(legend.position = "none"),
    	B=plt1[[2]]+theme_void()+ theme(legend.position = "none"),
    	C=plt1[[3]]+theme_void()+ theme(legend.position = "none"),
    	D=plt1[[4]]+theme_void()+ theme(legend.position = "none"),
    	E=plt1[[5]]+theme_void()+ theme(legend.position = "none"),
    	F=plt1[[6]]+theme_void()+ theme(legend.position = "none"),
    	G=plt1[[7]]+theme_void()+ theme(legend.position = "none"),
    	Z=plt2+ theme(legend.position = "none"),
    	design=layout) + plot_layout(guides="collect")
    legend <- cowplot::get_legend(plt2)
		} else if (res_count==9) {
		layout <- "
		AZZ
		BZZ
		CZZ
		DZZ
		EZZ
		FZZ
		GZZ
		HZZ
		IZZ"
    plt<-wrap_plots(
    	A=plt1[[1]]+theme_void()+ theme(legend.position = "none"),
    	B=plt1[[2]]+theme_void()+ theme(legend.position = "none"),
    	C=plt1[[3]]+theme_void()+ theme(legend.position = "none"),
    	D=plt1[[4]]+theme_void()+ theme(legend.position = "none"),
    	E=plt1[[5]]+theme_void()+ theme(legend.position = "none"),
    	F=plt1[[6]]+theme_void()+ theme(legend.position = "none"),
    	G=plt1[[7]]+theme_void()+ theme(legend.position = "none"),
    	H=plt1[[8]]+theme_void()+ theme(legend.position = "none"),
    	I=plt1[[9]]+theme_void()+ theme(legend.position = "none"),
    	Z=plt2+ theme(legend.position = "none"),
    	design=layout) + plot_layout(guides="collect")
    legend <- cowplot::get_legend(plt2)
		}

		ggsave(plt,file=paste(output_name,"clustree.pdf",sep="."),width=8,height=8,units="in")
    system(paste0("slack -F ",paste(output_name,"clustree.pdf",sep=".")," ryan_todo"))
    ggsave(legend,file=paste(output_name,"clustree.legend.pdf",sep="."))
    system(paste0("slack -F ",paste(output_name,"clustree.legend.pdf",sep=".")," ryan_todo"))
    #plt<-clustree_overlay(dat, prefix = "peaks_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2",red_dim="umap")
    #ggsave(plt,file=paste(output_name,"clustree.overlay.pdf",sep="."),width=15,height=15)
    #system(paste0("slack -F ",paste(output_name,"clustree.overlay.pdf",sep=".")," ryan_todo"))

}

#Cell types for clustree

#Cell type object locations
oligo<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.oligodendrocytes/combined.olig.updated.SeuratObject.Rds"
opc<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.opc.updated.SeuratObject.Rds"
astro<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/combined.astro.updated.SeuratObject.Rds"
micro<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.microglia/combined.microglia.updated.SeuratObject.Rds"

#Processing

clustering_loop(in_path=oligo,output_name="oligodendrocytes",res_count=6)
clustering_loop(in_path=opc,output_name="opc",res_count=5)
clustering_loop(in_path=astro,output_name="astrocytes",res_count=9)
clustering_loop(in_path=micro,output_name="microglia",res_count=7)


```



# Batch correction analysis using LISI

```R
#devtools::install_github("immunogenomics/lisi")
library(lisi)
library(Signac)
library(ggplot2)
library(dplyr)

setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
#hgap_cortex<-readRDS("cortex.SeuratObject.Rds")
#hgap_hippo<-readRDS("hippocampus.SeuratObject.Rds")

astro="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/combined.astro.updated.SeuratObject.Rds"
micro="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.microglia/combined.microglia.updated.SeuratObject.Rds"
oligo="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.oligodendrocytes/combined.olig.updated.SeuratObject.Rds"
opc="/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.opc.updated.SeuratObject.Rds"

run_lisi<-function(matrix_in,name_in,metadat_in,label,cells){
	X<-matrix_in
	meta_data<-metadat_in
	res <- compute_lisi(X, meta_data, label)
	colnames(res)<-"lisi"
	res$label<-label
	res$dim<-name_in
	res$cell=cells
	return(res)
}

loop_lisi<-function(dat,run_harmony=TRUE,cell){
	dat<-readRDS(dat)
	if(run_harmony){
		list_in=setNames(nm=c("harmony","cistopic","umap"),
				list(as.data.frame(dat@reductions$harmony@cell.embeddings),
				as.data.frame(dat@reductions$cistopic@cell.embeddings),
				as.data.frame(dat@reductions$umap@cell.embeddings)))
	} else {
		list_in=setNames(nm=c("cistopic","umap"),
				list(
				as.data.frame(dat@reductions$cistopic@cell.embeddings),
				as.data.frame(dat@reductions$umap@cell.embeddings)))
	}
	labels<-c("Experiment","Region","Individual","Sample","seurat_clusters_marker")
	out<-lapply(labels, function(y)
			lapply(names(list_in),function(x) 
				run_lisi(matrix_in=list_in[[x]],name_in=x,metadat_in=dat@meta.data,label=y,cells=cell)))
	return(out)
}

for(i in c(astro,micro,oligo,opc)){
	dat<-readRDS(i)
	print(paste(i,length(unique(dat$seurat_clusters_marker))))
}

astro_out<-loop_lisi(dat=astro,run_harmony=TRUE,cell="astro")
astro_out<-rbind(do.call(rbind,lapply(c(1,2,3),function(x) do.call(rbind,lapply(astro_out,"[[",x)))))

micro_out<-loop_lisi(dat=micro,run_harmony=TRUE,cell="micro")
micro_out<-rbind(do.call(rbind,lapply(c(1,2,3),function(x) do.call(rbind,lapply(micro_out,"[[",x)))))

oligo_out<-loop_lisi(dat=oligo,run_harmony=TRUE,cell="oligo")
oligo_out<-rbind(do.call(rbind,lapply(c(1,2,3),function(x) do.call(rbind,lapply(oligo_out,"[[",x)))))

opc_out<-loop_lisi(dat=opc,run_harmony=FALSE,cell="opc")
opc_out<-rbind(do.call(rbind,lapply(c(1,2),function(x) do.call(rbind,lapply(opc_out,"[[",x)))))

out<-rbind(astro_out,micro_out,oligo_out,opc_out)
saveRDS(out,file="lisi_cluster.Rds")

write.table(df_out,file="lisi_values_flat.tsv",sep="\t",col.names=T,row.names=T)
system("slack -F lisi_values_flat.tsv ryan_todo")

out$cellid<-row.names(out)

df_out<- out %>% group_by(dim,label,cell) %>% summarize(
	mean=mean(lisi),median=median(lisi),count=n(),sd=sd(lisi)
	) %>% as.data.frame()

write.table(df_out,file="lisi_summary_statistics.tsv",sep="\t",col.names=T,row.names=F)
system("slack -F lisi_summary_statistics.tsv ryan_todo")

out<-out[out$dim %in% c("cistopic","harmony"),]
plt<-ggplot(out,aes(y=lisi,x=cell,fill=dim))+geom_violin()+geom_boxplot(outlier.shape=NA)+facet_wrap(~label,scales="free")+theme_minimal()
ggsave(plt,file="lisi_cluster.Sample.pdf")
system("slack -F lisi_cluster.Sample.pdf ryan_todo")
```