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
| Allen Brain Map | Postmortem brain (MTG, ACC, V1C, M1C, S1C) |  snATAC (SMART-seq) | 49495 | https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq |
| Human Hippo Axis | Hippocampus | snRNA | 129,908 | https://cells.ucsc.edu/?ds=human-hippo-axis |
| Human Hippocampus Lifespan | Hippocampus | snRNA| 224,464 | https://cells.ucsc.edu/?bp=brain&org=Human+(H.+sapiens)&ds=hippo-lifespan |

# Plan for public data set integrations: 
* Cortex: 
** Allen Brain Map: focus on M1C (Primary Motor Cortex) (RNA) #Done
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

### Allen Brain Span M1C Cortex (RNA)
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

### Integration of Cortex
```R
library(Seurat)
library(Signac) 
library(patchwork)
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#Public RNA data
cortex_brainspan<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")

#Ours
hgap_cortex<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/cortex.SeuratObject.Rds")


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


generate_transfer_anchor<-function(in_dat,ref_dat,feat,prefix){
	#generate LSI matrix for downstream normalization
	in_dat <- RunTFIDF(in_dat)
	in_dat <- FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)

	#generate gene activity for just variable genes, add to object and normalize
	ga_mat<-GeneActivity(in_dat,features = feat, process_n = 10000)
	in_dat[['RNA']] <- CreateAssayObject(counts = ga_mat)
	in_dat<- NormalizeData(
	  object = in_dat,
	  assay = 'RNA',
	  normalization.method = 'LogNormalize',
	  scale.factor = median(in_dat$nCount_RNA)
	)

	DefaultAssay(in_dat)<-"RNA"
	transfer.anchors <- FindTransferAnchors(
	  reference = ref_dat,
	  query = in_dat,
	  reduction = 'cca'
	)
	saveRDS(transfer.anchors,paste0(prefix,".transferanchors.rds"))
	saveRDS(in_dat,file=paste0(prefix,".SeuratObject.Rds"))
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

#Process HGAP cortex and brainspan integration
	features=VariableFeatures(cortex_brainspan,selection.method="vst",assay="RNA")
	#FindTopFeatures instead?
	hgap_cortex<-generate_transfer_anchor(in_dat=hgap_cortex,ref_dat=cortex_brainspan,prefix="cortex",feat=features)
	transfer_anchors=readRDS("cortex.transferanchors.rds")
	hgap_cortex<-sample_label_transfer(in_dat=hgap_cortex,prefix="brainspan_class_",transfer.anchors.=transfer_anchors,transfer_label="class_label")
	hgap_cortex<-sample_label_transfer(in_dat=hgap_cortex,prefix="brainspan_subclass_",transfer.anchors.=transfer_anchors,transfer_label="subclass_label")
	saveRDS(hgap_cortex,file="cortex.SeuratObject.Rds")

	#Plot Label Transfer
	plt1<-DimPlot(hgap_cortex,group.by="Putative_Celltype")#,reduction="umap.atac")
	plt2<-DimPlot(hgap_cortex,group.by="brainspan_class_predicted.id")
	plt3<-DimPlot(hgap_cortex,group.by="brainspan_subclass_predicted.id")
	plt<-plt1/plt2/plt3
	ggsave(plt1,file="hgap_cortex.predicted.brainspan.umap.png")
	system("slack -F hgap_cortex.predicted.brainspan.umap.png ryan_todo")

	#Generate coembedding
	cortex_coembed<-coembed_data(in_dat=hgap_cortex,ref_dat=cortex_brainspan,transfer.anchors.=transfer_anchors,feat=features,prefix="cortex",assay_name="brainspan_imputation")
	plt<-DimPlot(coembed, group.by = c("orig.ident", "class_label","subclass_label"))
	ggsave(plt,file="hgap_cortex.brainspan.coembed.umap.pdf")
	system("slack -F hgap_cortex.brainspan.coembed.umap.pdf ryan_todo")


```


### Integration of Hippocampus
```R
library(Seurat)
library(Signac) 
library(patchwork)
library(ggplot2)
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#Public RNA data
hippo_axis<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.rds")
hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.rds")

#Ours
hgap_hippo<-readRDS("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/hippocampus.SeuratObject.Rds")


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


generate_transfer_anchor<-function(in_dat,ref_dat,feat,prefix){
	#generate LSI matrix for downstream normalization
	in_dat <- RunTFIDF(in_dat)
	in_dat <- FindTopFeatures(in_dat, min.cutoff = 'q0')
	in_dat <- RunSVD(in_dat)

	#generate gene activity for just variable genes, add to object and normalize
	ga_mat<-GeneActivity(in_dat,features = feat, process_n = 10000)
	in_dat[['RNA']] <- CreateAssayObject(counts = ga_mat)
	in_dat<- NormalizeData(
	  object = in_dat,
	  assay = 'RNA',
	  normalization.method = 'LogNormalize',
	  scale.factor = median(in_dat$nCount_RNA)
	)

	DefaultAssay(in_dat)<-"RNA"
	transfer.anchors <- FindTransferAnchors(
	  reference = ref_dat,
	  query = in_dat,
	  reduction = 'cca'
	)
	saveRDS(transfer.anchors,paste0(prefix,".transferanchors.rds"))
	saveRDS(in_dat,file=paste0(prefix,".SeuratObject.Rds"))
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

#Process HGAP hippocampus and Hippo Axis integration
	features=VariableFeatures(hippo_axis)
	hgap_hippo<-generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_axis,prefix="hippo",feat=features)
	transfer_anchors=readRDS("hippo.transferanchors.rds")
	#hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,prefix="brainspan_class_",transfer.anchors.=transfer_anchors,transfer_label="class_label")
	#hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,prefix="brainspan_subclass_",transfer.anchors.=transfer_anchors,transfer_label="subclass_label")
	saveRDS(hgap_hippo,file="hippo.SeuratObject.Rds")

	#Plot Label Transfer
	plt1<-DimPlot(hgap_hippo,group.by="Putative_Celltype")#,reduction="umap.atac")
	#plt2<-DimPlot(hgap_hippo,group.by="brainspan_class_predicted.id")
	#plt3<-DimPlot(hgap_hippo,group.by="brainspan_subclass_predicted.id")
	plt<-plt1/plt2/plt3
	ggsave(plt1,file="hgap_hippo.predicted.axis.umap.png")
	system("slack -F hgap_hippo.predicted.axis.umap.png ryan_todo")

	#Generate coembedding
	hippo_coembed<-coembed_data(in_dat=hgap_hippo,ref_dat=hippo_axis,transfer.anchors.=transfer_anchors,feat=features,prefix="hippo",assay_name="axis_imputation")
	plt<-DimPlot(coembed, group.by = c("orig.ident", "class_label","subclass_label"))
	ggsave(plt,file="hgap_hippo.axis.coembed.umap.pdf")
	system("slack -F hgap_hippo.axis.coembed.umap.pdf ryan_todo")

#Process HGAP hippocampus and Hippo Lifespan integration
	features=VariableFeatures(hippo_lifespan)
	hgap_hippo<-generate_transfer_anchor(in_dat=hgap_hippo,ref_dat=hippo_lifespan,prefix="hippo",feat=features)
	transfer_anchors=readRDS("hippo.transferanchors.rds")
	#hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,prefix="brainspan_class_",transfer.anchors.=transfer_anchors,transfer_label="class_label")
	#hgap_hippo<-sample_label_transfer(in_dat=hgap_hippo,prefix="brainspan_subclass_",transfer.anchors.=transfer_anchors,transfer_label="subclass_label")
	saveRDS(hgap_hippo,file="hippo.SeuratObject.Rds")

	#Plot Label Transfer
	plt1<-DimPlot(hgap_hippo,group.by="Putative_Celltype")#,reduction="umap.atac")
	#plt2<-DimPlot(hgap_hippo,group.by="brainspan_class_predicted.id")
	#plt3<-DimPlot(hgap_hippo,group.by="brainspan_subclass_predicted.id")
	plt<-plt1/plt2/plt3
	ggsave(plt1,file="hgap_hippo.predicted.lifespan.umap.png")
	system("slack -F hgap_hippo.predicted.lifespan.umap.png ryan_todo")

	#Generate coembedding
	hippo_coembed<-coembed_data(in_dat=hgap_hippo,ref_dat=hippo_lifespan,transfer.anchors.=transfer_anchors,feat=features,prefix="hippo",assay_name="lifespan_imputation")
	plt<-DimPlot(coembed, group.by = c("orig.ident", "class_label","subclass_label"))
	ggsave(plt,file="hgap_hippo.lifespan.coembed.umap.pdf")
	system("slack -F hgap_hippo.lifespan.coembed.umap.pdf ryan_todo")

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

OPC Subset. Integration using cisTopic Modelling.

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
  # cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  # row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  # sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  # print("made cistopic object")
  # sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,topic=c(25,28,30,33,35,40),nCores=6,addModels=FALSE)
  # saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
  
  # saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.Rds"))
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

setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
hgap_hippo<-readRDS("hippocampus.SeuratObject.Rds");hgap_hippo$orig.ident<-"HGAP_hippo"
hgap_cortex<-readRDS("cortex.SeuratObject.Rds");hgap_cortex$orig.ident<-"HGAP_cortex"

hgap_hippo<-subset(hgap_hippo,Putative_Celltype=="OPCs")
hgap_cortex<-subset(hgap_cortex,Putative_Celltype=="Putative_OPCs")

#seurat clusters 0: OPC_BCL11B; 1:OPC_MAG
opc_cellstate<-read.table("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.OPCs.cistopic.seurat_clusters.annot",header=F)
opc_cellstate<-setNames(opc_cellstate$V2,nm=opc_cellstate$V1)
opc_cellstate[opc_cellstate==0]<-"OPC_BCL11B"
opc_cellstate[opc_cellstate==1]<-"OPC_MAG"

hgap_hippo<-AddMetaData(hgap_hippo,opc_cellstate,col.name="opc_subtype")
hgap_cortex<-AddMetaData(hgap_cortex,opc_cellstate,col.name="opc_subtype")

#Read in Corces et al Data
setwd("/home/groups/CEDAR/mulqueen/human_brain_ref/corces_gwas/")
mfg_atac<-readRDS("mfg.SeuratObject.Rds");mfg_atac$orig.ident<-"Corces_MFG"
hippo1_atac<-readRDS("hippo1.SeuratObject.Rds");hippo1_atac$orig.ident<-"Corces_Hippo1"
hippo2_atac<-readRDS("hippo2.SeuratObject.Rds");hippo2_atac$orig.ident<-"Corces_Hippo2"
mfg_atac<-subset(mfg_atac,Cluster_Description=="OPCs")
hippo1_atac<-subset(hippo1_atac,Cluster_Description=="OPCs")
hippo2_atac<-subset(hippo2_atac,Cluster_Description=="OPCs")

opc<-merge(hgap_hippo,y=c(hgap_cortex,mfg_atac,hippo1_atac,hippo2_atac),add.cell.ids=c("HGAP_hip","HGAP_cortex","Corces_mfg","Corces_hip1","Corces_hip2"))

opc$Region<-"Cortex"
opc@meta.data[grepl("hip",row.names(opc@meta.data)),]$Region<-"Hippocampus"
opc@meta.data[grepl("HGAP",row.names(opc@meta.data)),]$orig.ident<-"HGAP"
opc$assay<-"HGAP"
opc@meta.data[grepl("Corces",row.names(opc@meta.data)),]$assay<-"Corces"
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
opc<-normalize_and_umap(opc) #standard TFIDF and LSI normalization pipeline
opc<-cistopic_generation(opc) #cistopic dimensionality reduction
#Run projection before integration
opc<-readRDS("opc.integrated.SeuratObject.rds")

opc<-RunHarmony(opc,group.by.vars=c("assay"),reduction.save="harmony_atac",assay.use="peaks",reduction="cistopic",dims.use=1:dim(opc@reductions$cistopic)[2],project.dim=F,max.iter.harmony=20)
opc<-RunUMAP(opc,reduction.name="harmonyumap_cistopic",reduction = "harmony_atac",dims=1:dim(opc@reductions$harmony_atac)[2],n.components=3) 
opc$assigned_celltype<-opc$Putative_Celltype
opc@meta.data[opc$orig.ident!="SeuratProject",]$assigned_celltype<-opc@meta.data[opc$orig.ident!="SeuratProject",]$Cluster_Description

#Plot
plt1<-DimPlot(opc,dims=c(2,3),group.by="opc_subtype",reduction="harmonyumap_cistopic",raster=F)
plt2<-DimPlot(opc,dims=c(2,3),group.by="orig.ident",reduction="harmonyumap_cistopic",raster=F)
plt3<-DimPlot(opc,dims=c(2,3),group.by="Putative_Celltype",reduction="harmonyumap_cistopic",raster=F)
plt4<-DimPlot(opc,dims=c(2,3),group.by="Cluster",reduction="harmonyumap_cistopic",raster=F)
plt5<-DimPlot(opc,dims=c(2,3),group.by="Region",reduction="harmonyumap_cistopic",split.by="orig.ident",raster=FALSE)
ggsave((plt1|plt2)/(plt3|plt4)/(plt5),file="opc.atac.integrated.umap.png",width=15,height=20,limitsize=F)
system(paste0("slack -F ","opc.atac.integrated.umap.png"," ryan_todo"))

#Plot PDFs
plot_name="HGAP_defined_subtype"
ggsave(plt1,file=paste0("opc.atac.integrated.",plot_name,".umap.pdf"))
system(paste0("slack -F ",paste0("opc.atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

plot_name="dataset"
ggsave(plt2,file=paste0("opc.atac.integrated.",plot_name,".umap.pdf"))
system(paste0("slack -F ",paste0("opc.atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

plot_name="HGAP_defined_celltype"
ggsave(plt3,file=paste0("opc.atac.integrated.",plot_name,".umap.pdf"))
system(paste0("slack -F ",paste0("opc.atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

plot_name="Corces_cluster"
ggsave(plt4,file=paste0("opc.atac.integrated.",plot_name,".umap.pdf"))
system(paste0("slack -F ",paste0("opc.atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))

plot_name="dataset_by_region"
ggsave(plt5,file=paste0("opc.atac.integrated.",plot_name,".umap.pdf"),width=20)
system(paste0("slack -F ",paste0("opc.atac.integrated.",plot_name,".umap.pdf") ," ryan_todo"))


saveRDS(opc,"opc.integrated.SeuratObject.rds")

opc<-subset(opc,assay=="HGAP")
names_out<-colnames(opc)
names_out<-names_out[grepl(names_out,pattern="HGAP")]
names_out_cellid<-unlist(lapply(strsplit(names_out,"_"),"[[",3))
name_conversion_table<-as.data.frame(cbind(names_out,names_out_cellid))
write.table(name_conversion_table,file="opc.cellIDs.tsv",col.names=F,quote=F,row.names=F,sep="\t")

```
## Perform Gene Activity and OPC integration on RNA sets

First subseting the atac fragments files to just these cell IDs. This is going to greatly increase the speed for gene activity calculation.


```python
import gzip
import pandas as pd

cellid_in=pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.cellIDs.tsv",header=None,sep=" ") #read in cellIDs

read_lines=[0,1e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9]

for i in range(0,len(read_lines)-1):
	df = pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/coverage.plots.BAM_DataFreeze_Cortex_and_Hippocampus/BAM_DataFreeze_Cortex_and_Hippocampus.bbrd.q10.filt.bam.fragments.tsv.gz",sep="\t",header=None,nrows=1e9,skiprows=i) #read in fragments file
	df_filt=df[df[3].isin(cellid_in[1])] #filter to cell IDs in OPC
	print(df_filt)
	df_filt.to_csv('/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.q10.filt.bam.fragments.'+str(read_lines[i])+'tsv.gz',compression="gzip",sep="\t") #write out compressed opc fragment file

```

```bash
cat opc.q10.*tsv.gz > opc.q10.filt.bam.fragments.tsv.gz &

bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
zcat opc.q10.filt.bam.fragments.tsv.gz | awk 'OFS="\t" {print $2,$3,$4,$5,$6}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > opc.q10.filt.bam.fragments.sorted.tsv.gz; wait ;
zcat opc.q10.filt.bam.fragments.sorted.tsv.gz | tail -n +9 | $bgzip > opc.q10.filt.bam.fragments.sorted2.tsv.gz; wait ; #this is to remove header from concatenation
```

```python
import gzip
import pandas as pd
import numpy as np

cellid_in=pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.cellIDs.tsv",header=None,sep=" ") #read in cellIDs

df = pd.read_csv("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.q10.filt.bam.fragments.sorted2.tsv.gz",sep="\t",header=None) #read in fragments file

df[3]=[cellid_in.loc[cellid_in[1]==x,0].values[0] for x in df[3][1:1000]] #replace cellIDs with those consistent with seurat object

df.to_csv('/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.q10.filt.bam.fragments.sorted2.tsv.gz',compression="gzip",sep="\t") #write out compressed opc fragment file

```

```bash
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
zcat opc.q10.filt.bam.fragments.sorted2.tsv.gz | tail -n +9 | $bgzip > opc.q10.filt.bam.fragments.sorted.tsv.gz; wait ; 
$tabix -p bed opc.q10.filt.bam.fragments.sorted2.tsv.gz&
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
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
opc<-readRDS("opc.integrated.SeuratObject.rds")
opc<-subset(opc,assay=="HGAP")
DefaultAssay(opc)<-"peaks"
opc<-SeuratObject::RenameCells(opc, new.names = unlist(lapply(strsplit(colnames(opc),"_"),"[",3)),old.names=colnames(opc)) #rename back to just cellID


fragments <- CreateFragmentObject(
   path = "/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses/opc.q10.filt.bam.fragments.sorted2.tsv.gz",
   cells = colnames(opc),
   validate.fragments = FALSE
 )

Fragments(opc[["peaks"]])<-NULL
Fragments(opc[["peaks"]]) <- fragments

#generate gene activity for just variable genes, add to object and normalize
ga_mat<-GeneActivity(opc)
saveRDS(ga_mat,file="opc.geneactivity.rds")
ga_mat<-readRDS(file="opc.geneactivity.rds")

opc[['RNA']] <- CreateAssayObject(counts = ga_mat)
opc<- NormalizeData(
  object = opc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(opc$nCount_RNA)
)

saveRDS(opc,"opc.ga.SeuratObject.rds")

```

Integrate Gene activity with Public RNA data
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
setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")
opc<-readRDS("opc.ga.SeuratObject.rds")
DefaultAssay(opc)<-"RNA"
opc <- ScaleData(opc, features = rownames(opc))
opc<-FindTopFeatures(object = opc)
opc <- RunPCA(opc, npcs = 30,features=VariableFeatures(opc))
opc<-RunUMAP(opc,reduction.name="geneactivity_umap",reduction = "pca",dims=1:dim(opc@reductions$pca)[2],n.components=3) 

#Public RNA data
#hippo_axis<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.rds")
#hippo_axis$orig.ident<-"hippo_axis"
#hippo_axis<-subset(hippo_axis,Cluster %in% c("OPC1","OPC2","OPC3","OPC4"))
#saveRDS(hippo_axis,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.opc.rds")

#hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.rds")
#hippo_lifespan$orig.ident<-"hippo_lifespan"
#hippo_lifespan<-subset(hippo_lifespan,MajorCellTypes=="OPC")
#saveRDS(hippo_lifespan,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.opc.rds")
hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.opc.rds")

#cortex_brainspan<-readRDS("/home/groups/oroaklab/adey_lab/projects/sciDROP/public_data/allen_brainspan_humancortex/allen_brainspan_humancortex.rds")
#cortex_brainspan$orig.ident<-"cortex_brainspan"
#cortex_brainspan<-subset(cortex_brainspan,cluster_label=="OPC L1-6 PDGFRA COL20A1")
#saveRDS(cortex_brainspan,file="/home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex/allen_brainspan_humancortex.opc.rds")
cortex_brainspan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex/allen_brainspan_humancortex.opc.rds")

integrate_atac_rna<-function(ours=opc,theirs=hippo_axis){
	DefaultAssay(ours)<-"RNA"

	transfer.anchors <- FindTransferAnchors(
	  reference = ours,
	  query = theirs,
	  features=VariableFeatures(object = ours),
	  reduction = 'cca'
	)

	predicted.labels <- TransferData(
	  anchorset = transfer.anchors,
	  refdata = ours$opc_subtype,
	  weight.reduction = theirs[['pca']],
	  dims = 2:30
	)

	theirs <- AddMetaData(object = theirs, metadata = predicted.labels)
return(theirs)}

plot_output<-function(theirs=hippo_axis,outname="hippo_axis",orig.cluster="Cluster"){
	plot1 <- DimPlot(
	object = theirs,
	group.by = orig.cluster,
	label = TRUE,
	repel = TRUE) + NoLegend() + ggtitle('Their Clusters')

	plot2 <- DimPlot(
	object = theirs,
	group.by = 'predicted.id',
	label = TRUE,
	repel = TRUE) + NoLegend() + ggtitle('Our Clusters')

	plot3 <- FeaturePlot(
	object = theirs,
	features= "BCL11B") + NoLegend() + ggtitle('BCL11B')

	plot4 <- FeaturePlot(
	object = theirs,
	features= "MAG") + NoLegend() + ggtitle('MAG')

	plot5<-VlnPlot(object=theirs,
		features=c("BCL11B","MAG"),
		group.by="predicted.id")

	plt<-plot1 + plot2 + plot3 + plot4 + plot5
	ggsave(plt,filename=paste0(outname,".integrated.GA.OPCsubtype.pdf"),width=15)
	system(paste0("slack -F ",paste0(outname,".integrated.GA.OPCsubtype.pdf")," ryan_todo"))
}

hippo_axis<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.opc.rds")
hippo_axis<-integrate_atac_rna(ours=opc,theirs=hippo_axis,outname="hippo_axis",orig.cluster="Cluster")
saveRDS(hippo_axis,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_axis/hippo_axis.opc.rds")
plot_output(theirs=hippo_axis,outname="hippo_axis",orig.cluster="Cluster")

hippo_lifespan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.opc.rds")
hippo_lifespan<-integrate_atac_rna(ours=opc,theirs=hippo_lifespan,outname="hippo_lifespan",orig.cluster="seurat_clusters")
saveRDS(hippo_lifespan,file="/home/groups/CEDAR/mulqueen/human_brain_ref/human_hippo_lifespan/hippo_lifespan.opc.rds")
plot_output(theirs=hippo_lifespan,outname="hippo_lifespan",orig.cluster="seurat_clusters")

cortex_brainspan<-readRDS("/home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex/allen_brainspan_humancortex.opc.rds")
cortex_brainspan<-integrate_atac_rna(ours=opc,theirs=cortex_brainspan,outname="cortex_brainspan",orig.cluster="cell_type_alias_label")
saveRDS(cortex_brainspan,file="/home/groups/CEDAR/mulqueen/human_brain_ref/allen_brainspan_humancortex/allen_brainspan_humancortex.opc.rds")
plot_output(theirs=cortex_brainspan,outname="cortex_brainspan",orig.cluster="cell_type_alias_label")



# #Use function MapQuery?
# # find transfer anchors
# transfer.anchors <- FindTransferAnchors(
#   reference = pbmc.multi,
#   query = pbmc.atac,
#   reference.reduction = "lsi",
#   reduction = "lsiproject",
#   dims = 2:30
# )

# # map query onto the reference dataset
# pbmc.atac <- MapQuery(
#   anchorset = transfer.anchors,
#   reference = pbmc.multi,
#   query = pbmc.atac,
#   refdata = pbmc.multi$predicted.id,
#   reference.reduction = "lsi",
#   new.reduction.name = "ref.lsi",
#   reduction.model = 'umap'
# )
DefaultAssay(opc)<-"RNA"
plt<-CoveragePlot(opc,group.by="Cluster_Description",features=c("BCL11B","MAG"))
ggsave(plt,file="opc_marker_cov.png",width=15,height=20,limitsize=F)
system(paste0("slack -F ","opc_marker_cov.png"," ryan_todo"))

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
set.seed(1234)
library(dplyr)
library(ggrepel)
library(clustree)

setwd("/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/rm_integration_reviewerresponses")

#UMAP Projection and clustering on selected cistopic model
clustering_loop<-function(in_path,output_name){
		dat<-readRDS(in_path)
		res_list<-colnames(dat@meta.data)[startsWith(prefix="peaks_snn_res",colnames(dat@meta.data))]
    plt<-DimPlot(dat,group.by=res_list)
    ggsave(plt,file=paste(output_name,"clustering.dimplot.pdf",sep="."),width=15,height=15)
    system(paste0("slack -F ",paste(output_name,"clustering.dimplot.pdf",sep=".")," ryan_todo"))

    plt<-clustree(dat, prefix = "peaks_snn_res.")
    ggsave(plt,file=paste(output_name,"clustree.pdf",sep="."),width=15,height=15)
    system(paste0("slack -F ",paste(output_name,"clustree.pdf",sep=".")," ryan_todo"))

    plt<-clustree_overlay(dat, prefix = "peaks_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2",red_dim="umap")
    ggsave(plt,file=paste(output_name,"clustree.overlay.pdf",sep="."),width=15,height=15)
    system(paste0("slack -F ",paste(output_name,"clustree.overlay.pdf",sep=".")," ryan_todo"))

}

#Cell types for clustree

#Cell type object locations
oligo<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.oligodendrocytes/combined.olig.updated.SeuratObject.Rds"
opc<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.OPCs/combined.opc.updated.SeuratObject.Rds"
astro<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.astroependymal/combined.astro.updated.SeuratObject.Rds"
micro<-"/home/groups/oroaklab/adey_lab/projects/maga/00_DataFreeze_Cortex_and_Hippocampus/glialatlas.celltype.analyses/combined.microglia/combined.microglia.updated.SeuratObject.Rds"

#Processing

clustering_loop(in_path=oligo,output_name="oligodendrocytes")
clustering_loop(in_path=opc,output_name="opc")
clustering_loop(in_path=astro,output_name="astrocytes")
clustering_loop(in_path=micro,output_name="microglia")


```