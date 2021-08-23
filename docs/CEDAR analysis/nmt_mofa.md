---
title: MOFA2 on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /mofa_nmt/
category: CEDAR
---

http://www.bioconductor.org/packages/release/bioc/html/MOFA2.html

# RNA Processing

## Check RNA data from Aaron Doe Processing

```bash
mkdir /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna
multiqc -o /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples
```

I think his filtering is way to aggressive given the quality of the data. So I'm going to realign and generate my own counts matrix.

## Raw data is located here: 
```bash
/home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples/raw
```

### Make a symbolic link to RNA file directories

```bash
ln -s /home/groups/CEDAR/doe/projects/my_NMT/MCF7_T47D/scRNA_SMARTseq2/samples/raw /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw
```

### Align raw files with batch scripting
Using the genome reference and version of STAR packages with cellranger. Trimming first 15 bp and last 3 bp as per fastqc output.
Trimming first 15 bases and last 3 bases, only aligning read 1
Perform counting (I think this accounts for duplicate reads, but I'm not entirely sure)
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-672
#SBATCH --tasks-per-node=30 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=1 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=2gb ## request gigabyte per cpu
#SBATCH --time=1:00:00 ## ask for 1 hour on the node
#SBATCH --

array_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/*_R1.fq.gz | wc -l`

fq_in=`ls /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/*_R1.fq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
name="$(basename -- $fq_in)";
name=${name%.fq.gz} 
echo $name;
/home/groups/CEDAR/mulqueen/src/cellranger/cellranger-6.0.1/lib/bin/STAR --runMode alignReads \
--genomeDir /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star \
--readFilesIn /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna_raw/${name}.fq.gz \
--outFileNamePrefix /home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam/${name}. \
--readFilesCommand zcat \
--clip3pNbases 3 \
--clip5pNbases 15 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts; 
```

## Generate a Seurat Object for our RNA

Read in ReadsPerGene.out.tab files and format into a data frame. Then make a seurat object for filtering and processing cells.

Following this: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

```R
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(EnsDb.Hsapiens.v86)
library(dplyr)

setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam")
ff <- list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 ) #first 4 lines are mapping QC, might be good to add as a metadata table
counts_in <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "_R1[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/", "", ff )
colnames(counts_in) <- ff
row.names(counts_in) <- counts.files[[1]]$V1


##Decided to just keep them all in ENSG notation for now##
#convert gene names to symbols in counts data
#genenames<-AnnotationDbi::select(org.Hs.eg.db, keys = row.names(counts), keytype = 'ENSEMBL', columns = 'SYMBOL',multiVals=first)
#genenames<-genenames[!is.na(genenames[,2]),] #remove rows that dont match ENSG and gene symbol
#genenames<-genenames[!duplicated(genenames[1]),] #remove duplicate names
#genenames<-genenames[!duplicated(genenames[2]),] #remove duplicate names
#row.names(counts[genenames[,1],])<-genenames[,2]

#make seurat object
dat <- CreateSeuratObject(counts = counts_in, project = "cellline", min.cells = 3, min.features = 1000)

#add mito reads to feature data
mito.features= genes(EnsDb.Hsapiens.v86, filter = ~ seq_name == "MT")
#dat[["percent.mt"]] <- PercentageFeatureSet(dat, assay="RNA",features=mito.features$gene_id)

#qc plot
pdf("rna_qc.pdf")
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
dev.off()
system("slack -F rna_qc.pdf ryan_todo")

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("rna_topvar.pdf")
plot1 / plot2
dev.off()
system("slack -F rna_topvar.pdf ryan_todo")

#dim reduction
all.genes <- rownames(dat)
dat <- ScaleData(dat, features = all.genes)
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

pdf("rna_pcaloadings.pdf")
VizDimLoadings(dat, dims = 1:2, reduction = "pca")
dev.off()
system("slack -F rna_pcaloadings.pdf ryan_todo")

pdf("rna_pca.pdf")
DimPlot(dat, reduction = "pca")
dev.off()
system("slack -F rna_pca.pdf ryan_todo")

pdf("rna_pcaheatmap.pdf")
DimHeatmap(dat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
system("slack -F rna_pcaheatmap.pdf ryan_todo")

#determine dimensionality
dat <- JackStraw(dat, num.replicate = 100)
dat <- ScoreJackStraw(dat, dims = 1:20)

pdf("rna_elbowplot.pdf")
JackStrawPlot(dat, dims = 1:15)
ElbowPlot(dat)
dev.off()
system("slack -F rna_elbowplot.pdf ryan_todo")

#Cluster and UMAP
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.5)
dat <- RunUMAP(dat, dims = 1:10)

pdf("rna_umap.pdf")
DimPlot(dat, reduction = "umap",group.by=c("cluster","orig.ident"))
dev.off()
system("slack -F rna_umap.pdf ryan_todo")

dat.markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dat.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)

saveRDS("SeuratObject.rds")
#Clustering looks pretty good to me just based on tutorial defaults. Using this seurat object for MOFA analysis

```
## Following a walkthrough with our data
http://www.bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/getting_started_R.html

Samples are stored in columns and features in rows.

## Installation
```R
BiocManager::install("MOFA2")
```

## Create a MOFA object

Copying input of cistopic object data from [here.](https://mulqueenr.github.io/cistopic_nmt/#combine-multiple-binarized-matrices-prior-to-running-cistopic)

Following tutorial here: https://raw.githack.com/bioFAM/MEFISTO_tutorials/master/scnmt_mefisto_vignette.html

```R
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")
library(MOFA2)
library(data.table)
library(ggplot2)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(cisTopic)
library(Seurat)

#args <- commandArgs(trailingOnly=TRUE)
args<-list()
args[1]<-"CpG.bcEnhance.cistopic_object.Rds"
args[2]<-"GpC.promoter.cistopic_object.Rds"
args[3]<-"CpG.promoter.cistopic_object.Rds"
args[4]<-"GpC.bcEnhance.cistopic_object.Rds"
args<-unlist(args)

#preprocessed RNA data from AD
rna<-readRDS("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/rna/rna_bam/SeuratObject.rds")

#add sample cell line names as metadata
rna@meta.data$cellLine<-unlist(lapply(strsplit(row.names(rna@meta.data),"_"),"[",1))

#set up treatment conditions
rna@meta.data$treatment<-substr(unlist(lapply(strsplit(row.names(rna@meta.data),"_"),"[",1)),3,3)
rna@meta.data$batch<-paste(rna@meta.data$cellLine,rna@meta.data$treatment,sep="_")
rna@meta.data[startsWith(rna@meta.data$cellLine,"M7"),]$cellLine<-"MCF7"
rna@meta.data[startsWith(rna@meta.data$cellLine,"TD"),]$cellLine<-"T47D"



pdf("rna_umap.pdf")
DimPlot(rna,group.by=c("orig.ident","Phase","seurat_clusters"),ncol=2)
dev.off()
system("slack -F rna_umap.pdf ryan_todo")

dat<-lapply(args,readRDS) #read in data as list of cistopic objects
dat<-lapply(dat,function(x){
        annotateRegions(x, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')}
        )

c_type<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][1])) #extract c type
regions<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][2])) #extract region info

#rename cell names for consistency
dat<-lapply(dat,function(x){
        colnames(x@count.matrix)<-unlist(lapply(strsplit(colnames(x@count.matrix),"[.]"),"[",1))
        row.names(x@cell.data)<-unlist(lapply(strsplit(row.names(x@cell.data),"[.]"),"[",1))
        return(x)
        })
########TO DO############
#set row names to be the same
row.names(x@cell.data)<-unlist(lapply(strsplit(row.names(x@cell.data),"[.]"),"[",1))

#determine shared cells across cistopic objects
cells_to_keep<-Reduce(intersect,lapply(dat,function(x){row.names(x@cell.data)}))

#rename region names for interpretability
dat<-lapply(1:length(dat),function(x){
        row.names(dat[[x]]@count.matrix)<-paste(c_type[x],regions[x],dat[[x]]@region.data$SYMBOL,row.names(dat[[x]]@count.matrix),sep="_")
        return(dat[[x]])
        })

X<-lapply(dat,function(x){x@count.matrix[,cells_to_keep]}) # List of data matrix
names(X)<-unlist(lapply(1:length(c_type), function(i) paste(c_type[i],regions[i],sep="_"))) #set up named list

labels<-as.data.frame(dat[[1]]@cell.data[which(cells_to_keep %in% row.names(dat[[1]]@cell.data)),]$cellLine_treatment)  # the collected cell line and treatment of cells, which is used for validation
row.names(labels)<-cells_to_keep

#remove rows with 0 variance
X<-lapply(X,function(y) {y[apply(y, 1, var)>0,]})
MOFAobject <- create_mofa(X)

pdf("test.out.pdf")
plot_data_overview(MOFAobject)
dev.off()

system("slack -F test.out.pdf ryan_todo")

#Set up options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views=TRUE 
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter<-100
#perform scaling and high variable feature selection?
#also should return to methylation rate and not binarized methylation state (gaussian and not bernoulli distributions)
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("/home/groups/CEDAR/mulqueen/temp/model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)

model <- load_model("/home/groups/CEDAR/mulqueen/temp/model.hdf5")

#add metadata
Nsamples = sum(model@dimensions$N)

sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  treatment = dat[[1]]@cell.data[samples_names(model)[[1]],]$treatment,
  cellline = dat[[1]]@cell.data[samples_names(model)[[1]],]$cellLine
)

sample_metadata$cellLine_treatment=paste(sample_metadata$treatment,sample_metadata$cellline,sep="_")

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)

#variance explained per modality
head(model@cache$variance_explained$r2_total[[1]]) # group 1

pdf("test.pdf")
plot_variance_explained(model, x="view", y="factor")
plot_factor(model, 
  factors = c(1:10),
  color_by = "cellLine_treatment",
  group_by="cellLine_treatment",
  dot_size = 1,        # change dot size
  dodge = T,           # dodge points with different colors
  legend = T,          # remove legend
  add_violin = T,      # add violin plots,
  violin_alpha = 0.25  # transparency of violin plots
)
plot_factors(model, 
  factors = c(1:10),
  color_by = "cellLine_treatment"
)
dev.off()

system("slack -F test.pdf ryan_todo")

pdf("test.pdf")
plot_weights(model,
        view = "CpG_promoter",
  factor = 1,
  nfeatures = 10,     # Number of features to highlight
  scale = T,          # Scale weights from -1 to 1
  abs = F             # Take the absolute value?
)

plot_top_weights(model,
  view = "CpG_promoter",
  factor = 1,
  nfeatures = 10
)
dev.off()

system("slack -F test.pdf ryan_todo")


pdf("test.pdf")
plot_data_heatmap(model,
  view = "CpG_promoter",         # view of interest
  factor = 1,             # factor of interest
  features = 20,          # number of features to plot (they are selected by weight)
  
  # extra arguments that are passed to the `pheatmap` function
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE
)
dev.off()
system("slack -F test.pdf ryan_todo")

model <- run_umap(model)

pdf("test.pdf")
plot_dimred(model,
  method = "UMAP",  # method can be either "TSNE" or "UMAP"
  color_by = "cellLine_treatment"
)
dev.off()
system("slack -F test.pdf ryan_todo")

```



