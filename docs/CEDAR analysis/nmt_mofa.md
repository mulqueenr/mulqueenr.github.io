---
title: MOFA2 on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /mofa_nmt/
category: CEDAR
---

http://www.bioconductor.org/packages/release/bioc/html/MOFA2.html


## Installation
```R
BiocManager::install("MOFA2")
```

## Following a walkthrough with our data
http://www.bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/getting_started_R.html

Samples are stored in columns and features in rows.

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

#args <- commandArgs(trailingOnly=TRUE)
args<-list()
args[1]<-"CpG.bcEnhance.cistopic_object.Rds"
args[2]<-"GpC.promoter.cistopic_object.Rds"
args[3]<-"CpG.promoter.cistopic_object.Rds"
args[4]<-"GpC.bcEnhance.cistopic_object.Rds"
args<-unlist(args)

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



