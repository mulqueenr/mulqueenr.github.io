---
title: scAI on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /scai_nmt/
category: CEDAR
---

Going to try to use the recently published scAI for multiomic integration (https://github.com/sqjin/scAI). 
Ended up giving up on it because it only takes 2 inputs at a time.

## Installation
```R
devtools::install_github("sqjin/scAI")
```

## Following a walkthrough with our data
https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_mESC_dataset.html


## Load data
```R
library(scAI)
library(dplyr)
library(cowplot)
library(ggplot2)
```

## Look at the provided data for formatting
```R
load("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/data_mESC.rda") #example format of input data

X <- data_mESC$data # List of data matrix X$RNA  X$DNA both sparse format
labels <- data_mESC$labels # the collected time of cells, which is used for validation, rowname is cellID, labels column with labels
```

## Create an scAI Object

Copying input of cistopic object data from [here.](https://mulqueenr.github.io/cistopic_nmt/#combine-multiple-binarized-matrices-prior-to-running-cistopic)

```R
setwd("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/methylation_regions")
library(scAI)
library(dplyr)
library(cowplot)
library(ggplot2)
library(doParallel)

#args <- commandArgs(trailingOnly=TRUE)
args<-list()
args[1]<-"CpG.promoter.cistopic_object.Rds"
args[2]<-"CpG.gene.cistopic_object.Rds"
args[3]<-"GpC.gene.cistopic_object.Rds"
args<-unlist(args)

dat<-lapply(args,readRDS) #read in data as list of cistopic objects

c_type<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][1])) #extract c type
regions<-unlist(lapply(args,function(x) strsplit(x,"[.]")[[1]][2])) #extract region info

#rename cell names for consistency
dat<-lapply(dat,function(x){
        colnames(x@binary.count.matrix)<-unlist(lapply(strsplit(colnames(x@binary.count.matrix),"[.]"),"[",1))
        row.names(x@cell.data)<-unlist(lapply(strsplit(row.names(x@cell.data),"[.]"),"[",1))
        return(x)
        })

#determine shared cells across cistopic objects
cells_to_keep<-Reduce(intersect,lapply(dat,function(x){row.names(x@cell.data)}))
X<-lapply(dat,function(x){x@binary.count.matrix[,cells_to_keep]}) # List of data matrix
names(X)<-unlist(lapply(1:length(c_type), function(i) paste(c_type[i],regions[i],sep="_"))) #set up named list

labels<-as.data.frame(dat[[1]]@cell.data[which(cells_to_keep %in% row.names(dat[[1]]@cell.data)),]$cellLine_treatment)  # the collected cell line and treatment of cells, which is used for validation
row.names(labels)<-cells_to_keep

#Generate scAI object and preprocess
scAI_outs <- create_scAIobject(raw.data = X)
scAI_outs <- preprocessing(scAI_outs, assay = NULL, minFeatures = 200, minCells = 1,
                           libararyflag = F, logNormalize = F)
scAI_outs <- addpData(scAI_outs, pdata = labels, pdata.name = "Conditions")

#run the model
scAI_outs <- run_scAI(scAI_outs, K = 4, nrun = 5, do.fast=T)

#identify cell clusters
scAI_outs <- identifyClusters(scAI_outs, resolution = 1)
levels(scAI_outs@identity)

#visualize the cells
scAI_outs <- reducedDims(scAI_outs, method = "umap",do.scale = F)
gg1 <- cellVisualization(scAI_outs, scAI_outs@embed$umap, color.by = "Conditions",show.legend = T, title = "Conditions")
gg2 <- cellVisualization(scAI_outs, scAI_outs@embed$umap, color.by = "cluster", ylabel = NULL, title = "scAI clusters")
cowplot::plot_grid(gg1, gg2)

cell_coords.RNA <- reducedDims(scAI_outs, data.use = scAI_outs@norm.data$RNA, do.scale = T, method = "pca", return.object = F)

cell_coords.DNA <- reducedDims(scAI_outs, data.use = scAI_outs@norm.data$DNA, do.scale = T, method = "pca", return.object = F)

cell_coords.DNAagg <- reducedDims(scAI_outs, data.use = scAI_outs@agg.data, do.scale = T, method = "pca", return.object = F)


gg1 <- cellVisualization(scAI_outs, cell_coords.RNA, color.by = "cluster",  show.legend = F, title = "scRNA-seq",xlabel = "PCA1", ylabel = "PCA2")
gg2 <- cellVisualization(scAI_outs, cell_coords.DNA, color.by = "cluster",show.legend = F, xlabel = "PCA1",ylabel = NULL,title = "scDNA-seq")
gg3 <- cellVisualization(scAI_outs, cell_coords.DNAagg, color.by = "cluster", xlabel = "PCA1",ylabel = NULL, title = "Aggregated scDNA-seq")
cowplot::plot_grid(gg1, gg2, gg3, ncol = 3)

featureScoreVisualization(scAI_outs, feature.scores = t(scAI_outs@fit$H), feature.use = c('factor1','factor2','factor3'),  method = "umap", nCol = 3, cell.size = 0.1, show.legend = T, show.legend.combined = F)

#Rank features
# show the markers of GR activation
feature_genes = c('Zfp42','Esrrb','Morc1','Fbxo15','Jam2','Klf4','Tcl1','Tbx3',
                  'Tex19.1','Krt8','Cald1','Anxa5','Tagln','Ahnak','Dsp','Anxa3','Krt19','Fgf5');

featureRankingPlot(scAI_outs, assay = 'RNA', feature.show = feature_genes, top.p = 0.5, ylabel = "Gene score")

#Identify markers per clusters
markers.RNA.cluster <- identifyClusterMarkers(scAI_outs, assay = "RNA")
markers.DNA.cluster <- identifyClusterMarkers(scAI_outs, assay = 'DNA')

n.top = 10
# RNA top 10
markers.RNA.clusterTop <- markers.RNA.cluster %>% group_by(clusters) %>% top_n(n.top, logFC) %>% slice(1:n.top)
featureHeatmap(scAI_outs, assay = "RNA", feature.use = markers.RNA.clusterTop$features, group.by = "cluster")
# methylation top 10
markers.DNA.clusterTop <- markers.DNA.cluster %>% group_by(clusters) %>% top_n(n.top, logFC) %>% slice(1:n.top)
featureHeatmap(scAI_outs, assay = "DNA", feature.use = markers.DNA.clusterTop$features, group.by = "cluster")

## Embed cells, genes loci and factors into 2D space

scAI_outs <- getEmbeddings(scAI_outs)
VscAIplot(scAI_outs, gene.use = feature_genes, loci.use = NULL, loci.use.names = NULL, color.by = "cluster")

#use a feature plot to pick out specifics in cell embeddings
featureVisualization(scAI_outs, assay = "RNA", feature.use = c('Tcl1','Krt19','Dsp','Fgf5'),  method = "umap", nCol = 4, cell.size = 0.1, show.legend = F, show.legend.combined = F)

```



