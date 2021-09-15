---
title: CisTopic on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /cistopic_nmt/
category: CEDAR
---


Note: this takes up a lot of memory

```R
library(readr)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(EnsDb.Hsapiens.v86)

read_cov<-function(infile){
	print(paste("Reading in",infile))
	d<- as.data.frame(read_tsv(infile,col_names=c("chr","start","end","perc","met","unmet")))
	d<-makeGRangesFromDataFrame(d,keep.extra.columns=T,ignore.strand=T)
	return(d)}

summarize_met_over_tiles<-function(x){
	hits <- GenomicRanges::findOverlaps(dat[[x]], reg_tiled)
    xhits <- dat[[x]][queryHits(hits)]
    yhits <- reg_tiled[subjectHits(hits)]
    dat_tmp<-data.frame(feature=yhits$feat_name,tile=yhits$tile,met=xhits$met,unmet=xhits$unmet)
    dat_tmp<- as.data.frame(dat_tmp %>% group_by(tile) %>% summarize(met=sum(met),unmet=sum(unmet)))
    dat_tmp$cell_name<-names(dat)[x]
    return(dat_tmp)}

summarize_met_over_tiles_feat<-function(x){
	hits <- GenomicRanges::findOverlaps(dat[[x]], reg_tiled)
    xhits <- dat[[x]][queryHits(hits)]
    yhits <- reg_tiled[subjectHits(hits)]
    dat_tmp<-data.frame(feature=yhits$feat_name,tile=yhits$tile,met=xhits$met,unmet=xhits$unmet)
    dat_tmp<- as.data.frame(dat_tmp %>% group_by(tile,features) %>% summarize(met=sum(met),unmet=sum(unmet)))
    dat_tmp$cell_name<-names(dat)[x]
    return(dat_tmp)}


#variables
in_dir="/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/bismark_cov_out"
reg_in="/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed"
tile_n=200
c_context="CpG" #GpC or CpG
out_reg_name="gene"

#set up cell names
setwd(in_dir)
infile_list<-list.files(in_dir, pattern=paste0(".NOMe.",c_context,".cov.gz$"))
infile_names<-unlist(lapply(1:length(infile_list),function(x){strsplit(infile_list,"[.]")[[x]][1]}))

#add sample cell line names as metadata
annot<-data.frame(cell_name=infile_names)
annot$cellLine<-unlist(lapply(strsplit(annot$cell_name,"_"),"[",1))
annot$treatment<-substr(unlist(lapply(strsplit(annot$cell_name,"_"),"[",2)),1,1) #set up T47D cell treatment first
annot[startsWith(annot$cellLine,"M7"),]$treatment<-substr(annot[startsWith(annot$cellLine,"M7"),]$cellLine,3,3)
annot[startsWith(annot$cellLine,"BSM7"),]$treatment<-substr(annot[startsWith(annot$cellLine,"BSM7"),]$cellLine,5,5)
annot$cellLine_treatment<-paste(annot$cellLine,annot$treatment,sep="_")
annot[annot$cellLine %in% c("BSM7E6","M7C1A","M7C2B","M7E4C"),]$cellLine<-"MCF7" #clean up metadata a bit more
annot[annot$cellLine %in% c("T"),]$cellLine<-"T47D"
annot[annot$treatment %in% c("C"),]$treatment<-"control"
annot[annot$treatment %in% c("E"),]$treatment<-"estrogen"

####set up regulatory sites#######
reg_bed<-read_tsv(reg_in,col_names=c("chr","start","end","name")) #add upstream and downstream expansion (TODO)
reg_bed<-reg_bed[reg_bed$chr %in% (paste0("chr",c(1:22,"X"))),] #limit to autosomes and X chr
reg<-makeGRangesFromDataFrame(reg_bed,keep.extra.columns=T,ignore.strand=T) #make granges
#reg<-reduce(reg) #this line removes the gene names, causes a bug

#reg<-reg+expand_bp #expand gene ranges

#for gene body
	reg<-reg[width(reg@ranges)>=1000]
	reg_tiled_list<-tile(reg, n=tile_n) #tile genes into 100 equally spaced sites (can also provide width for bp steps)
	names(reg_tiled_list)<-reg$name
	reg_tiled<-unlist(reg_tiled_list) #unlist back into single granges
	reg_tiled$tile<-c(1:tile_n)
	reg_tiled$feat_name<-names(reg_tiled) #add gene name to tile

#for TFS
	#reg_tiled_list<-resize(reg, width = 50+(50*2), fix = "center") #resize centered
	#reg_tiled_list<-tile(reg_tiled_list, width=1) #tile genes by 100 bp steps
	#names(reg_tiled_list)<-reg$name
	#reg_tiled<-unlist(reg_tiled_list) #unlist back into single granges
	#reg_tiled$tile<-c(1:150)
	#reg_tiled$feat_name<-names(reg_tiled) #add gene name to tile

####set up data#######
dat<-mclapply(infile_list,read_cov,mc.cores=20) #read in cov data and make Grangeslist
passing_list<-unlist(lapply(dat,function(x) is(x,"GRanges")))
dat<-dat[passing_list] #NOTE there is an error reading in GpC data for some reason. For now just removing the failing cells
dat<-GRangesList(dat)
names(dat)<-infile_names[passing_list] #set list names
dat_tiled<-mclapply(1:length(dat),summarize_met_over_tiles,mc.cores=5) #sum methylation data over tiles (can use multicore also)
names(dat_tiled)<-names(dat)
dat_tiled<-do.call("rbind",dat_tiled)
dat_tiled<-merge(dat_tiled,annot,by="cell_name")

#summarize by cell line and treatment
dat<-as.data.frame(dat_tiled %>% group_by(cellLine,treatment,tile) %>% summarize(met=sum(met),unmet=sum(unmet)))
dat$met_perc<-(dat$met/(dat$met+dat$unmet))*100

#plot along the tiled feature
plt<-ggplot(dat,aes(x=as.numeric(tile),y=met_perc,color=paste(cellLine,treatment)))+geom_smooth(span=0.1)+theme_minimal()+ylim(c(0,100))
ggsave(plt,file=paste(out_reg_name,c_context,"metsummary.pdf",sep="."))
system(paste0("slack -F ",paste(out_reg_name,c_context,"metsummary.pdf",sep=".")," ryan_todo"))

# #summarize by cell line and treatment
# dat_feat<-as.data.frame(dat_tiled  %>% group_by(cellLine,treatment,tile,feature) %>% summarize(met=sum(met),unmet=sum(unmet)))

# #convert gene names to symbols in counts data
# genenames<-AnnotationDbi::select(org.Hs.eg.db, keys = dat_feat$feature, keytype = 'ENSEMBL', columns = 'SYMBOL',multiVals=first)
# dat_feat<-merge(dat_feat,genenames,by.x="feature",by.y="ENSEMBL")
# dat_feat$met_perc<-(dat_feat$met/(dat_feat$met+dat_feat$unmet))*100

# #plot along the tiled feature
# plt<-ggplot(dat_feat[dat_feat$feature=="ENSG00000136826",],aes(x=as.numeric(tile),y=met_perc,color=paste(cellLine,treatment)))+geom_smooth(span=0.1)+theme_minimal()

# ggsave(plt,file=paste(out_reg_name,c_context,"metsummary.pdf",sep="."))
# system(paste0("slack -F ",paste(out_reg_name,c_context,"metsummary.pdf",sep=".")," ryan_todo"))
```

Add tf List from chromvar
```R
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
    tfList <- getMatrixByID(JASPAR2020, ID=row.names(atac_sub@assays$chromvar@data)) 
    tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
```