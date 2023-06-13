---
title: WGD Analysis
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /wgd/
category: mdanderson
---

## Flow Cytometry data of WGD induction on 230606
fcs files were transfered to
```bash
/volumes/seq/projects/wgd/230606_MCF10A_+-P53_WGD/fcs_data
```

Installing ggcyto and flowCore for R analysis
Using a new conda env for this.



```R
#install.packages("BiocManager")
#BiocManager::install(c("flowCore","ggcyto","flowStats"))

library(flowCore)
library(ggcyto)
library(flowStats)
library(patchwork)
setwd("/Volumes/seq/projects/wgd/230606_MCF10A_+-P53_WGD/fcs_data")
fcs.dir<-getwd()
frames <- lapply(dir(fcs.dir, pattern="*.fcs" ,full.names=TRUE), read.FCS)
frame_names<-basename(dir(fcs.dir, pattern="*.fcs" ,full.names=TRUE))
names(frames)<-frame_names
frames <- as(frames, "flowSet")



rectGate <- rectangleGate(filterId="cell_filter","FSC-A"=c(50000, 250000), "SSC-A"=c(20000, 200000))

plt<-autoplot(frames,"FSC-A","SSC-A",bins=100)+theme_minimal()+geom_gate(rectGate)
ggsave(plt,file="fscA_sscA_cellfilter.pdf")

frames_subset<-Subset(frames,rectGate)
plt<-autoplot(frames_subset,"FSC-A","SSC-A",bins=100)+theme_minimal()
ggsave(plt,file="fscA_sscA_cellfilter_postfilter.pdf")

cells<-filter(frames,rectGate)
summary(cells)

myTrans<-transformList(from=c("*DAPI*-A","*DAPI*-H"),to=c("*DAPI*-A","*DAPI*-H"),tfun=c(log,log))
cells_trans<-transform(frames_subset, myTrans)



autoplot(cells_trans,"*DAPI*-A","*PE-Texas Red*-A",bins=1000)+theme_minimal() +theme(strip.text.x = element_text(size = 3))+coord_cartesian(ylim=c(0,1000))
ggsave(plt,file="cells_dapi_gem.pdf")

cells_trans<-transform(frames_subset,transformList(c("*DAPI*-A"),c(log)))
gate_2n <- rectangleGate(filterId="2N_filter","*DAPI*-A"=c(10.3, 10.7), "*PE-Texas Red*-A"=c(0, 750))
autoplot(cells_trans,"*DAPI*-A","*PE-Texas Red*-A",bins=500)+theme_minimal()+geom_gate(gate_2n,colour="#CCCCFF")+coord_cartesian(xlim=c(10,12),ylim=c(0,1000))
gate_4n <- rectangleGate(filterId="4N_filter","*DAPI*-A"=c(11, 11.35), "*PE-Texas Red*-A"=c(0, 750))
autoplot(cells_trans,"*DAPI*-A","*PE-Texas Red*-A",bins=500)+theme_minimal()+geom_gate(gate_4n,colour="#9FE2BF")+coord_cartesian(xlim=c(10,12),ylim=c(0,1000))
gate_8n <- rectangleGate(filterId="8N_filter","*DAPI*-A"=c(11.5, 11.7), "*PE-Texas Red*-A"=c(0, 750))
autoplot(cells_trans,"*DAPI*-A","*PE-Texas Red*-A",bins=500)+theme_minimal()+geom_gate(gate_8n,colour="#6495ED")+coord_cartesian(xlim=c(10,12),ylim=c(0,1000))
plt<-autoplot(cells_trans,"*DAPI*-A","*PE-Texas Red*-A",bins=500)+coord_cartesian(xlim=c(10,12),ylim=c(0,1000))+theme_minimal()+geom_gate(gate_2n,colour="#CCCCFF")+geom_gate(gate_4n,colour="#9FE2BF")+geom_gate(gate_8n,colour="#6495ED")
ggsave(plt,file="cells_dapi_gem_gated.pdf")


#apply gates based on discrete numbers and plot.

lapply(1:length(cells_trans),function(x) {
		dip=Subset(as(cells_trans[[x]],"flowFrame"),gate_2n)
		gate_wgd <- rectangleGate(filterId="wgd_filter","*DAPI*-A"=c(10.3, 11.35), "*PE-Texas Red*-A"=c(0, 250))
		gate_g2m <- rectangleGate(filterId="wgd_filter","*DAPI*-A"=c(10.3, 11.35), "*PE-Texas Red*-A"=c(250, 750))
		tet_and_dip=Subset(as(cells_trans[[x]],"flowFrame"),gate_2n|gate_4n)
		tet=Subset(as(cells_trans[[x]],"flowFrame"),gate_4n)

		gate_tet<-as.ggplot(autoplot(tet_and_dip,"*DAPI*-A","*PE-Texas Red*-A",bins=100)+theme_minimal()+geom_gate(gate_2n,colour="#CCCCFF")+geom_gate(gate_4n,colour="#9FE2BF")+geom_gate(gate_wgd,colour="red")+geom_gate(gate_g2m,colour="purple")+coord_cartesian(xlim=c(10,12),ylim=c(0,1000)))
		plt_dip=as.ggplot(autoplot(dip,"*PE-Texas Red*-A"))+geom_vline(xintercept=250)+coord_flip()+ggtitle("2N")+theme_minimal()+coord_cartesian(xlim=c(0,1000),ylim=c(0,1000))
		plt_tet=as.ggplot(autoplot(tet,"*PE-Texas Red*-A"))+geom_vline(xintercept=250)+coord_flip()+ggtitle("4N")+theme_minimal()+coord_cartesian(xlim=c(0,1000),ylim=c(0,1000))
		plt<-(gate_tet)/(plt_dip|plt_tet)
		ggsave(plt,file=paste0(cells_trans[[x]]@description$`$SMNO`,"_gem.pdf"))
})


lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)


```
