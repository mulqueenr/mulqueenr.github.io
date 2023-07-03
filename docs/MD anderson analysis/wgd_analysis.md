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
library(ComplexHeatmap)
library(circlize)

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

out<-lapply(1:length(cells_trans),function(x) {
		dip=Subset(as(cells_trans[[x]],"flowFrame"),gate_2n)
		gate_wgd <- rectangleGate(filterId="wgd_filter","*DAPI*-A"=c(10.3, 11.35), "*PE-Texas Red*-A"=c(0, 250))
		gate_g2m <- rectangleGate(filterId="wgd_filter","*DAPI*-A"=c(10.3, 11.35), "*PE-Texas Red*-A"=c(250, 750))
		tet_and_dip=Subset(as(cells_trans[[x]],"flowFrame"),gate_2n|gate_4n)
		tet=Subset(as(cells_trans[[x]],"flowFrame"),gate_4n)
		if(nrow(tet)>1){
		dip_wgd<-(sum(dip$"*PE-Texas Red*-A"<250)/nrow(dip)*100)
		tet_wgd<-(sum(tet$"*PE-Texas Red*-A"<250)/nrow(tet)*100)
		gate_tet<-as.ggplot(autoplot(tet_and_dip,"*DAPI*-A","*PE-Texas Red*-A",bins=100)+theme_minimal()+geom_gate(gate_2n,colour="#CCCCFF")+geom_gate(gate_4n,colour="#9FE2BF")+geom_gate(gate_wgd,colour="red")+geom_gate(gate_g2m,colour="purple")+coord_cartesian(xlim=c(10,12),ylim=c(0,1000)))
		plt_dip=as.ggplot(autoplot(dip,"*PE-Texas Red*-A"))+geom_vline(xintercept=250,color="purple")+coord_flip()+ggtitle("2N")+theme_minimal()+coord_cartesian(xlim=c(0,1000))+geom_text(aes(label = as.character(dip_wgd,length=4),y=0.005,x=400), vjust = "inward", hjust = "inward")
		plt_tet=as.ggplot(autoplot(tet,"*PE-Texas Red*-A"))+geom_vline(xintercept=250,color="purple")+coord_flip()+ggtitle("4N")+theme_minimal()+coord_cartesian(xlim=c(0,1000))+geom_text(aes(label = as.character(tet_wgd),y=0.005,x=400), vjust = "inward", hjust = "inward")

		plt<-(gate_tet)/(plt_dip|plt_tet)
		print(plt)
		ggsave(plt,file=paste0(cells_trans[[x]]@description$`$SMNO`,"_gem.pdf"))
		} else {
			dip_wgd=0
			tet_wgd=0
		}
		return(c(cells_trans[[x]]@description$`$SMNO`,sum(dip$"*PE-Texas Red*-A"<250),sum(tet$"*PE-Texas Red*-A"<250),nrow(dip),nrow(tet)))
})

out<-as.data.frame(do.call("rbind",out))
colnames(out)<-c("sample","diploid_gem_neg","tetraploid_gem_neg","dip_all","tet_all")
out$diploid_gem_neg<-as.numeric(out$diploid_gem_neg)
out$tetraploid_gem_neg<-as.numeric(out$tetraploid_gem_neg)
out$dip_all<-as.numeric(out$dip_all)
out$tet_all<-as.numeric(out$tet_all)

out$diploid_gem_pos<-100-((out$diploid_gem_neg/(out$dip_all+out$tet_all))*100)
out$tetraploid_gem_pos<-100-((out$tetraploid_gem_neg/(out$dip_all+out$tet_all))*100)
out$wgd_perc<-(out$tetraploid_gem_neg/(out$dip_all+out$tet_all))*100
out$g2m_perc<-((out$tet_all-out$tetraploid_gem_neg)/(out$dip_all+out$tet_all)*100)
out<-out[out$sample!="MCF10A_Veh_nodapi",]
row.names(out)<-out$sample
out$cell_line<-c("MCF10A_P53KO","MCF10A_P53KO","MCF10A_P53KO","MCF10A_P53KO","MCF10A_P53KO","MCF10A_WT","MCF10A_WT","MCF10A_WT","MCF10A_WT","MCF10A_P53KO","MCF10A_P53KO","MCF10A_WT","MCF10A_WT")#just hard coding this for now, because we didnt name fcs files in a single format
out$treatment<-c("DCB","DCB","Veh_no2nd","Veh","SP600125","DCB","DCB","monastrol+MPI","RO-3306","RO-3306","monastrol+MPI","SP600125","Veh")
out$time<-c("24h","2h","2h","2h","2h","24h","2h","2h","24h","24h","2h","2h","2h")
saveRDS(out,file="summarized_data.rds")
row_ha = rowAnnotation(cell_line = out$cell_line, treatment= out$treatment,time=out$time,
	col = list(cell_line=c("MCF10A_P53KO" = "red", "MCF10A_WT"="blue"),
		treatment=c("DCB"="blue","Veh"="gray","SP600125"="green","monastrol+MPI"="orange","Veh_no2nd"="gray","RO-3306"="purple"),
		time=c("2h"="white","24h"="black")))

col_fun = colorRamp2(c(0, max(out$wgd_perc*1.2)), c("white", "red"))

plt<-Heatmap(out[c("wgd_perc","g2m_perc")],left_annotation=row_ha,show_row_names = FALSE,row_split=out$cell_line,col=col_fun)

pdf("summarized_wgd_percentage_of_tetcells.pdf")
plt
dev.off()

)#just hard coding this for now, because we didnt name fcs files in a single format
```
