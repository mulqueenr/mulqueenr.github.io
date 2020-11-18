---
title: Tbr1 Mutations
layout: page
author: Ryan Mulqueen
permalink: /tbr1/
category: sciATAC
---

# Processing for sciATAC portion for Tbr1 Patient specific mutation mouse models.

## BCL File Locations

For a test run we got about 10% of the run pool after PCR amplification of two plates.

Full breakdown of experimental setup is located [here.](https://docs.google.com/spreadsheets/d/1Px1OAE8vIi3GUXPny7OaVvYJHgGESCp4fyZKnLWW0UE/edit#gid=823628902)

{% include text-expand.html %}
```bash
  #First Test Run (Two PCR Plates)
  /home/groups/oroaklab/seq/madbum/201116_NS500556_0437_AH72CMAFX2

```

## Initial Processing of Files
Includes conversion of bcl to fastq files, barcode assignment, fastq splitting, alignment, removal of duplicate reads, calling peaks and looking at TSS enrichment.

### BCL to FASTQ Conversion
For this we use a wrapper function NextSeq2fastq which wraps around bcl2fastq (v2.19). The wrapper is just to make it easier, since it infers where the run folder and output folders are based on our directory structure on the clusters. 

This is read in bcl files from the raw run folder in:
/home/groups/oroaklab/seq/madbum

And output fastq files in:
/home/groups/oroaklab/fastq

```bash
  NextSeq2fastq -R 201116_NS500556_0437_AH72CMAFX2

  #This assumes the following variables for bcl2fastq:
  # DEFAULT VARIABLES
  #$run_path = "/home/groups/oroaklab/seq/madbum";
  #$fastq_paths = "/home/groups/oroaklab/fastq,/home/groups/oroakdata/fastq";
  #$bcl2fastq_version = "bcl2fastq/2.19.0";
  #$bcl_opts = "with-failed-reads,no-lane-splitting,fastq-compression-level=9,create-fastq-for-index-reads";
  #$bcl_ignore_opts = "with-failed-reads,no-lane-splitting,fastq-compression-level=9,create-fastq-for-index-reads,ignore-missing-bcls,ignore-missing-filter,ignore-missing-positions,ignore-missing-controls";
  #$run_processing_log_file = "/home/groups/oroaklab/fastq/run_processing.log";
  #@POSSIBLE_OUTS = ("Undetermined_S0_R1_001.fastq.gz", "Undetermined_S0_R2_001.fastq.gz", "Undetermined_S0_I1_001.fastq.gz", "Undetermined_S0_I2_001.fastq.gz");

```

### Demultiplexing the fastq files
After fastq files are generated we can then demultiplex them. By this, I mean that we are going to assign our index sequences based on the index reads from the run. 

Index cycles on the Nextseq are substantially more error prone than read cycles, so we account for an error rate. For our 8 and 10 bp indexes, we allow 2 base mismatches for each (Hamming distance of 2). This is enough to still unambiguously assign the proper original primer for any index.

We use a scitools function which is a perl script to do this. This demultiplexer also assumes the same directory structure as NextSeq2fastq, meaning only the run name must be specified on our clusters.

*The script does the following:*
  1. Reads in the supplied index files. By default the index file is located here: /home/groups/oroaklab/src/scitools/scitools-dev/SCI_Indexes.txt. 
  These are of format: <IDX_NAME><\t><IDX_NUMBER><\t><IDX_SEQUENCE>.
  2. For each index, it then creates a hash (think python dictionary) for possible base changes for each index.
  3. It then reads in the fastq data, and splits index reads to appropriate lengths for sci-chemistry.
  4. It assigns the proper index sequence to all indexes based on the hash lookup table.
  5. It writes out properly assigned reads (all four indexes have a proper match in the hash) in the sci-format, where the read name becomes the corrected list of indexes (referred to as a library barcode).

```bash
  scitools fastq-dump -R 201116_NS500556_0437_AH72CMAFX2
```

This will output to /home/groups/oroaklab/fastq/201116_NS500556_0437_AH72CMAFX2.

I then set up a working directory a moved the properly assigned reads.

```bash
  mkdir /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates

  mv /home/groups/oroaklab/fastq/201116_NS500556_0437_AH72CMAFX2/201116_NS500556_0437_AH72CMAFX2.1.fq.gz \
  /home/groups/oroaklab/fastq/201116_NS500556_0437_AH72CMAFX2/201116_NS500556_0437_AH72CMAFX2.2.fq.gz \
  /home/groups/oroaklab/adey_lab/projects/tbr1_mus/201117_firstplates

```
### Generation of thorough annotation file and all meta data per cell 

Now that we have reads that assign to known scitools indexes, we have to get more specific. We are going to generate a proper annotation for our experiment based on our PCR and Tn5 primers used. We will do this for all possible index combinations as a ".annot" file.
scitools assumes annot files are in the following format: BARCODE   ANNOTATION

BARCODE structure is based on the sequencer. Because we use a Nextseq most commonly, we set them up as:

| 8bp_Tn5_i7_idx | 10bp_PCR_i7_idx | 8bp_Tn5_i5_idx | 10bp_PCR_i5_idx|

We used six tn5 plates: AB,CB,CC,BB,AC,CA (listed as Tn5_i5_idx,Tn5_i7_idx) so we will also limit barcodes to just those attainable from these tn5 combinations.

For this, I am once again looking at our experimental design from  [here.](https://docs.google.com/spreadsheets/d/1Px1OAE8vIi3GUXPny7OaVvYJHgGESCp4fyZKnLWW0UE/edit#gid=823628902)

I'm going to use a scitools helper function to do this, but a simple for loop through the index master list would work as well.

``` bash
#Since all plates are a random assortment of all Tn5 tagmentation, we can generate a simplified annotation schematic for PCR plates.

#simplified experiment annot, I just made this in the same directory with nano text editor.
#firstplates_annot.txt
Plate   Plate_SDSBSA_Condition  PCR_Index_i5    PCR_Index_i7    PCR_Cycles
1       Fresh   A       A       17
10      Old     A       B       17

scitools make-annot \
Plate_1+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,AA=ALL \
Plate_10+NEX,AB=ALL+NEX,CB=ALL+NEX,CC=ALL+NEX,BB=ALL+NEX,AC=ALL+NEX,CA=ALL+PCR,AB=ALL \
> firstplates.annot
```
### Splitting out our reads from the demultiplexed fastqs

Now that we know which barcodes belong to our reads, we can split them out from the full pool.

To do this we will use a scitools function that looks at fastq read 1 and read 2 and splits it into new files based on matches to our annotation.

```bash
scitools fastq-split -X -A firstplates.annot \
201116_NS500556_0437_AH72CMAFX2.1.fq.gz \
201116_NS500556_0437_AH72CMAFX2.2.fq.gz &

#The -X flag tells it to not write out barcodes which don't match. Those would be other sci formatted experiments on the same run
```
## Alignment

We have our reads, so now we can align them to the mouse reference genome. ATAC data is count based. 

We use another scitools function for convenience. It wraps bwa mem. We will use -t 10 threads for alignment and -r 10 threads for samtools sort afterwards.

```bash
  #For plate 1
  scitools fastq-align -t 10 -r 10 mm10 plate1 firstplates.Plate_1.1.fq.gz firstplates.Plate_1.2.fq.gz &
  #For plate 10
  scitools fastq-align -t 10 -r 10 mm10 plate10 firstplates.Plate_10.1.fq.gz firstplates.Plate_10.2.fq.gz &
```

### Deduplicate

Once we have aligned reads, we can mark PCR duplicates. Because we are sampling across the genome, it is highly unlikely that we capture the same exact start and end region twice. So we can use a combination of our barcode, and the start and end positions of a read to mark duplication rates.


```bash
  #For plate 1
  scitools bam-rmdup plate1.bam &
  #For plate 10
  scitools bam-rmdup plate10.bam &

  #Once these finish, plot the complexity per cell
  scitools plot-complexity plate10.complexity.txt &
  scitools plot-complexity plate1.complexity.txt &
```

## Looks good! Need more sequencing to get a better sense of it.

<!---
### Tabix fragment file generation


- Column Number  Name    Description

- 1 chrom   Reference genome chromosome of fragment
- 2 chromStart  Adjusted start position of fragment on chromosome.
- 3 chromEnd    Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval.
- 4 barcode The 10x cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment.
- 5 duplicateCount  The number of PCR duplicate read pairs observed for this fragment. Sequencer-created duplicates, such as Exclusion Amp duplicates created by the NovaSeqT instrument are excluded from this count.



```bash
  #Organoid processing
  input_bam="orgo.ID.bam"
  output_name="orgo"
  tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
  bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
  samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
  $tabix -p bed $output_name.fragments.tsv.gz &
```

# sciATAC Full Processing in R

## Generating Seurat Objects

Using R v4.0 and Signac v1.0 for processing.


```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  # make counts matrix from sparse matrix
  IN<-as.matrix(read.table("orgo.500.counts.sparseMatrix.values.gz"))
  IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
  COLS<-read.table("orgo.500.counts.sparseMatrix.cols.gz")
  colnames(IN)<-COLS$V1
  ROWS<-read.table("orgo.500.counts.sparseMatrix.rows.gz")
  row.names(IN)<-ROWS$V1

  #Read in fragment path for coverage plots
  orgo_fragment.path="./orgo.fragments.tsv.gz"

  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

  # change to UCSC style since the data was mapped to hg38
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

  #Generate ChromatinAssay Objects
  orgo_chromatinassay <- CreateChromatinAssay(
    counts = IN,
    genome="hg38",
    min.cells = 1,
    annotation=annotations,
    sep=c("_","_"),
    fragments=orgo_fragment.path
  )

  #Create Seurat Object
  orgo_atac <- CreateSeuratObject(
    counts = orgo_chromatinassay,
    assay = "peaks",
  )

  #Meta.data to be updated after clustering


  #saving unprocessed SeuratObject
  saveRDS(orgo_atac,file="orgo_SeuratObject.Rds")
```

### Plotting and updating metadata

```R
  #renaming annot for simplified annotation file making
  #rename processing_ processing. *annot
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  #Set up annotation summaries to contain the same information, in same column order
  first_annot_append<-read.table("first_prep_summary_statistics_per_cell.tsv",header=T)
  first_annot_append<-first_annot_append[c("cellID",
                                           "tn5_plate",
                                           "column","row",
                                           "tn5_i5","tn5_i5_idx_name","tn5_i5_idx_seq",
                                           "tn5_i7","tn5_i7_idx_name","tn5_i7_idx_seq",
                                           "pcr_i7_idx_seq","pcr_i5_idx_seq",
                                           "total_reads","uniq_reads","perc_uniq",
                                           "prep","orgID","cell_line","differentiation_exp","DIV",
                                           "freezing_protocol","sort_gate","treatment","organoid")]

  second_annot_append<-read.table("second_prep_summary_statistics_per_cell.tsv",header=T)
  second_annot_append$freezing_protocol<-"Flash_Frozen" #change this for the DIV90 cirm 43 diff exp 5 organoids
  second_annot_append[(second_annot_append$DIV=="90" & second_annot_append$differentiation_exp=="5"),]$freezing_protocol<-"Slow_Freeze"
  second_annot_append$sort_gate<-"NA"
  second_annot_append$treatment<-"No"
  second_annot_append$organoid<-second_annot_append$orgID

  second_annot_append<-second_annot_append[c("cellID",
                                           "tn5_plate",
                                           "column","row",
                                           "tn5_i5","tn5_i5_idx_name","tn5_i5_idx_seq",
                                           "tn5_i7","tn5_i7_idx_name","tn5_i7_idx_seq",
                                           "pcr_i7_idx_seq","pcr_i5_idx_seq",
                                           "total_reads","uniq_reads","perc_uniq",
                                           "prep","orgID","cell_line","differentiation_exp","DIV",
                                           "freezing_protocol","sort_gate","treatment","organoid")]

  annot_append<-rbind(first_annot_append,second_annot_append)
  #orgID and prep need to be accounted for to get unique organoids (there are duplicates in orgID)

  original_cluster<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/180703_sciATAC_Organoids/NatureLetterData/cellID_fullcharacterization.txt",header=T)
  original_cluster$cellID<-paste0(original_cluster$cellID,"_1")
  original_cluster<-original_cluster[c("cellID","Phenograph_Cluster")]
  colnames(original_cluster)<-c("cellID","original_cluster")

  orgo_atac<-readRDS(file="orgo_SeuratObject.Rds")
  orgo_atac@meta.data$cellID<-row.names(orgo_atac@meta.data)

  annot<-as.data.frame(orgo_atac@meta.data)
  annot<-merge(annot,annot_append,by="cellID",all.x=T)
  annot<-merge(annot,original_cluster,by="cellID",all.x=T)
  row.names(annot)<-annot$cellID

  orgo_atac@meta.data<-annot
  saveRDS(orgo_atac,file="orgo_SeuratObject.Rds")
  write.table(annot,file="merged_summary_statistics_per_cell.tsv",col.names=T,row.names=T,sep="\t",quote=F)

  #excluding differentiation experiment 4
  orgo_atac<-subset(orgo_atac,differentiation_exp %in% c("5","7"))
  orgo_cirm43<-subset(orgo_atac,cell_line=="CIRM43")

  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")
```

## Performing cisTopic and UMAP

```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(cisTopic)
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  cistopic_processing<-function(seurat_input,prefix){
      cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
      row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
      atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
      #Run warp LDA on objects
      atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(5,10,20:30,40,50,55),nCores=15,addModels=FALSE)
      print("Saving cistopic models.")
      saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
  }
          

  cistopic_processing(seurat_input=orgo_cirm43,prefix="orgo_cirm43")

  cirm43_cistopic_models<-readRDS("orgo_cirm43.CisTopicObject.Rds")


  #Setting up topic count selection
  pdf("cirm43_model_selection.pdf")
  par(mfrow=c(1,3))
  cirm43_cistopic_models <- selectModel(cirm43_cistopic_models, type='derivative')
  dev.off()
  system("slack -F cirm43_model_selection.pdf ryan_todo")


  ###############################################
  #Loop through cistopic models
  cistopic_loop<-function(topic_number,object_input,models_input){
      models_input<-selectModel(models_input,select=topic_number)
      #perform UMAP on topics
      topic_df<-as.data.frame(models_input@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      print("Performed UMAP.")
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")
     
      #combine with seurat object    
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$umap<-cistopic_umap
      
      #finally plot
      plt<-DimPlot(object_input,group.by=c('DIV','cell_line'),size=0.1)+ggtitle(as.character(topic_number))
      return(plt)
  }

  library(patchwork)

  plt_list<-lapply(cirm43_cistopic_models@calc.params$runWarpLDAModels$topic,
                     FUN=cistopic_loop,
                     object_input=orgo_cirm43,
                     models_input=cirm43_cistopic_models)
  plt_list<-wrap_plots(plt_list)
  ggsave(plt_list,file="cirm43.umap_multipleTopicModels_clustering.png",height=20,width=60,limitsize=FALSE)

  ###############################################



  #set topics based on derivative
  cirm43_selected_topic=27
  cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,select=cirm43_selected_topic,keepModels=T)

  #saving model selected RDS
  saveRDS(cirm43_cisTopicObject,file="orgo_cirm43.CisTopicObject.Rds")

  ####Function to include topics and umap in seurat object
  cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


      #run UMAP on topics
      topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
      row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
      dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
      row.names(dims)<-colnames(topic_df)
      colnames(dims)<-c("x","y")
      dims$cellID<-row.names(dims)
      dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


      #Add cell embeddings into seurat
      cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
      colnames(cell_embeddings)<-cisTopicObject@cell.names
      n_topics<-nrow(cell_embeddings)
      row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
      cell_embeddings<-as.data.frame(t(cell_embeddings))

      #Add feature loadings into seurat
      feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
      row.names(feature_loadings)<-paste0("topic_",1:n_topics)
      feature_loadings<-as.data.frame(t(feature_loadings))

      #combined cistopic results (cistopic loadings and umap with seurat object)
      cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
      umap_dims<-as.data.frame(as.matrix(dims[2:3]))
      colnames(umap_dims)<-c("UMAP_1","UMAP_2")
      row.names(umap_dims)<-dims$cellID
      cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
      object_input@reductions$cistopic<-cistopic_obj
      object_input@reductions$umap<-cistopic_umap

      n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))

      object_input <- FindNeighbors(
        object = object_input,
        reduction = 'cistopic',
        dims = 1:n_topics
      )
      object_input <- FindClusters(
        object = object_input,
        verbose = TRUE,
        resolution=resolution
      )

  return(object_input)}

  orgo_cirm43<-cistopic_wrapper(object_input=orgo_cirm43,cisTopicObject=cirm43_cisTopicObject,resolution=0.5)

  plt<-DimPlot(orgo_cirm43,group.by=c('DIV','cell_line','prep','orgID','differentiation_exp','seurat_clusters','original_cluster'),size=0.1)
  ggsave(plt,file="cirm43.umap.png",width=20)
  ggsave(plt,file="cirm43.umap.pdf",width=20)

  i="cirm43.umap.png"
  system(paste0("slack -F ",i," ryan_todo"))#post to ryan_todo
         
  ###save Seurat file
  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")
```

### Statistics on cell reads

```R
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(dplyr)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  #Add FRIP to meta data
  frip<-read.table("orgo.500.fracOnTarget.values")
  colnames(frip)<-c("cellID","frip")
  orgo_atac$FRIP<-frip[match(orgo_atac$cellID,frip$cellID,),]$frip
  orgo_cirm43$FRIP<-frip[match(orgo_cirm43$cellID,frip$cellID,),]$frip
  orgo_cirm87$FRIP<-frip[match(orgo_cirm87$cellID,frip$cellID,),]$frip
  orgo_atac<-saveRDS(orgo_atac,"orgo_SeuratObject.Rds")
  orgo_cirm43<-saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.Rds")

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  #Cluster summaries
  dat<-orgo_cirm43@meta.data
  dat_sum<-as.data.frame(dat %>% 
  group_by(differentiation_exp,DIV,treatment,seurat_clusters) %>% 
  summarize(mean=mean(uniq_reads),sd=sd(uniq_reads),median=median(uniq_reads),mean_FRIP=mean(FRIP),cell_count=n(),organoid_count=length(unique(orgID))))
  write.table(dat_sum,"cirm43_cluster_summary_statistics.tsv",col.names=T,row.names=T,quote=F,sep="\t")

  system("slack -F cirm43_cluster_summary_statistics.tsv ryan_todo")
```

### Differential Accessibillity on Clusters

```R
  
```

### Performing GREAT on DA peaks


```R
  #mkdir GREAT_analysis

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  #To perform GREAT on peaks for enrichment per cluster
  write("Performing GREAT on all enriched sites per annotation group", stderr())
  library(rGREAT)

  #format data as bed file all seurat objects have the same peak list
  write("Preparing Background Set as all called peaks.", stderr())
  orgo_bg_bed<-do.call("rbind",strsplit(unlist(orgo_cirm43@assays$peaks@counts@Dimnames[1]),"[-]"))
  orgo_bg_bed<-as.data.frame(orgo_bg_bed)
  colnames(orgo_bg_bed)<-c("chr","start","end")
  orgo_bg_bed$start<-as.numeric(as.character(orgo_bg_bed$start))
  orgo_bg_bed$end<-as.numeric(as.character(orgo_bg_bed$end))

  cirm43_da_peaks<-read.table("cirm43.onevone.da_peaks.txt",header=T)

  write("Beginning loop through all annotation groups.", stderr())

  great_processing<-function(enriched_group_input,peak_dataframe,prefix){
      #subset bed file to peaks enriched in input group
      orgo_bed<-as.data.frame(do.call("rbind",strsplit(orgo_da_peaks[orgo_da_peaks$enriched_group==enriched_group_input,]$da_region,"-")))
      colnames(orgo_bed)<-c("chr","start","end")
      orgo_bed$start<-as.numeric(as.character(orgo_bed$start))
      orgo_bed$end<-as.numeric(as.character(orgo_bed$end))
      
      #run GREAT using all peaks as background
      write(paste("Using",nrow(orgo_bed), "DA peaks from",enriched_group_input), stderr())
      job = submitGreatJob(orgo_bed,orgo_bg_bed,species="hg38",request_interval=30)
      tb = getEnrichmentTables(job, ontology = c("GO Molecular Function", "GO Biological Process","GO Cellular Component"))
      tb = getEnrichmentTables(job, category = c("GO","Phenotype","Genes"))
      #Plot gene association
      pdf(paste0("./GREAT_analysis/",prefix,"_DApeaks_",enriched_group_input,".GeneAssociation.pdf"))
      plotRegionGeneAssociationGraphs(job)
      dev.off()

      for (j in 1:length(names(tb))){
            write(paste("Outputting DA GREAT Analysis for", enriched_group_input, as.character(names(tb))[j]), stderr())
            tabl_name<-gsub(" ","",as.character(names(tb))[j])
            write.table(as.data.frame(tb[[j]]),file=paste0("./GREAT_analysis/",prefix,"_DApeaks_",enriched_group_input,".",tabl_name,".txt"),sep="\t",col.names=T,row.names=T,quote=F)
        }
  }

  library(parallel)
  mclapply(unique(cirm43_da_peaks$enriched_group), FUN=great_processing, peak_dataframe=cirm43_da_peaks,prefix="cirm43",mc.cores=10)
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

  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  #Read in data and modify to monocle CDS file
  #read in RDS file.

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species =9606, all_versions = FALSE))

  # Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
  motif.matrix <- CreateMotifMatrix(
    features = granges(orgo_atac),
    pwm = pfm,
    genome = 'hg38',
    use.counts = FALSE)

  # Create a new Mofif object to store the results
  motif <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm)


  # Add the Motif object to the assays and run ChromVar
  ###CIRM43###
  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  orgo_cirm43 <- SetAssayData(
    object = orgo_cirm43,
    assay = 'peaks',
    slot = 'motifs',
    new.data = motif)
  orgo_cirm43 <- RegionStats(object = orgo_cirm43, genome = BSgenome.Hsapiens.UCSC.hg38)
  orgo_cirm43 <- RunChromVAR( object = orgo_cirm43,genome = BSgenome.Hsapiens.UCSC.hg38)
  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")
```
### Differential Motif Accessibility

```R
  ###Differential TF Accessibility by cluster###
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")

  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(parallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")


  #Perform One vs. rest DA enrichment

  write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

  DefaultAssay(orgo_cirm43) <- 'chromvar'

  #set up an empty list for looping through
  cirm43_tf<-list()

  #define DA functions for parallelization
  #Use LR test for atac data
  da_one_v_rest<-function(i,obj,group){
      da_peaks_tmp <- FindMarkers(
          object = obj,
          ident.1 = i,
          group.by = group,
          test.use = 'LR',
          latent.vars = 'nCount_peaks',
          only.pos=T
          )
      da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
      da_peaks_tmp$enriched_group<-c(i)
      da_peaks_tmp$compared_group<-c("all_other_cells")
      return(da_peaks_tmp)
    }

  da_one_v_one<-function(i,obj,group,j_list){
      i<-as.character(i)
      da_tmp_2<-list()
      for (j in j_list){
          if ( i != j){
          da_peaks_tmp <- FindMarkers(
              object = obj,
              ident.1 = i,
              ident.2 = j,
              group.by = group,
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              only.pos=T
              )
          da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
          da_peaks_tmp$enriched_group<-c(i)
          da_peaks_tmp$compared_group<-c(j)
          da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
          }
      }
      return(da_tmp_2)
    }

  #Perform parallel application of DA test
  library(parallel)
  n.cores=length(unique(orgo_cirm43@meta.data$seurat_clusters))
  cirm43_tf<-mclapply(
      unique(orgo_cirm43@meta.data$seurat_clusters),
      FUN=da_one_v_rest,
      obj=orgo_cirm43,
      group="seurat_clusters",
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1vrest DA
  cirm43_tf<-do.call("rbind",cirm43_tf)

  write("Outputting One v Rest DA Table.", stderr())
  write.table(cirm43_tf,file="cirm43.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

  dat<-read.table("cirm43.onevrest.da_tf.txt",header=T,sep="\t")
  #To convert JASPAR ID TO TF NAME
  dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(dat,file="cirm43.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)
  dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% slice(1:2) #grabbing top 2 most significant peaks to label
  plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(dat=dat_select,aes(label=da_tf),force=3)+theme_bw()
  ggsave(plt,file="cirm43_oncevrest.da_tf.pdf")

  #Empty list to rerun for 1v1 comparisons
  cirm43_tf<-list()
      
  n.cores=length(unique(orgo_cirm43@meta.data$seurat_clusters))
  cirm43_tf<-mclapply(
      unique(orgo_cirm43@meta.data$seurat_clusters),
      FUN=da_one_v_one,
      obj=orgo_cirm43,
      group="seurat_clusters",
      j_list=do.call("as.character",list(unique(orgo_cirm43@meta.data$seurat_clusters))),
      mc.cores=n.cores)

  #Merge the final data frame from the list for 1v1 DA
  cirm43_tf<-do.call("rbind",do.call("rbind",cirm43_tf))

  write("Outputting One v One DA Table.", stderr())
  write.table(cirm43_tf,file="cirm43.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

  dat<-read.table("cirm43.onevone.da_tf.txt",header=T,sep="\t")
  #To convert JASPAR ID TO TF NAME
  dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
  write.table(dat,file="cirm43.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)
```

## Cicero for Coaccessible Networks


```R
  library(Signac)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(monocle3)
  library(cicero)
  library(EnsDb.Hsapiens.v86)
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")


  #Cicero processing function
  cicero_processing<-function(object_input=orgo_atac,prefix="orgo_atac"){

      #Generate CDS format from Seurat object
      atac.cds <- as.cell_data_set(object_input,group_by="seurat_clusters")

      # convert to CellDataSet format and make the cicero object
      print("Making Cicero format CDS file")
      atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)
      saveRDS(atac.cicero,paste(prefix,"atac_cicero_cds.Rds",sep="_"))
      
      genome <- seqlengths(object_input) # get the chromosome sizes from the Seurat object
      genome.df <- data.frame("chr" = names(genome), "length" = genome) # convert chromosome sizes to a dataframe
      
      print("Running Cicero to generate connections.")
      conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 10) # run cicero
      saveRDS(conns,paste(prefix,"atac_cicero_conns.Rds",sep="_"))
      
      print("Generating CCANs")
      ccans <- generate_ccans(conns) # generate ccans
      saveRDS(ccans,paste(prefix,"atac_cicero_ccans.Rds",sep="_"))
      
      print("Adding CCAN links into Seurat Object and Returning.")
      links <- ConnectionsToLinks(conns = conns, ccans = ccans) #Add connections back to Seurat object as links
      Links(object_input) <- links
      return(object_input)
  }


  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

  orgo_cirm43<-cicero_processing(object_input=orgo_cirm43,prefix="orgo_cirm43")

  saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.Rds")

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

  conns<-as.data.frame(readRDS("orgo_cirm43_atac_cicero_conns.Rds"))
  orgo_cirm43.cicero<-readRDS("orgo_cirm43_atac_cicero_cds.Rds")
  geneactivity_processing(cds_input=as.cell_data_set(orgo_cirm43,group_by="seurat_clusters"),conns_input=conns,prefix="cirm43_atac")

  #These can be added to the seurat object as a new assay later

  #Read in unnormalized GA
  cicero_gene_activities<-readRDS("cirm43_atac.unnorm_GA.Rds")
  orgo_cirm43[['GeneActivity']]<- CreateAssayObject(counts = cicero_gene_activities) 

  # normalize
  orgo_cirm43 <- NormalizeData(
    object = orgo_cirm43,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(orgo_cirm43$nCount_peaks)
  )
  saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.Rds")
```
# TO BE ADDED scRNA PREPROCESSING SECTION


## Celltype Assignment of Clusters

Cell Type Assignment of Organoid Clusters
Doing this in three parts.
1. Using bulk sorted RG, IPC, eN and iN RNA markers compared to our ATAC cluster gene activity scores
2. Using bulk sorted RG, IPC, eN and iN ATAC motifs compared to our ATAC cluster motifs
3. Using single-cell Primary Cortex RG, IPC, eN and iN annotated cells to define signatures and perform CCA for label transfer


```R
  #https://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  set.seed(1234)
  library(reshape2)
  library(dplyr)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)
  # Load the pre-processed scRNA-seq and scATAC-seq data

  #Public RNA
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data")
  pubprimary<-readRDS("PublicPrimary.SeuratObject.rds")
  #perform random subsampling, since cross data integration is robust to cell count and its taking forever
  #using 1/10th the cells (~10k)
  pubprimary <- subset(pubprimary, cells = sample(x = colnames(pubprimary@assays$RNA@data), size = length(colnames(pubprimary@assays$RNA@data))/10) )
  pubprimary<-SetIdent(pubprimary,value="Type")
  #subset to cell types expected to occur in organoids
  pubprimary<-subset(pubprimary,idents=c("Excitatory Neuron","Inhibitory Neuron","IPC","Radial Glia"))

  #Our RNA
  #setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/rna_processing")
  #orgo_rna<-readRDS("orgo_rna.SeuratObject.rds")

  #Our ATAC
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS(file="orgo_cirm43.SeuratObject.Rds")


  #1. Using bulk sorted RG, IPC, eN and iN RNA markers compared to our ATAC cluster gene activity scores
  #Corticogenic data on basic cell types.
  #data from http://data.nemoarchive.org/5923ca16c51011e99da31f7757ebac1c/
  #described in https://www.nature.com/articles/s41586-020-2825-4?WT.ec_id=NATURE-202010&sap-outbound-id=60313C942AEB24BFE1AF8DD74FB2E05B7385E720#data-availability
  #Bulk ATAC Peaks for marker sorted cell types located in /home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/Song_2020
  #For overlapping these narrowPeak features do:

    markers<-c("CTCF","EMX1","EMX2","LHX2","PAX6","RFX4","SOX2",
               "TBR1","EOMES","NEUROD1","NEUROD2","NEUROG1","TGIF1","TGIF2",
               "DLX1","DLX2","DLX6","GSX2","LHX6",
               "POU3F3","POU3F2","TFAP4")
    #Setting up chromvar matrix from CIRM43
    tfList <- getMatrixByID(JASPAR2020, ID=row.names(orgo_cirm43@assays$chromvar@data)) 
    tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
    dat_tf<-orgo_cirm43@assays$chromvar@data
    row.names(dat_tf)<-tfList
    dat_tf<-data.frame(t(dat_tf))
    dat_tf$cellID<-row.names(dat_tf)

    #append cluster ID to column
    dat_tf$seurat_clusters<-orgo_cirm43@meta.data[match(orgo_cirm43@meta.data$cellID,dat_tf$cellID),]$seurat_clusters
    #subset to markers
    dat_tf<-dat_tf[colnames(dat_tf) %in% c("seurat_clusters",markers)]
    #reshape to long format
    dat_tf<-melt(dat_tf)
    #group by to summarize markers by column
    dat_tf<-as.data.frame(dat_tf %>% group_by(seurat_clusters,variable) %>% summarize(mean_chromvar=mean(value)))
    #plot as heatmap
    dat_tf<-dcast(dat_tf,seurat_clusters~variable)
    row.names(dat_tf)<-dat_tf$seurat_clusters
    dat_tf<-dat_tf[colnames(dat_tf) %in% markers]
    #set na values to 0 for clustering
    dat_tf[which(is.na(dat_tf),arr.ind=T)]<-0
    dat_tf<-as.data.frame(t(dat_tf))
    clus_order<-c("12","4","1","0","5","3","9","2","8","6","7","10","11","13")
    dat_tf<-dat_tf[colnames(dat_tf) %in% clus_order]

    plt<-Heatmap(dat_tf,
                row_order=match(markers,row.names(dat_tf))[!is.na(match(markers,row.names(dat_tf)))],
                column_order=clus_order
                               )
    pdf("cirm43_celltype_tfHeatmap.pdf")
    plt
    dev.off()
    system("slack -F cirm43_celltype_tfHeatmap.pdf ryan_todo")

  #2. Using bulk sorted RG, IPC, eN and iN ATAC motifs compared to our ATAC cluster motifs
    #Setting up gene activity matrix
    markers<-c("SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3",
               "NEUROG1","SSTR2","EOMES","PPP1R17","NEUROD4",
               "SLC17A7","NEUROD6","SATB2","TBR1","SLA",
               "DLX2","DLX1","LHX6","GAD1")
    dat_ga<-orgo_cirm43@assays$GeneActivity@data
    dat_ga<-data.frame(t(dat_ga))
    dat_ga$cellID<-row.names(dat_ga)
    #append cluster ID to column
    dat_ga$seurat_clusters<-orgo_cirm43@meta.data[match(orgo_cirm43@meta.data$cellID,dat_ga$cellID),]$seurat_clusters
    #subset to markers
    dat_ga<-dat_ga[colnames(dat_ga) %in% c("seurat_clusters",markers)]
    #reshape to long format
    dat_ga<-melt(dat_ga)
    #group by to summarize markers by column
    dat_ga<-as.data.frame(dat_ga %>% group_by(seurat_clusters,variable) %>% summarize(mean_ga=mean(value)))
    #plot as heatmap
    dat_ga<-dcast(dat_ga,seurat_clusters~variable)
    row.names(dat_ga)<-dat_ga$seurat_clusters
    dat_ga<-dat_ga[colnames(dat_ga) %in% markers]
    #zscore values
    dat_ga<-scale(dat_ga)
    #set na values to 0 for clustering
    dat_ga<-data.frame(t(dat_ga))
    colnames(dat_ga)<-as.character(0:(ncol(dat_ga)-1))
    clus_order<-c("12","4","1","0","5","3","9","2","8","6","7","10","11","13")
    dat_ga<-dat_ga[colnames(dat_ga) %in% clus_order]
    
    plt<-Heatmap(dat_ga,
    row_order=match(markers,row.names(dat_ga))[!is.na(match(markers,row.names(dat_ga)))],
    column_order=clus_order
                )
    pdf("cirm43_celltype_gaHeatmap.pdf")
    plt
    dev.off()
    system("slack -F cirm43_celltype_gaHeatmap.pdf ryan_todo")

  #3. Using single-cell Primary Cortex RG, IPC, eN and iN annotated cells to define signatures and perform CCA for label transfer
    markers<-read.table("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/pubprimary.markers.txt",header=T,sep="\t")

    transfer.anchors <- FindTransferAnchors(
      reference = pubprimary,
      reference.assay="RNA",
      query = orgo_cirm43,
      query.assay="GeneActivity",
      reduction   = 'cca',
      features=markers$gene,
      verbose=T
    )
    saveRDS(transfer.anchors,"orgo_cirm43.PublicPrimary.transferanchors.rds")

    predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      refdata = pubprimary$Type,
      weight.reduction = "cca",
      dims = 1:10
    )

    orgo_cirm43 <- AddMetaData(object = orgo_cirm43, metadata = predicted.labels)
    saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")

    plt1<-DimPlot(orgo_cirm43,group.by=c('predicted.id'),size=0.1)
    plt2<-FeaturePlot(orgo_cirm43,features=c('prediction.score.Radial.Glia'),pt.size=0.1)
    plt3<-FeaturePlot(orgo_cirm43,features=c('prediction.score.Excitatory.Neuron'),pt.size=0.1)
    plt4<-FeaturePlot(orgo_cirm43,features=c('prediction.score.Inhibitory.Neuron'),pt.size=0.1)
    plt5<-FeaturePlot(orgo_cirm43,features=c('prediction.score.IPC'),pt.size=0.1)

    plt<-plt1/plt2/plt3/plt4/plt5
    ggsave(plt,file="cirm43.predictedid.umap.png",width=10,height=30,limitsize=F)
    ggsave(plt,file="cirm43.predictedid.umap.pdf",width=10,height=30,limitsize=F)
    system("slack -F cirm43.predictedid.umap.png ryan_todo")

    predictdat<-orgo_cirm43@meta.data
    predictdat<-predictdat[startsWith(colnames(predictdat),"prediction.score")| colnames(predictdat) %in% c("seurat_clusters")]
    predictdat<-predictdat[!(colnames(predictdat) %in% c("prediction.score.max","predicted.id"))]

    predictdat<-melt(predictdat)
    predictdat<-as.data.frame(predictdat %>% group_by(seurat_clusters,variable) %>% summarize(average=mean(value)))

    predictdat$variable<-substr(predictdat$variable,18,length(predictdat$variable))
    predictdat<-predictdat[predictdat$variable %in% c("Excitatory.Neuron","Inhibitory.Neuron","IPC","Radial.Glia"),]
    predictdat<-dcast(predictdat,seurat_clusters~variable)
    row.names(predictdat)<-predictdat$seurat_clusters
    predictdat<-predictdat[!(colnames(predictdat) %in% c("seurat_clusters"))]
    predictdat<-as.data.frame(t(scale(predictdat,scale=F)))
    clus_order<-c("12","4","1","0","5","3","9","2","8","6","7","10","11","13")
    predictdat<-predictdat[colnames(predictdat) %in% clus_order]
    plt<-Heatmap(predictdat,
    row_order=c("Radial.Glia","IPC","Excitatory.Neuron","Inhibitory.Neuron"),
    column_order=clus_order)
    pdf("predictedid.heatmap.pdf")
    plt
    dev.off()
    system("slack -F predictedid.heatmap.pdf ryan_todo")

  #Based on the three separate measures assigning the following cell type ids.
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  orgo_cirm43<-readRDS(file="orgo_cirm43.SeuratObject.Rds")
  orgo_cirm43@meta.data$celltype<-"unknown"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$seurat_clusters %in% c("10"),]$celltype<-"neuroepithelial"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$seurat_clusters %in% c("12","4","1","0"),]$celltype<-"radial_glia"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$seurat_clusters %in% c("5","3"),]$celltype<-"intermediate_progenitor"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$seurat_clusters %in% c("9","2","8"),]$celltype<-"excitatory_neuron"
  orgo_cirm43@meta.data[orgo_cirm43@meta.data$seurat_clusters %in% c("6"),]$celltype<-"inhibitory_neuron"
  saveRDS(orgo_cirm43,file="orgo_cirm43.SeuratObject.Rds")
```

## Monocle


```R
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(ggplot2)
  set.seed(1234)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(cicero)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(patchwork)

  orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
  cirm43_subset<-subset(orgo_cirm43,cells=which(orgo_cirm43$celltype %in% c("radial_glia","intermediate_progenitor","excitatory_neuron")))

  monocle_processing<-function(prefix, seurat_input){
      atac.cds <- as.cell_data_set(seurat_input)
      atac.cds <- cluster_cells(cds = atac.cds, reduction_method = "UMAP") 
      #Read in cds from cicero processing earlier and continue processing
      atac.cds<- learn_graph(atac.cds, 
                             use_partition = F, 
                             learn_graph_control=list(
                                 minimal_branch_len=10,
                                 orthogonal_proj_tip=F,
                                 prune_graph=T))
      #plot to see nodes for anchoring
      plt1<-plot_cells(
                      cds = atac.cds,
                      show_trajectory_graph = TRUE,
                      color_cells_by="DIV",
                      label_leaves=T,
                      label_branch_points=F,
                      label_roots=T)
      plt2<-plot_cells(
                      cds = atac.cds,
                      show_trajectory_graph = TRUE,
                      label_leaves=T,
                      label_branch_points=F,
                      label_roots=T)    
      #Also make a plot of just node names for easier identification
      root_nodes<-as.data.frame(t(atac.cds@principal_graph_aux$UMAP$dp_mst))
      root_nodes$label<-row.names(root_nodes)
      plt3<-ggplot(
          root_nodes,
          aes(x=UMAP_1,y=UMAP_2))+
          geom_text(aes(label=label),size=3)+
          theme_bw()
      plt<-(plt1+plt2)/plt3
      ggsave(plt,file=paste(prefix,"DIV_trajectory.pdf",sep="_"),width=20)
      system(paste0("slack -F ",paste(prefix,"DIV_trajectory.pdf",sep="_")," ryan_todo"))
      return(atac.cds)
  }

  cirm43_subset.cicero<-monocle_processing(seurat_input=cirm43_subset,prefix="cirm43")

  #Then determine root nodes via plots and assign by order cells function.
  cirm43.cds <- order_cells(cirm43_subset.cicero, reduction_method = "UMAP", root_pr_nodes = c("Y_253")) #Chose youngest cells as root

  #Now replotting with pseudotime
  pdf("orgo_cirm43_trajectory.pseudotime.pdf")
  plot_cells(
    cds = cirm43.cds,
    show_trajectory_graph = TRUE,
  color_cells_by = "pseudotime"
  )
  dev.off()
  system("slack -F orgo_cirm43_trajectory.pseudotime.pdf ryan_todo")


  #Append pseudotime to meta data of seurat object
  cirm43_pseudotime<-as.data.frame(cirm43.cds@principal_graph_aux@listData$UMAP$pseudotime)
  colnames(cirm43_pseudotime)<-c("pseudotime")
  cirm43_pseudotime$cellID<-row.names(cirm43_pseudotime)
  pseudotime<-cirm43_pseudotime$pseudotime
  names(pseudotime)<-cirm43_pseudotime$cellID
  orgo_cirm43 <- AddMetaData(object = orgo_cirm43, metadata = pseudotime,col.name="pseudotime")
  saveRDS(orgo_cirm43,"orgo_cirm43.SeuratObject.Rds")
```

## Plot interactive scatter plot

```R
#Generating a 3D Plot via Plotly of the umap projection.
#Loading in additional libraries.
  setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
  library(Seurat)
  library(Signac)
  library(plotly)
  library(htmlwidgets)
  library(RColorBrewer)

orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

dat<-merge(orgo_cirm43@reductions$umap@cell.embeddings,orgo_cirm43@meta.data,by="row.names") 

dat$DIV<-as.character(dat$DIV) 
  
#Generate 3D Plot and standalone HTML widget
  p<-plot_ly(dat, type="scattergl", mode="markers", size=I(2),
    x= ~UMAP_1, y= ~UMAP_2,
    color=~celltype)
  p <- p %>% add_markers(color=~DIV)

htmlwidgets::saveWidget(as_widget(toWebGL(p)), "cirm43_umap.html",selfcontained=TRUE)

system("slack -F cirm43_umap.html ryan_todo")
```
### 3D Plotting for better trajectory visualization

```R
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(cicero)
library(SeuratWrappers)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(plotly)

orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
#cirm43_subset<-subset(orgo_cirm43,cells=which(orgo_cirm43$celltype %in% c("radial_glia","intermediate_progenitor","excitatory_neuron")))


#reading model selected RDS
library(cisTopic)

cirm43_cistopic_models<-readRDS(file="orgo_cirm43.CisTopicObject.Rds")
#set topics based on derivative
cirm43_selected_topic=27
cirm43_cisTopicObject<-cisTopic::selectModel(cirm43_cistopic_models,select=cirm43_selected_topic,keepModels=T)

####Function to include topics and umap in seurat object
cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


    #run UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
    dims<-as.data.frame(uwot::umap(t(topic_df),n_components=3))
    row.names(dims)<-colnames(topic_df)
    colnames(dims)<-c("x","y","z")
    dims$cellID<-row.names(dims)
    dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")


    #Add cell embeddings into seurat
    cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
    colnames(cell_embeddings)<-cisTopicObject@cell.names
    n_topics<-nrow(cell_embeddings)
    row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
    cell_embeddings<-as.data.frame(t(cell_embeddings))

    #Add feature loadings into seurat
    feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
    row.names(feature_loadings)<-paste0("topic_",1:n_topics)
    feature_loadings<-as.data.frame(t(feature_loadings))

    #combined cistopic results (cistopic loadings and umap with seurat object)
    cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
    umap_dims<-as.data.frame(as.matrix(dims[2:4]))
    colnames(umap_dims)<-c("UMAP_1","UMAP_2","UMAP_3")
    row.names(umap_dims)<-dims$cellID
    cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
    object_input@reductions$cistopic<-cistopic_obj
    object_input@reductions$umap<-cistopic_umap

    n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))

    object_input <- FindNeighbors(
      object = object_input,
      reduction = 'cistopic',
      dims = 1:n_topics
    )
    object_input <- FindClusters(
      object = object_input,
      verbose = TRUE,
      resolution=resolution
    )

return(object_input)}

orgo_cirm43<-cistopic_wrapper(object_input=orgo_cirm43,
  cisTopicObject=cirm43_cisTopicObject,
  resolution=0.5)


monocle_processing<-function(prefix, seurat_input){
    atac.cds <- as.cell_data_set(seurat_input)
    atac.cds <- cluster_cells(cds = atac.cds, reduction_method = "UMAP") 
    #Read in cds from cicero processing earlier and continue processing
    atac.cds<- learn_graph(atac.cds, 
                           use_partition = F, 
                           learn_graph_control=list(
                               minimal_branch_len=10,
                               orthogonal_proj_tip=F,
                               prune_graph=T))
    return(atac.cds)
}

cirm43_subset.cicero<-monocle_processing(seurat_input=orgo_cirm43,prefix="cirm43")

#Then determine root nodes via plots and assign by order cells function.
#cirm43.cds <- order_cells(cirm43_subset.cicero, 
#reduction_method = "UMAP", 
#root_pr_nodes = c("Y_253")) #Chose youngest cells as root


#Generating a 3D Plot via Plotly of the umap projection.
#Loading in additional libraries.
  library(plotly)
  library(htmlwidgets)
  library(RColorBrewer)

#cirm43_subset.cicero@principal_graph$UMAP

dat<-merge(orgo_cirm43@reductions$umap@cell.embeddings,orgo_cirm43@meta.data,by="row.names") 
dat$DIV<-as.character(dat$DIV) 
dat$differentiation_exp<-as.character(dat$differentiation_exp) 

  #Generate 3D Plot and standalone HTML widget
  p<-plot_ly(dat,
    x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    type="scatter3d", mode="markers",
    size=I(4), hoverinfo="none",
    color=~differentiation_exp) %>%
  add_markers(color=~DIV) %>%
  add_markers(color=~seurat_clusters)
  add_markers(color=~pseudotime) %>%

  htmlwidgets::saveWidget(as_widget(partial_bundle(p)), "cirm43_umap.3d.html",selfcontained=TRUE)

system("slack -F cirm43_umap.3d.html ryan_todo")


```

### Plotting ChromVAR motifs through pseudotime

 The following code is exploratory but in the end wasn't included in analysis for the manuscript.
I mainly just wanted to play around with network analysis a bit.

```R
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(cicero)
library(SeuratWrappers)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")

#Setting up chromvar matrix
tfList <- getMatrixByID(JASPAR2020, ID=row.names(orgo_cirm43@assays$chromvar@data)) 
tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
dat_tf<-orgo_cirm43@assays$chromvar@data
row.names(dat_tf)<-tfList
dat_tf<-data.frame(t(dat_tf))
dat_tf$cellID<-row.names(dat_tf)
cirm43_tf<-dat_tf

cirm43_pseudotime<-orgo_cirm43@meta.data[orgo_cirm43@meta.data$cell_line=="CIRM43",c("cellID","pseudotime")]
cirm43_tf<-merge(cirm43_tf,cirm43_pseudotime,by="cellID")
cirm43_tf<-cirm43_tf[complete.cases(cirm43_tf),]
cirm43_cell_row_order=order(cirm43_tf$pseudotime,decreasing=T)


cirm43_row_ha = rowAnnotation(
    cluster_id=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$seurat_clusters),
    organoid=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$orgID),
    DIV=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$DIV),
    experiment=as.factor(orgo_cirm43@meta.data[match(cirm43_tf$cellID,orgo_cirm43@meta.data$cellID),]$differentiation_exp)
)

#what is the distribution of pseudotime?
plt<-ggplot()+geom_density(aes(x=cirm43_tf$pseudotime))+theme_bw()
ggsave(plt,file="pseudotime_distribution.pdf")
system("slack -F pseudotime_distribution.pdf ryan_todo")
#divide pseudotime into equally sized (by cell count) bins

cirm43_tf<-cirm43_tf[!(colnames(cirm43_tf) %in% c("cellID","pseudotime"))]

vars <- apply(cirm43_tf, 2, var)
cirm43_tf<-cirm43_tf[vars > quantile(vars, 0.7)]

pdf("cirm43_tfmotif_pseudotimeordering.pdf")
Heatmap(as.matrix(cirm43_tf),column_km=7,left_annotation = cirm43_row_ha,
        row_order=cirm43_cell_row_order,
        show_row_names=F,
        clustering_distance_columns="maximum",
        column_names_gp = gpar(fontsize = 5)
)
dev.off()
system("slack -F cirm43_tfmotif_pseudotimeordering.pdf ryan_todo")

library(igraph)
cor_mat<-cor(cirm43_tf,method="pearson")
cor_mat[which(cor_mat<0.3,arr.ind=T)]<-0
graph<-graph.adjacency(cor_mat,weighted=TRUE,"directed",diag=F)
pdf("tf_motif.igraph.pdf")
plot.igraph(graph,edge.size=E(graph)$weight,vertex.size=3,vertex.label.cex=0.8)
dev.off()
system("slack -F tf_motif.igraph.pdf ryan_todo")
                 

```

## Analysis of ChromVAR TF motifs and Gene Activity Through Pseudotime

Decided to use a binning strategy to assess pseudotime signal.

```R

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(cicero)
library(SeuratWrappers)
library(ComplexHeatmap)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(parallel) 
library(zoo)

orgo_cirm43<-readRDS("orgo_cirm43.SeuratObject.Rds")
orgo_cirm43<-subset(orgo_cirm43,cells=which(orgo_cirm43$celltype %in% c("radial_glia","intermediate_progenitor","excitatory_neuron")))

#Set up a sliding window
# 1% bins with 0.33% step
ord_cells=orgo_cirm43$cellID[order(orgo_cirm43$pseudotime,decreasing=F)]
binwidth=as.integer(length(ord_cells)/100)
stepsize=as.integer(length(ord_cells)/300)

#Set up a nested list of sliding windows, binwidth is 1% data (~300 cells) and step is 0.3% (~100 cells)
timebin<-lapply(0:round(length(ord_cells)/stepsize), 
  FUN= function(i) ord_cells[(i*stepsize)+1:(binwidth+(i*stepsize))])

saveRDS(timebin,"cirm43_pseudotime_bins.rds")


# Generate gene annotations for CCAN assignment of marker genes
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

#Assign CCAN to RNA cell type markers by promoter overlap
markers<-c("CTCF","EMX1","EMX2","LHX2","PAX6","RFX4","SOX2",
             "TBR1","EOMES","NEUROD1","NEUROD2","NEUROG1","TGIF1","TGIF2",
             "DLX1","DLX2","DLX6","GSX2","LHX6",
             "POU3F3","POU3F2","TFAP4")
gene_annotation<-gene_annotation[gene_annotation$gene %in% markers,]
gene_annotation<-do.call("rbind",
  lapply(unique(gene_annotation$gene), function(x) {
  tmp<-gene_annotation[gene_annotation$gene==x,]
  return(tmp[which.max(abs(tmp$end-tmp$start)),])})) #select longest transcript of each gene for ccan overlap

gene_annotation<-makeGRangesFromDataFrame(gene_annotation,keep.extra.columns=T)

#Set up TF markers
tf_markers<-c("SOX2","PAX6","HES1","HOPX","VIM","GFAP","TNC","GPX3",
           "NEUROG1","SSTR2","EOMES","PPP1R17","NEUROD4",
           "SLC17A7","NEUROD6","SATB2","TBR1","SLA",
           "DLX2","DLX1","LHX6","GAD1")

#CCANs over pseudotime
  #First generating list of peaks in CCANs
  ccan<-Links(orgo_cirm43)
  ccan<-ccan[!is.na(ccan$group),]
  ccan<-split(ccan,f=ccan$group)
  #Plotting amount of peak membership per CCAN
    ccan_membership<-unlist(lapply(ccan,length))
    plt<-ggplot()+geom_density(aes(x=ccan_membership))+theme_bw()+xlim(c(0,500))
    ggsave(plt,file="ccan_membership.pdf")
    system("slack -F ccan_membership.pdf ryan_todo")

  #Getting list of peaks
    peaks<-row.names(orgo_cirm43@assays$peaks@data)
    seqname=unlist(lapply(strsplit(peaks,"-"),"[",1))
    start=unlist(lapply(strsplit(peaks,"-"),"[",2))
    end=unlist(lapply(strsplit(peaks,"-"),"[",3))
  peak_info<-makeGRangesFromDataFrame(data.frame(seqname=seqname,start=start,end=end,peak_id=peaks),keep.extra.columns=T)

  #find overlaps between peaks and CCAN regions (cicero merges similar peaks, so take all within overlap)
  #then from those overlapped peak per CCAN, grab peak data from the seurat object
  #then perform colSums (per cell peak accessibility summarized by the CCAN)
  #then normalize the data based on the cell specific size factor
  cell_accessibility_per_ccan<-function(x){
    ccan_tmp<-ccan[x]
    peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
    dat_tmp<-as.data.frame(orgo_cirm43@assays$peaks@data[peak_subset,])
    norm_dat_tmp<-colSums(dat_tmp)/orgo_cirm43@meta.data[match(orgo_cirm43@meta.data$cellID,colnames(dat_tmp)),]$nCount_peaks
    norm_dat_tmp
  }

  ccan_sum<-mclapply(1:length(ccan),cell_accessibility_per_ccan,mc.cores=20)
  ccan_sum<-as.data.frame(do.call("rbind",ccan_sum))
  row.names(ccan_sum)<-names(ccan)
  dim(ccan_sum)
  #[1]  3927 30293

  marker_annote_per_ccan<-function(x){
    ccan_tmp<-ccan[x]
    markers<-ifelse(countOverlaps(query=ccan_tmp,subject=gene_annotation)>0,
      unique(gene_annotation[unique(findOverlaps(query=ccan_tmp,subject=gene_annotation)@to),]$gene),
      "NA")
    markers
  }

  ccan_markers<-mclapply(1:length(ccan),marker_annote_per_ccan,mc.cores=10)
  ccan_markers<-as.data.frame(do.call("rbind",ccan_markers))
  row.names(ccan_markers)<-names(ccan)

  #summarize ccans through pseudotime cell bins
  #using mean here as the summary statistic
  ccan_summarize_per_bin<-function(x){
    bin_tmp<-timebin[[x]]
    return(apply(ccan_sum[colnames(ccan_sum) %in% bin_tmp],1,mean))
  }

  ccan_mean<-lapply(1:length(timebin),ccan_summarize_per_bin)
  ccan_mean<-as.data.frame(do.call("rbind",ccan_mean))
  row.names(ccan_mean)<-names(timebin)
  colnames(ccan_mean)<-names(ccan)
  saveRDS(ccan_mean,file="pseudotime_cirm43_ccan.rds")


  #Transcription Factor activity (ChromVAR TF motifs) over pseudotime bins
  chromvar_motifs_per_bin<-function(x){
    bin_tmp<-timebin[[x]]
    chromvar<-orgo_cirm43@assays$chromvar@data
    return(
      apply(chromvar[,colnames(chromvar) %in% bin_tmp],1,mean))
  }

  motif_mean<-mclapply(1:length(timebin),chromvar_motifs_per_bin,mc.cores=10)
  motif_mean<-as.data.frame(do.call("rbind",motif_mean))
    #Assign human readable TF motif names
  tfList <- getMatrixByID(JASPAR2020, ID=colnames(motif_mean))
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  colnames(motif_mean)<-tfList
  motif_mean<-as.data.frame(t(motif_mean))
  motif_mean<-motif_mean[complete.cases(motif_mean),]

  dim(motif_mean)
  saveRDS(motif_mean,file="pseudotime_cirm43_chromvarMotifs.rds")

  #Question: do we see a bias in CCAN opening through pseudotime with TF motif presence?
  #Take out motif presence in overlapping peaks with CCANs
  #This is taking from the chromvar motif peak association
  #generate background set per TF motif from all peaks
  #going to use a hypergeometric score since number of peaks per ccan is variable

  #Building contingency table
  ##             TFmotif  non-TFmotif
  ## CCAN            a  b
  ## Non-CCAN        c  d
  ## Total    f (a+c) g(b+d)

  dat_tmp_bg_f.<-colSums(orgo_cirm43@assays$peaks@motifs@data>0)
  dat_tmp_bg_g.<-colSums(orgo_cirm43@assays$peaks@motifs@data==0)

  fisher_tfmotifs_per_ccan<-function(ccan_idx=x,dat_tmp_bg_f=dat_tmp_bg_f.,dat_tmp_bg_g=dat_tmp_bg_g.){
    ccan_tmp<-ccan[ccan_idx]
    peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
    dat_tmp<-as.data.frame(orgo_cirm43@assays$peaks@motifs@data[peak_subset,])
    a_list<-colSums(dat_tmp>0)
    b_list<-colSums(dat_tmp==0)
    c_list<-dat_tmp_bg_f - a_list
    d_list<-dat_tmp_bg_g - b_list
    tf_dat<-unlist(
      lapply(colnames(dat_tmp),FUN=function(y) {
        a<-a_list[y]
        b<-b_list[y]
        c<-c_list[y]
        d<-d_list[y]
        contingency_table<-matrix(c(a,b,c,d),nrow=2)
        fisher.test(contingency_table,alternative="greater")$p.value}
      ))
    names(tf_dat)<-colnames(dat_tmp)
    return(tf_dat)
  }
    
  ccan_tf_enrich<-mclapply(1:length(ccan),function(x) fisher_tfmotifs_per_ccan(ccan_idx=x),mc.cores=10)
  ccan_tf_enrich<-as.data.frame(do.call("rbind",ccan_tf_enrich))
  #Assign human readable TF motif names
  tfList <- getMatrixByID(JASPAR2020, ID=colnames(ccan_tf_enrich))
  tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
  colnames(ccan_tf_enrich)<-tfList
  sig_tfs<-as.data.frame(t(apply(ccan_tf_enrich,1,p.adjust,method="bonferroni")))
  #sig_tfs <- sig_tfs[colSums(sig_tfs<0.05)>length(ccan)/20] #limit to TF motifs that reach significance in at least 5% ccans
  sig_tfs <- -log10(sig_tfs)
  saveRDS(sig_tfs,file="pseudotime_cirm43_ccan_tfmotif_enrichment.rds")

  #Question: do we see bias in GWAS associated sites within CCANs through pseudotime?
  #Take out GWAS sites in overlapping peaks with CCANs
  #This is taking from https://www.ebi.ac.uk/gwas/docs/file-downloads
  #generate background set per condition from all peaks
  #going to use a hypergeometric score since number of peaks per ccan is variable

  #Building contingency table
  ##             GWAS   non-GWAS
  ## CCAN            a  b
  ## Non-CCAN        c  d
  ## Total    f (a+c) g(b+d)
  gwas<-read.csv("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/gwas_catalog_v1.0-associations_e100_r2020-10-20.tsv",header=T,sep="\t")
  gwas$seqname<-paste0("chr",gwas$CHR_ID)
  gwas$start<-as.numeric(gwas$CHR_POS)
  gwas<-gwas[colnames(gwas) %in% c("seqname","start","MAPPED_GENE","DISEASE.TRAIT")]
  gwas$end<-gwas$start+1
  gwas<-gwas[complete.cases(gwas),]
  gwas<-makeGRangesFromDataFrame(gwas,keep.extra.columns=T)
  gwas<-split(gwas,f=gwas$DISEASE.TRAIT)

#GWAS analysis to be changed
  #dat_tmp_bg_f<-colSums(orgo_cirm43@assays$peaks@motifs@data>0)
  #dat_tmp_bg_g<-colSums(orgo_cirm43@assays$peaks@motifs@data==0)

  #tfmotifs_per_ccan<-function(x){
#   ccan_tmp<-ccan[x]
#   peak_subset<-peak_info[unique(findOverlaps(query=ccan_tmp,subject=peak_info)@to),]$peak_id
#   dat_tmp<-as.data.frame(orgo_cirm43@assays$peaks@motifs@data[peak_subset,])
#   a_list<-colSums(dat_tmp>0)
#   b_list<-colSums(dat_tmp==0)
#   c_list<-dat_tmp_bg_f - a_list
#   d_list<-dat_tmp_bg_g - b_list
#   tf_dat<-unlist(lapply(1:ncol(dat_tmp),FUN=fisher_test_columnwise,a_list=a_list,b_list=b_list,c_list=c_list,d_list=d_list))
#   names(tf_dat)<-colnames(dat_tmp)
#   return(tf_dat)
# }

# ccan_tf_enrich<-mclapply(1:length(ccan),tfmotifs_per_ccan,mc.cores=10)
# ccan_tf_enrich<-as.data.frame(do.call("rbind",ccan_tf_enrich))
  #Assign human readable TF motif names
# tfList <- getMatrixByID(JASPAR2020, ID=colnames(ccan_tf_enrich))
# tfList <-unlist(lapply(names(tfList), function(x) name(tfList[[x]])))
# colnames(ccan_tf_enrich)<-tfList
# sig_tfs<-as.data.frame(t(apply(ccan_tf_enrich,1,p.adjust,method="bonferroni")))
# sig_tfs <- sig_tfs[colSums(sig_tfs<0.05)>length(ccan)/20] #limit to TF motifs that reach significance in at least 5% ccans
# sig_tfs <- -log10(sig_tfs)

  sig_tfs<-readRDS(file="pseudotime_cirm43_ccan_tfmotif_enrichment.rds")
  motif_mean<-readRDS(file="pseudotime_cirm43_chromvarMotifs.rds")
  ccan_mean<-readRDS(file="pseudotime_cirm43_ccan.rds")

  timebin<-readRDS("cirm43_pseudotime_bins.rds")
  timebin_df<-as.data.frame(do.call("rbind",lapply(1:length(timebin),function(x) cbind(names(timebin)[x],unlist(timebin[x])))))
  colnames(timebin_df)<-c("pseudotime_bin","cellID")
  timebin_df<-merge(orgo_cirm43@meta.data,timebin_df,by="cellID")
  
  library(dplyr)
  library(reshape2)
  
  library(ComplexHeatmap)
  library(circlize)

  #Set up different annotation stacked barplots
  annotation_bin_summary<-function(x){
  tmp_timebin<-as.data.frame(timebin_df %>% group_by_("pseudotime_bin",x) %>% summarize(count=n()))
  colnames(tmp_timebin)<-c("bin","var","count")
  tmp_timebin<-dcast(tmp_timebin,bin~var,value.var="count",fill=0)
  row.names(tmp_timebin)<-tmp_timebin$bin
  tmp_timebin<-as.data.frame(t(tmp_timebin[2:ncol(tmp_timebin)]))
  tmp_timebin = as.data.frame(scale(tmp_timebin, center = FALSE, 
               scale = colSums(tmp_timebin)))
  return(tmp_timebin)}
  
  div_timebin<-annotation_bin_summary(x="DIV")
  celltype_timebin<-annotation_bin_summary(x="celltype")
  cluster_timebin<-annotation_bin_summary(x="seurat_clusters")
  organoid_timebin<-annotation_bin_summary(x="orgID")

  ha_barplots = HeatmapAnnotation(
    DIV = anno_barplot(t(div_timebin), name="DIV",bar_width=1,gp = gpar(color = 1:nrow(div_timebin),fill = 1:nrow(div_timebin))),
    celltype = anno_barplot(t(celltype_timebin), name="cell type",bar_width=1,gp = gpar(color = 1:nrow(celltype_timebin),fill= 1:nrow(celltype_timebin))),
    cluster = anno_barplot(t(cluster_timebin), name="cluster",bar_width=1,gp = gpar(color = 1:nrow(cluster_timebin),fill = 1:nrow(cluster_timebin))),
    organoid = anno_barplot(t(organoid_timebin), name="organoid",bar_width=1,gp = gpar(color = 1:nrow(organoid_timebin),fill = 1:nrow(organoid_timebin)))
    )

    lgd_list = list(
    Legend(labels = row.names(div_timebin), title = "DIV", type = "points", pch = 16, legend_gp = gpar(col = 1:nrow(div_timebin))),
    Legend(labels = row.names(celltype_timebin), title = "celltype", type = "points", pch = 16, legend_gp = gpar(col = 1:nrow(celltype_timebin))),
    Legend(labels = row.names(cluster_timebin), title = "cluster", type = "points", pch = 16, legend_gp = gpar(col = 1:nrow(cluster_timebin))),
    Legend(labels = row.names(organoid_timebin), title = "organoid", type = "points", pch = 16, legend_gp = gpar(col = 1:nrow(organoid_timebin)))
)

  vars<-apply(motif_mean, 1, var)
  motif_mean<-motif_mean[vars > quantile(vars, 0.8),]

  col_fun = colorRamp2(c(-4, -2,0, 2,4), rev(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")))

  plt3<-Heatmap(motif_mean,
    column_order=1:ncol(motif_mean),
    row_names_gp = gpar(fontsize = 3),
    clustering_distance_rows="spearman",
    col=col_fun,
    show_column_names=F,
    bottom_annotation=ha_barplots,
    row_km=3,
    show_heatmap_legend=T)

  pdf("cirm43_pseudotime.motifTF.heatmap.pdf",width=30)
  draw(plt3,annotation_legend_list=lgd_list)
  dev.off()
  system("slack -F cirm43_pseudotime.motifTF.heatmap.pdf ryan_todo")
  
  col_fun = colorRamp2(c(-4, -2,0, 2,4), rev(c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0")))

  plt1<-Heatmap(t(scale(ccan_mean)),
    column_order=1:nrow(ccan_mean),
    show_row_names=F,
    clustering_distance_rows="spearman",
    row_km=6,
    bottom_annotation=ha_barplots,
    col=col_fun
    )
  plt2<-Heatmap(sig_tfs,
    column_km=3,
    col = colorRamp2(c("white", "red"),breaks=c(0,5)),
    column_names_gp = gpar(fontsize = 3),
    clustering_distance_columns="spearman")

  plt<-plt1+plt2
  pdf("cirm43_pseudotime.ccan.heatmap.pdf")
  draw(plt1,annotation_legend_list=lgd_list)
  dev.off()

  system("slack -F cirm43_pseudotime.ccan.heatmap.pdf ryan_todo")
    
```

# 10x scRNA Seq Analysis


```python
Files stored in: /home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids
```


```python
#Get refence for hg38:
curl -O http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz  
tar -xvf refdata-cellranger-GRCh38-3.0.0.tar.gz
#stored in /home/groups/oroaklab/src/cellranger/cellranger-3.0.1
```


```python
#Merging all sequencing runs for scRNA seq
#191113 sequencing
RUN_PATH="/home/groups/oroaklab/seq/madbum/191113_NS500556_0361_AHTVFLAFXY"  
SAMPLE_INDEX_PATH="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids/10xRNA_SampleSheet.csv"  
OUTPUT_DIR="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids"  
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger mkfastq --run=$RUN_PATH --sample-sheet=$SAMPLE_INDEX_PATH --maxjobs=20 --mempercore=2 --output-dir=$OUTPUT_DIR --delete-undetermined --qc --ignore-dual-index --localmem=100  
#191118 sequencing
RUN_PATH="/home/groups/oroaklab/seq/madbum/191118_NS500556_0362_AHVYV7AFXY"
SAMPLE_INDEX_PATH="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids/10xRNA_SampleSheet.csv"
OUTPUT_DIR="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids"
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger mkfastq --run=$RUN_PATH --sample-sheet=$SAMPLE_INDEX_PATH --maxjobs=20 --mempercore=2 --output-dir=$OUTPUT_DIR --delete-undetermined --qc --ignore-dual-index --localmem=100
#191119 sequencing
RUN_PATH="/home/groups/oroaklab/seq/madbum/191119_NS500556_0363_AHTVL7AFXY"
SAMPLE_INDEX_PATH="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids/10xRNA_SampleSheet.csv"
OUTPUT_DIR="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids"
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger mkfastq --run=$RUN_PATH --sample-sheet=$SAMPLE_INDEX_PATH --maxjobs=20 --mempercore=2 --output-dir=$OUTPUT_DIR --delete-undetermined --qc --ignore-dual-index --localmem=100
```


```python
#Running all fastqs generated at once for each DIV
FASTQ_PATH="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids"
#DIV30
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger count --id=RM_Organoid_MergedRuns_DIV30 --fastqs=$FASTQ_path/HM3HFBGXC,$FASTQ_path/HTVFLAFXY,$FASTQ_path/HTVL7AFXY,$FASTQ_path/HVYV7AFXY, --transcriptome=/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/refdata-cellranger-GRCh38-3.0.0 --localcores=15 --localmem=75 --mempercore=3 --sample=RM_Organoid_MergedRuns_DIV30 &  
  
#DIV60
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger count --id=RM_Organoid_MergedRuns_DIV60 --fastqs=$FASTQ_path/HM3HFBGXC,$FASTQ_path/HTVFLAFXY,$FASTQ_path/HTVL7AFXY,$FASTQ_path/HVYV7AFXY, --transcriptome=/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/refdata-cellranger-GRCh38-3.0.0 --localcores=15 --localmem=75 --mempercore=3 --sample=RM_Organoid_MergedRuns_DIV60 &  
  
#DIV90
/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/cellranger count --id=RM_Organoid_MergedRuns_DIV90 --fastqs=$FASTQ_path/HM3HFBGXC,$FASTQ_path/HTVFLAFXY,$FASTQ_path/HTVL7AFXY,$FASTQ_path/HVYV7AFXY, --transcriptome=/home/groups/oroaklab/src/cellranger/cellranger-3.0.1/refdata-cellranger-GRCh38-3.0.0 --localcores=15 --localmem=75 --mempercore=3 --sample=RM_Organoid_MergedRuns_DIV90 &  

#Move to final working directory
rna_dir="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/rna_processing"
mv $FASTQ_PATH/RM_Organoid_MergedRuns_DIV30 $rna_dir
mv $FASTQ_PATH/RM_Organoid_MergedRuns_DIV60 $rna_dir
mv $FASTQ_PATH/RM_Organoid_MergedRuns_DIV90 $rna_dir

```

# R Script for scRNA Processing
This is based on the Seurat guided tutorials for 10x scRNA analysis. Guided PBMC Tutorial
## load and process 10x rna data



```python
#Set up merged Seurat Object
#Set up working directory and load libraries.
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/rna_processing")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)

#Read in 10X format filtered feature barcode matrices. Make into Seurat objects and merge.
#DIV90 data is bad, and will be discarded
orgo_rna_30<-Read10X("./RM_Organoid_MergedRuns_DIV30/outs/filtered_feature_bc_matrix")
orgo_rna_60<-Read10X("./RM_Organoid_MergedRuns_DIV60/outs/filtered_feature_bc_matrix")
#orgo_rna_90<-Read10X("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/191114_scRNA_Organoids/RM_Organoid_MergedRuns_DIV90/outs/filtered_feature_bc_matrix")
orgo_rna_30 <- CreateSeuratObject(counts = orgo_rna_30, project = "DIV30", min.cells = 5, min.features = 500)
orgo_rna_60 <- CreateSeuratObject(counts = orgo_rna_60, project = "DIV60", min.cells = 5, min.features = 500)
#orgo_rna_90 <- CreateSeuratObject(counts = orgo_rna_90, project = "DIV90", min.cells = 5, min.features = 500)
orgo_rna<-merge(orgo_rna_30,y=c(orgo_rna_60),add.cell.ids=c("DIV30","DIV60"),project="orgo")
saveRDS(orgo_rna,"orgo_rna.SeuratObject.rds")
```

## Perform preprocessing steps.
Largely following https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html


```python
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/rna_processing")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(patchwork)

orgo_rna<-readRDS("orgo_rna.SeuratObject.rds")
  
#Check low quality and dying cells via mitochrondrial read count
orgo_rna[["percent.mt"]] <- PercentageFeatureSet(orgo_rna, pattern = "^MT-")
plt<-VlnPlot(orgo_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plt,file="mitochondrial_feature.pdf")
system("slack -F mitochondrial_feature.pdf ryan_todo")
#looks like a lot are dying, going to process without filter for now

#Normalize and feature selection on data.
#Normalize merged data frame.
orgo_rna <- NormalizeData(orgo_rna)  

#Find variable features to be used in clustering.
orgo_rna<- FindVariableFeatures(orgo_rna, selection.method = "vst",mean.cutoff=c(0.1,10),nfeatures=3000)

#Perform a sanity check on variable feature selection.   
# Identify the 10 most variable genes
top10 <- head(VariableFeatures(orgo_rna), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(orgo_rna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plt<-plot1+plot2

ggsave(plt,file="rna_variance_perfeature.pdf",width=20)
system("slack -F rna_variance_perfeature.pdf ryan_todo")

#Dimensionality Reduction
#Scale data prior to PCA dimensionality reduction.
#set up row names
all.genes <- rownames(orgo_rna)

#scale data on all features (determined through the row names)
orgo_rna <- ScaleData(orgo_rna, features = all.genes)
  
#run PCA on data
orgo_rna <- RunPCA(orgo_rna, features = VariableFeatures(object = orgo_rna))

#Correlation of PCA with mitocondrial percentage
unlist(lapply(1:ncol(orgo_rna@reductions$pca@cell.embeddings),FUN=function(x) cor(orgo_rna@reductions$pca@cell.embeddings[,x],orgo_rna$percent.mt)))
# [1] -0.5229113669 -0.0685008854  0.2652615611 -0.2038899637  0.0756491031
# [6]  0.3945594648  0.2382365691 -0.1976708987  0.2105167401  0.1840434315
#[11]  0.0131652769  0.0220330030  0.0179212504 -0.0754405852 -0.0285540666
#[16] -0.0035171083  0.0212756056  0.0061314751 -0.0081926300 -0.0107569148
#[21]  0.0087961995  0.0049021395 -0.0178669753 -0.0208020419 -0.0233018669
#[26]  0.0029644556 -0.0006758250  0.0117499771  0.0330717885  0.0212195089
#[31]  0.0388860188 -0.0183217812 -0.0048459335  0.0092408376  0.0177255669
#[36] -0.0072072212 -0.0015833953 -0.0143009906  0.0008082347 -0.0423683706
#[41]  0.0092484075  0.0087918485  0.0114118100 -0.0076350385 -0.0077646584
#[46] -0.0005258069 -0.0330634015 -0.0153099380  0.0051146937  0.0023854102

#Sanity checks on PCA loadings.
#Generate a list of top features for principal components.
print(orgo_rna[["pca"]], dims = 2:5, nfeatures = 5)

#Plot feature PCA loadings for first handful of dimensions
plt<-VizDimLoadings(orgo_rna, dims = 2:6, reduction = "pca")
ggsave(plt,file="rna_pcaloadings.pdf",width=30,height=20)
system("slack -F rna_pcaloadings.pdf ryan_todo")

#Plot feature heatmaps for first handful of dimensions
plt<-DimHeatmap(orgo_rna, dims = 2:30, cells = 500, balanced = TRUE)
ggsave(plt,file="rna_pcaheatmap.pdf",height=100,width=50,limitsize=F)
system("slack -F rna_pcaheatmap.pdf ryan_todo")


#Plot scatter plots for first handful of dimensions
plot1<-DimPlot(orgo_rna, reduction = "pca")
plot2<-DimPlot(orgo_rna,reduction="pca",dims=c(3,4))
plot3<-DimPlot(orgo_rna,reduction="pca",dims=c(5,6))
plt<-plot1+plot2+plot3
ggsave(plt,file="rna_pcadim1through6.pdf",width=40)
system("slack -F rna_pcadim1through6.pdf ryan_todo")



#Determine PC dimensions to use for clustering and 2/3D projection.
plt<-ElbowPlot(orgo_rna)
ggsave(plt,file="rna_pcaelbow.pdf")
system("slack -F rna_pcaelbow.pdf ryan_todo")

#Going with 9 dimensions based on the elbow plot.
#Then clustering cells based on those 9 dimensions.
#Then performing UMAP for 3D projection using those 9 dimensions.

pca_dims<-8
orgo_rna <- FindNeighbors(
    object = orgo_rna,
    reduction = 'pca',
    dims = 2:pca_dims
    )
orgo_rna <- FindClusters(
    object = orgo_rna,
    verbose = TRUE,
    resolution=0.5
    )

orgo_rna <- RunUMAP(
    object = orgo_rna,
    reduction = 'pca',
    dims = 2:pca_dims,
    n.components = 2
    )

plt1<-DimPlot(orgo_rna, reduction = "umap")
plt2<-DimPlot(orgo_rna, reduction = "umap",group.by="orig.ident")
plt<-plt1+plt2
ggsave(plt,file="orgo_rna.umap.pdf",width=30)
system("slack -F orgo_rna.umap.pdf ryan_todo")

            
#Finding cluster markers for cell type identification
#Using seurats default find marker function.
#Loading dplyr for shorthand grouping of output.
library(dplyr)

#Running marker discovery
orgo_rna.markers <- FindAllMarkers(orgo_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
#Use dplyr for generation of marker lists and select top 10 genes (features) for plotting.
library(dplyr)
orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

top10 <- orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#Write out marker table
write.table(orgo_rna.markers,"orgo_rna.markers.txt",col.names=T,quote=F,sep="\t")

#Generate Heatmap of marker features
pdf("191222_rna_clustermarkers.pdf",width=30)
DoHeatmap(orgo_rna, features = top10$gene)
dev.off()

#Generate 2d umap projection plot
pdf("191222_rna_umap.pdf")
DimPlot(orgo_rna, reduction = "umap")
dev.off()

#Generate a pdf of stitched feature value plots for markers per cluster.
pdf("191222_rna_umapfeatureplots.pdf",width=100,height=100)
feat_genes<-orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
FeaturePlot(
object = orgo_rna,
features = feat_genes$gene,
pt.size = 0.1,
max.cutoff = 'q95',
ncol = 5
)
dev.off()
    
#Save the Seurat object for future cell type identification and coembedding with ATAC
saveRDS(orgo_rna,file="orgo_rna.SeuratObject.rds")
```

## Processing Public Data Sets 


```python
R
library(data.table)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(patchwork)

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data")
mat <- fread("curl http://cells.ucsc.edu/organoidreportcard/primary10X/exprMatrix.tsv.gz | zcat")
meta <- data.frame(fread("http://cells.ucsc.edu/organoidreportcard/primary10X/meta.tsv"), row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
so <- CreateSeuratObject(counts = mat, project = "PublicPrimary", meta.data=meta)
saveRDS(so,"PublicPrimary.SeuratObject.rds")

mat <- fread("curl http://cells.ucsc.edu/organoidreportcard/organoids10X/exprMatrix.tsv.gz | zcat")
meta <- data.frame(fread("http://cells.ucsc.edu/organoidreportcard/organoids10X/meta.tsv"), row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
so <- CreateSeuratObject(counts = mat, project = "PublicOrganoids", meta.data=meta)
saveRDS(so,"PublicOrganoid.SeuratObject.rds")
```

Processing Public Data with Oroak lab data


```python
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(patchwork)

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data")

puborgo<-readRDS("PublicOrganoid.SeuratObject.rds")
pubprimary<-readRDS("PublicPrimary.SeuratObject.rds")
orgo_rna<-readRDS("../organoid_finalanalysis/rna_processing/orgo_rna.SeuratObject.rds")

cortical<-merge(puborgo, y = c(pubprimary,orgo_rna), add.cell.ids = c("pub_orgo","pub_primary","orgo_oroak"), project = "merge")
saveRDS(cortical,"Merged.SeuratObject.rds")

#Check low quality and dying cells via mitochrondrial read count
cortical[["percent.mt"]] <- PercentageFeatureSet(cortical, pattern = "^MT-")

plt<-VlnPlot(cortical, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plt,file="mitochondrial_feature.pdf")
system("slack -F mitochondrial_feature.pdf ryan_todo")
#looks like a lot are dying, going to process without filter for now

#Normalize and feature selection on data.
#Normalize merged data frame.
cortical <- NormalizeData(cortical)  

#Find variable features to be used in clustering.
cortical<- FindVariableFeatures(cortical, selection.method = "vst", mean.cutoff=c(0.1,10),nfeatures=3000)

#Perform a sanity check on variable feature selection.   
# Identify the 10 most variable genes
top10 <- head(VariableFeatures(cortical), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cortical)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plt<-plot1+plot2

ggsave(plt,file="rna_variance_perfeature.pdf",width=20)
system("slack -F rna_variance_perfeature.pdf ryan_todo")

#Dimensionality Reduction
#Scale data prior to PCA dimensionality reduction.
#set up row names
all.genes <- rownames(cortical)

#scale data on all features (determined through the row names)
cortical <- ScaleData(cortical, features = all.genes)
  
#run PCA on data
cortical <- RunPCA(cortical, features = VariableFeatures(object = cortical))

"""
PC_ 1
Positive:  STMN2, MEF2C, GPM6A, RTN1, SLA, MYT1L, NEUROD6, MAPT, ZBTB18, ARPP21
           DCX, NRXN1, SATB2, SHTN1, HMP19, SYT1, NEUROD2, MEG3, GRIA2, NFIB
           PCLO, ANK2, CEP170, NCAM1, CHL1, GRIN2B, KIDINS220, NFIX, TCF4, ADCY1
Negative:  VIM, CRABP1, HSPB1, HMGB2, ID3, PTTG1, PKM, CRABP2, SOX2, PFN1
           GNG5, TTYH1, MDK, TPI1, ENO1, HES1, KIAA0101, UBE2C, RANBP1, CTGF
           CKS1B, CKS2, TOP2A, CLU, BIRC5, PHGDH, RPS27L, DDIT4, SERF2, HES4
PC_ 2
Positive:  MALAT1, MT-CO1, XIST, MT-CO3, MT-CO2, FOS, MT-ND4L, MT-ND3, DLX6-AS1, AC090498.1
           MT-ATP6, HSPA1A, SPP1, SCRG1, CCL3, C1orf61, FABP7, RPL17, GFAP, HEPN1
           DLX1, LGALS1, EGR1, DLX5, AIF1, RP11-436D23.1, BCAN, FOSB, HOPX, RP11-76I14.1
Negative:  RPS2, FTH1, MARCKSL1, C4orf48, ZNF428, PCSK1N, YBX1, CKB, YWHAH, RPL13
           YWHAQ, TERF2IP, CXADR, TUBB2A, RANBP1, PFN1, RAC1, RPS4Y1, VAMP2, CCNI
           HMGN1, CENPV, MIF, CRMP1, TBCB, CALM1, BNIP3, TTC3, TTC9B, LSM4
PC_ 3
Positive:  BCAN, SLC1A3, HEPN1, HOPX, PTPRZ1, GFAP, AQP4, FOS, MT2A, ATP1A2
           TNC, PTN, GATM, TFPI, SCRG1, C8orf4, EGR1, LGALS1, FAM107A, PEA15
           ZFP36L1, SLCO1C1, SOX9, GPX3, B2M, HES1, PON2, PMP2, CLU, FABP5
Negative:  MALAT1, MLLT11, MT-CO2, MT-CO1, C4orf48, GAP43, STMN2, PCSK1N, FTH1, TTR
           MAB21L1, SOX4, NEFL, RPS2, YBX1, TTC9B, CRABP1, SNCG, MAB21L2, LHX9
           IGFBP2, NHLH2, TAC1, SOX11, EBF1, LHX1, OCIAD2, LHX5-AS1, NRN1, SEC61G
PC_ 4
Positive:  BCAN, PCSK1N, HOPX, BNIP3, HEPN1, AQP4, FTH1, TNC, GFAP, MIF
           C4orf48, TTC9B, LGALS3, MT3, CPE, ATP1A2, DDIT4, TFPI, SLCO1C1, IGFBP2
           PEA15, SLC1A3, TERF2IP, PMP2, GATM, C8orf4, FAM107A, MT1E, ZFP36, IGFBP5
Negative:  TOP2A, BIRC5, CENPF, NUSAP1, UBE2C, CCNB2, AURKB, CCNB1, PBK, CDC20
           HIST1H4C, MAD2L1, SPC25, PTTG1, NUF2, CDK1, CKS1B, HMGB2, PLK1, KIAA0101
           CDKN3, CKS2, ASPM, TUBA1B, NUCKS1, TMSB15A, KPNA2, NCL, EEF1B2, CD24
PC_ 5
Positive:  S100A11, CRABP1, LHX5-AS1, AIF1, CRABP2, CCL3, SPP1, TYROBP, TPBG, CHCHD2
           NEAT1, RGS10, CCL3L3, APOC1, SERPINF1, SNCG, NEFM, NEFL, MT-ND1, WLS
           FCGRT, TTR, MAF, C1QC, APOE, CCL4, A2M, MT-ND5, MT-ND2, TRH
Negative:  C1orf61, NFIB, LINC01551, NEUROD6, EMX1, SFRP1, NFIA, MPPED2, CDK1, ID4
           IGFBP2, FOXG1, HEY1, DCK, PBK, TFAP2C, BCL11A, IFI44L, NUSAP1, NEUROD2
           SPC25, PHLDA1, CREB5, AMBN, NUF2, LMO3, TOP2A, BCL11B, CAPZA2, ASPM
"""            

#Correlation of PCA with mitocondrial percentage
unlist(lapply(1:ncol(cortical@reductions$pca@cell.embeddings),FUN=function(x) cor(cortical@reductions$pca@cell.embeddings[,x],cortical$percent.mt)))
""" [1]  0.321220450  0.589311081 -0.186444244 -0.040557701  0.094008965
 [6]  0.032307199 -0.075635766  0.039737821  0.025347671 -0.051287403
[11]  0.088106373  0.024331557 -0.228733358  0.101015557  0.170805763
[16] -0.030679608  0.089366212 -0.134942972 -0.169712539 -0.071244543
[21] -0.147319178  0.033716260  0.092703220  0.115814228  0.051682211
[26]  0.035515382  0.128533533  0.001860703 -0.087658975  0.086798654
[31] -0.025774687  0.076835677 -0.008687632 -0.039099727 -0.049125063
[36]  0.104619187  0.068383824 -0.018604314  0.030639035 -0.105359735
[41]  0.017081060  0.004682327 -0.032829498  0.020922410 -0.041354740
[46] -0.011782319  0.057609091  0.005762146  0.031202524 -0.009937891"""

#Sanity checks on PCA loadings.
#Plot feature PCA loadings for first handful of dimensions
plt<-VizDimLoadings(cortical, dims = 1:6, reduction = "pca")
ggsave(plt,file="rna_pcaloadings.pdf",width=30,height=20)
system("slack -F rna_pcaloadings.pdf ryan_todo")

#Plot feature heatmaps for first handful of dimensions
pdf("rna_pcaheatmap.pdf",height=100,width=50)
DimHeatmap(cortical, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
system("slack -F rna_pcaheatmap.pdf ryan_todo")

#Determine PC dimensions to use for clustering and 2/3D projection.
plt<-ElbowPlot(cortical)
ggsave(plt,file="rna_pcaelbow.pdf")
system("slack -F rna_pcaelbow.pdf ryan_todo")

#Going with 10 dimensions based on the elbow plot.
#Then clustering cells based on those 9 dimensions.
#Then performing UMAP for 3D projection using those 9 dimensions.

pca_dims<-10
cortical <- FindNeighbors(
    object = cortical,
    reduction = 'pca',
    dims = 1:pca_dims
    )
orgo_rna <- FindClusters(
    object = cortical,
    verbose = TRUE,
    resolution=0.5
    )

cortical <- RunUMAP(
    object = cortical,
    reduction = 'pca',
    dims = 1:pca_dims,
    n.components = 2
    )

plt1<-DimPlot(cortical, reduction = "umap")
plt2<-DimPlot(cortical, reduction = "umap",group.by="Cluster")
plt<-plt1+plt2
ggsave(plt,file="cortical.umap.png",width=30)
system("slack -F cortical.umap.png ryan_todo")

            
#Finding cluster markers for cell type identification
#Using seurats default find marker function.
#Loading dplyr for shorthand grouping of output.
library(dplyr)

#Running marker discovery
orgo_rna.markers <- FindAllMarkers(orgo_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
#Use dplyr for generation of marker lists and select top 10 genes (features) for plotting.
library(dplyr)
orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

top10 <- orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#Write out marker table
write.table(orgo_rna.markers,"orgo_rna.markers.txt",col.names=T,quote=F,sep="\t")

#Generate Heatmap of marker features
pdf("191222_rna_clustermarkers.pdf",width=30)
DoHeatmap(orgo_rna, features = top10$gene)
dev.off()

#Generate 2d umap projection plot
pdf("191222_rna_umap.pdf")
DimPlot(orgo_rna, reduction = "umap")
dev.off()

#Generate a pdf of stitched feature value plots for markers per cluster.
pdf("191222_rna_umapfeatureplots.pdf",width=100,height=100)
feat_genes<-orgo_rna.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
FeaturePlot(
object = orgo_rna,
features = feat_genes$gene,
pt.size = 0.1,
max.cutoff = 'q95',
ncol = 5
)
dev.off()
    
#Save the Seurat object for future cell type identification and coembedding with ATAC
saveRDS(orgo_rna,file="orgo_rna.Rds")
```

Further code, not yet implemented


```python

#Generating a 3D Plot via Plotly of the umap projection.
#Loading in additional libraries.
  library(plotly)
  library(htmlwidgets)
  library(RColorBrewer)

#Extracting UMAP dimensions from Seurat Object
  dims<-as.data.frame(orgo_rna[["umap"]]@cell.embeddings)
  dims$cellID<-row.names(dims)
  
#Extracting Meta Data from Seurat Object for coloring plot.
  annot<-as.data.frame(orgo_rna@meta.data)
  annot$cellID<-row.names(annot)
  
#Merging dimension data with annotations.
  dat<-merge(dims,annot,by="cellID",all.x=T)
  
  #Generate 3D Plot and standalone HTML widget
  p<-plot_ly(type="scatter3d")
  p<- p %>% add_trace(p,mode="markers",data=dat,x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color=~seurat_clusters,size=5,opacity=0.7,hovertext=~cellID)
  p<- p %>% add_trace(p,mode="markers",data=dat,x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color=~orig.ident,size=5,opacity=0.7,hovertext=~cellID)
  htmlwidgets::saveWidget(as_widget(p), "rna_umap.html",selfcontained=TRUE)

```

## Processing scATAC Data


```python

#metadata file placed into same folder
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")
dat<-read.table("E-MTAB-8089.sdrf.txt",header=T,sep="\t") #from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8089/files/

#unique(dat$Characteristics.organism.part.) #[1] "skin"                  "mixed embryo and skin"
# unique(dat$Characteristics.cell.line.)
#[1] "409B2"              "mixed 409B2 and H9"
# unique(dat$Characteristics.single.cell.well.quality.)
#[1] "multiple cells" "no cell"        "OK"             "2 cells"
#[5] "3 cells"        "debris"

#Filter to dat$Characteristics.single.cell.well.quality. == "OK"
dat<-dat[dat$Characteristics.single.cell.well.quality. == "OK",]

destdir="/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC/kanton_scATAC"
download_data<-function(x){
    out_name=unlist(lapply(strsplit(x,"/"),"[",9))
    dest_file=paste(destdir,out_name,sep="/")
    download.file(x,destfile=dest_file)
}
library(parallel)
mclapply(unique(dat$Comment.FASTQ_URI.), download_data,mc.cores=10)

#> length(unique(dat$Source.Name))
#[1] 525
#So 525 cells, all with PE sequencing (1050 files total)

cat ./kanton_scATAC/*1.fastq.gz > kanton.1.fastq.gz &
cat ./kanton_scATAC/*2.fastq.gz > kanton.2.fastq.gz &

#deleting files because they take up a ton of space (would take about an hour or so to redownload)

ref="/home/groups/oroaklab/refs/hg38/hg38.fa"

bwa mem -t 30 $ref kanton.1.fastq.gz kanton.2.fastq.gz | samtools view -@ 5 -b - > kanton.scATAC.bam

#modify bam file to follow scitools naming convention (readID#readnumber#0/1 example: TTGGAACTATCCTTCAACGCTCCAACTAATTCCGGT_3:77791856#0)
((samtools view -H kanton.scATAC.bam) && \
 (samtools view kanton.scATAC.bam | \
  awk 'OFS="\t" {split($1,a,"."); if(and($2,0x40)) $1=a[1]":"a[2]"#0"; else $1=a[1]":"a[2]"#1" ;print $0}')) | \
samtools view -b - > kanton.scATAC.sciformat.bam &

#scitools barcode based rmdup
scitools bam-rmdup kanton.scATAC.sciformat.bam &

#scitools make counts matrix based off our data
scitools atac-counts -O kanton_oroakpeaks kanton.scATAC.sciformat.bbrd.q10.bam ../../organoid_finalanalysis/orgo.500.bed &


#Tabix processing
input_bam="kanton.scATAC.sciformat.bbrd.q10.bam"
output_name="kanton"
tabix="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/tabix"
bgzip="/home/groups/oroaklab/src/cellranger-atac/cellranger-atac-1.1.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bgzip"
samtools view --threads 10 $input_bam | awk 'OFS="\t" {split($1,a,":"); print $3,$4,$8,a[1],1}' | sort -S 2G -T . --parallel=30 -k1,1 -k2,2n -k3,3n | $bgzip > $output_name.fragments.tsv.gz
$tabix -p bed $output_name.fragments.tsv.gz &
```


      File "<ipython-input-2-84d39b9df021>", line 1
        wget -r "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR345/*"   #download fastq files for 689 cells
                                                            ^
    SyntaxError: invalid syntax



## Generating Seurat Object For Kanton Data Set


```python
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")

# make counts matrix from sparse matrix
IN<-as.matrix(read.table("kanton_oroakpeaks.counts.sparseMatrix.values.gz"))
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table("kanton_oroakpeaks.counts.sparseMatrix.cols.gz")
colnames(IN)<-COLS$V1
ROWS<-read.table("kanton_oroakpeaks.counts.sparseMatrix.rows.gz")
row.names(IN)<-ROWS$V1

#Read in fragment path for coverage plots
kanton_fragment.path="./kanton.fragments.tsv.gz"

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

#Generate ChromatinAssay Objects
kanton_chromatinassay <- CreateChromatinAssay(
  counts = IN,
  genome="hg38",
  min.cells = 1,
  annotation=annotations,
  sep=c("_","_"),
  fragments=kanton_fragment.path
)

#Create Seurat Object
kanton_atac <- CreateSeuratObject(
  counts = kanton_chromatinassay,
  assay = "peaks",
)

#Meta.data updating
metadat<-read.table("E-MTAB-8089.sdrf.txt",header=T,sep="\t") #from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8089/files/
metadat<-metadat[grepl("1.fastq.gz",metadat$Comment.SUBMITTED_FILE_NAME.),] #remove read 2 duplicate meta data
metadat$cellID<-metadat$Comment.ENA_RUN.
kanton_atac$cellID<-row.names(kanton_atac@meta.data)
kanton_atac@meta.data<-merge(kanton_atac@meta.data,metadat,by="cellID")
row.names(kanton_atac@meta.data)<-kanton_atac@meta.data$cellID
#saving unprocessed SeuratObject
saveRDS(kanton_atac,file="kanton_SeuratObject.Rds")
```

## Performing cisTopic and UMAP On Kanton Data Set 


```python
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(Matrix)
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")

library(cisTopic)
kanton_atac<-readRDS(file="kanton_SeuratObject.Rds")
orgo_atac<-readRDS("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/orgo_SeuratObject.Rds")
atac_combined<- merge(kanton_atac, y = orgo_atac, add.cell.ids = c("kanton", "oroak"), project = "all_cells")

cistopic_processing<-function(seurat_input,prefix){
    cistopic_counts_frmt<-seurat_input$peaks@counts #grabbing counts matrices
    row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt)) #renaming row names to fit granges expectation of format
    atac_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt) #set up CisTopicObjects
    #Run warp LDA on objects
    atac_cistopic_models<-cisTopic::runWarpLDAModels(atac_cistopic,topic=c(5,10,20:30,40,50,55),nCores=15,addModels=FALSE)
    print("Saving cistopic models.")
    saveRDS(atac_cistopic_models,file=paste(prefix,"CisTopicObject.Rds",sep=".")) 
}
        
cistopic_processing(seurat_input=kanton_atac,prefix="kanton_atac")
cistopic_processing(seurat_input=atac_combined,prefix="all_combined")


kanton_atac_cistopic_models<-readRDS("kanton_atac.CisTopicObject.Rds")
atac_combined_cistopic_models<-readRDS("all_combined.CisTopicObject.Rds")



#Setting up topic count selection
pdf("kanton_atac_model_selection.pdf")
par(mfrow=c(1,3))
kanton_atac_cistopic_models <- selectModel(kanton_atac_cistopic_models, type='derivative')
dev.off()
system("slack -F kanton_atac_model_selection.pdf ryan_todo")

pdf("all_combined_model_selection.pdf")
par(mfrow=c(1,3))
kanton_atac_cistopic_models <- selectModel(atac_combined_cistopic_models, type='derivative')
dev.off()
system("slack -F all_combined_model_selection.pdf ryan_todo")

#set topics based on derivative
kanton_selected_topic=29
kanton_cisTopicObject<-cisTopic::selectModel(kanton_atac_cistopic_models,select=kanton_selected_topic,keepModels=T)
combined_selected_topic=21
combined_cisTopicObject<-cisTopic::selectModel(atac_combined_cistopic_models,select=combined_selected_topic,keepModels=T)

#saving model selected RDS
saveRDS(kanton_cisTopicObject,file="kanton_CisTopicObject.Rds")
saveRDS(combined_cisTopicObject,file="allcombined_CisTopicObject.Rds")
saveRDS(atac_combined,file="allcombined_SeuratObject.Rds")

####Function to include topics and umap in seurat object
cistopic_wrapper<-function(object_input=orgo_atac,cisTopicObject=orgo_cisTopicObject,resolution=0.8){   


    #run UMAP on topics
    topic_df<-as.data.frame(cisTopicObject@selected.model$document_expects)
    row.names(topic_df)<-paste0("Topic_",row.names(topic_df))
    dims<-as.data.frame(uwot::umap(t(topic_df),n_components=2))
    row.names(dims)<-colnames(topic_df)
    colnames(dims)<-c("x","y")
    dims$cellID<-row.names(dims)
    dims<-merge(dims,object_input@meta.data,by.x="cellID",by.y="row.names")

    #Add cell embeddings into seurat
    cell_embeddings<-as.data.frame(cisTopicObject@selected.model$document_expects)
    colnames(cell_embeddings)<-cisTopicObject@cell.names
    n_topics<-nrow(cell_embeddings)
    row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
    cell_embeddings<-as.data.frame(t(cell_embeddings))

    #Add feature loadings into seurat
    feature_loadings<-as.data.frame(cisTopicObject@selected.model$topics)
    row.names(feature_loadings)<-paste0("topic_",1:n_topics)
    feature_loadings<-as.data.frame(t(feature_loadings))

    #combined cistopic results (cistopic loadings and umap with seurat object)
    cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
    umap_dims<-as.data.frame(as.matrix(dims[2:3]))
    colnames(umap_dims)<-c("UMAP_1","UMAP_2")
    row.names(umap_dims)<-dims$cellID
    cistopic_umap<-CreateDimReducObject(embeddings=as.matrix(umap_dims),assay="peaks",key="UMAP_")
    object_input@reductions$cistopic<-cistopic_obj
    object_input@reductions$umap<-cistopic_umap

    n_topics<-ncol(Embeddings(object_input,reduction="cistopic"))

    object_input <- FindNeighbors(
      object = object_input,
      reduction = 'cistopic',
      dims = 1:n_topics
    )
    object_input <- FindClusters(
      object = object_input,
      verbose = TRUE,
      resolution=resolution
    )
}


kanton_atac<-cistopic_wrapper(object_input=kanton_atac,cisTopicObject=kanton_cisTopicObject,resolution=0.5)   
atac_combined<-cistopic_wrapper(object_input=atac_combined,cisTopicObject=combined_cisTopicObject,resolution=0.5)   

plt<-DimPlot(kanton_atac,group.by=c("Factor.Value.time.",
                                    "Characteristics.developmental.stage.",
                                    "Characteristics.cell.line.",
                                   "Characteristics.cell.type.",
                                   "Factor.Value.phenotype.",
                                   "seurat_clusters"),
                                    size=0.1)
ggsave(plt,file="kanton_atac.umap.png",width=20)
ggsave(plt,file="kanton_atac.umap.pdf",width=20)


plt<-DimPlot(atac_combined,group.by=c("Factor.Value.time.",
                                    "Characteristics.developmental.stage.",
                                    "Characteristics.cell.line.",
                                   "Characteristics.cell.type.",
                                   "Factor.Value.phenotype."),
                                    size=0.1)
ggsave(plt,file="atac_combined.umap.png",width=20)
ggsave(plt,file="atac_combined.umap.pdf",width=20)


i="kanton_atac.umap.png"
system(paste0("slack -F ",i," ryan_todo"))#post to ryan_todo

###save Seurat file
saveRDS(kanton_atac,file="kanton_SeuratObject.Rds")
saveRDS(atac_combined,file="allcombined_SeuratObject.Rds")
```


```python
#DA on sites by cell type
#Continued from above
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(ggplot2)
library(ggrepel)
library(dplyr)

kanton_atac<-readRDS(file="kanton_SeuratObject.Rds")

#Perform One vs. rest DA enrichment

write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

#set up an empty list for looping through
#orgo_da_peaks<-list()
da_peaks<-list()


#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group){
    da_peaks_tmp <- FindMarkers(
        object = obj,
        ident.1 = i,
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T
        )
    da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
    closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
    da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
    da_peaks_tmp$enriched_group<-c(i)
    da_peaks_tmp$compared_group<-c("all_other_cells")
    return(da_peaks_tmp)
  }

da_one_v_one<-function(i,obj,group,j_list){
    i<-as.character(i)
    da_tmp_2<-list()
    for (j in j_list){
        if ( i != j){
        da_peaks_tmp <- FindMarkers(
            object = obj,
            ident.1 = i,
            ident.2 = j,
            group.by = group,
            test.use = 'LR',
            latent.vars = 'nCount_peaks',
            only.pos=T
            )
        da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
        closest_genes <- ClosestFeature(obj,da_peaks_tmp$da_region)
        da_peaks_tmp<-cbind(da_peaks_tmp,closest_genes)
        da_peaks_tmp$enriched_group<-c(i)
        da_peaks_tmp$compared_group<-c(j)
        da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
        }
    }
    return(da_tmp_2)
  }

#Perform parallel application of DA test
library(parallel)
n.cores=length(unique(kanton_atac@meta.data$Factor.Value.time.))
da_peaks<-mclapply(
    unique(kanton_atac@meta.data$Factor.Value.time.),
    FUN=da_one_v_rest,
    obj=kanton_atac,
    group="Factor.Value.time.",
    mc.cores=n.cores)

#Merge the final data frame from the list for 1vrest DA
da_peaks<-do.call("rbind",da_peaks)

write("Outputting One v Rest DA Table.", stderr())
write.table(da_peaks,file="kanton.onevrest.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

#Plot out top peaks and associated gene name for each cluster
dat<-read.table("kanton.onevrest.da_peaks.txt",header=T,sep="\t")
dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% slice(1:2) #grabbing top 2 most significant peaks to label
plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(dat=dat_select,aes(label=gene_name,size=-distance),force=3)+theme_bw()
ggsave(plt,file="kanton_da_peaks.pdf")


#Empty list to rerun for 1v1 comparisons
da_peaks<-list()

n.cores=length(unique(kanton_atac@meta.data$Factor.Value.time.))
da_peaks<-mclapply(
    unique(kanton_atac@meta.data$Factor.Value.time.),
    FUN=da_one_v_one,
    obj=kanton_atac,
    group="Factor.Value.time.",
    j_list=do.call("as.character",list(unique(kanton_atac@meta.data$Factor.Value.time.))),
    mc.cores=n.cores)

#Merge the final data frame from the list for 1v1 DA
da_peaks<-do.call("rbind",do.call("rbind",da_peaks))


write("Outputting One v One DA Table.", stderr())
write.table(da_peaks,file="kanton.onevone.da_peaks.txt",sep="\t",col.names=T,row.names=T,quote=F)

```


```python
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

#Read in data and modify to monocle CDS file
#read in RDS file.

setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")
kanton_atac<-readRDS(file="kanton_SeuratObject.Rds")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
motif.matrix <- CreateMotifMatrix(
  features = granges(kanton_atac),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm)

# Add the Motif object to the assays and run ChromVar
###Full Datset###
kanton_atac <- SetAssayData(
  object = kanton_atac,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)
kanton_atac <- RegionStats(object = kanton_atac, genome = BSgenome.Hsapiens.UCSC.hg38)
kanton_atac<- RunChromVAR( object = kanton_atac,genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(kanton_atac,file="kanton_SeuratObject.Rds")
```


```python
###Differential TF Accessibility by cluster###
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
set.seed(1234)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")
kanton_atac<-readRDS(file="kanton_SeuratObject.Rds")


#Perform One vs. rest DA enrichment

write("Performing one vs. rest DA enrichment per annotation grouping supplied.", stderr())

DefaultAssay(kanton_atac) <- 'chromvar'

#set up an empty list for looping through
kanton_tf<-list()

#define DA functions for parallelization
#Use LR test for atac data
da_one_v_rest<-function(i,obj,group){
    da_peaks_tmp <- FindMarkers(
        object = obj,
        ident.1 = as.character(i),
        group.by = group,
        test.use = 'LR',
        latent.vars = 'nCount_peaks',
        only.pos=T
        )
    da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
    da_peaks_tmp$enriched_group<-c(as.character(i))
    da_peaks_tmp$compared_group<-c("all_other_cells")
    return(da_peaks_tmp)
  }

da_one_v_one<-function(i,obj,group,j_list){
    i<-as.character(i)
    da_tmp_2<-list()
    for (j in j_list){
        if ( i != j){
        da_peaks_tmp <- FindMarkers(
            object = obj,
            ident.1 = i,
            ident.2 = j,
            group.by = group,
            test.use = 'LR',
            latent.vars = 'nCount_peaks',
            only.pos=T
            )
        da_peaks_tmp$da_region<-row.names(da_peaks_tmp)
        da_peaks_tmp$enriched_group<-c(i)
        da_peaks_tmp$compared_group<-c(j)
        da_tmp_2[[paste(i,j)]]<-da_peaks_tmp
        }
    }
    return(da_tmp_2)
  }

#Perform parallel application of DA test
library(parallel)
n.cores<-length(unique(kanton_atac@meta.data$Factor.Value.time.))
kanton_tf<-mclapply(
    sapply(unique(kanton_atac@meta.data$Factor.Value.time.),as.character),
    FUN=da_one_v_rest,
    obj=kanton_atac,
    group="Factor.Value.time.",
    mc.cores=n.cores)

#Merge the final data frame from the list for 1vrest DA
kanton_tf<-do.call("rbind",kanton_tf)

write("Outputting One v Rest DA Table.", stderr())
write.table(kanton_tf,file="kanton.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

dat<-read.table("kanton.onevrest.da_tf.txt",header=T,sep="\t")
#To convert JASPAR ID TO TF NAME
dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
dat$enriched_group<-as.character(dat$enriched_group)
write.table(dat,file="kanton.onevrest.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)
dat_select<-dat %>% arrange(rev(desc(p_val_adj))) %>% group_by(enriched_group) %>% dplyr::slice(1:2) #grabbing top 2 most significant peaks to label
plt<-ggplot(dat,aes(x=avg_logFC,y=(-log(p_val)),color=as.factor(enriched_group)))+geom_point(aes(alpha=0.1))+geom_label_repel(dat=dat_select,aes(label=da_tf),force=3)+theme_bw()
ggsave(plt,file="kanton_oncevrest.da_tf.pdf")
system("slack -F kanton_oncevrest.da_tf.pdf ryan_todo")
#Empty list to rerun for 1v1 comparisons
kanton_tf<-list()
    
n.cores=length(unique(kanton_atac@meta.data$Factor.Value.time.))
kanton_tf<-mclapply(
    unique(kanton_atac@meta.data$Factor.Value.time.),
    FUN=da_one_v_one,
    obj=kanton_atac,
    group="Factor.Value.time.",
    j_list=do.call("as.character",list(unique(kanton_atac@meta.data$Factor.Value.time.))),
    mc.cores=n.cores)

#Merge the final data frame from the list for 1v1 DA
kanton_tf<-do.call("rbind",do.call("rbind",kanton_tf))

write("Outputting One v One DA Table.", stderr())
write.table(kanton_tf,file="kanton.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

dat<-read.table("kanton.onevone.da_tf.txt",header=T,sep="\t")
#To convert JASPAR ID TO TF NAME
dat$da_tf <- unlist(lapply(unlist(lapply(dat$da_region, function(x) getMatrixByID(JASPAR2020,ID=x))),function(y) name(y)))
write.table(dat,file="kanton.onevone.da_tf.txt",sep="\t",col.names=T,row.names=T,quote=F)

```


```python
#################################################################
#Coembedding of ATAC GA and RNA
#################################################################
#https://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html://satijalab.org/seurat/v3.1/atacseq_integration_vignette.html
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
set.seed(1234)
# Load the pre-processed scRNA-seq and scATAC-seq data

#Public RNA
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data")
puborgo<-readRDS("PublicOrganoid.SeuratObject.rds")
pubprimary<-readRDS("PublicPrimary.SeuratObject.rds")

#Public ATAC #needs to go through gene activity analysis
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/Public_Data/scATAC")
kanton_atac<-readRDS(file="kanton_SeuratObject.Rds")

#Our RNA
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis/rna_processing")
orgo_rna<-readRDS("orgo_rna.SeuratObject.rds")

#Our ATAC
setwd("/home/groups/oroaklab/adey_lab/projects/BRAINS_Oroak_Collab/organoid_finalanalysis")
orgo_atac<-readRDS(file="orgo_atac.SeuratObject.GeneActivity.RDS")

transfer.anchors <- FindTransferAnchors(
  reference = orgo_rna,
  query = orgo_atac,
  reduction = 'cca',
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = orgo_rna$RNA_snn_res.0.5,
  weight.reduction = orgo_atac[['lsi']],
)

orgo_atac <- AddMetaData(object = orgo_atac, metadata = predicted.labels)
DefaultAssay(orgo_atac) <- 'peaks'

plot1 <- DimPlot(
  object = orgo_rna,
  group.by = 'RNA_snn_res.0.5',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = orgo_atac,
  group.by = 'predicted.id',
  reduction="umap",
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("integrated.data.pdf")
CombinePlots(list(plot1,plot2), ncol = 2)
dev.off()

```
--->