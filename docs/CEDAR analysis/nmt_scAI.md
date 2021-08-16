---
title: scAI on NMTseq
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /scai_nmt/
category: CEDAR
---

Going to try to use the recently published scAI for multiomic integration (https://github.com/sqjin/scAI). 

## Installation
```R
devtools::install_github("sqjin/scAI")
```

## Following a walkthrough with our data
https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_simulation.html


## Load data
```R
library(scAI)
library(dplyr)
library(cowplot)
library(ggplot2)
```
## Create a scAI Object

```R
load("/home/groups/CEDAR/mulqueen/projects/nmt/nmt_test/data_mESC.rda")
X <- data_mESC$data # List of data matrix
labels <- data_mESC$labels # the collected time of cells, which is used for validation
```

## Preprocess data
```R
scAI_outs <- preprocessing(scAI_outs, assay = NULL, minFeatures = 200, minCells = 1,
                           libararyflag = F, logNormalize = F)
scAI_outs <- addpData(scAI_outs, pdata = labels, pdata.name = "Conditions")
```

## Run scAI model
```R

```