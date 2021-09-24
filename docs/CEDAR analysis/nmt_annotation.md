---
title: NMTseq Annotations
layout: cedar_analysis
author: Ryan Mulqueen
permalink: /nmt_annot/
category: CEDAR
---

## Set up genomic annotations

Generate bed files for annotation regions to sum methylation data over.

Gene location information stored here:

```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes
```

Generation of bed files:
```bash
awk 'OFS="\t" {split($10,a,"\;"); split(a[1],b,"\""); if($3=="gene") print $1,$4,$5,b[2]}' genes.gtf > genes.bed
```

Regulatory information stored here:
```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds
```

Some of these features have redundancies (i.e. same start and end sites but different names). Have to filter that out.

Generation of bed files:
```bash
#Regulatory information from ensembl regulatory build
wget http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz

zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="enhancer") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > enhancers.bed
zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="promoter") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > promoters.bed
zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz | awk 'OFS="\t" {if ($3=="CTCF_binding_site") split($9,a,"=");split(a[2],b,";"); if (b[1] != "") print "chr"$1,$4,$5,b[1]}' | sort -u -k1,1 -k2,2n -k3,3n > ctcf.bed

awk 'OFS="\t" {print $1,1,$2}' /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt | bedtools makewindows -b - -w 100000 | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}'> 100kb.bed
awk 'OFS="\t" {print $1,1,$2}' /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt | bedtools makewindows -b - -w 10000 | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}'> 10kb.bed

``` 

Also given a list of enhancer regions that are breast cancer specific from Hisham
````bash
#Raw file stored here:
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed

#Converted to same bed format with metadata bin names
grep -v "_" /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/Enhancer_hg38_SI_RI.bed | awk 'OFS="\t" {print $1,$2,$3,"bin_"NR}' > /home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed

````

ENCODE ChIP Peaks
```bash
cd /home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata
#FOXA1 
#MCF7 https://www.encodeproject.org/files/ENCFF112JVK/
wget https://www.encodeproject.org/files/ENCFF112JVK/@@download/ENCFF112JVK.bed.gz
#T47D https://www.encodeproject.org/files/ENCFF758GJL/
wget https://www.encodeproject.org/files/ENCFF758GJL/@@download/ENCFF758GJL.bed.gz

zcat ENCFF112JVK.bed.gz | awk 'OFS="\t" {print $1,$2,$3,"MCF7_"NR}' > FOXA1.bed
zcat ENCFF758GJL.bed.gz | awk 'OFS="\t" {print $1,$2,$3,"T47D_"NR}' >> FOXA1.bed
```
List of genomic annotation bed files:
```bash
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/promoters.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/enhancers.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/ctcf.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/genes/genes.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/100kb.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/10kb.bed
/home/groups/CEDAR/mulqueen/ref/refdata-gex-GRCh38-2020-A/regulatory_beds/breastcancer_enhancers_SI_RI.bed
/home/groups/CEDAR/mulqueen/ref/public_cellline_chipdata/FOXA1.bed
```
