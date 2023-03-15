---
title: Auto-ideogram project
layout: side_projects
author: Ryan Mulqueen
permalink: /auto_ideogram/
category: side_projects
---

## Workshopping some names.

```bash
#Cytogenetic
#Ideogram
#Karyotype
#Structural
#Variant
```

Download wikimedia chromosomes, download NCBI genome assembly table with bands, apply locations to pixel heights (by percentage),
cut and paste chromosomes around according to karyotype, ooo maybe change the color too

# send a HTTP request to the server and save
# the HTTP response in a response object called r

Set up environment
```bash
mkdir ~/autoideo
cd ~/autoideo
conda create -n autoideo
conda activate autoideo
pip3 install wget svgpathtools pandas cairosvg

```

brew install imagemagick -i;
./configure --disable-osx-universal-binary --prefix=/opt/homebrew/Cellar/imagemagick/7.1.1-3 --disable-silent-rules --with-x11
make install

Get to scriptin! This is mainly to practice my python scripting.
```python
import wget
import svgpathtools as svg
import os
import pandas as pd
import cairosvg 
import time
import seaborn as sns


chr_loc_dict_subdir={
1:"https://upload.wikimedia.org/wikipedia/commons/0/07/Human_chromosome_1_ideogram_vertical.svg",
2:"https://upload.wikimedia.org/wikipedia/commons/1/11/Human_chromosome_2_ideogram_vertical.svg",
3:"https://upload.wikimedia.org/wikipedia/commons/5/58/Human_chromosome_3_ideogram_vertical.svg",
4:"https://upload.wikimedia.org/wikipedia/commons/1/11/Human_chromosome_4_ideogram_vertical.svg",
5:"https://upload.wikimedia.org/wikipedia/commons/4/47/Human_chromosome_5_ideogram_vertical.svg",
6:"https://upload.wikimedia.org/wikipedia/commons/2/2c/Human_chromosome_6_ideogram_vertical.svg",
7:"https://upload.wikimedia.org/wikipedia/commons/7/7f/Human_chromosome_7_ideogram_vertical.svg",
8:"https://upload.wikimedia.org/wikipedia/commons/7/7d/Human_chromosome_8_ideogram_vertical.svg",
9:"https://upload.wikimedia.org/wikipedia/commons/0/01/Human_chromosome_9_ideogram_vertical.svg",
10:"https://upload.wikimedia.org/wikipedia/commons/8/89/Human_chromosome_10_ideogram_vertical.svg",
11:"https://upload.wikimedia.org/wikipedia/commons/d/dd/Human_chromosome_11_ideogram_vertical.svg",
12:'https://upload.wikimedia.org/wikipedia/commons/d/d2/Human_chromosome_12_ideogram_vertical.svg',
13:"https://upload.wikimedia.org/wikipedia/commons/5/5c/Human_chromosome_13_ideogram_vertical.svg",
14:"https://upload.wikimedia.org/wikipedia/commons/f/f1/Human_chromosome_14_ideogram_vertical.svg",
15:"https://upload.wikimedia.org/wikipedia/commons/0/0f/Human_chromosome_15_ideogram_vertical.svg",
16:"https://upload.wikimedia.org/wikipedia/commons/9/96/Human_chromosome_16_ideogram_vertical.svg",
17:"https://upload.wikimedia.org/wikipedia/commons/b/bf/Human_chromosome_17_ideogram_vertical.svg",
18:"https://upload.wikimedia.org/wikipedia/commons/8/86/Human_chromosome_18_ideogram_vertical.svg",
19:"https://upload.wikimedia.org/wikipedia/commons/d/df/Human_chromosome_19_ideogram_vertical.svg",
20:"https://upload.wikimedia.org/wikipedia/commons/4/48/Human_chromosome_20_ideogram_vertical.svg",
21:"https://upload.wikimedia.org/wikipedia/commons/2/2c/Human_chromosome_21_ideogram_vertical.svg",
22:"https://upload.wikimedia.org/wikipedia/commons/f/fa/Human_chromosome_22_ideogram_vertical.svg",
"X":"https://upload.wikimedia.org/wikipedia/commons/d/d2/Human_chromosome_X_ideogram_vertical.svg",
"Y":"https://upload.wikimedia.org/wikipedia/commons/b/b8/Human_chromosome_Y_ideogram_vertical.svg"}


def set_chr_colors():
	""" Function to set colors per chromosome """
	chr_list=[*(range(1,22,1))]+["X","Y"]
	pal = sns.color_palette('Spectral', len(chr_list))
	pal=pal.as_hex()
	chr_color = dict(zip(chr_list, pal))
	return(chr_color)

#Banded chromosomes from: https://commons.wikimedia.org #check what band resolution those are at
def download_chr_svgs():
	""" Function to download chromosome public domain SVGs from wikimedia commons"""
	chr_list=[*(range(1,22,1))]+["X","Y"]
	for i in chr_list:
		out_name="chr"+str(i)+".svg"
		print("Checking for: "+out_name+"\n")
		if out_name not in os.listdir():
			image_url=chr_loc_dict_subdir[i]
			null=wget.download(image_url,out=out_name)
			time.sleep(15)
		else:
			print(out_name+" Found! Proceeding..."+"\n")

#Band to genomic loci data set from NCBI: ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/ideogram_9606_GCF_000001305.14_850_V1
def download_geneloc_to_band():
	"""Function to download public domain NCBI chromosome band info"""
	table_ftp="ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/ideogram_9606_GCF_000001305.14_850_V1"
	table_name=table_ftp.split("/")[-1]
	if table_name not in os.listdir():
		wget.download(table_ftp)
	else:
		print(table_name+" Found! Proceeding...")

#Recolor attribute function (to make the recolor_chr_paths function more readable)
def recolor_attribute(i,attribute):
	chr_col=chr_color[i]
	fill='fill:'+chr_col
	stroke="stroke:"+chr_col
	attribute_style=attribute['style'].split(";")
	attribute_style=[fill if i.startswith('fill:') else i for i in attribute_style]
	attribute_style=[fill if i.startswith('stroke:') else i for i in attribute_style]
	attribute['style']=";".join(attribute_style)
	return(attribute)

#Supply i as chromosome to define color
def recolor_chr(i,attributes):
	attributes=[recolor_attribute(i,attribute) for attribute in attributes]
	return(attributes)

def translate_chr(i,paths,x,y):
	[path.translated(eval(x)+eval(y)j) for path in paths]


pdf = cairosvg.svg2pdf(image_name)

download_chr_svgs()
download_geneloc_to_band()
chr_color=set_chr_colors()

for i in [*(range(1,22,1))]+["X","Y"]:
	image_name="chr"+str(i)+".svg"
	paths, attributes, svg_attributes = svg.svg2paths2(image_name)
	attributes=recolor_chr(i,attributes)

#def cut_chr():
#this opens illustrator with the svg
svg.disvg(paths[300:600])
i=1
image_name="chr"+str(i)+".svg"
test_out="out.svg"
paths, attributes, svg_attributes = svg.svg2paths2(image_name)
svg.wsvg(paths, attributes=attributes, svg_attributes=svg_attributes, filename='output1.svg')


#maybe switch to pdfs or something. svgs are being annoying
