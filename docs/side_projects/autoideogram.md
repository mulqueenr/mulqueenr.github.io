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

<!--
Get to scriptin! This is mainly to practice my python scripting.
```python
import wget
import os
import pandas as pd
from cairosvg import svg2png
import time
import seaborn as sns
import skimage
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from skimage import color

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
	chr_list=[*(range(1,23,1))]+["X","Y"]
	pal = sns.color_palette('Spectral', len(chr_list))
	chr_color = dict(zip(chr_list, pal))
	return(chr_color)

#Banded chromosomes from: https://commons.wikimedia.org #check what band resolution those are at
def download_chr_svgs():
	""" Function to download chromosome public domain SVGs from wikimedia commons"""
	chr_list=[*(range(1,23,1))]+["X","Y"]
	for i in chr_list:
		out_name="chr"+str(i)+".svg"
		print("Checking for: "+out_name+"\n")
		if out_name not in os.listdir():
			image_url=chr_loc_dict_subdir[i]
			null=wget.download(image_url,out=out_name)
			time.sleep(10)#to not overload the server with requests
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


def color_chroms(chr_color,chr_length,crop_text=False):
	"""Function to color chromosomes by chr_color and output pngs
	Input: chr_color is dict of RGB color codes to chromosome (keys)
	chr_length is dict of bp length to chromosomes (keys)
	crop_text boolean if to crop the band info or not"""
	max_length=chr_length['1']
	for key in chr_color:
		i=key
		image_name="chr"+str(i)+".svg"
		svg2png(url=image_name, write_to="chr"+str(i)+".png",scale=2.0) #all images are sized the same so cropping is consistent
		image=skimage.io.imread("chr"+str(i)+".png",as_gray=True)
		if crop_text:
			image2=image[13:2370,172:275]#crop to chromosomes
			image2=skimage.transform.resize(image2,(image2.shape[0],image2.shape[1] // (chr_length[str(i)]/max_length)), anti_aliasing=True)
		else:
			image2=image
		col_multiplier=list(chr_color[key]) #take color for chromosome
		col_multiplier=col_multiplier
		image2=skimage.color.label2rgb(image2!=1, image=image2, colors=[col_multiplier], alpha=0.4, image_alpha=0.9, bg_color=[1,1,1],kind='overlay',saturation=0.5)
		skimage.io.imsave(arr=image2,fname="chr"+str(i)+".png")#save chr



download_geneloc_to_band()
iscn=pd.read_table("ideogram_9606_GCF_000001305.14_850_V1") #read in iscn band data
max_length=iscn[iscn['bp_stop']==iscn.groupby('#chromosome')['bp_stop'].transform('max')] #get chromosome length for resizing png files to appropriate size
chr_length=dict(zip(max_length['#chromosome'],max_length['bp_stop']))
download_chr_svgs()
chr_color=set_chr_colors()
color_chroms(chr_color,chr_length,crop_text=True)

```
--->

```python
import pandas as pd
import wget
import seaborn as sns
import os
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import colorsys
os.chdir('/Users/rmulqueen/autoideo')

"""Color set up"""
def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def lighten_bands(band_color,lighten=False):
	"""Function to apply staining saturation to chromosome defined colors or maintain saturated value for outline"""
	for i in range(0,len(band_color)):
		band_tmp=list(band_color[i])
		if lighten==True:
			band_tmp[1]=1-(band_tmp[1] * iscn.loc[iscn.index[i],'density']/100)
		band_tmp_rgb=colorsys.hls_to_rgb(band_tmp[0],band_tmp[1],band_tmp[2])
		band_tmp_hex=rgb_to_hex(int(band_tmp_rgb[0]*255),int(band_tmp_rgb[1]*255),int(band_tmp_rgb[2]*255))
		band_color[i]=band_tmp_hex
	return(band_color)

def set_chr_colors(iscn,pallet="husl"):
	""" Function to set colors per chromosome """
	chr_list=iscn['#chromosome'].unique()
	chr_list=[str(i) for i in chr_list]
	pal = sns.color_palette(pallet, len(chr_list))
	chr_color = dict(zip(chr_list, [colorsys.rgb_to_hls(i[0],i[1],i[2]) for i in pal]))
	band_color=[chr_color[i] for i in list(iscn["#chromosome"])] #add band color based on chr_color
	iscn["chr_color"]=lighten_bands(band_color,lighten=False)
	band_color=[chr_color[i] for i in list(iscn["#chromosome"])] #add band color based on band
	iscn["band_color"]=lighten_bands(band_color,lighten=True)
	return(iscn)

"""Input reference data"""
#Band to genomic loci data set from NCBI: ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/ideogram_9606_GCF_000001305.14_850_V1
def download_geneloc_to_band():
	"""Function to download public domain NCBI chromosome band info"""
	table_ftp="ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/ideogram_9606_GCF_000001305.14_850_V1"
	table_name=table_ftp.split("/")[-1]
	if table_name not in os.listdir():
		wget.download(table_ftp)
	else:
		print(table_name+" Found! Proceeding...")

"""Karyotype Building Functions"""
def get_sub_chr(chrom,arm,band,subband):
	"""Function to return subchr """
	subchr=iscn[iscn["#chromosome"]==chrom].copy() #subset by chr
	if arm == "q":
		p_arm=subchr[subchr["arm"]=="p"].copy() #grab all the p arm
		q_arm=subchr[subchr["arm"]=="q"].copy()
		q_arm=q_arm[(q_arm["large_band"].astype(int)<int(band)) | ((q_arm["large_band"].astype(int)==int(band)) & (q_arm["sub_band"].astype(int)<=int(subband)))].copy()
		merged_subchr=pd.concat([p_arm,q_arm])
	else:
		p_arm=subchr[subchr["arm"]=="p"].copy() #grab all the p arm
		p_arm=p_arm[(p_arm["large_band"].astype(int)<int(band)) | ((p_arm["large_band"].astype(int)==int(band)) & (p_arm["sub_band"].astype(int)<=int(subband)))].copy()
		merged_subchr=p_arm
	return(merged_subchr)
	####Think this through to grab the proper sections

def get_translocation_chr(N):
	"""Function to generate translocation chr, of format t(9;22)(q34.1;q11.2)"""
	out=pd.DataFrame()
	chr_in=N.split("t(")[1].split(")")[0].split(";")
	chr_in_boundaries=N.split("(")[2].split(")")[0].split(";")
	chr_in_boundary_arms=[i[0] for i in chr_in_boundaries] #isolate chr arms
	chr_in_boundary_bands=[i[1:3] for i in chr_in_boundaries] #isolate chr bands
	chr_in_boundary_subbands=[i[3] for i in chr_in_boundaries] #isolate chr subbands
	chr_in_boundary_subbands=[str(i).ljust(3,'0') for i in chr_in_boundary_subbands] #add 0 to expand precision if not included
	for n in range(0,len(chr_in)):
		chrom=chr_in[n]
		arm=chr_in_boundary_arms[n]
		band=chr_in_boundary_bands[n]
		subband=chr_in_boundary_subbands[n]
		out=pd.concat([out,get_sub_chr(chrom,arm,band,subband)])
	return(out)

#def get_derivative_chr(N):
#"""Function to generate translocation chr, of format +der(22)t(9;22)(q34.1;q11.2)"""

#def get_i_chr isocentric mirror chr

def set_up_chr(N):
	"""Return reference data on chrN"""
	if N in iscn['#chromosome'].unique():
		out=iscn[iscn['#chromosome'] == N].copy()
	elif N.startswith("t"):
		out=get_translocation_chr(N)
	return(out)

def set_up_karyotype(kary):
	i=0
	out_karyotype=pd.DataFrame()
	for N in karyotype:
		uniq_N=str(N)+"_"+str(i)
		out_N=set_up_chr(str(N))
		out_N['uniq_chr']=[uniq_N]*out_N.shape[0]
		i+=1
		out_karyotype=pd.concat([out_karyotype,out_N])
	return(out_karyotype)

### Set up reference Data ###
download_geneloc_to_band() #download geneloc table
iscn=pd.read_table("ideogram_9606_GCF_000001305.14_850_V1") #read in iscn band data
#ISCN description of bands includes stain information (including density) meaning we can scale chromosome bands on greyscale by density and then recolor
#Convert stain +/- and density to greyscale color
iscn.loc[iscn['density'].isnull(),['density']] = 0 #set gneg rows to 0 density
iscn['size']=iscn['iscn_stop']-iscn['iscn_start']
iscn['large_band']=[str(i).split(".")[0] for i in iscn['band']]
iscn['sub_band']=[str(str(i).split(".")[1].ljust(3,'0')) for i in iscn['band']] #add 0 to expand precision if not included
iscn["band_name"]=iscn["#chromosome"]+iscn["arm"]+iscn["band"].astype('string') #set band name for hover display

### Set band colors based on user pallet ###
iscn=set_chr_colors(iscn,pallet="Spectral") #supply chr color palette you want

### Build Karyotype ###
#Input is standard ISCN karyotype
#Example from cydas 47<2n>,XY,-7,+8,t(9;22)(q34.1;q11.2),i(17)(q10),+der(22)t(9;22)(q34.1;q11.2)
#Convert shorthand to long format
karyotype=["t(9;22)(q341;q112)",1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,'X','Y']
karyotype=[str(i) for i in karyotype]
kary=set_up_karyotype(kary=karyotype)

kary_labels=dict(zip(kary['uniq_chr'].unique(),karyotype))

###

fig = px.bar(kary, x='uniq_chr', y='size', color="band_name",color_discrete_sequence=list(kary["band_color"]),hover_name="band_name",hover_data=['#chromosome','bp_start','bp_stop'])
fig=fig.update_yaxes(showgrid=False, zeroline=False, autorange="reversed")
fig=fig.update_xaxes(ticktext=[kary_labels[i] for i in kary_labels],
		tickvals=[i for i in kary_labels])
fig=fig.update_traces(marker_line_color=list(kary["chr_color"]),marker_line_width=1.5, opacity=0.4)
fig=fig.update_layout(showlegend=False)
fig.show()


#for separation of chr, might have to make chr groups and make subplots per chr group.



###Insert space between translocations, insert i and der and read karytype into long format. also reciprocal t, fix karyo order
```