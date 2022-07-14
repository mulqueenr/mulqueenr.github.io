---
title: Autowall Battlemaps
layout: side_projects
author: Ryan Mulqueen
permalink: /battlemaps_auto/
category: side_projects
---



https://pytorch.org/get-started/locally/

## Install pytorch

```python
conda create -n battlemap
conda activate battlemap
conda install opencv pytorch torchvision torchaudio cpuonly -c pytorch
conda install scipy
```


## Testing our ability to cut out chromosomes

Basing of this to start: https://www.quora.com/How-can-I-detect-an-object-from-static-image-and-crop-it-from-the-image-using-openCV
Starting with easy-peasy pre-"grammed" karyogram.

```python
import cv2 
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
image_in="/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/images.jfif"

#reading the image  
image = cv2.imread(image_in) 
paper = cv2.resize(image, (1000,1000))

cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/" +"resize" + '.png', paper) 
ret, thresh_gray = cv2.threshold(cv2.cvtColor(paper, cv2.COLOR_BGR2GRAY), 245, 255, cv2.THRESH_OTSU) #otsu thresholding
#ret, thresh_gray = cv2.threshold(cv2.cvtColor(paper, cv2.COLOR_BGR2GRAY), 247, 255, cv2.THRESH_BINARY) #binary thresholding
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/" +"graythreshold" + '.png', thresh_gray) 

canny = cv2.Canny(thresh_gray, 20, 300) 
contours, hier = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

canny = paper
canny = cv2.drawContours(canny, contours, -1, (0,255,0), 1)
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/" +"contoured" + '.png', canny ) 


idx = 0 
for c in contours: 
    x,y,w,h = cv2.boundingRect(c) 
    area = cv2.contourArea(c)
    if area > 20: # Fill very small contours with zero (erase small contours).
        idx+=1 
        new_img=paper[y:y+h+int(0.25*h),x:x+w+int(0.25*w)] #adding 25% height and width to each output
        cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/chr/" +str(idx) + '.png', new_img) 

chrom_in = "/mnt/c/Documents and Settings/mulqueen/Documents/Dnd/maps/chr/10.png"
image = cv2.imread(chrom_in) 

Tcsr = minimum_spanning_tree(image.data)
Tcsr.toarray().astype(int)
```


#Trying depth detection for fun?

conda create -n monodepth2 python=3.6.6 anaconda
conda activate monodepth2
conda install pytorch=0.4.1 torchvision=0.2.1 -c pytorch
pip install tensorboardX==1.4
conda install opencv=3.3.1   # just needed for evaluation

https://github.com/nianticlabs/monodepth2

```python

```