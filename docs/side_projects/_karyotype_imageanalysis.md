---
title: Molecular renders
layout: side_projects
author: Ryan Mulqueen
permalink: /molecular_render/
category: side_projects
---



https://pytorch.org/get-started/locally/

## Install pytorch

```python
conda create -n kary
conda activate kary
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
image_in="/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/1280px-NHGRI_human_male_karyotype.png"

#reading the image  
image = cv2.imread(image_in) 

ret, thresh_gray = cv2.threshold(cv2.cvtColor(image, cv2.COLOR_BGR2GRAY), 245, 255, cv2.THRESH_OTSU) #otsu thresholding

cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"graythreshold" + '.png', thresh_gray) 

canny = cv2.Canny(thresh_gray, 20, 300) 
contours, hier = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

canny = cv2.drawContours(canny, contours, -1, (0,255,0), 1)
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"contoured" + '.png', canny ) 

idx = 0 
for c in contours: 
    x,y,w,h = cv2.boundingRect(c) 
    area = cv2.contourArea(c)
    if area > 20: # Fill very small contours with zero (erase small contours).
        idx+=1 
        new_img=image[y:y+h+int(0.25*h),x:x+w+int(0.25*w)] #adding 25% height and width to each output
        cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/chr/" +str(idx) + '.png', new_img) 

chrom_in = "/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/chr/10.png"
image = cv2.imread(chrom_in) 

Tcsr = minimum_spanning_tree(image.data)
Tcsr.toarray().astype(int)
```

This works pretty well but it looks like we get more countours than we need (some that are within p or q arms). Going to try to merge nearby ones. This will probably not work on a real case scenario because chromosomes will overlap quite a bit on a regular squash.

```python
import cv2

image = cv2.imread("example.png")

gray = cvw.bgr2gray(image)

thresh = cvw.threshold_otsu(gray, inverse=True)

# dilation
img_dilation = cvw.dilate(thresh, 5)

# Find contours
contours = cvw.find_external_contours(img_dilation)
# Map contours to bounding rectangles, using bounding_rect property
rects = map(lambda c: c.bounding_rect, contours)
# Sort rects by top-left x (rect.x == rect.tl.x)
sorted_rects = sorted(rects, key=lambda r: r.x)

# Distance threshold
dt = 5

# List of final, joined rectangles
final_rects = [sorted_rects[0]]

for rect in sorted_rects[1:]:
    prev_rect = final_rects[-1]

    # Shift rectangle `dt` back, to find out if they overlap
    shifted_rect = cvw.Rect(rect.tl.x - dt, rect.tl.y, rect.width, rect.height)
    intersection = cvw.rect_intersection(prev_rect, shifted_rect)
    if intersection is not None:
        # Join the two rectangles
        min_y = min((prev_rect.tl.y, rect.tl.y))
        max_y = max((prev_rect.bl.y, rect.bl.y))
        max_x = max((prev_rect.br.x, rect.br.x))
        width = max_x - prev_rect.tl.x
        height = max_y - min_y
        new_rect = cvw.Rect(prev_rect.tl.x, min_y, width, height)
        # Add new rectangle to final list, making it the new prev_rect
        # in the next iteration
        final_rects[-1] = new_rect
    else:
        # If no intersection, add the box
        final_rects.append(rect)

for rect in sorted_rects:
    cvw.rectangle(image, rect, cvw.Color.MAGENTA, line_style=cvw.LineStyle.DASHED)

for rect in final_rects:
    cvw.rectangle(image, rect, cvw.Color.GREEN, thickness=2)

cv2.imwrite("result.png", image)
```


#image segmentation using kmeans?
https://www.youtube.com/watch?v=6CqRnx6Ic48

```python
import cv2 
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
image_in="/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/1280px-NHGRI_human_male_karyotype.png"

#reading the image  
img = cv2.imread(image_in) 
img2 = img.reshape((-1,3))
img2 = np.float32(img2)

criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10,1.0)

#Clusters
k=2
attempts=10
ret,label,center=cv2.kmeans(img2,k,None,criteria,attempts,cv2.KMEANS_PP_CENTERS)

center= np.uint8(center)
res = center[label.flatten()]
res2 = res.reshape((img.shape))
cv2.imwrite(image_in[:-3]+"seg.png",res2)
res2 = cv2.cvtColor(res2, cv2.COLOR_BGR2GRAY)

contours, hierarchy = cv2.findContours(res2, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

idx = 0 
for c in contours: 
    x,y,w,h = cv2.boundingRect(c) 
    area = cv2.contourArea(c)
    if area > 20: # Fill very small contours with zero (erase small contours).
        idx+=1 
        new_img=img[y:y+h+int(0.25*h),x:x+w+int(0.25*w)] #adding 25% height and width to each output
        cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/chr/" +str(idx) + '.png', new_img) 


```


https://github.com/suryavb95/Image-Segmentation-spanning-trees/blob/master/segmentation.py

```bash
pip install networkx matplotlib Pillow
```

```python

from networkx import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image
from PIL import ImageFilter
import sys


# ###################################################### #
# Definitions
# ###################################################### #

image = "/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/1280px-NHGRI_human_male_karyotype.png"

k = 3000
gaussianfilter = True
showeachstep = True
intDiff = 1
G = Graph()
channel = 0
nshape = 's' # shape of nodes when drawn: 'o' circles, 's' squares
nsize = 2500 # divided by image width


# ###################################################### #
# Algorithm functions
# ###################################################### #


def Int(C):
    if len(C) == 1:
        return intDiff
    T = minimum_spanning_tree(G.subgraph(C))
    maxEdge = max(T.edges(), key=lambda e:T[e[0]][e[1]]['weight'])
    return T[maxEdge[0]][maxEdge[1]]['weight']


def Tau(C):
    return k/len(C)


def MInt(C1, C2):
    return min(Int(C1) + Tau(C1), Int(C2) + Tau(C2))


def Dif(C1, C2):
    min = float('inf')
    for v1 in C1:
        if v1 not in G:
            continue
        edges1 = G[v1]
        for v2 in C2:
            if v2 not in edges1:
                continue
            if edges1[v2]['weight'] < min:
                min = edges1[v2]['weight']
    return min


def D(C1, C2):
    return Dif(C1,C2) > MInt(C1,C2)


# ###################################################### #
# Helping functions
# ###################################################### #


def getComponentsIndices(S, e):
    indices = []
    for i,c in enumerate(S):
        if e[0] in c:
            indices.append(i)
        if e[1] in c:
            indices.append(i)
    return indices


def drawComponents(S, mode, extra_edge=None):
    node_colors = [0] * len(G.nodes())
    edge_colors = []
    edges = []
    for i,C in enumerate(S):
        if mode == 'random':
            color = i
        else:
            color = (min((Int(C) + Tau(C)), 255))   # max e MST weight
        for node in C:
            node_colors[node] = color
        sub_g = G.subgraph(C)
        sub_g_edges = sub_g.edges()
        if sub_g_edges != []:
            edges += sub_g_edges
            edge_colors += [color]*len(sub_g_edges)
    map = 'Paired'
    max = len(S)
    shape = nshape
    if mode != 'random':
        map = 'RdYlGn'
        max = 255
        shape = 'o'
    draw_networkx_nodes(G, positions, node_color=node_colors, cmap=plt.get_cmap(map), vmin=0 , vmax=max, node_shape=shape, node_size=nsize)
    draw_networkx_edges(G, positions, edgelist=edges, width=2, edge_color=edge_colors, edge_cmap=plt.get_cmap(map), edge_vmin=0 , edge_vmax=max)
    if extra_edge != None:
        draw_networkx_edges(G, positions, edgelist=[extra_edge], width=10, edge_color=[G[extra_edge[0]][extra_edge[1]]['weight']], edge_cmap=plt.get_cmap(map), edge_vmin=0 , edge_vmax=max)


# ###################################################### #
# Graph initialization from image
# ###################################################### #


print("Creating model...")

img = Image.open(image)
img = img.filter(ImageFilter.GaussianBlur(radius=0.8))
width, height = img.size
nsize /= width
pixels = img.load()


# create nodes
for i in range(0, width):
    for j in range(0, height):
        G.add_node(i+j*width)

# create edges
for i in range(1, width-1):
    for j in range(1, height):
        n = i+j*width
        w = abs(pixels[i,j][channel] - pixels[i-1,j][channel])
        G.add_edge(n, n-1, weight=w)
        w = abs(pixels[i,j][channel] - pixels[i-1,j-1][channel])
        G.add_edge(n, n-width-1, weight=w)
        w = abs(pixels[i,j][channel] - pixels[i,j-1][channel])
        G.add_edge(n, n-width, weight=w)
        w = abs(pixels[i,j][channel] - pixels[i+1,j-1][channel])
        G.add_edge(n, n-width+1, weight=w)

# create edges for border pixels
# top
for i in range(1, width):
    w = abs(pixels[i,0][channel] - pixels[i-1,0][channel])
    G.add_edge(i, i-1, weight=w)
# left
for j in range(1, height):
    w = abs(pixels[0,j][channel] - pixels[0,j-1][channel])
    G.add_edge(j*width, (j-1)*width, weight=w)
    w = abs(pixels[0,j][channel] - pixels[1,j-1][channel])
    G.add_edge(j*width, (j-1)*width+1, weight=w)
# right
for j in range(1, height):
    w = abs(pixels[width-1,j][channel] - pixels[width-1,j-1][channel])
    G.add_edge(j*width+(width-1), (j-1)*width+(width-1), weight=w)
    w = abs(pixels[width-1,j][channel] - pixels[width-2,j-1][channel])
    G.add_edge(j*width+(width-1), (j-1)*width+(width-2), weight=w)
    w = abs(pixels[width-1,j][channel] - pixels[width-2,j][channel])
    G.add_edge(j*width+(width-1), j*width+(width-2), weight=w)


# ###################################################### #
# Graph drawing
# ###################################################### #


print("Drawing model...")

plt.imshow(img, interpolation='nearest')
positions = dict()
for i in range(0,width):
    for j in range(0,height):
        positions[i+j*width] = [i,j]

colors = []
for e in G.edges():
    colors.append(G[e[0]][e[1]]['weight'])

draw_networkx_nodes(G, positions, node_size=nsize)
draw_networkx_edges(G, positions, width=2, edge_color=colors, edge_cmap=plt.get_cmap('RdYlGn'), edge_vmin=0 , edge_vmax=255)
#plt.show()
plt.imsave(image[:-3]+"seg.png",img)


# ###################################################### #
# Algorithm
# ###################################################### #


print("Running segmentation algorithm...")

S = []
for v in G:
    S.append([v])



i = 1
n = 0
for e in sorted(G.edges(data=True), key = lambda (a,b,att) :att['weight']):
    print(str(i)+"/"+str(len(G.edges())))
    i += 1
    c1,c2 = getComponentsIndices(S,e)
    if n > 0:
        n -= 1
    if (showeachstep and c1 != c2 and n==0) or (showeachstep==None and i==len(G.edges())):
        plt.close('all')
        plt.imshow(img, interpolation='nearest')
        drawComponents(S, 'weights', e)
        plt.show(block=False)
        inp = str(raw_input())
        if inp == 'c':
            showeachstep = None
        try:
            n = int(inp)
        except:
            pass
    if c1 == c2 or e[2]['weight'] > MInt(S[c1],S[c2]):
        continue
    else:
        S[c1] += S[c2]
        S.remove(S[c2])


# ###################################################### #
# Components drawing
# ###################################################### #


print("Drawing "+str(len(S))+" segments...")
plt.close('all')
plt.imshow(img, interpolation='nearest')
drawComponents(S, 'random')
plt.show()
plt.imsave(image[:-3]+"seg2.png",img)

```