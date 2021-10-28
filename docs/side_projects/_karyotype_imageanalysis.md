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

```

## Testing our ability to cut out chromosomes

Basing of this to start: https://www.quora.com/How-can-I-detect-an-object-from-static-image-and-crop-it-from-the-image-using-openCV
Starting with easy-peasy pre-"grammed" karyogram.

```python
import cv2 
import numpy as np
image_in="/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/1280px-NHGRI_human_male_karyotype.png"

#reading the image  
image = cv2.imread(image_in) 
paper = cv2.resize(image, (500,500))

cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"resize" + '.png', paper) 

ret, thresh_gray = cv2.threshold(cv2.cvtColor(paper, cv2.COLOR_BGR2GRAY), 245, 255, cv2.THRESH_BINARY)
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"graythreshold" + '.png', thresh_gray) 

canny = cv2.Canny(thresh_gray, 20, 300) 
contours, hier = cv2.findContours(canny, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

canny = paper
canny = cv2.drawContours(canny, contours, -1, (0,255,0), 1)
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"contoured" + '.png', canny ) 

 # Erase small contours, and contours which small aspect ratio (close to a square)
for c in contours:
     area = cv2.contourArea(c)
     if area < 50: # Fill very small contours with zero (erase small contours).
         cv2.fillPoly(thresh_gray, pts=[c], color=0)
         continue

canny=paper
canny=cv2.drawContours(canny, contours, -1, (255,0,0), 1)
cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/" +"contoured_filter" + '.png', canny) 


idx = 0 
for c in contours: 
    x,y,w,h = cv2.boundingRect(c) 
    area = cv2.contourArea(c)
    if area > 5: # Fill very small contours with zero (erase small contours).
        idx+=1 
        new_img=paper[y:y+h,x:x+w] 
        cv2.imwrite("/mnt/c/Documents and Settings/mulqueen/Documents/karyotype/chr/" +str(idx) + '.png', new_img) 


cv2.waitKey(0) 


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