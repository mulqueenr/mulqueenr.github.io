---
title: 3D Scatter Plots
layout: side_projects
author: Ryan Mulqueen
permalink: /3d_scatterplots/
category: side_projects
---


<style>

h1 {
  text-transform: uppercase;
  font-size: clamp(5rem, 5vw + 0.5rem, 5rem);
  line-height: 0.9;
  margin: 0;
  color: antiquewhite;
  mix-blend-mode: hard-light;
  position:absolute; left:2rem; top: 4em;

  code {
    display: block;
    width: max-content;
    background: white;
    font-size: 0.5em;
    color: #355f08;
    padding: 0.1em;
    border-radius: 0.125em;
    margin-bottom: 0.05em;
    font-family: 'PT Sans';
    font-weight: 1000;
  }
}

.background_img {
  min-height: 75vh;
  background-image: url("{{site.baseurl}}/assets/images/OPC_scatterplot.png");
  color: white;
  display: flex;
  isolation: isolate;
  width: 100%;

  /* Create the parallax scrolling effect */
  background-attachment: fixed;
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;

}

.hero-intro {
  width: 100%;
  padding: 2em;
}

</style>

<div class="background_img">
  <div class="hero-intro">
    <h1><code style="font-weight: 1000; font-family: 'PT Sans'" >3D Scatter Plots</code></h1>
  </div>
</div>

### Introduction

Here's how to generate 3D scatterplots of single cell experiments in Blender. To do this, we are going to use a similar shader to the one I use in the [molecular render]((https://mulqueenr.github.io/3d_scatterplots/)) project. Additionally, I will be using blender's python console to read in a simple text file to populate a UMAP projection. 

For the text file, generate a header, tab-separated table with the following columns:

Cell_ID | X | Y | Z | cluster | cluster_color_hexvalue |
:--:|:--:|:--: | :--:|:--:|:--: | 
String, unique cell identifier| float of x coordinate | float of y coordinate | float of z coordinate | String, of cluster grouping | Hex value of fill color | 

Here's an example lines used for the render at the top of the page.

```bash
CellID  UMAP_X  UMAP_Y  UMAP_Z  Annot Hex
AACGCAATAAGAAGTCTAACGTCCTGGAACTTGCGG  -1.844343639  -2.362912544  -0.114864404  0 #EA8C55
AACGCAATAAGAAGTCTAGTCAGTACGTTCTCTCCT  -0.460297069  -1.17707301 1.235511487 0 #EA8C55
AACGCAATAAGAAGTCTAGTTAATTCCTGCATTACG  1.809149885 1.921362272 1.354163831 0 #EA8C55
AACGCAATAAGAAGTCTATAAGATTATCTCCTCCTG  0.518666887 0.555381766 -0.726500417  0 #EA8C55
AACGCAATAAGAAGTCTATCGTTGCGCGTTCTATCA  -1.022839314  0.342407963 -0.949501479  0 #EA8C55
```

Note: the annot column is not necessary for the render, and can be filled with a dummy variable. In a later render, I'll show a way to highlight select groups if you do want to use the annot column. 

Required software
- Blender v 2.93 or higher ([link](https://www.blender.org/))

Here I generated a gallery of GIFs showing the step by step process from start to finish to generate the render at the top of this page. 

<blockquote class="imgur-embed-pub" lang="en" data-id="a/S4WyQl3"  ><a href="//imgur.com/a/S4WyQl3">3D Scatterplots in blender</a></blockquote><script async src="//s.imgur.com/min/embed.js" charset="utf-8"></script>

Gif numbers (1-16) are listed as separate steps below with more information. If you want to pause them, rewind, or zoom in, you can click on the widget and it should bring you to the Imgur gallery. 

## Generate a cell shader
### 1. Open Blender and use the default cube to generate a cell shader.

- Select the cube and then go to the "Shader Editor" by clicking the clock looking button on the bottom left and expand the window.

- Add nodes by using "Shift+A" then use the search function to find the nodes by their name.
  - Note that "Multiply" node is actually called the "Math" node.

- Below is my final node set up if you want to copy the same material shader I use. We will be changing the input RGB value for each cell to color it based on the respective hex code in the table.

- Once you are happy with the shader, click on the name (i.e. "Material") and rename it to "mymaterial". This is important since the script uses "mymaterial" as the name for copying. Alternatively, you can just change the script to the material name of your choosing.

<img src="{{site.baseurl}}/assets/images/scatterplot_cellshader.png">

## Set up the scene

### 2. Check the shader configuration

- In the top right corner of the viewport, click the shaded looking ball icon on the furthest right to begin a real-time rendered view. 
- Change the rendering engine by clicking the "Render Properties" tab on the right panel. 
  - Change the render engine from "Eevee" to "Cycles"
  - In "Sampling" select "Denoising" and make sure the "Render NLM" box is checked.

- Now select the light and move it closer to the cube to see the lighting effects. To reposition the light use "G", use "X" "Y" or "Z" to move it along just one axis at a time. With the light selected, you can go to the Object Data Properties tab on the right (it looks like a green lightbulb). Here you can adjust the Wattage of your light bulb, the color that it emits (I keep with white) and the radius of its effect. 

### 3. Make the stage and hide the cube
Now to make a backdrop for our render.
- Hide the cube from the viewport and the render by turning off the eye and camera icons in the "Outliner" tab.
- On the top left of the viewport window select "Add" > "Mesh" > "Plane"
- Scale the plane to a large enough size that it will take up the whole camera angle.
- With the plane selected, go to "Edit" mode and select the three verteces and use "E" to extrude them further.
- Drag the extruded vertexes back and upwards.
- Select the modifier tab for the plane (the blue wrench) and add a "Subdivision Surface" modifier. 
- Set the render levels to 6 to really smooth it for our final render.
- Finally change the viewport back to object mode. This is just to make the python script run more quickly as the computer won't try to render the cells as they come in.

## Run the python script
This script is copied and edited from [Mattheu Gibert](http://www.gibert.biz/downloads/3dscatterplotswithblender). In short, it reads in the table line by line, parses out the X,Y, and Z coordinates and the hexcode to determine the color. For each line it generates a 3D sphere. It converts the hex into RGB for the RGB input we included in the shader. 

This line (13) in the script should be modified to point correctly to your own table.

```python
file_xyz=open("C:/Users/mulqueen/Downloads/3D.OPC.table","r")
```

Some other things of note:
- The loop skips the first line of the table, assuming it is the headers
- Variable "n" can be changed to modify the size of the spheres. 

### 4. Use the python console to run the script
- Make a new collection to hold all the cell objects you are about to generate. To do this, go to the outliner panel, right click, and press "New Collection"
- Then go to the python console window, by changing from the Shader editor at the bottom panel.
- Copy and paste this code below to populate blender with your 3D scatter plot.
  - I recommend first testing this with the first 100 or so cells, by limiting the "for" loop from [1:100]. This will give you a change to adjust the lighting, move the camera, adjust the stage height and generally set up the render, before clearing out the collection and repopulating with all the cells.

```python
import bpy
import math

#define a function to create a sphere at each cell location
def newSph(diam):
    tempSph = bpy.ops.mesh.primitive_uv_sphere_add(radius=diam,segments=32, ring_count=16) #higher segments and ring_counts will make a smoother sphere, but I dont think its necessary
    return tempSph

# Program begins here.
scn = bpy.context.scene
#removeObjects( scn ) 

file_xyz=open("C:/Users/mulqueen/Downloads/3D.OPC.tsv","r")
tabraw=file_xyz.readlines()
file_xyz.close()

#I read in points 1:1000 at a time, but that is because my computer is bad
for ligne in tabraw[1:1000]:
    ligne=ligne.replace('\n','')
    l=ligne.split('\t')
    print(ligne)
    x=float(l[1]) #location of spheres
    y=float(l[2])
    z=float(l[3])
    hexcode=l[5].lstrip("#")
    rgb=[int(hexcode[i:i+2], 16) for i in (0, 2, 4)]
    r=float(rgb[0])/255 #color of spheres, blender uses 0-1 scale
    g=float(rgb[1])/255
    b=float(rgb[2])/255
    name=str(l[0])
    n=0.1 #size of spheres
    newSph(n) #make new sphere
    ob = bpy.context.active_object
    ob.name = name
    ob.location=(x,y,z)
    me = ob.data
    mat_temp = bpy.data.materials.new(str(l[0])+".mat")
    mat_temp=bpy.data.materials["mymaterial"].copy()
    ob.data.materials.append(mat_temp) #add material we made, called mymaterial
    ob.active_material.node_tree.nodes["RGB"].outputs[0].default_value=(r,g,b,0.8) #change the color of the RGB input


```


## Final rendering
Finally lets snap a pic!
- I turned the Resolution X and Y under the "Output Properties" tab up to 3840 px X 2160 px for a 4K resolution render.
- I also turned the render sampling in the "Render Properties" tab up to 256 for Render
- Finally press "F12" key to render and watch a beautiful plot be born.


<img src="{{site.baseurl}}/assets/images/OPC_scatterplot.png">


## Bonus: Highlight a selected group of cells

Since our cells are translucent spheres, we can highlight the cells of a specific annotation by generating a small point light in each sphere.
To do this, we will use pretty much the same python "for" loop, but instance lights at the cells with annotation "1".

```python
import bpy
import math

#define a function to create a sphere at each cell location
def newSph(diam):
    tempSph = bpy.ops.mesh.primitive_uv_sphere_add(radius=diam,segments=32, ring_count=16) #higher segments and ring_counts will make a smoother sphere, but I dont think its necessary
    return tempSph

def newLight(name_in):
    light_data = bpy.data.lights.new(name=name_in, type='POINT') #make a new point light
    light_data = bpy.data.lights.new(name=name_in, type='POINT')
    lamp_object = bpy.data.objects.new(name=name_in, object_data=lamp_data)
    lamp_object.data.energy = 100 #set power to 100 watts
    return lamp_object

# Program begins here.
scn = bpy.context.scene
#removeObjects( scn ) 

file_xyz=open("C:/Users/mulqueen/Downloads/3D.OPC.tsv","r")
tabraw=file_xyz.readlines()
file_xyz.close()

#I read in points 1:1000 at a time, but that is because my computer is bad
for ligne in tabraw[1:]:
  ligne=ligne.replace('\n','')
  l=ligne.split('\t')
  if ligne[4]=="1": #select only cells in annotation 1
    print(ligne)
    x=float(l[1]) #location of spheres
    y=float(l[2])
    z=float(l[3])
    name=str(l[0])+".light"
    newLight(name_in) #make a light
    ob = bpy.context.active_object
    ob.location=(x,y,z) # Place lamp to a specified location
    ob.select = True # And finally select it make active
    scn.objects.active = lamp_object

```

