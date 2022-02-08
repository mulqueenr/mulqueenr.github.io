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
    top: 80vh;
    left: 4vw;
  }
}

.background_img {
  min-height: 100vh;
  background-image: url("{{site.baseurl}}/assets/images/cortex_scatterplot.png");
  color: white;
  display: flex;
  isolation: isolate;
  width: 70vw;

  /* Create the parallax scrolling effect */
  background-attachment: fixed;
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;

}


.background_img2 {
  min-height: 120vh;
  background-image: url("{{site.baseurl}}/assets/images/s3human_scatterplot.png");
  color: white;
  display: flex;
  isolation: isolate;
  width: 100vw;

  /* Create the parallax scrolling effect */
  background-attachment: fixed;
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;

}

.hero-intro {
  width: 100vw;
  padding: 2em;
}

</style>

<div class="background_img">
  <div class="hero-intro">
    <h1><code style="font-weight: 1000; font-family: 'PT Sans'" >3D Scatter Plots</code></h1>
  </div>
</div>

## Introduction

Here's how to generate 3D scatterplots of single cell experiments in Blender. To do this, we are going to use a similar shader to the one I use in the ["molecular render"](https://mulqueenr.github.io/molecular_render/) project. Additionally, I will be using blender's python console to read in a simple text file to populate a UMAP projection. 

For the text file, generate a header, tab-separated table with the following columns:

Cluster | cellID | X | Y | Z | cluster_color_hexvalue |
:--:|:--:|:--: | :--:|:--:|:--: | 
String, of cluster grouping | String, unique cell identifier| float of x coordinate | float of y coordinate | float of z coordinate | Hex value of fill color | 

Here's an example lines used for the render at the top of the page.

```bash
Astrocytes    TAGGTCCGACGTACTAGGGCCTCGGTCTATGGCCTA    4.24424248742567    -1.74691044949975    -6.48374510684418    #1C7D54
Astrocytes    ATTCAGAAGCATCGCGCAGCCAGACTCTATGGCCTA    3.60301401455387    -1.96493138894082    -6.47136162049336    #1C7D54
Astrocytes    TCAACGAGTTCGCGATGGTCAGAGCCCGCCGATATC    5.51775913941571    -1.87741656898663    -6.76243310557264    #1C7D54
```

Required software
- Blender v 2.93 or higher ([link](https://www.blender.org/))

## Description of the script and how to run

I set up a script to automate the full process. I'll list what is happening within the script so it is a bit easier to tweak. The blender-python (bpy) API is kind of in a state right now because there is a lot of active development. Sorry it isn't as readable as a good python script should be! To run the script, just open blender, change a window panel to "Python Console" or just press Shift+F4. Then run. The script makes a shader similar to what we use in the molecular render side project. It does this in python but the end result (when made nice to look at) is this:

<img src="{{site.baseurl}}/assets/images/scatterplot_cellshader.png">

This script uses linked objects (basically a way to save memory). To do this, it makes a "master" data point for each cluster, hides it, and then duplicates it for each data point. That means each data point isn't a unique object, but a copy with a unique location. Saves on a ton of time and memory.
In the future, the geometry or particle nodes might make this even easier, but this is what we have to work with now.
I've run this script on up to 200k cells. It seems to slow down the more data-blocks there are in the blender file. From my computer, I run the plotting in 1k cell batches. I find that for each 1K cells I plot, the next 1K cells take about 0.2-0.4 more seconds. This adds up with bigger files, but hey, at least it is linear. For my trash computer, 200k cells took about 2-3 hours.

Change these lines:
```python
#Input
file_xyz=open("C:/Users/Chad/Desktop/cortex.3D.table","r") #change the filepath to your tsv will the line format described above.
#Output
bpy.ops.wm.save_as_mainfile(filepath="C:/Users/Chad/Desktop/cortex.blend") #change filepath to the output blender file you would like
bpy.context.scene.render.filepath = 'cortex.test.png' #change file path to the output png you would like.

```

Here is the script as used to generate the hero image at the top of the page.

## The script

```python
#1. import modules
import bpy
import math
import time
import bmesh

#Read in file and store it in memory (this doesn't take up much memory)
file_xyz=open("C:/Users/mulqueen/Desktop/hg38_3d.umap.tsv","r") #change path to whatever filepath you want. I got my computer refurbished and it was named Chad. I swear it wasn't me.
tabraw=file_xyz.readlines()[1:]
data_count=len(tabraw)
file_xyz.close()

#initialize an object, a sphere, for our data points.
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.05,segments=64, ring_count=32) #higher segments and ring_counts will make a smoother sphere, but I dont think its necessary
obj=bpy.context.active_object #select the sphere we just made

#set up a master shader material
mat = bpy.data.materials.new(name='mymat')
mat.use_nodes = True #use node trees, these can be seen by switching a panel to the shader editor if you want. It will look like the above shader, just not nicely placed.
mat_nodes = mat.node_tree.nodes
mat_links = mat.node_tree.links
mat = bpy.data.materials['mymat'] #Get the material you want 
node_to_delete =  mat.node_tree.nodes['Principled BSDF'] #Get the node in its node tree (replace the name below)
mat.node_tree.nodes.remove( node_to_delete ) #Remove it
#add all the nodes, using col_node as variable of each node as it is being made. then using that to modify default value fields
col_node=mat_nodes.new('ShaderNodeRGB')
col_node=mat_nodes.new('ShaderNodeFresnel')
bpy.data.materials["mymat"].node_tree.nodes['Fresnel'].inputs[0].default_value = 1.33

col_node=mat_nodes.new('ShaderNodeHueSaturation')
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[0].default_value = 1
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[1].default_value = 0.7
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[2].default_value = 2
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[3].default_value = 0

col_node=mat_nodes.new('ShaderNodeMath')
bpy.data.materials["mymat"].node_tree.nodes["Math"].operation = 'MULTIPLY'

col_node=mat_nodes.new('ShaderNodeBsdfRefraction')
bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].inputs[1].default_value = 1

col_node=mat_nodes.new('ShaderNodeBsdfGlossy')
bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].inputs[1].default_value = 1

col_node=mat_nodes.new('ShaderNodeHueSaturation')
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[0].default_value = 1
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[1].default_value = 0.4
bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[2].default_value = 2

col_node=mat_nodes.new('ShaderNodeMixShader')
col_node=mat_nodes.new('ShaderNodeVolumeAbsorption')
bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].inputs[1].default_value = 0.3

col_node=mat_nodes.new('ShaderNodeBsdfTranslucent')
col_node=mat_nodes.new('ShaderNodeLightPath')
col_node=mat_nodes.new('ShaderNodeMixShader')

#build node tree links (going from left most inputs)
#sorry this is a monstrosity
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value"].inputs[4])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].inputs[4])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['RGB'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Fresnel'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Math"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Hue Saturation Value'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].inputs[0])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes['Hue Saturation Value'].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Math"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[0])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Refraction BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Glossy BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].inputs[2])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Hue Saturation Value.001"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Translucent BSDF"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Volume Absorption"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Material Output"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Translucent BSDF"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[2])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Mix Shader"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[1])
mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Light Path"].outputs[1], bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].inputs[0])

mat_links.new(bpy.data.materials["mymat"].node_tree.nodes["Mix Shader.001"].outputs[0], bpy.data.materials["mymat"].node_tree.nodes["Material Output"].inputs[0])


#set up render engine and scene
bpy.context.scene.render.engine="CYCLES" #set render engine to CYCLES
bpy.data.scenes["Scene"].cycles.denoiser="NLM" #set denoiser for render
bpy.data.scenes["Scene"].cycles.samples=512 #this is a whole lotta sampling
bpy.context.scene.render.image_settings.color_depth = '16' #more color channels!
bpy.context.scene.render.resolution_x = 3840 #up the resolution
bpy.context.scene.render.resolution_y = 2160
bpy.data.objects["Sphere"].hide_render = True # hide sphere in render
bpy.data.objects["Sphere"].hide_viewport=True
bpy.data.lights["Light"].energy = 100000 # increase light wattage
bpy.data.lights["Light"].shadow_soft_size= 1
bpy.data.objects["Light"].location=(5,-5,10) #location and rotation i deteremined manually and just set up here for convenience

#set up stage by cutting up the default cube vertices and smoothing it
obj_cube=bpy.data.objects["Cube"]
obj_cube.scale=(30,30,30) #scale up the cube
#this is to cut out a vertex to make an open box
bpy.context.view_layer.objects.active = obj_cube
bpy.ops.object.mode_set(mode='EDIT')
bpy.ops.mesh.select_mode(type="VERT")  # Switch to edge select mode
bm = bmesh.from_edit_mesh(obj_cube.data)  # Create bmesh object for easy mesh evaluation
bm.verts.ensure_lookup_table()
bm.verts.remove(bm.verts[2]) # Write the mesh back
bmesh.update_edit_mesh(obj_cube.data)  # Update the mesh in edit mode
bpy.ops.object.mode_set(mode='OBJECT') #switch back to object mode when done
bpy.ops.object.modifier_add(type='SUBSURF') #make it smooth
bpy.data.objects["Cube"].modifiers["Subdivision"].render_levels=6
bpy.data.objects["Cube"].location=(-4,4.3,17.725) #change the location for more dramatic shadows

#move the camera and rotate
bpy.data.objects["Camera"].location=(34.61997604370117, -40.53969955444336, 25.66326904296875)
bpy.data.objects["Camera"].rotation_euler=(1.1093189716339111, 0.0, 0.8149281740188599)

#finally ready to start reading in our data
scene=bpy.context.scene

#set up a material per hex color, name as annotation
#this is looping through the file, grabbing the unique clusters and there color codes, then making a dictionary for look up later
start = time.time()
annot={}
for line in tabraw[1:]:
  line=line.replace('\n','')
  l=line.split('\t')
  if l[0] not in annot:
    hexcode=l[5].lstrip("#")
    rgb=[int(hexcode[i:i+2], 16) for i in (0, 2, 4)]
    r=float(rgb[0])/255 #color of spheres, blender uses 0-1 scale
    g=float(rgb[1])/255
    b=float(rgb[2])/255
    clust=str(l[0])
    annot[clust]=[r,g,b]

end = time.time()
print(end - start)

#make a custom material shader for each annotation (just changing color)
#this copies the material shader we set up earlier, and then changes the input color
master_mat=source_mat = bpy.data.materials["mymat"]
for i in annot.keys():
  copied_mat = master_mat.copy()
  copied_mat.name=i
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[0]=annot[i][0]
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[1]=annot[i][1]
  bpy.data.materials[i].node_tree.nodes["RGB"].outputs[0].default_value[2]=annot[i][2]


#make a custom collection for each annotation. this makes a "master sphere" to link for each cluster also
for i in annot.keys():
  collection = bpy.data.collections.new(i) #make new collection
  bpy.context.scene.collection.children.link(collection) #link new collection
  mat = bpy.data.materials.get(i) #set material properties of collection
  name=str(i)+"_master" #make name of master sphere
  new_obj = bpy.data.objects.new(name, scene.objects.get("Sphere").data) #make a new copy
  new_obj.data = scene.objects.get("Sphere").data.copy()
  bpy.data.collections[i].objects.link(new_obj) #link new object to collection
  new_obj.data.materials.append(mat) #add material
  bpy.data.objects[name].hide_render = True # hide masters
  bpy.data.objects[name].hide_viewport=True

#make a dictionary look up for copying master spheres
master_sphere={}
for i in annot.keys():
  master_sphere[i]=scene.objects.get(i+"_master").data

#define a nice function to copy data points and link them to the master spheres. also places the copies into nice cluster named collections for easier navigation.
def add_data_point(input_dat):
    line=input_dat
    line=line.replace('\n','')
    l=line.split('\t')
    #print(line)
    x=float(l[2]) #location of spheres
    y=float(l[3])
    z=float(l[4])
    name=str(l[1])
    clust=str(l[0])
    my_new_obj = bpy.data.objects.new(name,master_sphere[clust])
    my_new_obj.location = (x,y,z)       
    my_new_obj.hide_viewport=False
    my_new_obj.hide_render=False
    bpy.data.collections[clust].objects.link(my_new_obj)

n=1000
in_list = [tabraw[i * n:(i + 1) * n] for i in range((len(tabraw) + n + 1) // n )] 
for in_dat_list in in_list:
  start = time.time()
  out=[add_data_point(in_dat) for in_dat in in_dat_list] 
  end = time.time()
  print(end - start)


#some last minute tweaks, here are some convenience functions if you want to change things. I also encourage you to play around with lighting and camera positioning to get some interesting views of your data.
#to adjust size of points
for clust in annot.keys():
  for i in bpy.data.collections[clust].objects:
    i.scale=(0.5,0.5,0.5)


#to adjust alpha value and translucence of material properties
for clust in annot.keys():
  bpy.data.materials[clust].node_tree.nodes["RGB"].outputs[0].default_value[3] = 0.3
  bpy.data.materials[clust].node_tree.nodes["Volume Absorption"].inputs[1].default_value = 0.1


bpy.ops.wm.save_as_mainfile(filepath="C:/Users/mulqueen/Desktop/cortex.blend") #save blender file

bpy.context.scene.render.filepath = 'cortex.test.png'
bpy.ops.render.render(write_still=True) #render and save file

```

## Bonus
Here I took the human s3-ATAC brain data and colored by cell type. I then cranked up the emission strength and removed all external lights.


<div class="background_img2">
</div>