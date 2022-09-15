---
title: DNA Builder
layout: side_projects
author: Ryan Mulqueen
permalink: /dna_builder/
category: side_projects
---

## Introduction

Note: This is a big ol' WIP and I'm currently cleaning it up and adding more features. Going to release it as a git repo to work as a blender addon. (Also going to break up the file, so its not absolutely illegible.)

```python
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

####TO DO###########
#Fix protein input
#Adjust hybridization and protein import data to more correctly align with dna
#Update dna.blend to be unique

bl_info = {
    "name": "DNA Segmentator",
    "author": "Ryan Mulqueen <RMulqueen@mdanderson.org>",
    "version": (0, 4),
    "blender": (3, 2, 2),
    "location": "View3D",
    "description": "Adds a menu to build DNA molecules",
    "warning": "",
    "wiki_url": "",
    "category": "Mesh",
}


import bpy
import os, sys
from mathutils import Vector
from math import pi

def setup_initial_curve(
    curve_name="gDNA_path",
    collection_name="gDNA"):
    
    curve=bpy.ops.curve.primitive_bezier_curve_add(
        radius=1, 
        align='WORLD', 
        location=(0, 0, 0), 
        scale=(1, 1, 1),
        enter_editmode=True)
    bpy.context.active_object.name = curve_name
    #this is to set handles straight so there isn't a curve
    bpy.ops.transform.resize(value=(1, 0, 1), 
        orient_type='GLOBAL', 
        orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), 
        orient_matrix_type='GLOBAL') #set all handles to 0 so the curve is straight
    bpy.ops.object.mode_set(mode='OBJECT')#ensure we are back in object mode
    scene = bpy.context.scene
    if collection_name not in bpy.data.collections: #add collection if name isn't in collections yet
        bpy.data.collections.new(name=collection_name) #make new collection
    if bpy.data.collections[collection_name] not in list(scene.collection.children): #link collection to scene if it isn't in yet
        scene.collection.children.link(bpy.data.collections[collection_name]) #link new collection to scene
    obj=bpy.context.active_object
    if obj not in list(bpy.data.collections[collection_name].objects): #add collection if name isn't in collections yet
        bpy.data.collections[collection_name].objects.link(obj) #link curve to new collection
    if obj in list(bpy.context.scene.collection.objects):
        bpy.context.scene.collection.objects.unlink(obj) #remove link from scene collection
    if obj in list(bpy.data.collections["Collection"].objects):
        bpy.data.collections["Collection"].objects.unlink(obj) #remove link from default collection

def load_in_structure(
    file_path="", #change this to the directory for dna.blend
    import_object_name="DNA_Surface_real",
    collection_name="gDNA",
    object_name="gDNA",
    curve_name="gDNA_path",
    dna_length=100,
    dna_color=(0.3, 0.3, 0.8, 1)):
    dna_count=int(dna_length/11)#set up count of 11bp turns to use in DNA
    inner_path="Object"
    file="dna.blend"
    #in case the imported file name overlaps with existing data using an old and new set of names for omparison
    old_set = set(bpy.data.objects[:])
    bpy.ops.wm.append(
        filepath=os.path.join(file_path,inner_path, import_object_name),
        directory=os.path.join(file_path,inner_path),
        filename=import_object_name
        )
    bpy.ops.object.select_all(action='DESELECT')
    new_set = set(bpy.data.objects[:]) - old_set
    #Set up modifiers to make generative DNA
    #first set up array to duplicate DNA enough to fit the curve
    obj=list(new_set)[0]
    obj.name= object_name #change name of object we just made
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='ARRAY')#Add array modifier to dna and assign to path
    bpy.context.object.modifiers["Array"].fit_type ="FIXED_COUNT" #fit array of objects to curve
    bpy.context.object.modifiers["Array"].count = dna_count
    obj.modifiers["Array"].curve = bpy.data.objects[curve_name] #assign the curve as the one we created
    obj.modifiers["Array"].relative_offset_displace[0] = 0 #don't go by x value
    obj.modifiers["Array"].relative_offset_displace[2] = 0.88 #array by z value (height) and squich a bit to close model gap
    obj.modifiers["Array"].use_merge_vertices = True #merge to close any small gaps
    #now use curve modifier to deform dna to fit more snuggly
    bpy.ops.object.modifier_add(type='CURVE')
    obj.modifiers["Curve"].object = bpy.data.objects[curve_name]
    obj.modifiers["Curve"].deform_axis = 'POS_Z'
    #add a bit of wiggle to the dna in the case of animation
    bpy.ops.object.modifier_add(type='WAVE')
    obj.modifiers["Wave"].use_y = False
    obj.modifiers["Wave"].height = 0.01
    obj.modifiers["Wave"].speed = 0.1
    #make material and add to DNA
    mat = bpy.data.materials.new(name=object_name)
    mat.use_nodes = True #use node trees, these can be seen by switching a panel to the shader editor if you want. It will look like the above shader, just not nicely placed.
    mat_nodes = mat.node_tree.nodes
    mat_links = mat.node_tree.links
    mat = bpy.data.materials[object_name] #Get the material you want 
    bpy.data.materials[object_name].node_tree.nodes["Principled BSDF"].inputs[0].default_value = dna_color
    bpy.data.materials[object_name].node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0.2
    #Assign material
    obj.select_set(True)
    obj.data.materials.clear()#clear material from new dna
    mat = bpy.data.materials.get(object_name) #set material properties of collection
    obj.data.materials.append(mat) #add material
    #Add to collection and clean from other collections
    if obj not in list(bpy.data.collections[collection_name].objects): #add collection if name isn't in collections yet
        bpy.data.collections[collection_name].objects.link(obj) #link curve to new collection
    if obj in list(bpy.context.scene.collection.objects):
        bpy.context.scene.collection.objects.unlink(obj) #remove link from scene collection
    if obj in list(bpy.data.collections["Collection"].objects):
        bpy.data.collections["Collection"].objects.unlink(obj) #remove link from default collection

def trim_and_move_curve(
    curve_name="gDNA_path",
    collection_name="gDNA",
    object_name="gDNA",
    attach_to=None,
    attach_left_or_right=None,
    hybridize_to=None,
    hybridize_anchor_left_or_right=None):
    #Cut curve to size of dna bounding box
    ob=bpy.data.objects[object_name]
    ob.rotation_mode="XYZ"
    spline=bpy.data.objects[curve_name].data.splines[0]#spline defining curve, since all paths are a simple single spline, we can use 0 index
    spline.bezier_points[0].co.x=[ob.matrix_world @ Vector(corner) for corner in ob.bound_box][0].x #set left most x value of DNA
    spline.bezier_points[1].co.x=[ob.matrix_world @ Vector(corner) for corner in ob.bound_box][7].x #set right most x value of DNA
    if attach_to is not None and attach_left_or_right is not None: #get target object to attach to
        attach_ob=attach_to
        length_curve=abs(spline.bezier_points[0].co.x-spline.bezier_points[1].co.x) #gets DNA repeat defined length
        if attach_left_or_right=="left": #attach at right point (so target DNA is left of the attach to DNA)
            spline.bezier_points[1].co.x=[attach_ob.matrix_world @ Vector(corner) for corner in attach_ob.bound_box][0].x+0.02 #set right most x value of DNA to left most object point
            spline.bezier_points[0].co.x=spline.bezier_points[1].co.x-length_curve+0.02 #set left most DNA point to length distance from right point
        if attach_left_or_right=="right":
            spline.bezier_points[0].co.x=[attach_ob.matrix_world @ Vector(corner) for corner in attach_ob.bound_box][7].x-0.02#set right most x value of DNA to left most object point
            spline.bezier_points[1].co.x=spline.bezier_points[0].co.x+length_curve-0.02 #set left most DNA point to length distance from right point
    if hybridize_to is not None and hybridize_anchor_left_or_right is not None: #get target segment to hybridize to
        hyb_ob=hybridize_to
        length_curve=abs(spline.bezier_points[0].co.x-spline.bezier_points[1].co.x) #gets DNA repeat defined length
        if hybridize_anchor_left_or_right=="left": #attach at right point (so target DNA is left of the attach to DNA), this helps with anchoring, since annealed oligoes can hang off edge
            spline.bezier_points[1].co.x=[hyb_ob.matrix_world @ Vector(corner) for corner in hyb_ob.bound_box][0].x+0.02 #set right most x value of DNA to left most object point
            spline.bezier_points[0].co.x=spline.bezier_points[1].co.x-length_curve+0.02 #set left most DNA point to length distance from right point
            rotation=hyb_ob.rotation_euler.z+(270*pi/180)
            loc_change=(0.025,0.025,0.025)
        if hybridize_anchor_left_or_right=="right":
            rotation=hyb_ob.rotation_euler.z+(270*pi/180)
            loc_change=(0.025,0.025,-0.025)
            spline.bezier_points[0].co.x=[hyb_ob.matrix_world @ Vector(corner) for corner in hyb_ob.bound_box][7].x-0.02#set right most x value of DNA to left most object point
            spline.bezier_points[1].co.x=spline.bezier_points[0].co.x+length_curve-0.02 #set left most DNA point to length distance from right point
        ob=bpy.data.objects[object_name]
        ob.delta_rotation_euler.z=rotation # rotate the object #just euler stuff
        ob.matrix_world.translation=loc_change # nudge a bit to properly align
    spline.bezier_points[0].handle_left[0]=spline.bezier_points[0].co[0]-0.1 #Also change the bezier handles to be a bit smaller
    spline.bezier_points[0].handle_right[0]=spline.bezier_points[0].co[0]+0.1
    spline.bezier_points[1].handle_left[0]=spline.bezier_points[1].co[0]-0.1
    spline.bezier_points[1].handle_right[0]=spline.bezier_points[1].co[0]+0.1

def add_text_to_curve(
    curve_name="gDNA_path",
    object_name="gDNA",
    collection_name="gDNA",
    text_body="gDNA",
    yloc=0.1):
    #Align text to path and offset
    text_name=curve_name+"_text"
    bpy.ops.object.text_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
    bpy.context.active_object.name = text_name
    text_ob=bpy.data.objects[text_name]
    text_ob.data.extrude = 0.01 #solidify text
    text_ob.rotation_euler.z = pi/2 #justeulerthings rotate 180
    text_ob.scale=(0.4,0.4,0.4)#rescale
    text_ob.data.size = 0.5 #resize text
    spline=bpy.data.objects[curve_name].data.splines[0]#spline defining curve, since all paths are a simple single spline, we can use 0 index
    location_x=(spline.bezier_points[1].co.x+spline.bezier_points[0].co.x)/2 #get midpoint on path
    text_ob.location.z=0.1 #set locations
    text_ob.location.x=location_x
    text_ob.data.body=text_body #finally add in the text    
    text_ob.location.y=yloc
    mat = bpy.data.materials.get(object_name) #material name is same as object name by default
    text_ob.data.materials.append(mat)#set material to match object
    #add a bit of wiggle to the dna in the case of animation
    bpy.ops.object.modifier_add(type='WAVE')
    text_ob.modifiers["Wave"].use_y = False
    text_ob.modifiers["Wave"].height = 0.01
    text_ob.modifiers["Wave"].speed = 0.1
    #Add to collection and clean from other collections
    if text_ob not in list(bpy.data.collections[collection_name].objects): #add collection if name isn't in collections yet
        bpy.data.collections[collection_name].objects.link(text_ob) #link curve to new collection
    if text_ob in list(bpy.context.scene.collection.objects):
        bpy.context.scene.collection.objects.unlink(text_ob) #remove link from scene collection
    if text_ob in list(bpy.data.collections["Collection"].objects):
        bpy.data.collections["Collection"].objects.unlink(text_ob) #remove link from default collection

def add_segment(
    segment_name="gDNA",
    dna_color="#808080",
    dna_length=50,
    import_object_name="DNA_Surface_real",
    attach_to=None,
    attach_left_or_right=None,
    hybridize_to=None,
    hybridize_anchor_left_or_right=None,
    text_body=None,
    file_path="",
    yloc=0.1):
    if segment_name not in bpy.data.objects: #make sure we aren't trying to make duplicate objects
        #dna_color=hex_color_to_blender_rgb(dna_color)
        setup_initial_curve(
            curve_name=segment_name+"_path",
            collection_name=segment_name)
        load_in_structure(
            collection_name=segment_name,
            object_name=segment_name,
            curve_name=segment_name+"_path",
            dna_length=dna_length,
            dna_color=dna_color,
            import_object_name=import_object_name,
            file_path=file_path)
        trim_and_move_curve(
            curve_name=segment_name+"_path",
            collection_name=segment_name,
            object_name=segment_name,
            attach_to=attach_to,
            attach_left_or_right=attach_left_or_right)
        if text_body is not None:
            add_text_to_curve(curve_name=segment_name+"_path",
                object_name=segment_name,
                collection_name=segment_name,
                text_body=text_body,
                yloc=yloc)
    else:
        print("Segment name already in scene.")

#using a simple poll to decrease list of dna segments
def choseable_dna_segment(self, object):
    return not object.name.endswith("text") and not object.name.endswith("path")

def add_tn5(
    file_path="",
    collection_name="Tn5_gDNA_left",
    color=(0.3, 0.3, 0.8, 1),
    attach_to="gDNA",
    attach_left_or_right="left"
    ):
    #tn5_loaded_oligo_list=["ME_right","tn5_idx_left","truseq_i7_left"]
    inner_path="Object"
    scene = bpy.context.scene
    attach_to=attach_to.name+"_path" #select the material, but we actually want to just deform it to the path
    if collection_name not in bpy.data.collections: #add collection if name isn't in collections yet
        bpy.data.collections.new(name=collection_name) #make new collection
    if bpy.data.collections[collection_name] not in list(scene.collection.children): #link collection to scene if it isn't in yet
        scene.collection.children.link(bpy.data.collections[collection_name]) #link new collection to scene
    object_name=collection_name #collection name and object name the same
    #make material and add to tn5
    mat = bpy.data.materials.new(name=collection_name)
    mat.use_nodes = True #use node trees, these can be seen by switching a panel to the shader editor if you want. It will look like the above shader, just not nicely placed.
    mat_nodes = mat.node_tree.nodes
    mat_links = mat.node_tree.links
    mat = bpy.data.materials[object_name] #Get the material you want 
    bpy.data.materials[object_name].node_tree.nodes["Principled BSDF"].inputs[0].default_value = color
    bpy.data.materials[object_name].node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0.2
    #tn5 1
    import_object_name_1="Tn5_1"
    old_set = set(bpy.data.objects[:])
    bpy.ops.wm.append(
        filepath=os.path.join(file_path,inner_path, import_object_name_1),
        directory=os.path.join(file_path, inner_path),
        filename=import_object_name_1, use_recursive=True
        )
    new_set = set(bpy.data.objects[:]) - old_set
    obj=list(new_set)[0]
    obj.name= object_name+"_1" #change name of object we just made
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='CURVE')#Add array modifier to dna and assign to path
    obj.modifiers["Curve"].object = bpy.data.objects[attach_to]
    obj.modifiers["Curve"].deform_axis = 'POS_Z'
    #Assign material
    obj.select_set(True)
    obj.data.materials.clear()#clear material from new dna
    mat = bpy.data.materials.get(object_name) #set material properties of collection
    obj.data.materials.append(mat) #add material
    #tn5 2
    import_object_name_1="Tn5_2"
    old_set = set(bpy.data.objects[:])
    bpy.ops.wm.append(
        filepath=os.path.join(file_path,inner_path, import_object_name_1),
        directory=os.path.join(file_path,inner_path),
        filename=import_object_name_1, use_recursive=True
        )
    new_set = set(bpy.data.objects[:]) - old_set
    obj=list(new_set)[0]
    obj.name= object_name+"_2" #change name of object we just made
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='CURVE')#Add array modifier to dna and assign to path
    obj.modifiers["Curve"].object = bpy.data.objects[attach_to]
    obj.modifiers["Curve"].deform_axis = 'POS_Z'
    #Assign material
    obj.select_set(True)
    obj.data.materials.clear()#clear material from new dna
    mat = bpy.data.materials.get(object_name) #set material properties of collection
    obj.data.materials.append(mat) #add material
    for i in [object_name+"_1",object_name+"_2"]:#Add to collection and clean from other collections
        obj=bpy.data.objects[i]
        if obj not in list(bpy.data.collections[collection_name].objects): #add collection if name isn't in collections yet
            bpy.data.collections[collection_name].objects.link(obj) #link curve to new collection
        if obj in list(bpy.context.scene.collection.objects):
            bpy.context.scene.collection.objects.unlink(obj) #remove link from scene collection
        if obj in list(bpy.data.collections["Collection"].objects):
            bpy.data.collections["Collection"].objects.unlink(obj) #remove link from default collection


#should probably update this to be relative, and use a proper addon name
DNA_PROPS = [
    ("dna_length", bpy.props.IntProperty(
        name="Length of DNA",
        description="Length of DNA (bp) of segment to generate.",
        default=25,
        soft_max=100,
        min=1
        )),
    ("segment_name", bpy.props.StringProperty(
        name="Segment Name",
        description="Name of DNA segment to generate.",
        default="gDNA"
        )),
    ("color_set", bpy.props.FloatVectorProperty(
         name = "Color Picker",
         subtype = "COLOR",
         default = (1.0,1.0,1.0,1.0),
         size = 4
         )),
    ("dna_type", bpy.props.EnumProperty(
        name="DNA/RNA Type",
        description="Select DNA or RNA type.",
        default="DNA_Surface_real",
        items=[("DNA_Surface_real","DNA","DNA Surface"),("RNA_Surface_real","RNA","RNA Surface")]
        )),
    ("text_body", bpy.props.StringProperty(
        name="Floating Text",
        description="Adds a floating text descriptor next to the segment.",
        default="gDNA"
        )),
    ("text_body_position", bpy.props.FloatProperty(
        name="Text Position",
        subtype = "DISTANCE",
        unit = "LENGTH",
        description="Adjusts floating text descriptor location",
        default=0.1,
        min=-1,
        max=1
        )),
    ("attach_to", bpy.props.PointerProperty(
        name="Attach to:",
        type=bpy.types.Object,
        poll=choseable_dna_segment,
        description="Select DNA segment to attach to.")),
    ("attach_left_or_right", bpy.props.EnumProperty(
        name="Attach Position for Target",
        description="Attach to the left (-X) or right (+X) of target molecule?",
        default="None",
        items=[("None","None","None"),("left","left","left"),("right","right","right")]
        )),
    ("hybridize_to", bpy.props.PointerProperty(
        name="Hybridize to:",
        type=bpy.types.Object,
        poll=choseable_dna_segment,
        description="Select DNA segment to hybridize to.")),
    ("hybridize_anchor_left_or_right", bpy.props.EnumProperty(
        name="Anchoring Position for Generated Segment",
        description="Anchor hybrid DNA to the left (-X) or right (+X) of target molecule?",
        default="None",
        items=[("None","None","None"),("left","left","left"),("right","right","right")]
        )),
    ("dna_modifier_color", bpy.props.FloatVectorProperty(
        name="New Color",
        subtype='COLOR',
        default = (1.0,1.0,1.0,1.0),
        size = 4))
]

TN5_PROPS = [
    ("tn5_attach_to", bpy.props.PointerProperty(
        name="Attach to:",
        type=bpy.types.Object,
        poll=choseable_dna_segment,
        description="Select DNA segment to attach to.")),
    ("tn5_attach_left_or_right", bpy.props.EnumProperty(
        name="Attach Position for Target",
        description="Attach to the left (-X) or right (+X) of target molecule?",
        default="left",
        items=[("left","left","left"),("right","right","right")]
        )),
    ("tn5_collection_name", bpy.props.StringProperty(
        name="Collection Name",
        description="Name of Tn5 Collection Name",
        default="Tn5_gDNA_left"
        )),
    ("tn5_color_set", bpy.props.FloatVectorProperty(
         name = "Color Picker",
         subtype = "COLOR",
         default = (1.0,1.0,1.0,1.0),
         size = 4
         )),
]
class Preferences_dna_maker(bpy.types.AddonPreferences):
    
    bl_idname= __name__ #change this to __package__ when module is packaged together

    filepath_blendfile: bpy.props.StringProperty(
        name="Path to Packaged Blender File",
        subtype='FILE_PATH',
    )


    def draw(self, context):
        layout = self.layout
        layout.label(text="Please add path to packaged blend file containing PDB structures (dna.blend).")
        layout.prop(self, "filepath_blendfile")

class MESH_OT_dna_maker(bpy.types.Operator):
    """Let's make some DNA!"""
    bl_idname = "mesh.dna_maker"
    bl_label = "DNA Maker"
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props=context.scene
        add_segment(
        file_path=context.preferences.addons['dna_built'].preferences.filepath_blendfile,
        segment_name=props.segment_name,
        dna_length=props.dna_length,
        dna_color=props.color_set,
        import_object_name=props.dna_type,
        text_body=props.text_body,
        yloc=props.text_body_position,
        attach_to=props.attach_to,
        attach_left_or_right=props.attach_left_or_right,
        hybridize_to=props.hybridize_to,
        hybridize_anchor_left_or_right=props.hybridize_anchor_left_or_right)
        return {'FINISHED'}

class VIEW3D_PT_dna_maker(bpy.types.Panel):
   """Panel for DNA Molecule Builder"""
   bl_space_type="VIEW_3D"
   bl_region_type="UI"
   bl_category="DNA Maker"
   bl_label="Add Segments"

   def draw(self,context):
        layout = self.layout
        obj =  context.active_object
        #Initial default option for messing with gDNA
        col = layout.row()
        col.operator('mesh.dna_maker',
            text="Make DNA",
            icon="RNA")
        layout.label(text="DNA Layout:")
        for prop_name in ["segment_name","dna_length","color_set","dna_type"]:
            col = layout.row()
            col.prop(context.scene, prop_name)
        layout.label(text="Text:")
        for prop_name in ["text_body","text_body_position"]:
            col = layout.row()
            col.prop(context.scene, prop_name)
        layout.label(text="Attach to:")
        col=layout.row()
        col.prop(context.scene,"attach_to")
        if context.scene.attach_to != None:
            col=layout.row()
            col.prop(context.scene,"attach_left_or_right")
        col=layout.row()
        layout.label(text="Hybridize to:")
        col=layout.row()
        col.prop(context.scene,"hybridize_to")
        if context.scene.hybridize_to != None:
            col=layout.row()
            col.prop(context.scene,"hybridize_anchor_left_or_right")

class MESH_OT_tn5_addition(bpy.types.Operator):
    """Let's add some proteins!"""
    bl_idname = "mesh.tn5_addition"
    bl_label = "DNA Maker"
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props=context.scene
        add_tn5(
        file_path=context.preferences.addons['dna_built'].preferences.filepath_blendfile,
        collection_name=props.tn5_collection_name,
        color=props.tn5_color_set,
        attach_to=props.tn5_attach_to,
        attach_left_or_right=props.tn5_attach_left_or_right)
        return {'FINISHED'}

class VIEW3D_PT_tn5_addition(bpy.types.Panel):
   """Panel for DNA Molecule Builder"""
   bl_space_type="VIEW_3D"
   bl_region_type="UI"
   bl_category="DNA Maker"
   bl_label="Add Tn5"

   def draw(self,context):
        layout = self.layout
        obj =  context.object
        #Modify color of DNA and Text through manipulating the material properties
        col = layout.row()
        layout.label(text="Select a DNA segment to tagment.")
        layout.prop_search(context.scene, "tn5_attach_to", context.scene, "objects", text="DNA")
        for prop_name in ["tn5_attach_left_or_right","tn5_collection_name","tn5_color_set"]:
            col = layout.row()
            col.prop(context.scene, prop_name)
        col = layout.row()
        if context.scene.tn5_attach_to != None:
            col.operator('mesh.tn5_addition',
                text="Tagment DNA",
                icon="RNA")


CLASSES = [
    Preferences_dna_maker,
    MESH_OT_dna_maker,
    VIEW3D_PT_dna_maker,
    MESH_OT_tn5_addition,
    VIEW3D_PT_tn5_addition
]


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    for (prop_name, prop_value) in DNA_PROPS:
        setattr(bpy.types.Scene, prop_name, prop_value)
    for (prop_name, prop_value) in TN5_PROPS:
        setattr(bpy.types.Scene, prop_name, prop_value)
    

def unregister():
    for cls in CLASSES:
        bpy.utils.unregister_class(cls)
    for (prop_name, prop_value) in DNA_PROPS:
        delattr(bpy.types.Scene, prop_name)
    for (prop_name, prop_value) in TN5_PROPS:
        delattr(bpy.types.Scene, prop_name)




#hybridization looks like it has a weird logic flow or something, attaches correctly, but doesn't rotate or translate
#To Do
#Add proteins as additional function
#Adjust Tn5 angle to correct for DNA
#Add DNA copy for second Tn5
```
