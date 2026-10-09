"""Specific auxiliary fabrication finishes, authored entirely as shader nodes."""
import bpy

def finish(name,base,rough,metal=0,brush=False,emission=0):
    # Reuse finishes by specification; no per-fitting random tarnish.
    existing=bpy.data.materials.get(name)
    if existing:return existing
    mat=bpy.data.materials.new(name);mat.use_nodes=True
    n=mat.node_tree.nodes;l=mat.node_tree.links;p=n.get('Principled BSDF')
    p.inputs['Base Color'].default_value=(*base,1);p.inputs['Roughness'].default_value=rough;p.inputs['Metallic'].default_value=metal
    if emission:
        p.inputs['Emission Color'].default_value=(*base,1);p.inputs['Emission Strength'].default_value=emission
    if brush:
        p.inputs['Anisotropic'].default_value=.36
        coords=n.new('ShaderNodeTexCoord');mapping=n.new('ShaderNodeVectorMath');mapping.operation='MULTIPLY'
        mapping.inputs[1].default_value=(1800,1800,18);l.new(coords.outputs['Object'],mapping.inputs[0])
        noise=n.new('ShaderNodeTexNoise');noise.inputs['Scale'].default_value=1;noise.inputs['Detail'].default_value=1
        l.new(mapping.outputs[0],noise.inputs['Vector'])
        bump=n.new('ShaderNodeBump');bump.inputs['Distance'].default_value=.000003;bump.inputs['Strength'].default_value=.5
        l.new(noise.outputs['Fac'],bump.inputs['Height']);l.new(bump.outputs[0],p.inputs['Normal'])
    mat['finish_specification']=name;mat['source_image_inputs']=0
    return mat

def paint(name,base):
    mat=finish(name,base,.79)
    if mat.get('paint_structure'):return mat
    n=mat.node_tree.nodes;l=mat.node_tree.links;p=n.get('Principled BSDF')
    coords=n.new('ShaderNodeTexCoord');noise=n.new('ShaderNodeTexNoise');noise.inputs['Scale'].default_value=95;noise.inputs['Detail'].default_value=2
    l.new(coords.outputs['Object'],noise.inputs['Vector'])
    bump=n.new('ShaderNodeBump');bump.inputs['Distance'].default_value=.000035;bump.inputs['Strength'].default_value=.35
    l.new(noise.outputs['Fac'],bump.inputs['Height']);l.new(bump.outputs[0],p.inputs['Normal'])
    mat['paint_structure']='Fine rolled mineral paint; intact matte finish'
    return mat

def box(r,name,pos,size,mat,bevel=.002):
    ob=r.a.atelier.bevel_cube(name,pos,size,mat,min(bevel,min(size)*.4))
    for p in ob.data.polygons:p.use_smooth=True
    next(m for m in ob.modifiers if m.type=='BEVEL').harden_normals=True
    return ob

def cylinder(r,name,pos,radius,height,mat):
    spec={'tile_m':1,'height_scale_m':0}
    return r.a.lathe(name,[(.0001,0),(radius,0),(radius,height),(.0001,height)],mat,spec,pos,False,64)
