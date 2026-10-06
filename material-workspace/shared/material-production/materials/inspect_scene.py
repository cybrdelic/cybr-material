"""Blender-only, read-only source inspection. No save or render side effects."""
import bpy,hashlib,json,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from materials.core import read,sha,write_new
from materials.path_config import remap_scene_images,configuration_state
args=sys.argv[sys.argv.index('--')+1:];source,expected,dest=args
if tuple(bpy.app.version)!=(4,3,2):raise RuntimeError('Selected contracts require pinned Blender4.3.2')
if sha(source)!=expected:raise RuntimeError('Source changed before inspection')
bpy.ops.wm.open_mainfile(filepath=source)
relocation_state=configuration_state()
relocations=remap_scene_images(bpy,relocation_state)
s=bpy.context.scene
if abs(s.unit_settings.scale_length-1)>1e-8:raise RuntimeError('Expected SI world units')
images=[]
for im in bpy.data.images:
 if im.source!='FILE' or not im.filepath:continue
 p=Path(bpy.path.abspath(im.filepath));packed=bool(im.packed_file)
 if not packed and not p.is_file():raise RuntimeError('Missing external image: '+str(p))
 digest=hashlib.sha256(bytes(im.packed_file.data)).hexdigest() if packed else sha(p)
 images.append({'name':im.name,'path':str(p),'packed':packed,'sha256':digest,'resolution':list(im.size),'colorspace':im.colorspace_settings.name})
materials=[]
for m in bpy.data.materials:
 if not m.use_nodes:continue
 nodes=m.node_tree.nodes
 bindings=[]
 for n in nodes:
  if n.type=='TEX_IMAGE' and n.image:
   data=any(word in n.image.name.lower() for word in ['normal','roughness','height','metallic','opacity','mask','ior'])
   if data and n.image.colorspace_settings.name!='Non-Color':raise RuntimeError('Data texture incorrectly color transformed: '+n.image.name)
   bindings.append({'image':n.image.name,'outputs':[{'node':l.to_node.name,'socket':l.to_socket.name} for o in n.outputs for l in o.links]})
 materials.append({'name':m.name,'bindings':bindings,'normal_nodes':[{'name':n.name,'space':n.space,'strength':n.inputs['Strength'].default_value} for n in nodes if n.type=='NORMAL_MAP']})
objects=[]
mesh_face_counts={}
for o in s.objects:
 if o.type not in {'MESH','CURVE','CURVES','SURFACE'} or o.hide_render:continue
 if o.type=='MESH' and o.data.as_pointer() not in mesh_face_counts:mesh_face_counts[o.data.as_pointer()]=dict(Counter(str(p.material_index) for p in o.data.polygons))
 objects.append({'name':o.name,'type':o.type,'dimensions_m':list(o.dimensions),'scale':list(o.scale),'vertices':len(o.data.vertices) if o.type=='MESH' else None,'polygons':len(o.data.polygons) if o.type=='MESH' else None,'uv_layers':[u.name for u in o.data.uv_layers] if o.type=='MESH' else [],'materials':[m.name if m else None for m in o.data.materials], 'material_face_counts':mesh_face_counts[o.data.as_pointer()] if o.type=='MESH' else {}, 'modifiers':[{'type':v.type,**({'material':v.material} if v.type=='BEVEL' else {})} for v in o.modifiers]})
if sha(source)!=expected:raise RuntimeError('Source mutated during inspection')
if configuration_state()!=relocation_state:raise RuntimeError('Resource relocation configuration/helper changed during inspection')
write_new(dest,{'blender_version':bpy.app.version_string,'unit_scale':s.unit_settings.scale_length,'images':images,'materials':materials,'objects':objects,'camera':s.camera.name if s.camera else None,'source_unchanged':True,'geometry_and_shader_bindings_inspected':True,'resource_relocations':relocations,'path_configuration':relocation_state,'rendered':False})
print('SCENE_INSPECTED',dest)
