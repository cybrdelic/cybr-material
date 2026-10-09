"""Render-relevant static scene identity, independent of file paths and blend bytes."""
import hashlib,json

def scene_signature():
 import bpy,numpy as np
 def digest(b):return hashlib.sha256(b).hexdigest()
 s=bpy.context.scene;bpy.context.view_layer.update();r={'objects':{},'materials':{},'world':{},'camera':s.camera.name,'view':[s.view_settings.view_transform,s.view_settings.look,s.view_settings.exposure,s.view_settings.gamma]}
 for o in s.objects:
  v={'matrix':[list(row) for row in o.matrix_world],'type':o.type,'hide_render':o.hide_render}
  if o.type=='MESH':
   me=o.data
   for label,data,prop,width,dtype in [('vertices',me.vertices,'co',3,'f4'),('loops',me.loops,'vertex_index',1,'i4'),('materials',me.polygons,'material_index',1,'i4'),('smooth',me.polygons,'use_smooth',1,'i4')]:
    a=np.empty(len(data)*width,dtype);data.foreach_get(prop,a);v[label]=digest(a.tobytes())
   if me.uv_layers.active:
    a=np.empty(len(me.loops)*2,'f4');me.uv_layers.active.data.foreach_get('uv',a);v['uv']=digest(a.tobytes())
   v['material_slots']=[x.name for x in me.materials]
  elif o.type=='CAMERA':v['data']={k:getattr(o.data,k) for k in ['type','lens','ortho_scale','clip_start','clip_end','sensor_width','sensor_height','sensor_fit','shift_x','shift_y']};v['dof']={'use_dof':o.data.dof.use_dof,'focus_distance':o.data.dof.focus_distance,'aperture_fstop':o.data.dof.aperture_fstop}
  elif o.type=='LIGHT':v['data']={k:list(getattr(o.data,k)) if k=='color' else getattr(o.data,k) for k in ['type','energy','color','shape','size','size_y','spread','use_shadow'] if hasattr(o.data,k)}
  r['objects'][o.name]=v
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  nodes={}
  for n in m.node_tree.nodes:
   nr={'type':n.type,'inputs':{str(i):list(x.default_value) if x.type in {'RGBA','VECTOR'} else x.default_value for i,x in enumerate(n.inputs) if hasattr(x,'default_value')}}
   if n.type=='TEX_IMAGE':nr.update(image_sha=digest(n.image.packed_file.data),space=n.image.colorspace_settings.name,interpolation=n.interpolation,extension=n.extension,projection=n.projection)
   if n.type=='NORMAL_MAP':nr.update(space=n.space,uv_map=n.uv_map)
   nodes[n.name]=nr
  r['materials'][m.name]={'nodes':nodes,'links':sorted((l.from_node.name,list(l.from_node.outputs).index(l.from_socket),l.to_node.name,list(l.to_node.inputs).index(l.to_socket)) for l in m.node_tree.links)}
 bg=s.world.node_tree.nodes['Background'];r['world']={'color':list(bg.inputs[0].default_value),'strength':bg.inputs[1].default_value}
 return {'sha256':digest(json.dumps(r,sort_keys=True,separators=(',',':')).encode()),'sections':{k:digest(json.dumps(v,sort_keys=True,separators=(',',':')).encode()) for k,v in r.items()}}
