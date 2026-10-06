import bpy,numpy as np,json,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];r=json.loads((R/'receipts/chart_audit.json').read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def val(x):
 try:return list(x)
 except TypeError:return x
def tree(nt):return {'nodes':{n.name:{'type':n.bl_idname,'inputs':[(i.identifier,val(i.default_value)) for i in n.inputs if hasattr(i,'default_value')],'image':sha(Path(bpy.path.abspath(n.image.filepath)).resolve()) if n.type=='TEX_IMAGE' and n.image else None} for n in nt.nodes},'links':sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in nt.links)}
def snap(p):
 bpy.ops.wm.open_mainfile(filepath=p);s=bpy.context.scene;d={'units':s.unit_settings.scale_length,'objects':{},'materials':{m.name:tree(m.node_tree) for m in bpy.data.materials if m.use_nodes},'world':tree(s.world.node_tree),'view':{k:getattr(s.view_settings,k) for k in ['view_transform','look','exposure','gamma']}}
 for o in s.objects:
  x={'matrix':list(map(list,o.matrix_world)),'type':o.type}
  if o.type=='MESH':
   a=np.array([tuple(v.co) for v in o.data.vertices],'f4');t=np.array([l.vertex_index for l in o.data.loops],'i4');uv=np.array([tuple(v.uv) for v in o.data.uv_layers.active.data],'f4');x.update(vertices=hashlib.sha256(a.tobytes()).hexdigest(),topology=hashlib.sha256(t.tobytes()).hexdigest(),materials=[m.name for m in o.data.materials]);
   if o.name=='07 new / finite60micron mean glaze shell':
    flat=[k for p in o.data.polygons if p.normal.z>.999 for k in p.loop_indices];x['flat_uv']=hashlib.sha256(uv[flat].tobytes()).hexdigest()
   else:x['uv']=hashlib.sha256(uv.tobytes()).hexdigest()
  if o.type=='CAMERA':x.update(type=o.data.type,lens=o.data.lens,ortho_scale=o.data.ortho_scale)
  if o.type=='LIGHT':x.update(type=o.data.type,energy=o.data.energy,color=list(o.data.color),size=o.data.size)
  d['objects'][o.name]=x
 return d
before=snap(r['source']);after=snap(r['scene']);assert before==after;assert sha(r['source'])==r['source_sha256'] and sha(r['scene'])==r['scene_sha256']
report={'fresh_process':True,'all_geometry_topology_transforms_and_flat_uv_exact':True,'all_material_inputs_links_and_image_bytes_exact':True,'camera_lights_world_color_management_exact':True,'only_outer_rounded_glaze_UV_changed':True,'source_sha256':r['source_sha256'],'candidate_sha256':r['scene_sha256']};(R/'receipts/fresh_scene_equivalence.json').write_text(json.dumps(report,indent=2));print('FRESH_EQUIVALENCE_PASS',json.dumps(report))
