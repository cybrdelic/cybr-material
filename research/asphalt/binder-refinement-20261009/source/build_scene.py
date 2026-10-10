"""Rebuild the binder refinement coupon and matched studio from numeric source only."""
from pathlib import Path
import bpy,json,numpy as np,sys,hashlib
R=Path(__file__).resolve().parents[1]
assert bpy.app.version==(4,3,2),'Pinned scene contract requires Blender 4.3.2'
exec(compile((R/'source/verify_native.py').read_text(),str(R/'source/verify_native.py'),'exec'))
r=json.loads((R/'scene_recipe.json').read_text());bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene
for k,v in r['units'].items():setattr(s.unit_settings,k,v)
for name,d in r['materials'].items():
 m=bpy.data.materials.new(name);m.use_nodes=True;m.node_tree.nodes.clear()
 for name1,v in d['nodes'].items():
  n=m.node_tree.nodes.new(v['type']);n.name=name1
  for i,x in v['inputs'].items():n.inputs[int(i)].default_value=x
  for k in ['space','uv_map','interpolation','extension','projection']:
   if k in v:setattr(n,k,v[k])
  if 'image' in v:
   n.image=bpy.data.images.load(str(R/'native4096'/v['image']),check_existing=True);n.image.colorspace_settings.name=v['colorspace'];n.image.pack()
 for a,i,b,j in d['links']:m.node_tree.links.new(m.node_tree.nodes[a].outputs[i],m.node_tree.nodes[b].inputs[j])
for name,d in r['objects'].items():
 if d['type']=='MESH':
  if d.get('generator')=='asphalt_coupon':
   n=d['grid'];w=d['width_m'];g=np.load(R/'native4096/GeometryHeight.npy')[::-1];h=(np.float32(.008)+(g-np.float32(.5))*np.float32(.004)).ravel();verts=[((i/n-.5)*w,(j/n-.5)*w,float(h[j*(n+1)+i])) for j in range(n+1) for i in range(n+1)];uv=[(i/n,j/n) for j in range(n+1) for i in range(n+1)];faces=[]
   for j in range(n):
    for i in range(n):q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
   top=len(faces);boundary=list(range(n+1))+[j*(n+1)+n for j in range(1,n+1)]+[n*(n+1)+i for i in range(n-1,-1,-1)]+[j*(n+1) for j in range(n-1,0,-1)];bottom=[]
   for idx in boundary:bottom.append(len(verts));x,y,z=verts[idx];verts.append((x,y,d['bottom_m']));uv.append(uv[idx])
   for i,idx in enumerate(boundary):k=(i+1)%len(boundary);faces.append((idx,bottom[i],bottom[k],boundary[k]))
   faces.append(tuple(reversed(bottom)))
  else:verts=d['vertices'];faces=d['faces']
  data=bpy.data.meshes.new(name);data.from_pydata(verts,[],faces);data.update();o=bpy.data.objects.new(name,data)
  for mat in d['materials']:data.materials.append(bpy.data.materials[mat])
  layer=data.uv_layers.new(name='UVMap')
  for p in data.polygons:
   p.use_smooth=p.index<top if d.get('generator') else d['smooth'][p.index];p.material_index=(0 if p.index<top else 1) if d.get('generator') else d['indices'][p.index]
   for li in p.loop_indices:layer.data[li].uv=uv[data.loops[li].vertex_index] if d.get('generator') else d['uv'][li]
 elif d['type']=='CAMERA':
  data=bpy.data.cameras.new(name);o=bpy.data.objects.new(name,data)
  for k,v in d['data'].items():setattr(data,k,v)
  for k,v in d['dof'].items():setattr(data.dof,k,v)
 elif d['type']=='LIGHT':
  data=bpy.data.lights.new(name,d['data']['type']);o=bpy.data.objects.new(name,data)
  for k,v in d['data'].items():setattr(data,k,v)
 else:raise ValueError(d['type'])
 s.collection.objects.link(o)
 for k in ['location','rotation_euler','scale','hide_render']:setattr(o,k,d[k])
s.camera=bpy.data.objects[r['camera']];s.world=bpy.data.worlds.new('World');s.world.use_nodes=True;bg=s.world.node_tree.nodes['Background'];bg.inputs[0].default_value=r['world']['color'];bg.inputs[1].default_value=r['world']['strength']
for k,v in r['view'].items():setattr(s.view_settings,k,v)
s.render.engine='CYCLES';s.cycles.device='CPU';s.cycles.samples=128;s.cycles.use_adaptive_sampling=False;s.cycles.use_denoising=False;s.render.resolution_x=s.render.resolution_y=960;s.render.resolution_percentage=100;s.render.threads_mode='FIXED';s.render.threads=2
bpy.context.view_layer.update();sys.path.insert(0,str(R/'source'));from scene_contract import scene_signature
signature=scene_signature();expected=json.loads((R/'expected_scene.json').read_text());assert signature==expected, 'Generated scene differs from exact binder candidate render contract'
s['selected']=False;s['physical_quality_accepted']=False;s['CYBR_revision']='2026-10-09 / binder refinement numeric reconstruction';bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(R/'CYBR_Asphalt_Binder_Refinement_Native4096_Candidate.blend'),compress=True)
s['selected']=False
out=R/'CYBR_Asphalt_Binder_Refinement_Native4096_Candidate.blend'
job=json.loads((R/'capture_profile.json').read_text());target=R/'render';target.mkdir(exist_ok=True)
job.update(source=str(out),source_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),output=str(target/'INTERNAL_PANEL.png'))
(target/'job.json').write_text(json.dumps(job,indent=2))
(R/'build_receipt.json').write_text(json.dumps({'render_signature':signature,'exact_candidate_render_signature':True,'numeric_scene_reconstruction':True,'no_input_blend_required':True,'selected':False},indent=2))
print('Binder candidate geometry, UVs, shaders, maps, camera, lights, world and display reproduced exactly; matched render/job.json written.')
