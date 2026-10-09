import bpy,bmesh,numpy as np,json,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R.parents[1];WORK=P.parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rec=json.loads((R/'receipts/reconstruction.json').read_text());src=Path(rec['scene']);assert sha(src)==rec['scene_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene;ob=s.objects['23_petg_transparent / true1.2mm printed panel'];assert len(ob.data.vertices)==268578 and len(ob.data.polygons)==268576
bm=bmesh.new();bm.from_mesh(ob.data);assert all(len(e.link_faces)==2 for e in bm.edges);volume=bm.calc_volume(signed=True);bm.free();assert volume>0
v=np.empty(len(ob.data.vertices)*3,'f4');ob.data.vertices.foreach_get('co',v);v=v.reshape(-1,3);dims=v.max(0)-v.min(0);assert np.max(abs(dims-np.array([.02800000086426735,.026000000536441803,.001304594217799604])))<1e-9
uv=ob.data.uv_layers.active;collapsed=[];topface_count=2*128*1040
for f in ob.data.polygons:
 if f.index<topface_count:continue
 q=np.array([uv.data[i].uv for i in f.loop_indices]);area=abs(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)))/2
 if area<1e-16:collapsed.append(f.index)
def val(a):
 try:return list(a)
 except TypeError:return a
def signature(allowed=None):
 out={}
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  names=set(allowed[m.name]) if allowed else {n.name for n in m.node_tree.nodes}
  out[m.name]={'nodes':{n.name:{'type':n.bl_idname,'inputs':[(i.identifier,val(i.default_value)) for i in n.inputs if hasattr(i,'default_value')]} for n in m.node_tree.nodes if n.name in names},'links':sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links if l.from_node.name in names and l.to_node.name in names)}
 return out
before=signature();allowed={k:list(v['nodes']) for k,v in before.items()};guide=[]
for m in bpy.data.materials:
 if not m.use_nodes or not m.users:continue
 nd=m.node_tree.nodes;lk=m.node_tree.links;output=next(n for n in nd if n.type=='OUTPUT_MATERIAL' and n.is_active_output);p=output.inputs['Surface'].links[0].from_node;assert p.type=='BSDF_PRINCIPLED';assert not p.inputs['Alpha'].links and p.inputs['Alpha'].default_value==1
 transmission=p.inputs['Transmission Weight'];assert not transmission.is_linked;tw=transmission.default_value
 if output.inputs['Volume'].is_linked:
  vol=output.inputs['Volume'].links[0].from_node;assert vol.type=='VOLUME_ABSORPTION';assert not any(i.is_linked for i in vol.inputs);assert vol.inputs['Density'].default_value==8;assert tw>.95
 assert not p.inputs['Normal'].links
 al=nd.new('ShaderNodeOutputAOV');al.name='OIDN_SurfaceReflectivity';al.aov_name=al.name
 if tw>0:al.inputs['Color'].default_value=(1,1,1,1)
 elif p.inputs['Base Color'].links:lk.new(p.inputs['Base Color'].links[0].from_socket,al.inputs['Color'])
 else:al.inputs['Color'].default_value=p.inputs['Base Color'].default_value
 no=nd.new('ShaderNodeOutputAOV');no.name='OIDN_SurfaceNormal';no.aov_name=no.name;geo=nd.new('ShaderNodeNewGeometry');lk.new(geo.outputs['Normal'],no.inputs['Color'])
 guide.append({'material':m.name,'first_hit_feature':'unit dielectric feature' if tw else 'intrinsic matte color','transmission':tw,'homogeneous_absorption_retained':output.inputs['Volume'].is_linked})
for vl in s.view_layers:
 for name in ['OIDN_SurfaceReflectivity','OIDN_SurfaceNormal']:a=vl.aovs.add();a.name=name;a.type='COLOR'
assert signature(allowed)==before
dst=R/'scenes/23_petg_reconstructed_first_hit_guides.blend';bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(dst),compress=True);assert sha(src)==rec['scene_sha256']
bpy.ops.wm.open_mainfile(filepath=str(dst));assert signature(allowed)==before
report={'source':str(src),'source_sha256':sha(src),'scene':str(dst),'scene_sha256':sha(dst),'dimensions_m':dims.tolist(),'vertices':len(v),'faces':len(ob.data.polygons) if False else 268576,'closed_manifold':True,'positive_volume_m3':volume,'all_original_surface_and_volume_nodes_equal':True,'guides':guide,'existing_collapsed_cut_face_UVs':len(collapsed),'guide_reference':'https://www.openimagedenoise.org/documentation.html#albedos','guide_scope':'Documented first-hit unit feature for dielectric; actual shading normal. Homogeneous volume absorption is retained in beauty, not encoded as a view-dependent guide. Behind-surface detail is not protected by first-hit features.','physical_transparency_accepted':False,'noise_or_detail_acceptance':False}
(R/'receipts/source_inspection.json').write_text(json.dumps(report,indent=2)+'\n');print('CLEAR_GUIDE_SOURCE_CHECKED',json.dumps(report),flush=True)
