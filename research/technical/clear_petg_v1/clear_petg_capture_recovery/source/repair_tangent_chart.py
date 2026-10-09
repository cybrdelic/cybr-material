"""Repair only the undefined anisotropic frame on cut faces.
Texture sampling remains on the original active-render chart.
"""
import bpy,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R.parents[1];WORK=P.parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();report=json.loads((R/'receipts/source_inspection.json').read_text());src=Path(report['scene']);assert sha(src)==report['scene_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene;ob=s.objects['23_petg_transparent / true1.2mm printed panel'];me=ob.data;nl=len(me.loops);ntop=2*128*1040
v=np.empty(len(me.vertices)*3,'f4');me.vertices.foreach_get('co',v);old_geometry=hashlib.sha256(v.tobytes()).hexdigest();original=me.uv_layers['UVMap'];old_uv=np.empty(nl*2,'f4');original.data.foreach_get('uv',old_uv);old_uv=old_uv.reshape(-1,2);assert original.active_render
uv=me.uv_layers.new(name='PhysicalAnisotropyFrame');uv.data.foreach_set('uv',old_uv.ravel());a=old_uv.copy();area=[]
for f in me.polygons:
 if f.index<ntop:continue
 # Cut faces are planar vertical walls, with a longitudinal surface coordinate
 # and physical Z as the second coordinate. The cross-sectional texture itself
 # retains the old chart; only the anisotropic basis consumes this new chart.
 coord=0 if abs(f.normal.y)>abs(f.normal.x) else 1
 for k in f.loop_indices:
  q=me.vertices[me.loops[k].vertex_index].co;a[k]=(q[coord]/.02+.5,q.z/.02+.5)
 q=a[list(f.loop_indices)].astype('f8');ar=abs(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)))/2;area.append(ar)
assert min(area)>1e-10 and len(area)==2336
uv.data.foreach_set('uv',a.ravel());me.uv_layers.active_index=0;me.uv_layers['UVMap'].active_render=True
assert me.uv_layers.active.name=='UVMap' and me.uv_layers['UVMap'].active_render and not me.uv_layers['PhysicalAnisotropyFrame'].active_render
toploops=ntop*4;assert np.array_equal(a[:toploops],old_uv[:toploops]);m=me.materials[0];tangents=[n for n in m.node_tree.nodes if n.type=='TANGENT'];assert len(tangents)==1 and tangents[0].direction_type=='UV_MAP' and tangents[0].uv_map=='UVMap';tangents[0].uv_map=uv.name
assert not m.node_tree.nodes['Principled BSDF'].inputs['Normal'].links
bpy.context.view_layer.update();me.calc_tangents(uvmap=uv.name);t=np.empty(nl*3,'f4');me.loops.foreach_get('tangent',t);t=t.reshape(-1,3);length=np.linalg.norm(t,axis=1);assert np.isfinite(t).all() and length.min()>.999
# Sign/orientation must remain the exact inherited basis on both printed faces.
new_top=t[:toploops].copy();me.free_tangents();me.calc_tangents(uvmap='UVMap');old=np.empty(nl*3,'f4');me.loops.foreach_get('tangent',old);old=old.reshape(-1,3);difference=float(np.max(abs(new_top-old[:toploops])));assert difference<1e-6;me.free_tangents()
me.uv_layers.active_index=0;me.uv_layers['UVMap'].active_render=True
out=R/'scenes/23_petg_native4k_valid_tangent_guides_r2.blend';bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(src)==report['scene_sha256']
bpy.ops.wm.open_mainfile(filepath=str(out));me=bpy.context.scene.objects['23_petg_transparent / true1.2mm printed panel'].data;check=np.empty(len(me.vertices)*3,'f4');me.vertices.foreach_get('co',check);assert hashlib.sha256(check.tobytes()).hexdigest()==old_geometry;checkuv=np.empty(nl*2,'f4');me.uv_layers['UVMap'].data.foreach_get('uv',checkuv);assert np.array_equal(checkuv,old_uv.ravel());assert me.uv_layers['UVMap'].active_render and not me.uv_layers['PhysicalAnisotropyFrame'].active_render
d={'source':str(src),'source_sha256':sha(src),'scene':str(out),'scene_sha256':sha(out),'geometry_exact':True,'original_color_and_roughness_UVs_bit_exact':True,'printed_face_tangents_max_component_error':difference,'repaired_cut_faces':len(area),'minimum_cut_chart_area':min(area),'all_tangent_lengths_minimum':float(length.min()),'only_shader_binding_change':'Existing Tangent node reads PhysicalAnisotropyFrame; BaseColor/Roughness/Metallic retain the original active-render UVMap','all_optical_coefficients_and_volume_absorption_preserved':True,'physical_material_qualification':False}
(R/'receipts/tangent_repair.json').write_text(json.dumps(d,indent=2)+'\n');print('CLEAR_TANGENT_REPAIRED',json.dumps(d),flush=True)
