import bpy,json,hashlib
from pathlib import Path
R=Path('/workspace/shared/material-slabs');src=R/'scenes/05_champagne_brass_original_native4k_slab.blend';sha=hashlib.sha256(src.read_bytes()).hexdigest();bpy.ops.wm.open_mainfile(filepath=str(src));rows=[]
for mat in bpy.data.materials:
 if not mat.use_nodes:continue
 for p in mat.node_tree.nodes:
  if p.type!='BSDF_PRINCIPLED' or 'Anisotropic' not in p.inputs or p.inputs['Anisotropic'].default_value<=0:continue
  assert not p.inputs['Tangent'].is_linked
  n=mat.node_tree.nodes.new('ShaderNodeTangent');n.name='Manufacturing / explicit UV brush tangent';n.direction_type='UV_MAP';n.uv_map='UVMap';mat.node_tree.links.new(n.outputs['Tangent'],p.inputs['Tangent']);p.inputs['Anisotropic Rotation'].default_value=.25;rows.append({'material':mat.name,'anisotropic':p.inputs['Anisotropic'].default_value,'axis':'UV V','rotation':.25})
assert len(rows)==1
out=R/'scenes/05_champagne_brass_original_native4k_slab_r2_uv_tangent.blend';bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert hashlib.sha256(src.read_bytes()).hexdigest()==sha;r={'source':str(src),'source_sha256':sha,'scene':str(out),'scene_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'correction':rows,'native_height_gradient_rms_xy':[133.1051483154297,2.4721906185150146],'all_other_shader_geometry_camera_lighting_changes':False};(R/'receipts/05_brass_uv_tangent.json').write_text(json.dumps(r,indent=2));print('BRASS_TANGENT_READY',r['scene_sha256'],flush=True)
