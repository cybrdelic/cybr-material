"""Remove the legacy added height-variance broadening in one unselected control.
The source's assumed intrinsic roughness is retained exactly, never fitted to pixels.
"""
from pathlib import Path
import bpy,json,hashlib,numpy as np
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();old=json.loads((R/'receipts/scene_native.json').read_text());src=Path(old['scene']);assert sha(src)==old['scene_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
scope=json.loads((R/'receipts/roughness_overlap_audit.json').read_text());assert scope['intrinsic_roughness_retained']==.475;assert scope['maximum_height_displacement_pixels_at_current_view_upper_bound']<.05
ob=s.objects['21 / flat production interior / native'];m=ob.data.materials[0];p=m.node_tree.nodes['Principled BSDF'];assert len(p.inputs['Roughness'].links)==1
for l in list(p.inputs['Roughness'].links):m.node_tree.links.remove(l)
p.inputs['Roughness'].default_value=scope['intrinsic_roughness_retained']
for n in list(m.node_tree.nodes):
 if n.type not in {'OUTPUT_MATERIAL','BSDF_PRINCIPLED'}:m.node_tree.nodes.remove(n)
for im in list(bpy.data.images):
 if im.users==0:bpy.data.images.remove(im)
m['known_limitation']='Intrinsic optical roughness remains the original assumed0.475. Resolved native relief is no longer counted again in the old fixed-Gaussian variance closure. Subtexel physical polymer optics and smooth-triangle normal convergence remain unvalidated.'
s['study_case']='native_intrinsic_only';s['selected']=False
out=R/'scenes/21_plastic_native_intrinsic_only.blend';bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(src)==old['scene_sha256'];d={'source':str(src),'source_sha256':old['scene_sha256'],'scene':str(out),'scene_sha256':sha(out),'changed':'Only native relief material roughness: disconnect old fixed-Gaussian broadening and retain its source intrinsic0.475. Remove now-unused nodes/images.','intrinsic_coefficient_unchanged':.475,'full_model_qualified':False,'baseline_untouched':True};(R/'receipts/intrinsic_build.json').write_text(json.dumps(d,indent=2));print('INTRINSIC_CONTROL_SAVED',d['scene_sha256'])
