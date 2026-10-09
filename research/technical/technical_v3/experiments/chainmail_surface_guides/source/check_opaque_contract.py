import bpy,json,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];p=R/'scenes/14_chainmail_loop_uv_declared_guides.blend';expected='6fd0c36241bac42983571f10a6c1280dfab296e2bb485c6c176e816a41bf9d0b'
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
bpy.ops.wm.open_mainfile(filepath=str(p));s=bpy.context.scene;rows=[]
for m in bpy.data.materials:
 if not m.use_nodes:continue
 for n in m.node_tree.nodes:
  if n.type!='BSDF_PRINCIPLED':continue
  a=n.inputs['Alpha'];assert not a.is_linked and a.default_value==1,(m.name,'alpha')
  ns=n.inputs['Normal'];src=ns.links[0].from_node if ns.is_linked else None
  if src:assert src.type in {'NORMAL_MAP','BUMP'} and ns.links[0].from_socket.name=='Normal'
  out=m.node_tree.nodes['OIDN_SurfaceNormal'];aov=out.inputs['Color'].links[0].from_socket
  assert (aov==ns.links[0].from_socket) if src else aov.node.type=='NEW_GEOMETRY' and aov.name=='Normal'
  rows.append(dict(material=m.name,alpha=1,alpha_linked=False,normal_source=src.type if src else 'Geometry Normal',feature_matches_actual_shading_normal=True))
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
r=dict(source_sha256=expected,passed=True,materials=rows,scene_unchanged=True)
(R/'receipts/opaque_contract.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
