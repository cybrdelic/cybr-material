"""Repair Blender4.3 AOV.aov_name; preserve all BSDF inputs and original sources."""
import bpy,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];EXP=ROOT.parent;WORK=EXP.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def value(v):
 try:return list(v)
 except TypeError:return v
def beauty():
 out={}
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  nodes={n.name for n in m.node_tree.nodes if n.type!="OUTPUT_AOV"}
  out[m.name]=dict(nodes={n.name:dict(type=n.type,inputs=[(i.identifier,value(i.default_value)) for i in n.inputs if hasattr(i,"default_value")]) for n in m.node_tree.nodes if n.name in nodes},links=sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links if l.from_node.name in nodes and l.to_node.name in nodes))
 return out
pairs=[("chainmail_surface_guides","chainmail_guide_names_r2"),("chainmail_uv_normal16","chainmail_normal16_guides_r2"),("wire_uv_normal16","wire_normal16_guides_r2")]
reports=[]
for oldid,newid in pairs:
 old=json.loads((EXP/oldid/"experiment.json").read_text());source=WORK/old["source_scene"]
 assert sha(source)==old["source_sha256"]
 bpy.ops.wm.open_mainfile(filepath=str(source));bpy.context.preferences.filepaths.save_version=0;s=bpy.context.scene;before=beauty();records=[]
 expected={"OIDN_SurfaceReflectivity","OIDN_SurfaceNormal"}
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  for n in m.node_tree.nodes:
   if n.type=="OUTPUT_AOV":
    assert n.name in expected,(m.name,n.name)
    assert hasattr(n,"aov_name"),"Expected explicit4.3 AOV RNA"
    oldname=n.aov_name;n.aov_name=n.name
    assert n.aov_name in {a.name for v in s.view_layers for a in v.aovs}
    assert n.inputs["Color"].is_linked or max(n.inputs["Color"].default_value[:3])>0
    records.append(dict(material=m.name,node=n.name,old_aov_name=oldname,new_aov_name=n.aov_name))
 assert records and beauty()==before
 eroot=EXP/newid;(eroot/"scenes").mkdir(parents=True,exist_ok=True);(eroot/"receipts").mkdir(exist_ok=True)
 path=eroot/"scenes"/(newid+".blend");bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True)
 assert sha(source)==old["source_sha256"]
 manifest=old.copy();manifest.update(experiment_id=newid,source_scene=str(path.relative_to(WORK)),source_sha256=sha(path))
 manifest["purpose"]+="; explicit Blender4.3 AOV pass names fixed"
 helper=ROOT/"source/surface_guides_v2.py";manifest["dependencies"].append(dict(path=str(helper.relative_to(WORK)),sha256=sha(helper)))
 manifest["pass_gates"].update(AOV_names_match_view_layer=True,beauty_closures_unchanged=True,nonzero_input_connections=True,actual_feature_render_pending=True)
 mp=eroot/"experiment.json";mp.write_text(json.dumps(manifest,indent=2))
 report=dict(source=str(source),source_sha256=old["source_sha256"],scene=str(path),scene_sha256=sha(path),beauty_closures_unchanged=True,nodes=records,scope="AOV pass-name binding only",remaining="Rendered feature ranges, nonzero coverage, same-beauty tolerance and visual detail gates")
 (eroot/"receipts/aov_names.json").write_text(json.dumps(report,indent=2));reports.append(report)
 hand=dict(schema_version=1,ready_for_capture=True,experiment_id=newid,experiment_manifest=str(mp),experiment_manifest_sha256=sha(mp),source_sha256=sha(path),requested_views=[manifest["capture"]["camera"]],capture_limit=1,guide_passes=["OIDN_SurfaceReflectivity","OIDN_SurfaceNormal"],required_gates=["Actual nonzero surface-guide coverage, not range alone","Raw beauty tolerance comparison; disclose measured non-bit-identical roundoff","Fine detail visual review"])
 (eroot/"capture_handoff.json").write_text(json.dumps(hand,indent=2))
 # Fresh reopen confirms actual serialized RNA, including each pass name.
 bpy.ops.wm.open_mainfile(filepath=str(path))
 for m in bpy.data.materials:
  if m.use_nodes:
   for n in m.node_tree.nodes:
    if n.type=="OUTPUT_AOV":assert n.aov_name==n.name and n.aov_name in expected
 assert beauty()==before
 print("AOV_NAME_FIX_READY",newid,sha(path),flush=True)
(ROOT/"aov_names_fixed.json").write_text(json.dumps(reports,indent=2))

