"""Add declared first-hit features to the frozen UV-corrected chain source."""
import bpy,numpy as np,sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent/"chainmail_loop_uv_fix"
SRC=BASE/"scenes/14_chainmail_native4k_loop_uv_r2.blend"
SHA="efd85608811ea7ef2eda19eb4a945d4f5bc26fae8227a4f8d2d4e9adae50606f"
sys.path.insert(0,str(BASE/"source"))
from bounded_surface_features import configure,ALBEDO,NORMAL
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def value(a):
    try:return list(a)
    except TypeError:return a
def node_state(n):
    d=dict(type=n.bl_idname,inputs=[(i.identifier,value(i.default_value)) for i in n.inputs if hasattr(i,"default_value")])
    for prop in ["operation","blend_type","space","uv_map","interpolation","extension","distribution","subsurface_method"]:
        if hasattr(n,prop):d[prop]=getattr(n,prop)
    if n.type=="TEX_IMAGE":d["image"]=n.image.name
    return d
def state(allowed=None):
    result={}
    for m in bpy.data.materials:
        if not m.use_nodes:continue
        names=set(allowed[m.name]) if allowed else {n.name for n in m.node_tree.nodes}
        result[m.name]=dict(nodes={n.name:node_state(n) for n in m.node_tree.nodes if n.name in names},links=sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links if l.from_node.name in names and l.to_node.name in names))
    return result
assert sha(SRC)==SHA
bpy.ops.wm.open_mainfile(filepath=str(SRC));bpy.context.preferences.filepaths.save_version=0;s=bpy.context.scene
old=state();allowed={m:list(d["nodes"]) for m,d in old.items()}
transforms={o.name:list(sum(([float(x) for x in row] for row in o.matrix_world),[])) for o in s.objects}
records=configure(s)
assert state(allowed)==old
assert transforms=={o.name:list(sum(([float(x) for x in row] for row in o.matrix_world),[])) for o in s.objects}
out=ROOT/"scenes/14_chainmail_loop_uv_declared_guides.blend"
bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True)
assert sha(SRC)==SHA
report=dict(source=str(SRC),source_sha256=SHA,scene=str(out),scene_sha256=sha(out),old_material_graphs_equal=True,all_object_transforms_equal=True,albedo_pass=ALBEDO,normal_pass=NORMAL,materials=records,documentation="https://www.openimagedenoise.org/documentation.html#albedos",remaining_gate="Actual AOV bounds and pixel-exact beauty comparison against frozen128sample capture; no visual promotion")
(ROOT/"receipts/build.json").write_text(json.dumps(report,indent=2))
experiment=json.loads((BASE/"experiment.json").read_text())
experiment.update(experiment_id="chainmail_surface_guides",source_scene=str(out.relative_to(ROOT.parents[3])),source_sha256=sha(out),purpose="Documented first-hit reflectivity and actual shading-normal guides; exact same chainmail beauty closures, UV, geometry and rig")
experiment["dependencies"].extend([dict(path=str((BASE/"source/bounded_surface_features.py").relative_to(ROOT.parents[3])),sha256=sha(BASE/"source/bounded_surface_features.py"))])
experiment["capture"]["guide_albedo_pass"]=ALBEDO
experiment["capture"]["guide_normal_pass"]=NORMAL
experiment["pass_gates"].update(original_material_graphs_equal=True,object_transforms_equal=True,feature_range_gate_pending=True,retained_beauty_equality_pending=True)
manifest=ROOT/"experiment.json";manifest.write_text(json.dumps(experiment,indent=2))
handoff=dict(schema_version=1,ready_for_capture=True,experiment_id=experiment["experiment_id"],experiment_manifest=str(manifest),experiment_manifest_sha256=sha(manifest),source_sha256=sha(out),requested_views=["Chainmail / seam inspection"],capture_limit=1,retained_beauty="/workspace/scratch/b4387906eb93/recovered-production/workspace/shared/material-production/runs/chainmail_uv_capture_20261006_0448/INTERNAL_PANEL_guides/beauty_0001.exr",required_gates=["AOV reflectivity[0,1] and normal[-1,1]","Raw beauty pixel equality at same sampling before any OIDN","Visual seam/detail preservation"],guide_passes=[ALBEDO,NORMAL])
(ROOT/"capture_handoff.json").write_text(json.dumps(handoff,indent=2))
print("DECLARED_CHAIN_GUIDES_READY",json.dumps(handoff),flush=True)

