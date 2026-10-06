"""Official OIDN first-hit feature definition; optical closures are untouched."""
import bpy,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=EXP.parents[2]
old=json.loads((EXP/"lcd_unfiltered_sampling/experiment.json").read_text())
base=json.loads((EXP/"display_sampling_recovery/receipts/19_crt.json").read_text())
old.update(material_id="19_crt",source_scene=str(Path(base["scene"]).relative_to(WORK)),source_sha256=base["scene_sha256"],physical_dimensions_m=base["dimensions_m"])
src=WORK/old["source_scene"]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def value(x):
 try:return list(x)
 except TypeError:return x
def state(allowed=None):
 result={}
 for m in bpy.data.materials:
  if not m.use_nodes:continue
  names=set(allowed[m.name]) if allowed else {n.name for n in m.node_tree.nodes}
  result[m.name]=dict(nodes={n.name:dict(type=n.bl_idname,inputs=[(i.identifier,value(i.default_value)) for i in n.inputs if hasattr(i,"default_value")]) for n in m.node_tree.nodes if n.name in names},links=sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links if l.from_node.name in names and l.to_node.name in names))
 return result
assert sha(src)==old["source_sha256"]
bpy.ops.wm.open_mainfile(filepath=str(src));bpy.context.preferences.filepaths.save_version=0;scene=bpy.context.scene
before=state();allowed={k:list(v["nodes"]) for k,v in before.items()}
transforms={o.name:[list(row) for row in o.matrix_world] for o in scene.objects}
records=[]
for m in bpy.data.materials:
 if not m.use_nodes:continue
 nd=m.node_tree.nodes;lk=m.node_tree.links
 outputs=[n for n in nd if n.type=="OUTPUT_MATERIAL" and n.is_active_output];assert len(outputs)==1
 output=outputs[0];assert not output.inputs["Volume"].is_linked
 p=output.inputs["Surface"].links[0].from_node;assert p.type=="BSDF_PRINCIPLED"
 assert not p.inputs["Alpha"].is_linked and p.inputs["Alpha"].default_value==1
 tw=p.inputs["Transmission Weight"];assert not tw.is_linked;transmission=float(tw.default_value);assert 0<=transmission<=1
 base=p.inputs["Base Color"];assert not base.is_linked and min(base.default_value[:3])>=0 and max(base.default_value[:3])<=1
 albedo=nd.new("ShaderNodeOutputAOV");albedo.name="OIDN_SurfaceReflectivity";albedo.aov_name=albedo.name
 albedo.inputs["Color"].default_value=(1,1,1,1) if transmission>0 else base.default_value
 normal=nd.new("ShaderNodeOutputAOV");normal.name="OIDN_SurfaceNormal";normal.aov_name=normal.name
 if p.inputs["Normal"].is_linked:
  n=p.inputs["Normal"].links[0].from_socket;assert n.node.type in {"NORMAL_MAP","BUMP"} and n.name=="Normal";lk.new(n,normal.inputs["Color"]);normal_source=n.node.type
 else:
  g=nd.new("ShaderNodeNewGeometry");g.name="First-hit feature normal only";lk.new(g.outputs["Normal"],normal.inputs["Color"]);normal_source="Geometry shading normal"
 records.append(dict(material=m.name,transmission=transmission,albedo_definition="Unit glass feature prescribed by OIDN" if transmission>0 else "Existing intrinsic Base Color; emitted radiance excluded",albedo=list(albedo.inputs["Color"].default_value),normal=normal_source,emission_strength_unchanged=float(p.inputs["Emission Strength"].default_value)))
for vl in scene.view_layers:
 for name in ["OIDN_SurfaceReflectivity","OIDN_SurfaceNormal"]:
  a=vl.aovs.add();a.name=name;a.type="COLOR"
assert state(allowed)==before and transforms=={o.name:[list(row) for row in o.matrix_world] for o in scene.objects}
out=R/"scenes/crt_first_hit_guides_for_remainder.blend";bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True)
assert sha(src)==old["source_sha256"]
bpy.ops.wm.open_mainfile(filepath=str(out))
for m in bpy.data.materials:
 if m.use_nodes:
  for n in m.node_tree.nodes:
   if n.type=="OUTPUT_AOV":assert n.aov_name==n.name
assert state(allowed)==before
report=dict(source=str(src),source_sha256=old["source_sha256"],scene=str(out),scene_sha256=sha(out),all_original_optical_nodes_inputs_links_equal=True,all_object_transforms_equal=True,materials=records,documentation="https://www.openimagedenoise.org/documentation.html#albedos",definition="First-hit glass auxiliary albedo1; ordinary intrinsic color on opaque surfaces. This is a documented guide proxy, not the glass's physical reflectance. Actual shading normals accompany it.",limitation="Official documentation warns that first-hit features do not carry details seen through transparent surfaces. Bounds do not establish preservation of emitted RGB stripes; actual512-sample pixel comparison remains mandatory.",noise_or_realism_acceptance=False)
(R/"receipts/source_contract.json").write_text(json.dumps(report,indent=2))
manifest=old.copy();manifest.update(experiment_id="crt_transport_split",source_scene=str(out.relative_to(WORK)),source_sha256=sha(out),purpose="Preserve CRT transmission/emission and denoise only remaining surface transport; unchanged optical source")
manifest["capture"].update(guide_albedo_pass="OIDN_SurfaceReflectivity",guide_normal_pass="OIDN_SurfaceNormal")
manifest["pass_gates"].update(all_optical_nodes_equal=True,source_AOV_names_verified=True,actual_nonzero_range_check_pending=True,raw_beauty_comparison_pending=True,stripe_retention_pending=True)
manifest["print_on_image"]="CRT · 54.4mm panel · retained phosphor structure"
mp=R/"experiment.json";mp.write_text(json.dumps(manifest,indent=2))
handoff=dict(ready_for_source_review=True,experiment_id="crt_transport_split",source_sha256=sha(out),source_scene=str(out),optical_source_sha256=old["source_sha256"],geometry_and_optics_unchanged=True,guide_use="Remainder only. Never filter transmitted phosphor structure.",capture_released=False,physical_dimensions_m=old["physical_dimensions_m"],comparison_scope="Fresh source node/transform invariance plus identity pass reconstruction; no prior512 CRT control exists")
(R/"capture_handoff.json").write_text(json.dumps(handoff,indent=2));print("CRT_COMPONENT_SOURCE_READY",sha(out),flush=True)
