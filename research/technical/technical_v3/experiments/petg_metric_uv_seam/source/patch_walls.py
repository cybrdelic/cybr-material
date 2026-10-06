"""Fix the wrapping wall texture chart without changing printed geometry/optics."""
import bpy,numpy as np,json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=Path("/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review")
MANIFEST=json.loads((ARCHIVE/"SHA256_MANIFEST.json").read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
reports=[]
for cid,filename in [("22_petg","22_petg_r6_seamproof.blend"),("23_petg_transparent","23_petg_transparent_r7_backlit_control.blend")]:
 source=ARCHIVE/"scenes"/filename;assert sha(source)==MANIFEST[str(source.relative_to(ARCHIVE))]
 bpy.ops.wm.open_mainfile(filepath=str(source));bpy.context.preferences.filepaths.save_version=0
 ob=bpy.context.scene.objects["PETG / unified printed wall rim and solid floor"];me=ob.data
 nt=160;nz=1040;ni=992;side_faces=(nz+ni)*nt;assert len(me.vertices)==(nz+ni+2)*nt+2
 v=np.empty(len(me.vertices)*3,"f4");me.vertices.foreach_get("co",v)
 l=np.empty(len(me.loops),"i4");me.loops.foreach_get("vertex_index",l)
 before=np.empty(len(me.loops)*2,"f4");me.uv_layers.active.data.foreach_get("uv",before)
 changed=0;oldspan=[];newspan=[];circumference=math.tau*.014/.02
 for f in list(me.polygons)[:side_faces]:
  if f.index%nt!=nt-1:continue
  oldspan.append(max(me.uv_layers.active.data[k].uv.x for k in f.loop_indices)-min(me.uv_layers.active.data[k].uv.x for k in f.loop_indices))
  for k in f.loop_indices:
   vertex=me.loops[k].vertex_index
   if vertex%nt==0:me.uv_layers.active.data[k].uv.x=circumference;changed+=1
  newspan.append(max(me.uv_layers.active.data[k].uv.x for k in f.loop_indices)-min(me.uv_layers.active.data[k].uv.x for k in f.loop_indices))
 after=np.empty_like(before);me.uv_layers.active.data.foreach_get("uv",after)
 assert np.array_equal(before.reshape(-1,2)[:,1],after.reshape(-1,2)[:,1])
 assert np.allclose(newspan,circumference/nt,atol=1e-6)
 v2=np.empty_like(v);me.vertices.foreach_get("co",v2);l2=np.empty_like(l);me.loops.foreach_get("vertex_index",l2)
 assert np.array_equal(v,v2) and np.array_equal(l,l2)
 out=ROOT/"scenes"/(cid+"_wall_uv_fixed.blend");bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(source)==MANIFEST[str(source.relative_to(ARCHIVE))]
 row=dict(material_id=cid,source=str(source),source_sha256=sha(source),scene=str(out),scene_sha256=sha(out),vertex_topology_bytes_equal=True,shader_image_bytes_unchanged=True,wall_thickness_m=float(ob["wall_thickness_m"]),layer_pitch_m=float(ob["layer_height_m"]),changed_loop_u_count=changed,changed_face_count=len(oldspan),old_U_span=float(oldspan[0]),new_U_span=float(newspan[0]),desired_U_span=circumference/nt,all_V_coordinates_bit_exact=True,scope="Angular wrap interpolation on outer/inner walls only. Existing physical U scale retained.",not_changed="No IOR, transmission, absorption, roughness, wall, bead, seam or light/camera modification.",remaining=["Rim/cap tangent chart is a separate modeling contract","Old packed source maps retained for controlled patch; this is not the later native4K slab or a ready user panel","Rough printed resin still has no internal-void/interdiffusion scattering model"],ready_for_capture=False)
 reports.append(row);print("PETG_WALL_UV_FIXED",cid,row["new_U_span"],flush=True)
(ROOT/"receipts/wall_uv_fix.json").write_text(json.dumps(reports,indent=2))

