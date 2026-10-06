"""Resolve the already-authored circular wires; no new finish or centerline."""
import bpy,bmesh,numpy as np,json,hashlib,math
from pathlib import Path
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=EXP.parents[2]
old=json.loads((EXP/"wire_normal16_guides_r2/experiment.json").read_text());src=WORK/old["source_scene"]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(src)==old["source_sha256"]
bpy.ops.wm.open_mainfile(filepath=str(src));bpy.context.preferences.filepaths.save_version=0;s=bpy.context.scene
ob=s.objects["Lattice / actual strand network"];me=ob.data;coords=np.empty(len(me.vertices)*3,"f4");me.vertices.foreach_get("co",coords)
q=coords.reshape(28,157,8,3);C=q.astype(float).mean(2)
theta=np.arange(8)*math.tau/8;A=2*np.mean((q-C[:,:,None])*np.cos(theta)[None,None,:,None],axis=2);B=2*np.mean((q-C[:,:,None])*np.sin(theta)[None,None,:,None],axis=2)
# Least-squares first harmonic recovers the unchanged circular section frame.
pred=C[:,:,None]+A[:,:,None]*np.cos(theta)[None,None,:,None]+B[:,:,None]*np.sin(theta)[None,None,:,None]
fit_error=float(np.max(np.linalg.norm(pred-q,axis=-1)));assert fit_error<1e-8
assert np.max(abs(np.linalg.norm(A,axis=-1)-.00042))<1e-8
assert np.max(abs(np.linalg.norm(B,axis=-1)-.00042))<1e-8
NS=32;th=np.arange(NS)*math.tau/NS
points=(C[:,:,None]+A[:,:,None]*np.cos(th)[None,None,:,None]+B[:,:,None]*np.sin(th)[None,None,:,None]).astype("f4")
points[:,:,::4]=q  # Every pre-existing circle sample remains bit-exact.
verts=points.reshape(-1,3);faces=[];uv=[];smooth=[]
for c in range(28):
 old_face0=c*(156*8+2)
 vs=[float(me.uv_layers.active.data[me.polygons[old_face0+i*8].loop_indices[0]].uv.y) for i in range(156)]+[1.]
 start=c*157*NS
 for i in range(156):
  for j in range(NS):
   faces.append((start+i*NS+j,start+i*NS+(j+1)%NS,start+(i+1)*NS+(j+1)%NS,start+(i+1)*NS+j))
   uv.extend([(j/NS,vs[i]),((j+1)/NS,vs[i]),((j+1)/NS,vs[i+1]),(j/NS,vs[i+1])]);smooth.append(True)
 disk=np.stack([.5+.5*np.cos(th),.5+.5*np.sin(th)],axis=1)
 faces.append(tuple(start+j for j in reversed(range(NS))));uv.extend(disk[::-1]);smooth.append(bool(me.polygons[old_face0+156*8].use_smooth))
 faces.append(tuple(start+156*NS+j for j in range(NS)));uv.extend(disk);smooth.append(bool(me.polygons[old_face0+156*8+1].use_smooth))
data=bpy.data.meshes.new("Woven wire /32-sided circle convergence");data.from_pydata(verts,[],faces);data.update()
for m in me.materials:data.materials.append(m)
layer=data.uv_layers.new(name="UVMap");layer.active_render=True;layer.data.foreach_set("uv",np.asarray(uv,"f4").ravel())
for f,flag in zip(data.polygons,smooth):f.use_smooth=flag
ob.data=data
actual=np.empty(len(data.vertices)*3,"f4");data.vertices.foreach_get("co",actual);actual=actual.reshape(28,157,NS,3)
assert np.array_equal(actual[:,:,::4],q)
center_error=float(np.max(np.linalg.norm(actual.astype(float).mean(2)-C,axis=-1)));assert center_error<1e-8
bm=bmesh.new();bm.from_mesh(data);bm.normal_update();assert all(e.is_manifold and e.is_contiguous for e in bm.edges);assert bm.calc_volume(signed=True)>0;bm.free()
# Geometry target/error is specified before capture, not chosen from appearance.
rad=.00042;errors={str(n):rad*(1-math.cos(math.pi/n)) for n in [8,16,32,64]}
out=R/"scenes/wire_circle32_native16.blend";bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(src)==old["source_sha256"]
report=dict(source=str(src),source_sha256=sha(src),scene=str(out),scene_sha256=sha(out),source_circle_fit_error_m=fit_error,all_old_vertices_bit_exact=True,centerline_error_m=center_error,radius_m=rad,side_counts=[8,32],maximum_inscribed_polygon_radial_deficit_m=errors,macro_FOV_m=.018,old_deficit_pixels_at960=errors["8"]/.018*960,new_deficit_pixels_at960=errors["32"]/.018*960,closed_outward=True,material_and_maps_unchanged=True,all_existing_light_camera_object_transforms_unchanged=True,scope="Numerical representation of the documented circular wire; same authored crimp and material response",not_claimed="New manufacturing, force equilibrium, wear or measured reflectance",visual_acceptance=False)
(R/"receipts/build.json").write_text(json.dumps(report,indent=2))
manifest=old.copy();manifest.update(experiment_id="wire_circle_resolution",source_scene=str(out.relative_to(WORK)),source_sha256=sha(out),purpose="Reduce octagonal-wire shape error under the same18mm capture, retaining exact original samples/centerlines/radius and all optical parameters")
manifest["pass_gates"].update(old_vertices_preserved=True,centerline_error_m=center_error,geometric_error_target_pixels=.15,actual_geometric_error_pixels=report["new_deficit_pixels_at960"],visual_acceptance=False)
manifest["print_on_image"]="Woven wire ·18mm detail · circular-section refinement"
mp=R/"experiment.json";mp.write_text(json.dumps(manifest,indent=2))
(R/"capture_handoff.json").write_text(json.dumps(dict(schema_version=1,ready_for_capture=True,experiment_id="wire_circle_resolution",experiment_manifest=str(mp),experiment_manifest_sha256=sha(mp),source_sha256=sha(out),requested_views=[old["capture"]["camera"]],capture_limit=1,matched_control=old["source_scene"],hypothesis="8-sided circular approximation contributes visible flattening;32 sides reduce projected shape error below0.15px without optical changes"),indent=2))
print("WIRE_CIRCLE_REFINEMENT_READY",json.dumps(report),flush=True)

