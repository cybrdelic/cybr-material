"""Bind only lossless RGB16 normal encoding; retain optical inputs and geometry."""
import bpy,numpy as np,json,sys,hashlib,struct,zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];EXP=ROOT.parent
sys.path.insert(0,str(EXP/"opaque_surface_guide_contract/source"))
from surface_guides import configure
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def decode(p):
 data=p.read_bytes();pos=8;zs=[]
 while pos<len(data):
  n=struct.unpack(">I",data[pos:pos+4])[0];name=data[pos+4:pos+8];v=data[pos+8:pos+8+n];pos+=12+n
  if name==b"IDAT":zs.append(v)
 raw=zlib.decompress(b"".join(zs));n=4096;stride=n*6+1
 a=np.frombuffer(raw,np.uint8).reshape(n,stride);assert np.all(a[:,0]==0)
 return a[:,1:].copy().view(">u2").reshape(n,n,3)
def geom():
 d={}
 for ob in bpy.context.scene.objects:
  if ob.type=="MESH":
   v=np.empty(len(ob.data.vertices)*3,"f4");ob.data.vertices.foreach_get("co",v)
   l=np.empty(len(ob.data.loops),"i4");ob.data.loops.foreach_get("vertex_index",l)
   d[ob.name]=[hashlib.sha256(v.tobytes()+l.tobytes()).hexdigest(),[list(row) for row in ob.matrix_world]]
 return d
rows=json.loads((ROOT/"receipts/native16_maps.json").read_text())
for row in rows:
 cid=row["material_id"];is_chain=cid.startswith("14")
 base=EXP/("chainmail_surface_guides" if is_chain else "wire_loop_uv_fix")
 old=json.loads((base/"experiment.json").read_text())
 source=EXP.parents[2]/old["source_scene"];assert sha(source)==old["source_sha256"]
 path=Path(row["path"]);assert sha(path)==row["sha256"]
 bpy.ops.wm.open_mainfile(filepath=str(source));bpy.context.preferences.filepaths.save_version=0;s=bpy.context.scene
 previous=geom();bindings=[];newim=bpy.data.images.load(str(path),check_existing=False);newim.colorspace_settings.name="Non-Color"
 for mat in bpy.data.materials:
  if mat.use_nodes:
   for node in mat.node_tree.nodes:
    if node.type=="TEX_IMAGE" and node.image and Path(node.image.filepath).stem=="Normal_OpenGL":
     bindings.append(dict(material=mat.name,node=node.name,old_image=node.image.filepath));node.image=newim
 assert len(bindings)==1,bindings
 if not is_chain:configure(s)
 assert geom()==previous
 for im in list(bpy.data.images):
  if not im.users:bpy.data.images.remove(im)
 encoded=decode(path);rng=np.random.default_rng(64016);maximum=0.
 for y,x in rng.integers(0,4096,(96,2)):
  actual=np.array(newim.pixels[(int(y)*4096+int(x))*4:(int(y)*4096+int(x))*4+3]);expected=encoded[4095-y,x].astype(float)/65535
  maximum=max(maximum,float(np.max(abs(actual-expected))))
 assert maximum<1e-7,maximum
 out=ROOT/"scenes"/(cid+"_uv_native16_guides.blend")
 bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert sha(source)==old["source_sha256"]
 eid="chainmail_uv_normal16" if is_chain else "wire_uv_normal16";eroot=EXP/eid;eroot.mkdir(exist_ok=True)
 manifest=old.copy();manifest.update(experiment_id=eid,source_scene=str(out.relative_to(EXP.parents[2])),source_sha256=sha(out),purpose="UV chart correction plus unchanged analytic normal field encoded RGB16; opaque first-hit guides preserve actual shading-normal output")
 manifest["dependencies"]=[d for d in old["dependencies"] if not d["path"].endswith("/Normal_OpenGL.png")]
 manifest["dependencies"].append(dict(path=str(path.relative_to(EXP.parents[2])),sha256=sha(path)))
 manifest["capture"].update(guide_albedo_pass="OIDN_SurfaceReflectivity",guide_normal_pass="OIDN_SurfaceNormal")
 manifest["print_on_image"]=("Chainmail" if is_chain else "Woven wire")+" ·18 mm detail · UV /16-bit normal correction"
 manifest["pass_gates"].update(geometry_equal=True,normal_field_unchanged=True,native_RGB16_read_error=maximum,visual_acceptance=False)
 mp=eroot/"experiment.json";mp.write_text(json.dumps(manifest,indent=2))
 hand=dict(schema_version=1,ready_for_capture=True,experiment_id=eid,experiment_manifest=str(mp),experiment_manifest_sha256=sha(mp),source_sha256=sha(out),requested_views=[old["capture"]["camera"]],capture_limit=1,guide_passes=["OIDN_SurfaceReflectivity","OIDN_SurfaceNormal"],status="Verified encoding correction, visual review pending; no roughness or geometry change")
 (eroot/"capture_handoff.json").write_text(json.dumps(hand,indent=2))
 report=dict(source=str(source),source_sha256=old["source_sha256"],scene=str(out),scene_sha256=sha(out),geometry_equal=True,changed_normal_binding=bindings,native16_read_max_error=maximum,normal_field=row,source_frozen=True,handoff=str(eroot/"capture_handoff.json"))
 (ROOT/"receipts"/(cid+"_binding.json")).write_text(json.dumps(report,indent=2));print("NORMAL16_READY",json.dumps(report),flush=True)

