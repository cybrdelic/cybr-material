"""Append exact selected r6 links; only corrected face-corner UV differs."""
import bpy,numpy as np,json,hashlib,sys,shutil
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=Path("/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review")
SOURCE=ARCHIVE/"scenes/12_lattice_wire_r4_hero.blend"
EXPECTED="d3f30f120f8a0e7e8a70357df08ac04688868ff960cae27c7cac11fa01bf21cb"
STUDIO=Path("/workspace/scratch/b4387906eb93/material-publication/audit/checkpoints/2026-10-05/sources/material-slabs/source")
sys.path[:0]=[str(ROOT/"source"),str(STUDIO)]
from loop_chart import open_tube_chart
from studio import configure
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def arrays(me):
    v=np.empty(len(me.vertices)*3,"f4");me.vertices.foreach_get("co",v)
    l=np.empty(len(me.loops),"i4");me.loops.foreach_get("vertex_index",l)
    u=np.empty(len(me.loops)*2,"f4");me.uv_layers.active.data.foreach_get("uv",u)
    return v,l,u
def digest(x):return hashlib.sha256(x.tobytes()).hexdigest()
assert sha(SOURCE)==EXPECTED
maps=json.loads((ROOT/"receipts/native_maps.json").read_text())["maps"]
for row in maps:assert sha(row["path"])==row["sha256"]
bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.context.preferences.filepaths.save_version=0
s=bpy.context.scene
with bpy.data.libraries.load(str(SOURCE),link=False) as (a,b):
    b.objects=[name for name in a.objects if name.startswith(("Lattice /","Fixture /"))]
obs=[]
for ob in b.objects:s.collection.objects.link(ob);obs.append(ob)
chain=next(o for o in obs if o.name.startswith("Lattice /"));me=chain.data
v,l,uv=arrays(me)
assert len(me.vertices)==28*157*8 and len(me.polygons)==28*(156*8+2)
images={}
for ob in obs:
    for m in getattr(ob.data,"materials",[]):
        if m and m.use_nodes:
            for node in m.node_tree.nodes:
                if node.type=="TEX_IMAGE" and node.image:
                    name=Path(node.image.filepath).stem
                    if name not in images:
                        path=ROOT/"maps"/(name+".png");assert path.exists()
                        im=bpy.data.images.load(str(path));im.colorspace_settings.name=node.image.colorspace_settings.name
                        images[name]=im
                    node.image=images[name]
for im in list(bpy.data.images):
    if not im.users:bpy.data.images.remove(im)
report=configure("12_lattice_wire",obs)
for ob in s.objects:
    if ob.type=="LIGHT" and ob.name.startswith("Slab /"):ob.rotation_euler=(0,0,0)
s["slab_studio_revision"]="r3 restored: 1um ground clearance and downward area normals"
bpy.ops.object.camera_add(location=(.024,-.026,.030))
camera=bpy.context.object;camera.name="Lattice / seam inspection";camera.data.type="ORTHO";camera.data.ortho_scale=.018;camera.data.clip_start=.00001
target=Vector((.002,-.002,.003));camera.rotation_euler=(target-camera.location).to_track_quat("-Z","Y").to_euler()
s.camera=bpy.data.objects.get("Slab / shared camera")
control=ROOT/"scenes/12_lattice_wire_native4k_recovered_control.blend"
bpy.ops.wm.save_as_mainfile(filepath=str(control),compress=True)
centers=v.reshape(28,157,8,3).mean(axis=2)
charts=[open_tube_chart(points,8) for points in centers]
new=np.concatenate([a for a,length in charts]).reshape(-1)
me.uv_layers.active.data.foreach_set("uv",new);me.uv_layers.active.active_render=True;me.update()
v2,l2,u2=arrays(me)
assert np.array_equal(v,v2) and np.array_equal(l,l2) and np.array_equal(new,u2)
chain["UV_revision"]="r2: corner-domain circumference; normalized true arc length; planar cut caps. Exact original positions and finish recipe."
candidate=ROOT/"scenes/12_lattice_wire_native4k_loop_uv_r2.blend"
bpy.ops.wm.save_as_mainfile(filepath=str(candidate),compress=True)
assert sha(SOURCE)==EXPECTED
report.update(source_scene=str(SOURCE),source_sha256=EXPECTED,control=str(control),control_sha256=sha(control),candidate=str(candidate),candidate_sha256=sha(candidate),geometry_sha256=digest(v),topology_sha256=digest(l),old_uv_sha256=digest(uv),new_uv_sha256=digest(u2),geometry_bit_exact=True,topology_bit_exact=True,changed="Per-loop UV chart only between control and candidate",known_remaining="Absolute60mm finish scale on a curved wire is unqualified; these authored texture fields are not a wire-drawing process model.",source_recovery="Source r4 scene byte-exact; native4K scene reconstructed from frozen builder contract and same native generator. Not claimed byte-exact old slab recovery.")
(ROOT/"receipts/build.json").write_text(json.dumps(report,indent=2))
dep=[dict(path=str(ROOT/"maps"/(name+".png")),sha256=sha(ROOT/"maps"/(name+".png"))) for name in images]
for eid,scene in [("wire_loop_uv_control",control),("wire_loop_uv_fix",candidate)]:
    eroot=ROOT.parent/eid;eroot.mkdir(exist_ok=True)
    manifest=dict(schema_version=1,selected=False,quality_acceptance=False,experiment_id=eid,material_id="12_lattice_wire",source_scene=str(scene.relative_to(ROOT.parents[3])),source_sha256=sha(scene),dependencies=[dict(path=str(Path(d["path"]).relative_to(ROOT.parents[3])),sha256=d["sha256"]) for d in dep],capture=dict(camera="Lattice / seam inspection",engine="CYCLES",width=960,height=960,samples=128,light_tree=False,save_guides=True,color_depth=16,dither_intensity=0,rss_cap_mib=2048),contract=dict(map_resolution=[4096,4096]),purpose="Open-wire seam and cut-cap UV correction with exact selected geometry and unchanged finish recipe.",pass_gates=dict(geometry_bit_exact=True,topology_bit_exact=True,native_source512_quantized_equivalence=True,loop_chart_tests_pass=True,visual_acceptance=False),physical_dimensions_m=report["dimensions_m"],print_on_image="Woven wire · 18 mm detail · UV seam audit")
    (eroot/"experiment.json").write_text(json.dumps(manifest,indent=2))
hand=dict(schema_version=1,ready_for_capture=True,experiment_id="wire_loop_uv_fix",experiment_manifest=str(ROOT/"experiment.json"),experiment_manifest_sha256=sha(ROOT/"experiment.json"),source_sha256=sha(candidate),requested_views=["Lattice / seam inspection"],capture_limit=1,comparison_control=str(control),comparison_control_sha256=sha(control),preflight=str(ROOT/"receipts/build.json"),visual_claim="UV correctness candidate; finish-scale model remains unqualified")
(ROOT/"capture_handoff.json").write_text(json.dumps(hand,indent=2))
print("CHAINMAIL_READY",json.dumps(hand),flush=True)

