"""Reconstruct the selected analytic LCD/CRT panels without any optics change."""
import bpy,sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=Path("/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review")
STUDIO=Path("/workspace/scratch/b4387906eb93/material-publication/audit/checkpoints/2026-10-05/sources/material-slabs/source")
sys.path.insert(0,str(STUDIO))
from studio import configure
hashes=json.loads((ARCHIVE/"SHA256_MANIFEST.json").read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for cid in ["18_lcd","19_crt"]:
    source=ARCHIVE/"scenes"/(cid+"_r5_hero.blend")
    assert sha(source)==hashes[str(source.relative_to(ARCHIVE))]
    bpy.ops.wm.read_factory_settings(use_empty=True);bpy.context.preferences.filepaths.save_version=0;s=bpy.context.scene
    with bpy.data.libraries.load(str(source),link=False) as(a,b):
        b.objects=[name for name in a.objects if name.startswith(("Display /","LCD /","CRT /"))]
    obs=[]
    for ob in b.objects:
        if ob.type in {"MESH","CURVE","CURVES"}:s.collection.objects.link(ob);obs.append(ob)
    assert obs
    optical=[]
    for ob in obs:
        if ob.type=="MESH":
            for mat in ob.data.materials:
                if mat and mat.use_nodes:
                    for n in mat.node_tree.nodes:
                        assert not (n.type=="TEX_IMAGE" and n.image),"Unexpected proxy image"
                        if n.type=="BSDF_PRINCIPLED":
                            optical.append(dict(object=ob.name,material=mat.name,ior=float(n.inputs["IOR"].default_value),roughness=float(n.inputs["Roughness"].default_value),transmission=float(n.inputs["Transmission Weight"].default_value),emission_strength=float(n.inputs["Emission Strength"].default_value),base_color=list(n.inputs["Base Color"].default_value),emission_color=list(n.inputs["Emission Color"].default_value)))
    report=configure(cid,obs)
    for ob in s.objects:
        if ob.type=="LIGHT" and ob.name.startswith("Slab /"):ob.rotation_euler=(0,0,0)
    s["slab_studio_revision"]="r3 restored:1um ground clearance and downward area normals"
    s["capture_status"]="Unfiltered sampling route pending operator integration; emission guide failure not bypassed."
    out=ROOT/"scenes"/(cid+"_selected_optics_sampling.blend")
    bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True)
    assert sha(source)==hashes[str(source.relative_to(ARCHIVE))]
    report.update(source=str(source),source_sha256=sha(source),scene=str(out),scene_sha256=sha(out),optical_parameters=optical,geometry="Exact selected source objects appended without edits",shaders="Exact selected analytic emission and separate coverglass",image_dependencies=0,quality_acceptance=False,capture_ready=False,blocker="Operator must use explicit unfiltered capture mode; current guided OIDN albedo is out of contract.")
    (ROOT/"receipts"/(cid+".json")).write_text(json.dumps(report,indent=2))
    print("DISPLAY_SOURCE_PREPARED",cid,report["scene_sha256"],flush=True)

