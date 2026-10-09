"""One bounded actual-scene transport probe. Renderer owns execution."""
import bpy,numpy as np,json,hashlib,sys,time
from pathlib import Path
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=EXP.parents[2]
contract=json.loads((R/"experiment.json").read_text());source=WORK/contract["source_scene"]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(source)==contract["source_sha256"]
args=sys.argv[sys.argv.index("--")+1:];out=Path(args[0]).resolve();assert not out.exists(),"Fresh output directory required";out.mkdir(parents=True)
bpy.ops.wm.open_mainfile(filepath=str(source));s=bpy.context.scene;s.camera=s.objects["Slab / shared camera"]
s.render.engine="CYCLES";s.cycles.device="CPU";s.cycles.samples=512;s.cycles.use_adaptive_sampling=False;s.cycles.use_denoising=False;s.cycles.use_light_tree=False
s.render.resolution_x=s.render.resolution_y=960;s.render.resolution_percentage=100;s.render.threads_mode="FIXED";s.render.threads=4;s.render.dither_intensity=0
v=s.view_layers[0]
for kind in ["diffuse","glossy","transmission"]:
 for part in ["direct","indirect","color"]:setattr(v,"use_pass_"+kind+"_"+part,True)
v.use_pass_emit=True;v.use_pass_environment=True;v.cycles.use_pass_volume_direct=True;v.cycles.use_pass_volume_indirect=True
v.cycles.denoising_store_passes=True;v.update_render_passes()
s.use_nodes=True;s.render.use_compositing=True;nodes=s.node_tree.nodes;nodes.clear();rl=nodes.new("CompositorNodeRLayers");rl.scene=s;rl.layer=v.name
fn=nodes.new("CompositorNodeOutputFile");fn.base_path=str(out);fn.format.file_format="OPEN_EXR";fn.format.color_depth="32";fn.format.color_mode="RGB";fn.format.exr_codec="ZIP"
passes={"beauty":"Image","diffuse_direct":"DiffDir","diffuse_indirect":"DiffInd","diffuse_color":"DiffCol","glossy_direct":"GlossDir","glossy_indirect":"GlossInd","glossy_color":"GlossCol","transmission_direct":"TransDir","transmission_indirect":"TransInd","transmission_color":"TransCol","emission":"Emit","environment":"Env","volume_direct":"VolumeDir","volume_indirect":"VolumeInd","albedo":"OIDN_SurfaceReflectivity","normal":"OIDN_SurfaceNormal"}
paths={}
for i,(name,socket) in enumerate(passes.items()):
 assert socket in rl.outputs,(socket,list(rl.outputs.keys()))
 if i:fn.file_slots.new(name)
 fn.file_slots[i].path=name+"_";s.node_tree.links.new(rl.outputs[socket],fn.inputs[i]);paths[name]=out/(name+"_"+str(s.frame_current).zfill(4)+".exr")
s.render.image_settings.file_format="PNG";s.render.image_settings.color_depth="16";s.render.image_settings.color_mode="RGB";s.render.filepath=str(out/"INTERNAL_UNFILTERED_PROBE.png")
shaders_before={m.name:[(n.name,n.bl_idname,[(i.identifier,list(i.default_value) if hasattr(i.default_value,"__len__") else i.default_value) for i in n.inputs if hasattr(i,"default_value")]) for n in m.node_tree.nodes] for m in bpy.data.materials if m.use_nodes}
started=time.time();bpy.ops.render.render(write_still=True);elapsed=time.time()-started
assert sha(source)==contract["source_sha256"]
data={};stats={}
for name,path in paths.items():
 assert path.is_file()
 im=bpy.data.images.load(str(path),check_existing=False);w,h=im.size;assert [w,h]==[960,960]
 a=np.empty(w*h*4,"f4");im.pixels.foreach_get(a);a=a.reshape(h,w,4)[:,:,:3].copy();bpy.data.images.remove(im)
 assert np.isfinite(a).all(),name;data[name]=a.astype(float)
 stats[name]=dict(min=float(a.min()),max=float(a.max()),mean=float(a.mean()),sha256=sha(path))
parts={}
for kind in ["diffuse","glossy","transmission"]:parts[kind]=(data[kind+"_direct"]+data[kind+"_indirect"])*data[kind+"_color"]
parts.update(emission=data["emission"],environment=data["environment"],volume=data["volume_direct"]+data["volume_indirect"])
retained=parts["transmission"]+parts["emission"];remainder=parts["diffuse"]+parts["glossy"]+parts["environment"]+parts["volume"]
reconstructed=retained+remainder;beauty=data["beauty"];error=reconstructed-beauty;u=2.**-24
scale=np.maximum(np.sum([abs(a) for a in parts.values()],axis=0),1e-12)
metrics=dict(max_abs=float(abs(error).max()),rms=float(np.sqrt(np.mean(error*error))),relative_l2=float(np.linalg.norm(error)/np.linalg.norm(beauty)),float32_bit_equal=bool(np.array_equal(reconstructed.astype("f4"),beauty.astype("f4"))),maximum_error_units_of_float32_u=float(np.max(abs(error)/(u*scale))),nonnegative_components=bool(min(a.min() for a in parts.values())>=0),retained_signal_min=float(retained.min()),remainder_min=float(remainder.min()))
for name,a in [("retained_signal",retained),("stochastic_remainder",remainder),("reconstructed",reconstructed)]:
 im=bpy.data.images.new(name,width=960,height=960,float_buffer=True);rgba=np.ones((960,960,4),"f4");rgba[:,:,:3]=a;im.pixels.foreach_set(rgba.ravel());im.file_format="OPEN_EXR";im.filepath_raw=str(out/(name+".exr"));im.save();im.save_render(str(out/("INTERNAL_"+name+".png")),scene=s);bpy.data.images.remove(im)
report=dict(source=str(source),source_sha256=sha(source),probe_script_sha256=sha(__file__),resolution=[960,960],samples=512,adaptive_sampling=False,render_seconds=elapsed,pass_paths={k:str(v) for k,v in paths.items()},pass_statistics=stats,reconstruction=metrics,formula="Combined=(DiffDir+DiffInd)*DiffCol+(GlossDir+GlossInd)*GlossCol+(TransDir+TransInd)*TransCol+Emit+Env+VolumeDir+VolumeInd",retained="All transmission plus direct visible emission; retains the display signal through both actual glass surfaces",to_filter="Diffuse, glossy, environment and volume sum only",documentation="https://docs.blender.org/manual/en/4.3/render/layers/passes.html",status="Native matched512 pass export; no automatic filtering or noise acceptance",noise_acceptance=False,material_or_optics_changed=False)
(out/"probe_report.json").write_text(json.dumps(report,indent=2));print("LCD_TRANSPORT_PROBE",json.dumps(metrics),flush=True)

