"""Direct Cycles shader-output proof. No beauty or built-in Normal pass proxy."""
import bpy,sys,json,hashlib,math,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parents[1];root=R.parents[1];variant=sys.argv[sys.argv.index('--')+1]
assert variant=='native', 'Compact source proof supports the native flat coupon only'
sources={'native':R/'CYBR_Asphalt_Binder_Refinement_Native4096_Candidate.blend'}
out=R/'aov'/variant;out.mkdir(parents=True,exist_ok=True);bpy.ops.wm.open_mainfile(filepath=str(sources[variant]));s=bpy.context.scene
if variant in {'object_transform','instances'}:
 specimen=bpy.data.objects['Asphalt / real displaced face'];specimen.scale=(1.4,.65,1.1);specimen.rotation_euler=(.12,.08,.37)
 if variant=='instances':
  specimen.scale=(.6,.9,.85);specimen.location=(-.11,0,0);copy=specimen.copy();copy.data=specimen.data;s.collection.objects.link(copy);copy.scale=(.5,.8,.7);copy.location=(.11,0,0);copy.rotation_euler=(.1,-.1,-.3)
 bpy.context.view_layer.update()
m=bpy.data.materials['Asphalt / Rolled mineral aggregate'];n=m.node_tree.nodes;l=m.node_tree.links
tex=n.new('ShaderNodeTexImage');tex.image=bpy.data.images.load(str((R/'native4096/Normal_Object_RGB16.png' if variant=='native' else R/'Normal_Object_Metric_1536_RGB16.png')));tex.image.colorspace_settings.name='Non-Color';old=n['Image Texture.003'];tex.interpolation=old.interpolation;tex.extension=old.extension
for link in old.inputs['Vector'].links:l.new(link.from_socket,tex.inputs['Vector'])
def vec(op,a,b=None):
 v=n.new('ShaderNodeVectorMath');v.operation=op
 if isinstance(a,tuple):v.inputs[0].default_value=a
 else:l.new(a,v.inputs[0])
 if b is not None:
  if isinstance(b,tuple):v.inputs[1].default_value=b
  else:l.new(b,v.inputs[1])
 return v.outputs['Value' if op in {'LENGTH','DOT_PRODUCT'} else 'Vector']
decoded=vec('ADD',vec('MULTIPLY',tex.outputs['Color'],(2,2,2)),(-1,-1,-1))
tr=n.new('ShaderNodeVectorTransform');tr.vector_type='NORMAL';tr.convert_from='OBJECT';tr.convert_to='WORLD';l.new(decoded,tr.inputs[0]);target=vec('NORMALIZE',tr.outputs[0]);actual=n['Normal Map'].outputs['Normal'];delta=vec('LENGTH',vec('SUBTRACT',actual,target))
outputs={'Actual':actual,'Target':target,'ErrorLength':delta}
vl=s.view_layers[0]
for name,socket in outputs.items():
 a=vl.aovs.add();a.name=name;a.type='VALUE' if name=='ErrorLength' else 'COLOR';o=n.new('ShaderNodeOutputAOV');o.name=name;o.aov_name=name;l.new(socket,o.inputs['Value' if name=='ErrorLength' else 'Color'])
a=vl.aovs.add();a.name='MaterialMask';a.type='VALUE';o=n.new('ShaderNodeOutputAOV');o.aov_name='MaterialMask';o.inputs['Value'].default_value=1
s.use_nodes=True;s.node_tree.nodes.clear();nodes=s.node_tree.nodes;links=s.node_tree.links;layer=nodes.new('CompositorNodeRLayers');file=nodes.new('CompositorNodeOutputFile');file.base_path=str(out);file.format.file_format='OPEN_EXR';file.format.color_depth='32';file.format.color_mode='RGBA';file.format.exr_codec='ZIP';file.file_slots.clear()
for name in [*outputs,'MaterialMask']:
 file.file_slots.new(name);file.file_slots[name].path=name+'_';links.new(layer.outputs[name],file.inputs[name])
s.camera.data.type='ORTHO';s.camera.data.ortho_scale=.55 if variant=='instances' else .249;s.camera.location=(0,0,.5);s.camera.rotation_euler=(0,0,0)
s.render.engine='CYCLES';s.cycles.device='CPU';s.cycles.samples=1;s.cycles.use_adaptive_sampling=False;s.cycles.use_denoising=False;s.cycles.seed=217;s.render.threads_mode='FIXED';s.render.threads=2;s.render.resolution_x=s.render.resolution_y=320;s.render.resolution_percentage=100;s.render.film_transparent=True;s.render.image_settings.file_format='OPEN_EXR';s.render.image_settings.color_depth='32';s.render.filepath=str(out/'beauty.exr')
bpy.ops.render.render(write_still=True)
def load(name):
 p=out/(name+'_'+str(s.frame_current).zfill(4)+'.exr');im=bpy.data.images.load(str(p),check_existing=False);a=np.empty(im.size[0]*im.size[1]*4,'f4');im.pixels.foreach_get(a);return a.reshape(im.size[1],im.size[0],4)[...,:3]
a=load('Actual').astype('f8');t=load('Target').astype('f8');mask=load('MaterialMask')[...,0]>.999;error=load('ErrorLength')[...,0][mask];an=a[mask]/np.linalg.norm(a[mask],axis=-1,keepdims=True);tn=t[mask]/np.linalg.norm(t[mask],axis=-1,keepdims=True);angle=np.degrees(2*np.arcsin(np.clip(np.linalg.norm(an-tn,axis=-1)/2,0,1)))
r={'variant':variant,'source_sha256':hashlib.sha256(sources[variant].read_bytes()).hexdigest(),'cycles_version':bpy.app.version_string,'pixels':int(mask.sum()),'resolution':[320,320],'samples':1,'shader_error_length_max':float(error.max()),'shader_error_length_mean':float(error.mean()),'aov_angular_mean_deg':float(angle.mean()),'aov_angular_p95_deg':float(np.quantile(angle,.95)),'aov_angular_max_deg':float(angle.max()),'scope':'Direct Normal Map node output versus independently decoded full-normal texture at the same actual shading samples; before BSDF grazing-normal correction','native_map_resolution':4096 if variant=='native' else 1536,'pass':bool(error.max()<2e-6) if variant.startswith('object') or variant in {'instances','native'} else None}
(out/'receipt.json').write_text(json.dumps(r,indent=2));print(json.dumps(r));
if variant.startswith('object') or variant in {'instances','native'}:assert r['pass'],r
