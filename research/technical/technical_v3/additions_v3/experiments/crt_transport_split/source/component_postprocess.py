"""Explicit component filtering bridge. Does not render or launch OIDN."""
import bpy,numpy as np,json,sys,hashlib
from pathlib import Path
args=sys.argv[sys.argv.index('--')+1:];mode=args[0];run=Path(args[1]).resolve();prior=Path(args[2]).resolve() if len(args)>2 else None
rep=json.loads((run/'probe_report.json').read_text());source=Path(rep['source']);
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(source)==rep['source_sha256'];bpy.ops.wm.open_mainfile(filepath=str(source));scene=bpy.context.scene
scene.render.image_settings.file_format='PNG';scene.render.image_settings.color_depth='16';scene.render.image_settings.color_mode='RGB';scene.render.dither_intensity=0
pipeline=Path(__file__).resolve().parents[4]/'cloud-eevee';sys.path.insert(0,str(pipeline))
from png_precision_check import compare_opaque_pngs
from oidn_contract import validate_stats

def exr(p):
 im=bpy.data.images.load(str(p),check_existing=False);w,h=im.size;a=np.empty(w*h*4,'f4');im.pixels.foreach_get(a);bpy.data.images.remove(im);a=a.reshape(h,w,4)[:,:,:3].copy();assert np.isfinite(a).all();return a

def write_pfm(p,a):
 h,w,c=a.shape;assert c==3
 with p.open('xb') as f:f.write(f'PF\n{w} {h}\n-1.0\n'.encode());f.write(a.astype('<f4').tobytes())

def read_pfm(p):
 with p.open('rb') as f:
  assert f.readline().strip()==b'PF';w,h=map(int,f.readline().split());scale=float(f.readline());assert abs(scale)==1;a=np.frombuffer(f.read(),dtype='<f4' if scale<0 else '>f4').reshape(h,w,3).copy()
 assert np.isfinite(a).all();return a

def save(a,path,linear=None):
 h,w,c=a.shape;im=bpy.data.images.new(path.stem,width=w,height=h,float_buffer=True);rgba=np.ones((h,w,4),'f4');rgba[:,:,:3]=a;im.pixels.foreach_set(rgba.ravel());im.update()
 if linear:im.file_format='OPEN_EXR';im.filepath_raw=str(linear);im.save()
 im.save_render(str(path),scene=scene);bpy.data.images.remove(im)

def stats(a):return {'min':float(a.min()),'max':float(a.max()),'mean':float(a.mean()),'finite':bool(np.isfinite(a).all()),'nonzero_pixel_fraction':float(np.mean(np.any(a!=0,axis=2))),'size':[a.shape[1],a.shape[0]]}

out=run/'component_filter';accounting=run/'linear_component_accounting.npz'
if mode=='export':
 assert not out.exists(),'Fresh filter directory required';out.mkdir();assert prior is None or prior.is_dir()
 state=np.load(accounting);t=state['transmission_emission'];r=state['physical_remainder'];delta=state['roundoff_residual'];b=exr(rep['pass_paths']['beauty']);assert np.array_equal((t.astype('f8')+r.astype('f8')+delta).astype('f4'),b)
 assert r.min()>0,'Only strictly positive physical remainder may enter OIDN'
 if prior is not None:
  previous=json.loads((prior/'INTERNAL_PANEL.png.json').read_text());old=exr(previous['guide_paths'][0]);same=np.array_equal(b,old)
  comparison={'available':True,'equal_float32':same,'changed_components':int(np.count_nonzero(b!=old)),'max_abs':float(abs(b.astype('f8')-old).max()),'previous_report':str(prior/'INTERNAL_PANEL.png.json')}
  (out/'raw_beauty_comparison.json').write_text(json.dumps(comparison,indent=2));assert same,'Raw beauty differs; investigate before filtering'
 else:
  comparison={'available':False,'reason':'First native512 CRT component capture; no independent prior512 beauty exists. Source optical closure and transform equality were checked on fresh reopen; identity recombination uses this capture Combined.'}
  (out/'raw_beauty_comparison.json').write_text(json.dumps(comparison,indent=2))
 albedo=exr(rep['pass_paths']['albedo']);normal=exr(rep['pass_paths']['normal']);features={'beauty':stats(r),'albedo':stats(albedo),'normal':stats(normal)};contract=validate_stats(features,require_nonzero_features=True);assert contract['valid'],contract
 save(b,out/'identity_raw_passthrough.png');check=compare_opaque_pngs(run/'INTERNAL_UNFILTERED_PROBE.png',out/'identity_raw_passthrough.png');assert check['max_error_native']<=1,check
 for name,a in [('remainder',r),('albedo',albedo),('normal',normal)]:write_pfm(out/(name+'.pfm'),a)
 receipt={'mode':'component-recombined','source_sha256':sha(source),'script_sha256':sha(__file__),'component_state_sha256':sha(accounting),'pass_hashes':{k:sha(v) for k,v in rep['pass_paths'].items()},'filter_inputs':{n:sha(out/(n+'.pfm')) for n in ['remainder','albedo','normal']},'physical_remainder_stats':features['beauty'],'feature_contract':contract,'identity_exact_float32':True,'native16_identity_roundtrip':check,'raw_beauty_comparison':comparison,'roundoff_residual_max_abs':float(abs(delta).max()),'transmission_emission_untouched':True,'delta64_untouched':True,'clipping':False,'noise_acceptance':False}
 (out/'export_receipt.json').write_text(json.dumps(receipt,indent=2));print('COMPONENT_EXPORT_READY',json.dumps(receipt),flush=True)
elif mode=='import':
 receipt=json.loads((out/'export_receipt.json').read_text());assert receipt['script_sha256']==sha(__file__);assert receipt['source_sha256']==sha(source);assert receipt['component_state_sha256']==sha(accounting)
 for k,h in receipt['pass_hashes'].items():assert sha(rep['pass_paths'][k])==h
 for n,h in receipt['filter_inputs'].items():assert sha(out/(n+'.pfm'))==h
 d=read_pfm(out/'denoised_remainder.pfm');state=np.load(accounting);t=state['transmission_emission'];delta=state['roundoff_residual'];assert d.shape==t.shape and d.min()>=0
 final=(t.astype('f8')+d.astype('f8')+delta).astype('f4');assert np.isfinite(final).all() and final.min()>=0,'No signed output clipping is allowed'
 image=run/'COMPONENT_RECOMBINED_PANEL.png';assert not image.exists();save(final,image,run/'component_recombined_linear.exr');save(d,out/'INTERNAL_denoised_remainder.png')
 b=exr(rep['pass_paths']['beauty']);patches={}
 if False:  # LCD-only declared patches do not apply to this CRT framing.
  b=np.flipud(b);view=np.flipud(final)
  for name,(x0,y0,x1,y1) in {'orange':(407,323,443,343),'cyan':(537,371,567,391),'green':(587,397,611,415),'red':(666,425,690,443)}.items():
   def energy(a):return float(np.mean(np.diff(a,axis=0)**2)+np.mean(np.diff(a,axis=1)**2))
   control=b[y0:y1,x0:x1];candidate=view[y0:y1,x0:x1];patches[name]={'xyxy':[x0,y0,x1,y1],'gradient_energy_ratio':energy(candidate)/energy(control),'raw_mean_rgb':control.mean(axis=(0,1)).tolist(),'recombined_mean_rgb':candidate.mean(axis=(0,1)).tolist()}
 receipt.update(output=str(image),output_sha256=sha(image),linear_output_sha256=sha(run/'component_recombined_linear.exr'),denoised_remainder_sha256=sha(out/'denoised_remainder.pfm'),script_sha256=sha(__file__),stripe_patches=patches,classification='Component-recombined candidate; OIDN applied only to nonnegative surface remainder',noise_acceptance=False,visual_review_required=['retained transmission/emission','filtered remainder','final frame native pixels and RGB stripe contrast'],selected=False)
 (run/'component_recombined_receipt.json').write_text(json.dumps(receipt,indent=2));print('COMPONENT_RECOMBINED_READY',str(image),flush=True)
else:raise ValueError(mode)
assert sha(source)==rep['source_sha256']
