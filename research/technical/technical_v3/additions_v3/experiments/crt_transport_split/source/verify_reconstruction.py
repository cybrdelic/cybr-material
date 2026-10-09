"""Read retained EXRs only. No render, denoise, source mutation or clipping."""
import bpy,numpy as np,json,hashlib,sys
from pathlib import Path
root=Path(sys.argv[sys.argv.index('--')+1]);report=json.loads((root/'probe_report.json').read_text())
def read(p):
 im=bpy.data.images.load(str(p),check_existing=False);w,h=im.size;a=np.empty(w*h*4,'f4');im.pixels.foreach_get(a);bpy.data.images.remove(im);return a.reshape(h,w,4)[:,:,:3].astype('f8')
a={k:read(v) for k,v in report['pass_paths'].items()};parts={k:(a[k+'_direct']+a[k+'_indirect'])*a[k+'_color'] for k in ['diffuse','glossy','transmission']}
t=parts['transmission']+a['emission'];r=parts['diffuse']+parts['glossy']+a['environment']+a['volume_direct']+a['volume_indirect'];b=a['beauty']
# Persistent float32 component files require a separate exact accounting residual.
t32=t.astype('f4');r32=r.astype('f4');delta=b-t32.astype('f8')-r32.astype('f8')
rebuilt=(t32.astype('f8')+r32.astype('f8')+delta).astype('f4');identity=np.array_equal(rebuilt,b.astype('f4'));assert identity
assert r32.min()>0 and np.isfinite(delta).all()
np.savez_compressed(root/'linear_component_accounting.npz',transmission_emission=t32,physical_remainder=r32,roundoff_residual=delta)
d={
 'status':'IDENTITY_RECOMBINATION_PASSED; native component noise remains unqualified',
 'source_sha256':report['source_sha256'],'rendered_probe_sha256':report['pass_statistics']['beauty']['sha256'],
 'reconstruction':'float32(T32.astype(float64)+R32.astype(float64)+delta64) equals original beauty float32 exactly; delta64 is preserved, never denoised or clipped',
 'identity_float32_bit_exact':identity,'changed_float32_components':int(np.count_nonzero(rebuilt!=b.astype('f4'))),
 'transmission_emission_min':float(t32.min()),'physical_remainder_min':float(r32.min()),'roundoff_residual_max_abs':float(abs(delta).max()),'roundoff_residual_rms':float(np.sqrt(np.mean(delta**2))),
 'physical_pass_reconstruction':report['reconstruction'],
 'all_features_bounds_verified':bool(a['albedo'].min()>=0 and a['albedo'].max()<=1 and a['normal'].min()>=-1 and a['normal'].max()<=1),
 'feature_scope':'First-hit white glass and opaque surface intrinsic reflectivity; actual shading normals. Used only for remaining surface-lighting component. No claim they encode transmitted stripes.',
 'source_findings':[{'url':'https://raw.githubusercontent.com/blender/blender/v4.3.2/intern/cycles/kernel/film/light_passes.h','lines':[387,391,473,477],'finding':'Transmission weighting subtracts separately computed diffuse and glossy weights from one, permitting small signed cancellation residuals.'},{'url':'https://raw.githubusercontent.com/blender/blender/v4.3.2/intern/cycles/kernel/film/read.h','lines':[155,187],'finding':'Light-path export divides by color and applies exposure in float arithmetic.'}],
 'normalization_or_clipping':False,'material_or_optics_changed':False,'visual_review':{'retained_screen':'visibly noisy at64 samples','remainder':'floor/casing/reflection noise','full_frame_noise_acceptance':False},
 'next_gate':'Native960px512-sample decomposition; inspect both raw components, denoise only R32 with bounded surface features, carry delta64 unchanged; require actual noise and stripe preservation review before delivery'}
(root/'reconstruction_accounting.json').write_text(json.dumps(d,indent=2));print(json.dumps(d))
