from pathlib import Path
import sys,json,hashlib,ast
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parents[1];root=R.parents[1]
(R/'tests').mkdir(exist_ok=True)
sys.path[:0]=[str(R/'source'),str(R/'base_recipe/vendor')]
from graded_fragment_field import graded_fragment_field,contact_state
from fractured_construction import packed_grains as after
from reviewed_construction import packed_grains as before

# Coarse function is identical to the reviewed candidate, not just random calls.
def packed_ast(path):
 t=ast.parse(path.read_text());return ast.dump(next(x for x in t.body if isinstance(x,ast.FunctionDef) and x.name=='packed_grains'),include_attributes=False)
assert packed_ast(R/'source/fractured_construction.py')==packed_ast(R/'source/reviewed_construction.py')
bands=[(.0025,.0075,.40,6500),(.0011,.0035,.16,18000),(.00035,.0011,.13,22000)]
a=before(512,.25,1918,bands,.13);b=after(512,.25,1918,bands,.13)
for key in ['mask','height_m','rgb','roughness','ids']:assert np.array_equal(a[key],b[key]),key
assert a['primitive_sha256']==b['primitive_sha256']
N=384;tile=N*.25/4096;x=(np.arange(N)+.5)*tile/N;xx,yy=np.meshgrid(x,x)
mask=np.zeros((N,N),'f4')
for cx,cy,r in [(tile*.23,tile*.31,.0025),(tile*.55,tile*.30,.0031),(tile*.77,tile*.73,.0023),(tile*.22,tile*.79,.0019)]:
 mask=np.maximum(mask,np.clip((r-np.sqrt((xx-cx)**2+(yy-cy)**2))/(tile/N)+.5,0,1)).astype('f4')
f=graded_fragment_field(N,tile,2599,mask,strip_rows=128)
g=graded_fragment_field(N,tile,2599,mask,strip_rows=77)
for key in ['height_m','coverage','exposure','mineral_rgb','binder_bed_m','band_ids']:
 assert np.isfinite(f[key]).all(),key
 assert np.array_equal(f[key],g[key]),key+' deterministic and strip-invariant'
assert f['report']['coverage_pixels_without_owner']==0
assert f['report']['exposure_pixels_without_relief']==0
assert f['height_m'].min()>=0 and f['height_m'].max()<.00023
assert f['exposure'].min()>=0 and f['exposure'].max()<=1
assert np.all(f['exposure']<=f['coverage'])
assert f['report']['matrix_exposed_fraction_gt_half']<.25
assert all(v>0 for v in f['report']['matrix_owner_fraction_by_band'].values())
# Matrix on opposing borders gets its same contact state through a wrap.
rolled=np.roll(mask,64,axis=0)
for old,new in zip(contact_state(mask,tile),contact_state(rolled,tile)):
 assert np.allclose(np.roll(old,64,axis=0),new,atol=1e-7)
for args in [(2,tile,0,mask), (N,0,0,mask),(N,.001,0,mask)]:
 try:graded_fragment_field(*args)
 except ValueError:pass
 else:raise AssertionError('invalid input accepted')
binder=np.array([.100,.103,.100],'f4');fine=.055+f['mineral_rgb']*.42
rgb=(binder*(1-f['exposure'][...,None])+fine*f['exposure'][...,None])*(1-mask[...,None])+.23*mask[...,None]
Image.fromarray(np.rint(rgb*255).astype('u1')).save(R/'tests/fine_field_albedo_native_pitch.png')
Image.fromarray(np.rint(f['exposure']*255).astype('u1')).save(R/'tests/fine_field_exposure_native_pitch.png')
report={'pass':True,'exact_reviewed_coarse_function_ast':True,'exact_reviewed_coarse_arrays_at_512':True,'exact_reviewed_primitive_sha256':a['primitive_sha256'],'native_pitch_field_resolution':[N,N],'tile_m':tile,'native_pixel_pitch_m':tile/N,'deterministic_strip_invariance':True,'periodic_coarse_contact_state':True,'shared_surface_ownership':True,'invalid_inputs_rejected':True,'fine_field':f['report'],'scope':'Source/function preservation plus same-pitch fine-field tests; does not constitute beauty-render or physical validation'}
(R/'tests/refinement_test_receipt.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
