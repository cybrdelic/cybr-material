"""Exact selected-aggregate recovery primitives; outputs are always new directories."""
from pathlib import Path
import time,json,resource,numpy as np
from scipy.ndimage import gaussian_filter,map_coordinates
from .core import ROOT,ContractError,sha,write_new
from .generation import concrete_fields
EXPECTED_GEOMETRY='7cbbbea96e7268d94298dae2386738063c8ca8610449a73b5dfba93450781572'

def rebuild_concrete_scene(native,out):
 """Source reconstruction only; the historical binary identity is not asserted."""
 import os
 from .__main__ import destination,run
 out=destination(out);native=Path(native).resolve()
 meta=json.loads((native/'native_recovery.json').read_text())
 if not meta['all_pixels_exact'] or not meta['all_png_bytes_exact']:raise ContractError('Exact native recovery required')
 inputs=[ROOT/'vendor/aggregate-reconstruction'/n for n in ['geometry.py','fix_specimen_edges.py','rebind_slab.py','selected_recipe_1536.json']]
 inputs.extend([native/'native_recovery.json',native/'maps/GeometryHeight.npy',ROOT.parent/'material-slabs/source/studio.py',ROOT/'tests/fixtures/concrete_inspection.json',ROOT/'materials/aggregate_rebuild_stage.py'])
 cfg={'out':str(out),'native':str(native),'input_hashes':{str(p):sha(p) for p in inputs}}
 write_new(out/'rebuild_config.json',cfg)
 run([os.environ.get('BLENDER_BIN','blender'),'-b','--threads','1','--python-exit-code','1','--python',str(ROOT/'materials/aggregate_rebuild_stage.py'),'--',str(out/'rebuild_config.json')],out/'build.log')
 return json.loads((out/'scene_reconstruction.json').read_text())

def recover_concrete_geometry(out):
 out=Path(out)
 if out.exists():raise ContractError('Use a new recovery output directory')
 out.mkdir(parents=True);start=time.time();n=1536;values,meta=concrete_fields(n);h=values['Height'];macro=gaussian_filter(h,max(.5,n/512*1.15),mode='wrap');coord=np.linspace(-.5,n-.5,513,dtype=np.float32);yy,xx=np.meshgrid(coord,coord,indexing='ij');mesh=map_coordinates(macro,np.stack((yy,xx)),order=1,mode='grid-wrap').astype('f4');path=out/'GeometryHeight.npy';np.save(path,mesh);digest=sha(path);r={'id':'32_concrete_polished','role':'Exact selected1536-field to retained513²mesh reconstruction; no material variant','source':meta,'geometry_shape':list(mesh.shape),'normalized_height_range':[float(mesh.min()),float(mesh.max())],'height_scale_m':.002,'physical_width_m':.32,'path':str(path),'sha256':digest,'expected_sha256':EXPECTED_GEOMETRY,'exact_match':digest==EXPECTED_GEOMETRY,'seconds':time.time()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'selected_registry_changed':False};write_new(out/'geometry_recovery.json',r)
 if not r['exact_match']:raise ContractError('Geometry differs; preserved failed recovery, cannot use as selected geometry')
 return r
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();print(json.dumps(recover_concrete_geometry(a.out),indent=2))

class ConcreteRecoveryComparison:
 def __init__(self,cid,out,compression,geometry):
  from .core import read,material
  self.cid=cid;self.out=Path(out);self.root=self.out/'maps';self.root.mkdir();self.n=4096;self.compression=compression;self.start=time.time();self.results={};self.reference=read(ROOT/'runs/native_concrete_final_20261006_0233/parity.json');self.meta={'maps':self.reference['channels']};self.material=material(cid)
  import shutil
  shutil.copy2(geometry,self.root/'GeometryHeight.npy')
 def compare(self,name,array,bits,encoded=False):
  from PIL import Image
  import hashlib,gc
  expected=self.reference['channels'][name]
  if array.shape[:2]!=(4096,4096) or not np.isfinite(array).all():raise ContractError('Invalid native field '+name)
  q=array if encoded else np.rint(np.clip(array,0,1)*(65535 if bits==16 else 255)).astype('uint16' if bits==16 else 'uint8');pixels=hashlib.sha256(q.tobytes()).hexdigest()
  if pixels!=expected['canonical_pixel_sha256']:raise ContractError('Native integer pixels differ from preserved canonical hash: '+name)
  path=self.root/(name+'.png');Image.fromarray(q).save(path,compress_level=self.compression);digest=sha(path);row={'file':str(path),'sha256':digest,'bits_per_channel':bits,'resolution':[4096,4096],'native_integer_pixel_sha256':pixels,'expected_pixel_sha256':expected['canonical_pixel_sha256'],'pixels_exact':True,'png_bytes_exact':digest==expected['canonical_png_sha256'],'expected_png_sha256':expected['canonical_png_sha256']};self.results[name]=row;self.snapshot(False);print('RECOVERED_NATIVE32',name,'pixels_exact',True,'png_exact',row['png_bytes_exact'],flush=True)
  if not encoded:del q
  gc.collect()
 def snapshot(self,complete,extra=None):
  r={'id':self.cid,'complete':complete,'status':'exact_native_pixels_recovered' if complete else 'in_progress','channels':self.results,'all_pixels_exact':all(x['pixels_exact'] for x in self.results.values()),'all_png_bytes_exact':all(x['png_bytes_exact'] for x in self.results.values()),'seconds':time.time()-self.start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'selected_registry_changed':False,'original_binary_scene_recovered':False}
  if extra:r.update(extra)
  (self.out/'native_recovery.json').write_text(json.dumps(r,indent=2));return r

def recover_concrete_native(out,geometry):
 from . import native_parity
 from .core import check_hash
 out=Path(out);geometry=Path(geometry).resolve();check_hash(geometry,EXPECTED_GEOMETRY)
 if out.exists():raise ContractError('Use a new recovery output directory')
 out.mkdir(parents=True);old_comparison=native_parity.Comparison;old_resolve=native_parity.resolve
 try:
  native_parity.Comparison=lambda cid,dest,compression:ConcreteRecoveryComparison(cid,dest,compression,geometry)
  native_parity.resolve=lambda value:geometry if value=='shared/material-aggregate-rebuild/expansion/maps/32_concrete_polished/GeometryHeight.npy' else old_resolve(value)
  return native_parity.concrete(out)
 finally:native_parity.Comparison=old_comparison;native_parity.resolve=old_resolve
