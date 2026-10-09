"""Exact native4096 replay, bounded channel comparison, no selected-asset writes."""
import ast,gc,hashlib,io,json,math,os,platform,resource,sys,time,zlib
from pathlib import Path
import numpy as np
import PIL,scipy
from PIL import Image
from scipy.ndimage import gaussian_filter,map_coordinates
from .core import ROOT,material,validate_selection,read,resolve,sha,check_hash,png_header,write_new,ContractError
from .generation import imports_from,load,concrete_fields
from .png_fingerprint import PNGSink,fingerprint_file,classify
from PIL import features

def native_function(path,name,namespace):
 tree=ast.parse(path.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
 exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),namespace)
 return namespace[name]

class Comparison:
 def __init__(self,cid,out,compression):
  self.material=material(cid);self.cid=cid;self.out=Path(out);self.n=4096;self.compression=compression
  validate_selection(self.material);self.meta=read(resolve(self.material['map_contract']['source']));self.root=resolve(self.material['map_contract']['map_root']);self.results={};self.start=time.time()
  self.source_hashes={p.name:sha(p) for p in self.root.glob('*.png')};self.preflight={p.name:png_header(p) for p in self.root.glob('*.png')}
  if any(v['resolution']!=[self.n,self.n] for v in self.preflight.values()):raise ContractError('Every canonical map must actually be4096²')
 def compare(self,name,array,bits,encoded=False):
  path=self.root/(name+'.png');header=png_header(path)
  if header['bits']!=bits or header['resolution']!=[self.n,self.n]:raise ContractError('Canonical bit depth/resolution differs: '+name)
  if array.shape[:2]!=(self.n,self.n) or not np.isfinite(array).all():raise ContractError('Generated shape/nonfinite: '+name)
  if not encoded and (array.min()<-1e-5 or array.max()>1.00001):raise ContractError('Generated range: '+name)
  q=array if encoded else np.rint(np.clip(array,0,1)*(65535 if bits==16 else 255)).astype('uint16' if bits==16 else 'uint8')
  unequal=0;max_error=0;absolute_sum=0;pixel_unequal=0;raw_hash=hashlib.sha256();canonical_raw=hashlib.sha256()
  with Image.open(path) as image:
   if image.size!=(self.n,self.n):raise ContractError('Decoder size mismatch')
   for row in range(0,self.n,64):
    a=q[row:row+64];b=np.asarray(image.crop((0,row,self.n,min(row+64,self.n))))
    if a.shape!=b.shape:raise ContractError('Channel component shape mismatch: '+name)
    delta=np.abs(a.astype(np.int32)-b.astype(np.int32));mask=delta!=0
    unequal+=int(np.count_nonzero(mask));pixel_unequal+=int(np.count_nonzero(np.any(mask,axis=-1) if mask.ndim==3 else mask))
    max_error=max(max_error,int(delta.max()));absolute_sum+=int(delta.sum());raw_hash.update(a.tobytes());canonical_raw.update(b.astype(q.dtype).tobytes())
  sink=PNGSink();Image.fromarray(q).save(sink,format='PNG',compress_level=self.compression)
  generated_encoding=sink.fingerprint();canonical_encoding=fingerprint_file(path)
  record={'resolution':[self.n,self.n],'bits_per_channel':bits,'components':q.shape[2] if q.ndim==3 else 1,'source_colorspace':'sRGB' if name=='BaseColor' else 'Non-Color','decoded_pixels_equal':unequal==0,'unequal_components':unequal,'unequal_pixels':pixel_unequal,'total_components':int(q.size),'maximum_integer_error':max_error,'mean_integer_error':absolute_sum/q.size,'generated_pixel_sha256':raw_hash.hexdigest(),'canonical_pixel_sha256':canonical_raw.hexdigest(),'generated_png_sha256':sink.digest.hexdigest(),'canonical_png_sha256':self.source_hashes[path.name],'png_bytes_equal':sink.digest.hexdigest()==self.source_hashes[path.name],'generated_png_bytes':sink.bytes,'canonical_png_bytes':path.stat().st_size,'compression_level':self.compression}
  record.update(encoding_classification=classify(generated_encoding,canonical_encoding,unequal==0),generated_encoding=generated_encoding,canonical_encoding=canonical_encoding)
  self.results[name]=record;self.snapshot(False);print('PARITY_CHANNEL',self.cid,name,'pixels',record['decoded_pixels_equal'],'bytes',record['png_bytes_equal'],'different',unequal,flush=True)
  if not encoded:del q
  gc.collect()
 def snapshot(self,complete,extra=None):
  report={'material_id':self.cid,'resolution':[4096,4096],'complete':complete,'channels':self.results,'all_decoded_pixels_equal':all(x['decoded_pixels_equal'] for x in self.results.values()),'all_png_bytes_equal':all(x['png_bytes_equal'] for x in self.results.values()),'runtime_seconds':time.time()-self.start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'Pillow':PIL.__version__,'zlib':zlib.ZLIB_VERSION,'pillow_zlib':features.version('zlib'),'pillow_zlib_ng':features.version_feature('zlib_ng')},'canonical_files_unchanged':all(sha(self.root/k)==v for k,v in self.source_hashes.items()),'selected_scene_sha256':self.material['selected']['scene_sha256'],'units':'meters','source_resolution_contract':'Exact4096 native replay. Cross-resolution output is not claimed equal: rasterization/anti-alias/packing stopping rules depend on resolution.','visual_acceptance':False,'full_scene_regeneration':False,'generated_images_written_to_disk':False}
  if extra:report.update(extra)
  report['qualification']='pending' if not complete else 'exact_native_pixels_and_png_bytes' if report['all_png_bytes_equal'] else 'exact_native_pixels_only_png_byte_gate_failed' if report['all_decoded_pixels_equal'] else 'native_pixel_gate_failed'
  # Only the new run's progress receipt is mutable; selected inputs never are.
  (self.out/'parity.json').write_text(json.dumps(report,indent=2)+'\n')
  return report

def brass(out):
 c=Comparison('05_champagne_brass',out,6);root=ROOT/'vendor/original';lock=read(root/'SOURCE_LOCK.json')
 for e in lock['files']:check_hash(ROOT/e['path'],e['sha256'])
 with imports_from(root/'source',['surface_math','physical_features','finish_structures','wood_anatomy']):
  m=load(root/'source/generate_materials.py','native_original_materials');s=m.SPECS[4];n=4096
  assert s['tile_m']==c.meta['tile_m'] and s['height_scale_m']==c.meta['height_scale_m']
  x=((np.arange(n,dtype=np.float32)+.5)/n)[None,:];y=((np.arange(n,dtype=np.float32)+.5)/n)[:,None]
  values=m.GENERATORS[4](x,y,n);rgb,h,rough,ao,metal,wear=[m.clamp(v).astype(np.float32) for v in values[:6]];del values
  h=gaussian_filter(h,.35,mode='wrap');macro=gaussian_filter(h,max(.5,n/1024),mode='wrap')
  c.compare('BaseColor',rgb,8);del rgb;gc.collect()
  c.compare('Height',h,16);c.compare('Height_Macro',macro,16)
  micro=m.normal_from_height(h-macro,s,n);c.compare('Normal_Micro_OpenGL',micro,8);del micro,macro;gc.collect()
  orm=np.stack((ao,rough,metal),axis=-1);c.compare('ORM',orm,8);del orm;gc.collect()
  for key,a in [('Roughness',rough),('AO',ao),('Metallic',metal),('WearMask',wear)]:c.compare(key,a,8)
  del rough,ao,metal,wear,a;gc.collect()
  opacity=np.ones((n,n),np.float32);c.compare('Opacity',opacity,8);del opacity
  normal=m.normal_from_height(h,s,n);del h
  encoded=np.rint(m.clamp(normal)*255).astype(np.uint8);del normal;gc.collect()
  c.compare('Normal_OpenGL',encoded,8,True);encoded[...,1]=255-encoded[...,1];c.compare('Normal_DirectX',encoded,8,True);del encoded
 expected=set(p.stem for p in c.root.glob('*.png'))
 if set(c.results)!=expected:raise ContractError('Brass channel coverage incomplete')
 return c.snapshot(True,{'source_lock_sha256':sha(root/'SOURCE_LOCK.json'),'source_commit':lock['commit'],'tile_m':s['tile_m'],'height_scale_m':s['height_scale_m'],'coverage':f'All{len(c.results)} canonical exchange maps; source geometry is retained, not regenerated'})

def concrete(out):
 c=Comparison('32_concrete_polished',out,4);n=4096;root=ROOT/'vendor/aggregate';s=read(root/'recipe.json');fields,source=concrete_fields(n)
 normal=native_function(root/'source/maps.py','normals',{'np':np})
 h=fields.pop('Height');c.compare('Height',h,16)
 for key in ['BaseColor','Roughness','AggregateMask','Metallic','Opacity']:
  a=fields.pop(key);c.compare(key,a,16 if key=='Roughness' else 8);del a;gc.collect()
 fields.clear();gc.collect()
 nm=normal(h,s);c.compare('Normal_OpenGL',nm,8);nm[...,1]=1-nm[...,1];c.compare('Normal_DirectX',nm,8);del nm;gc.collect()
 geometry=resolve('shared/material-aggregate-rebuild/expansion/maps/32_concrete_polished/GeometryHeight.npy');canonical_geometry=c.root/'GeometryHeight.npy'
 check_hash(canonical_geometry,sha(geometry));mesh=np.load(geometry);grid=mesh.shape[0]-1;macro=np.empty((n,n),'f4');px=(np.arange(n,dtype='f4')+.5)*grid/n
 for row in range(0,n,128):
  py=(np.arange(row,min(row+128,n),dtype='f4')+.5)*grid/n;yy,xx=np.meshgrid(py,px,indexing='ij');macro[row:row+len(py)]=map_coordinates(mesh,np.stack((yy,xx)),order=1,mode='nearest')
 c.compare('Height_Macro',macro,16);residual=h-macro;micro=normal(residual,s);c.compare('Normal_Micro_OpenGL',micro,8);del residual,macro,micro,mesh;gc.collect()
 hm=(h-.5)*s['height_scale_m'];ao=np.zeros_like(h);steps=[.00025,.0007,.0017,.0045,.010]
 for angle in np.arange(8)*math.tau/8:
  horizon=np.zeros_like(h)
  for dist in steps:
   sx=int(round(math.cos(angle)*dist/s['tile_m']*n));sy=int(round(math.sin(angle)*dist/s['tile_m']*n))
   if sx==sy==0:continue
   distance=math.hypot(sx,sy)*s['tile_m']/n;shifted=np.roll(np.roll(hm,sy,0),sx,1);slope=(shifted-hm)/distance
   if s.get('tiling')=='finite specimen':
    if sy>0:slope[:sy]=0
    if sy<0:slope[sy:]=0
    if sx>0:slope[:,:sx]=0
    if sx<0:slope[:,sx:]=0
   np.maximum(horizon,slope,out=horizon)
  ao+=1/(1+horizon*horizon);del horizon,shifted,slope;gc.collect()
 ao/=8;c.compare('AO',ao,8);del ao,hm,h;gc.collect()
 if set(c.results)!=set(c.meta['maps']):raise ContractError('Concrete channel coverage incomplete')
 return c.snapshot(True,{'source_lock_sha256':sha(root/'SOURCE_LOCK.json'),'tile_m':s['tile_m'],'height_scale_m':s['height_scale_m'],'source_recipe_report':source['recipe'],'geometry_state':{'path':str(geometry),'sha256':sha(geometry),'contract':'Retained selected513² macro geometry state, explicitly required by native exporter; not regenerated'},'coverage':f'All{len(c.results)} canonical exchange maps plus exact retained geometry-state identity'})

def main():
 cid,out=sys.argv[1:];p=Path(out)
 if not p.is_dir() or (p/'parity.json').exists():raise ContractError('Require a new precreated run directory')
 result=brass(out) if cid=='05_champagne_brass' else concrete(out) if cid=='32_concrete_polished' else None
 if result is None:raise ContractError('Unsupported native parity material')
 print(json.dumps(result,indent=2));return 0 if result['all_decoded_pixels_equal'] and result['all_png_bytes_equal'] and result['canonical_files_unchanged'] else 2
if __name__=='__main__':raise SystemExit(main())
