"""Bounded adapters to pinned selected generators; no historical CLI side effects."""
import importlib.util,sys
from contextlib import contextmanager
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from .core import ROOT,resolve,read,sha,check_hash,ContractError

@contextmanager
def imports_from(path,names):
 """Legacy modules use bare imports; isolate those names instead of sharing state."""
 oldpath=sys.path[:];saved={n:sys.modules.get(n) for n in names}
 for n in names:sys.modules.pop(n,None)
 sys.path.insert(0,str(path))
 try:yield
 finally:
  sys.path[:]=oldpath
  for n,value in saved.items():
   sys.modules.pop(n,None)
   if value is not None:sys.modules[n]=value

def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def concrete_fields(n):
 root=ROOT/'vendor/aggregate'
 source_lock=read(root/'SOURCE_LOCK.json')
 for entry in source_lock['files']:check_hash(ROOT/entry['path'],entry['sha256'])
 lock={Path(e['path']).name:e['sha256'] for e in source_lock['files']}
 with imports_from(root/'source',['surface_math']):
  module=load(root/'source/mineral_construction_bounded.py','selected_concrete_construction')
  meta=read(root/'recipe.json')
  fields,report=module.construction(meta,n)
 return fields,{'mode':'procedural evaluation of selected authored recipe','source_sha256':lock,'tile_m':meta['tile_m'],'height_scale_m':meta['height_scale_m'],'recipe':report}

def brass_fields(n):
 root=ROOT/'vendor/original';lock=read(root/'SOURCE_LOCK.json')
 for entry in lock['files']:check_hash(ROOT/entry['path'],entry['sha256'])
 with imports_from(root/'source',['surface_math','physical_features','finish_structures','wood_anatomy']):
  module=load(root/'source/generate_materials.py','selected_original_materials')
  meta=module.SPECS[4];x=((np.arange(n,dtype=np.float32)+.5)/n)[None,:];y=((np.arange(n,dtype=np.float32)+.5)/n)[:,None]
  values=module.GENERATORS[4](x,y,n)
  fields={k:module.clamp(a).astype(np.float32) for k,a in zip(['BaseColor','Height','Roughness','AO','Metallic','WearMask'],values)}
  fields['Height']=gaussian_filter(fields['Height'],.35,mode='wrap')
  fields['Opacity']=np.ones((n,n),np.float32)
  fields['Normal_OpenGL']=module.normal_from_height(fields['Height'],meta,n)
 return fields,{'mode':'procedural evaluation of recovered exact selected brass source','original_commit':lock['commit'],'source_lock_sha256':sha(root/'SOURCE_LOCK.json'),'tile_m':meta['tile_m'],'height_scale_m':meta['height_scale_m'],'remaining_qualification':'Full4096 regeneration/byte comparison not run in constrained workspace'}

def asphalt_fields(n):
 """The explicitly reviewed fracture model, separate from later experiments."""
 root=ROOT/'vendor/asphalt-reviewed-20261009'
 lock=read(root/'SOURCE_LOCK.json')
 for row in lock['files']:check_hash(root/row['path'],row['sha256'])
 # The pinned sources use one legacy helper namespace; restore callers' modules.
 with imports_from(root/'base_recipe/vendor',['surface_math','fracture_field']):
  oldpath=sys.path[:]
  try:
   sys.path.insert(0,str(root/'source'))
   module=load(root/'source/fractured_construction.py','reviewed_asphalt_construction')
   meta=read(root/'base_recipe/material.json')
   fields,report=module.construction(meta,n)
  finally:sys.path[:]=oldpath
 return fields,{'mode':'procedural evaluation of user-reviewed asphalt fracture source',
  'source_lock_sha256':sha(root/'SOURCE_LOCK.json'),'tile_m':meta['tile_m'],
  'height_scale_m':meta['height_scale_m'],'recipe':report,
  'normal_space':'OBJECT: full metric height, not tangent residual',
  'physical_qualification':False}

def validate_fields(fields,n):
 required={'BaseColor','Height','Roughness','Metallic','Opacity'}
 if not required<=fields.keys():raise ContractError('Missing generated channel')
 for k,a in fields.items():
  if a.shape[:2]!=(n,n) or not np.isfinite(a).all():raise ContractError('Bad shape/nonfinite channel: '+k)
  if a.min()<-1e-6 or a.max()>1+1e-6:raise ContractError('Out-of-range channel: '+k)
 if fields['BaseColor'].shape!=(n,n,3):raise ContractError('BaseColor must be RGB')
 return {k:{'shape':list(a.shape),'min':float(a.min()),'max':float(a.max()),'sha256':__import__('hashlib').sha256(a.tobytes()).hexdigest()} for k,a in fields.items()}

def evaluate(identifier,n):
 if not 16<=n<=256:raise ContractError('Smoke evaluations require16..256; never production resolution')
 if identifier=='05_champagne_brass':fields,meta=brass_fields(n)
 elif identifier=='32_concrete_polished':fields,meta=concrete_fields(n)
 elif identifier=='28_asphalt':fields,meta=asphalt_fields(n)
 else:raise ContractError('No source-generation adapter qualified for '+identifier)
 meta.update(resolution=[n,n],tier='integration smoke only; cannot promote or replace selected4K assets',channels=validate_fields(fields,n))
 return fields,meta
