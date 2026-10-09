"""Recovery replay of pinned native sources into a new output root only."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from PIL import Image
from .core import ROOT,read,write_new,check_hash,sha,ContractError
from .generation import imports_from,load

def recover_brass(output):
 out=Path(output);original=ROOT/'vendor/original';lock=read(original/'SOURCE_LOCK.json')
 for e in lock['files']:check_hash(ROOT/e['path'],e['sha256'])
 reference=read(ROOT/'runs/native_brass_codec_verified_20261006_0230/parity.json')
 with imports_from(original/'source',['surface_math','physical_features','finish_structures','wood_anatomy']):
  m=load(original/'source/generate_materials.py','recovery_original_materials');oldargv=sys.argv[:]
  try:
   m.ROOT=out
   sys.argv=['generate_materials.py','--resolution','4096','--only','05_champagne_brass'];m.main()
  finally:sys.argv[:]=oldargv
 folder=out/'materials/05_champagne_brass';rows=[]
 for name,expected in reference['channels'].items():
  p=folder/(name+'.png');digest=hashlib.sha256()
  with Image.open(p) as im:
   if im.size!=(4096,4096):raise ContractError('Native resolution mismatch')
   dtype='<u2' if expected['bits_per_channel']==16 else 'u1'
   for y in range(0,4096,64):digest.update(np.asarray(im.crop((0,y,4096,min(y+64,4096)))).astype(dtype).tobytes())
  if digest.hexdigest()!=expected['canonical_pixel_sha256']:raise ContractError('Recovered native pixels differ: '+name)
  rows.append({'channel':name,'path':str(p),'sha256':sha(p),'native_integer_pixel_sha256':digest.hexdigest(),'canonical_integer_pixel_sha256':expected['canonical_pixel_sha256'],'pixels_exact':True,'png_bytes_match_historical':sha(p)==expected['canonical_png_sha256'],'resolution':[4096,4096],'bits':expected['bits_per_channel']})
 report={'id':'05_champagne_brass','status':'regenerated_native_pixels_exact','source_commit':lock['commit'],'source_lock_sha256':sha(original/'SOURCE_LOCK.json'),'material_json_sha256':sha(folder/'material.json'),'channels':rows,'selected_scene_recovered':False,'selected_registry_changed':False,'classification':'Exact native integer channel reconstruction from pinned original source; historical compressed PNG identity is not claimed','out':str(out)}
 write_new(out/'native_recovery.json',report);return report

def main():
 out=sys.argv[1];print(json.dumps(recover_brass(out),indent=2))
if __name__=='__main__':main()
