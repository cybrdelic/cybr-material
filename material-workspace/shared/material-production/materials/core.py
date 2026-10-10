"""Authoritative selection and fail-closed asset contracts, independent of Blender."""
from __future__ import annotations
import hashlib,json,os,struct
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
class ContractError(ValueError): pass

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())
def write_new(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')

def workspace():return Path(os.environ.get('MATERIAL_WORKSPACE',str(ROOT.parents[1]))).resolve()
def resolve(reference):
 reference=Path(reference)
 if reference.is_absolute() or '..' in reference.parts:raise ContractError('Registry references must be relative without parent traversal')
 path=(workspace()/reference).resolve()
 if not path.is_relative_to(workspace()):raise ContractError('Reference escapes MATERIAL_WORKSPACE')
 from .path_config import remap_path
 legacy=Path('/workspace')/reference;mapped=remap_path(legacy)
 return mapped if mapped!=legacy else path

def registry():
 r=read(ROOT/'registry.json');ids=[m['id'] for m in r['materials']]
 if len(ids)!=36 or len(set(ids))!=36 or sorted(int(i[:2]) for i in ids)!=list(range(1,37)):
  raise ContractError('Registry must contain each of 36 identities exactly once')
 return r

def material(identifier):
 rows=[m for m in registry()['materials'] if m['id']==identifier or m['id'][:2]==str(identifier).zfill(2)]
 if len(rows)!=1:raise ContractError('Unknown/ambiguous material: '+str(identifier))
 return rows[0]

def check_hash(path,expected):
 if not Path(path).is_file():raise ContractError('Missing dependency: '+str(path))
 actual=sha(path)
 if actual!=expected:raise ContractError('Changed dependency: '+str(path))
 return actual

def png_header(path):
 with Path(path).open('rb') as f:b=f.read(29)
 if len(b)<29 or b[:8]!=b'\x89PNG\r\n\x1a\n':raise ContractError('Not PNG: '+str(path))
 width,height=struct.unpack('>II',b[16:24]);return {'resolution':[width,height],'bits':b[24],'color_type':b[25]}

def validate_selection(m,deep=True):
 s=m.get('selected')
 if not s:raise ContractError(m['id']+': no selected panel; requires qualification')
 scene=resolve(s['scene'])
 if m['id']=='28_asphalt' and s.get('scene_is_generated'):
  if not scene.is_file():raise ContractError('Generate reviewed asphalt scene first: '+str(scene))
  proof=read(scene.parent/'scene_contract_receipt.json')
  if proof.get('render_contract_sha256')!=s.get('scene_contract_sha256'):raise ContractError('Generate reviewed asphalt with the current source contract')
  scene_hash=check_hash(scene,proof['scene_sha256'])
 else:scene_hash=check_hash(scene,s['scene_sha256'])
 result={'id':m['id'],'scene':str(scene),'scene_sha256':scene_hash,'map_dependencies':[],'quality':m['quality']}
 c=m.get('map_contract')
 if c:
  meta_path=resolve(c['source']);check_hash(meta_path,c['source_sha256']);meta=read(meta_path)
  for key,v in meta['maps'].items():
   if isinstance(v,dict):name=v['file'];expected=v['sha256'];resolution=v.get('resolution',[4096,4096])
   else:name=key;expected=c.get('map_hashes',{}).get(name);resolution=[meta['resolution']]*2
   p=resolve(c['map_root'])/name
   if not p.is_file():raise ContractError('Missing channel: '+str(p))
   hdr=png_header(p)
   if hdr['resolution']!=resolution:raise ContractError('Channel resolution mismatch: '+str(p))
   digest=sha(p) if deep else expected
   if expected and digest!=expected:raise ContractError('Channel hash mismatch: '+str(p))
   result['map_dependencies'].append({'path':str(p),'sha256':digest,**hdr})
 return result

def validate_inspection(m, inspection):
 if inspection['unit_scale'] != 1:raise ContractError('Non-SI scene units')
 if m['id']=='28_asphalt':
  if inspection.get('render_contract_sha256')!=m['selected'].get('scene_contract_sha256'):raise ContractError('Generated asphalt differs from reviewed render contract')
  images=inspection['images']
  if len(images)!=4 or any(x['resolution']!=[4096,4096] or not x['packed'] for x in images):raise ContractError('Reviewed asphalt requires four packed native4096 images')
  specimen=[o for o in inspection['objects'] if o['name']=='Asphalt / real displaced face']
  if len(specimen)!=1:raise ContractError('Reviewed asphalt specimen missing')
  ob=specimen[0]
  if any(abs(v-.25)>1e-6 for v in ob['dimensions_m'][:2]) or ob['scale']!=[1.,1.,1.] or 'UVMap' not in ob['uv_layers']:raise ContractError('Reviewed asphalt scale or chart changed')
  mats=[x for x in inspection['materials'] if x['name']=='Asphalt / Rolled mineral aggregate']
  normals=[n for mat in mats for n in mat['normal_nodes']]
  if len(normals)!=1 or normals[0]['space']!='OBJECT' or normals[0]['strength']!=1:raise ContractError('Reviewed asphalt requires full object normal')
  return
 if m['id'] not in ['05_champagne_brass','32_concrete_polished']:return
 images=inspection['images']
 if not images or any(x['resolution']!=[4096,4096] for x in images):raise ContractError('Vertical slice requires native4096 textures')
 objects=[o for o in inspection['objects'] if not o['name'].startswith('Slab /')]
 if len(objects)!=1:raise ContractError('Expected one selected specimen')
 ob=objects[0];extent=.2 if m['id']=='05_champagne_brass' else .32
 if any(abs(v-extent)>1e-6 for v in ob['dimensions_m'][:2]):raise ContractError('Physical specimen extent changed')
 if ob['scale']!=[1.,1.,1.] or 'UVMap' not in ob['uv_layers']:raise ContractError('Scale/UV chart changed')
 normals=[n for mat in inspection['materials'] for n in mat['normal_nodes']]
 spaces=sorted(n['space'] for n in normals)
 expected=['OBJECT','TANGENT'] if m['id']=='05_champagne_brass' else ['TANGENT']
 if spaces!=expected or any(n['strength']!=1 for n in normals):raise ContractError('Selected normal-space contract changed')
 if m['id']=='05_champagne_brass':
  counts=ob.get('material_face_counts',{})
  if counts.get('1')!=65536 or counts.get('0')!=1025:raise ContractError('Audited brass top/side partition changed')
  if not any(v['type']=='BEVEL' and v.get('material')==0 for v in ob.get('modifiers',[])):raise ContractError('Bevel must retain uncorrected control material')

def verify_build(build):
 m=material(build['id']);s=m.get('selected')
 if not s:raise ContractError('Build is no longer selected')
 if m['id']=='28_asphalt' and s.get('scene_is_generated'):
  if build['inspection'].get('render_contract_sha256')!=s.get('scene_contract_sha256'):raise ContractError('Reviewed asphalt contract changed')
  check_hash(resolve(s['scene']),build['source']['scene_sha256'])
 else:
  if build['source']['scene_sha256']!=s['scene_sha256']:raise ContractError('Build is no longer selected')
  check_hash(resolve(s['scene']),s['scene_sha256'])
 for dep in build['inspection']['images']:
  if not dep['packed']:check_hash(dep['path'],dep['sha256'])
 return m
