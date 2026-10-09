"""Explicit relocation map for restored resources. No source-file rewriting."""
import hashlib,json,os
from pathlib import Path

def configuration_state():
 config=os.environ.get('MATERIAL_PATH_MAP_FILE')
 raw=Path(config).read_bytes() if config else None
 data=json.loads(raw) if raw is not None else {'schema_version':1,'prefixes':{}}
 if data.get('schema_version')!=1 or not isinstance(data.get('prefixes'),dict):raise ValueError('Invalid material path-map schema')
 rows=[]
 for old,new in data['prefixes'].items():
  a,b=Path(old),Path(new)
  if not a.is_absolute() or not b.is_absolute() or '..' in a.parts or '..' in b.parts:raise ValueError('Path-map roots must be explicit absolute paths')
  rows.append((a,b))
 rows=sorted(rows,key=lambda x:len(x[0].parts),reverse=True)
 helper=Path(__file__).resolve()
 return {'config_path':str(Path(config).resolve()) if config else None,'config_sha256':hashlib.sha256(raw).hexdigest() if raw is not None else None,'prefixes':{str(a):str(b) for a,b in rows},'helper_path':str(helper),'helper_sha256':hashlib.sha256(helper.read_bytes()).hexdigest()}

def configured_mappings(state=None):
 state=configuration_state() if state is None else state
 return [(Path(a),Path(b)) for a,b in state['prefixes'].items()]

def remap_path(path,state=None):
 path=Path(path)
 if not path.is_absolute():raise ValueError('Only resolved absolute paths may be remapped')
 if '..' in path.parts:raise ValueError('Resolve parent segments before remapping')
 for old,new in configured_mappings(state):
  if path.is_relative_to(old):return new/path.relative_to(old)
 return path

def remap_scene_images(bpy,state=None):
 """Reload relocated external images in memory; packed bytes/source file stay intact."""
 state=configuration_state() if state is None else state
 changes=[]
 for im in bpy.data.images:
  if im.source!='FILE' or not im.filepath or im.packed_file:continue
  old=Path(bpy.path.abspath(im.filepath)).resolve();new=remap_path(old,state)
  if new!=old:
   if not new.is_file():raise FileNotFoundError('Relocated image is missing: '+str(new))
   im.filepath=str(new);im.reload();changes.append({'image':im.name,'from':str(old),'to':str(new)})
 return changes
