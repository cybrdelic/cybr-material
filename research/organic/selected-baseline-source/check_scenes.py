import bpy,sys,json,math,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];rows=[]
for path in sorted((ROOT/'scenes').glob('*.blend')):
 bpy.ops.wm.open_mainfile(filepath=str(path));scene=bpy.context.scene
 objects=[o for o in scene.objects if o.type=='MESH' and not o.name.startswith('Studio')]
 row={'scene':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'mesh_count':len(objects),'vertices':0,'faces':0,'finite_vertices':True,'degenerate_faces':0,'external_images':[], 'bounds_m':{},'render_claim':'not a visual acceptance test'}
 for ob in objects:
  a=np.empty(len(ob.data.vertices)*3,dtype='f4');ob.data.vertices.foreach_get('co',a);a=a.reshape(-1,3)
  row['vertices']+=len(a);row['faces']+=len(ob.data.polygons);row['finite_vertices'] &= bool(np.isfinite(a).all())
  row['bounds_m'][ob.name]={'min':a.min(0).tolist(),'max':a.max(0).tolist()}
  ar=np.empty(len(ob.data.polygons),dtype='f4');ob.data.polygons.foreach_get('area',ar);row['degenerate_faces']+=int((ar<1e-16).sum())
 row['external_images']=[im.filepath for im in bpy.data.images if im.source=='FILE' and not im.packed_file]
 row['pass']=row['finite_vertices'] and row['degenerate_faces']==0 and len(row['external_images'])==0
 rows.append(row)
(ROOT/'receipts/scene_geometry_tests.json').write_text(json.dumps(rows,indent=2));print(json.dumps([{k:r[k] for k in ('scene','pass','vertices','faces','degenerate_faces')} for r in rows]),flush=True)
assert all(r['pass'] for r in rows)
