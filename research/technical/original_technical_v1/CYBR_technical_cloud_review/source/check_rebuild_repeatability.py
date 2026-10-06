import bpy,sys,hashlib,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import technical_studies as t
root=Path(__file__).resolve().parents[1]
recipes=[('18_lcd',t.display),('19_crt',t.display),('20_plastic_smooth',t.plastic),('21_plastic_texture',t.plastic),('22_petg',t.petg),('23_petg_transparent',t.petg),('25_tile_ceramic',t.ceramic),('15_future_ceramic',t.engineered),('16_future_ribbed',t.engineered),('14_chainmail',t.chain)]
def digest():
    h=hashlib.sha256();nv=0;nf=0
    for o in bpy.context.scene.objects:
        if o.type!='MESH':continue
        v=np.empty(len(o.data.vertices)*3,np.float32);o.data.vertices.foreach_get('co',v);assert np.isfinite(v).all();h.update(v.tobytes())
        v=np.empty(len(o.data.loops),np.int32);o.data.loops.foreach_get('vertex_index',v);h.update(v.tobytes())
        h.update(np.asarray(o.matrix_world,np.float32).tobytes());nv+=len(o.data.vertices);nf+=len(o.data.polygons)
    return h.hexdigest(),nv,nf
out={}
for id,fn in recipes:
    bpy.ops.wm.read_factory_settings(use_empty=True);fn(id);a,nv,nf=digest()
    bpy.ops.wm.read_factory_settings(use_empty=True);fn(id);b,_,_=digest();assert a==b,id
    out[id]={'geometry_sha256':a,'repeat_equal':True,'vertices':nv,'faces':nf,'scope':'Unmodified procedural specimens, no lighting or render'}
    print('REPEAT_PASS',id,nv,nf,flush=True)
(root/'tests'/'rebuild_repeatability.json').write_text(json.dumps(out,indent=2)+'\n')
