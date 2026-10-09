import bpy,hashlib,json
import numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def signature():
    o=bpy.context.scene.objects['PETG / unified printed wall rim and solid floor'];h=hashlib.sha256()
    a=np.empty(len(o.data.vertices)*3,np.float32);o.data.vertices.foreach_get('co',a);h.update(a.tobytes());a=np.empty(len(o.data.loops),np.int32);o.data.loops.foreach_get('vertex_index',a);h.update(a.tobytes());h.update(np.asarray(o.matrix_world,np.float32).tobytes())
    m=o.data.materials[0];state=[]
    for n in m.node_tree.nodes:
        inp=[]
        for socket in n.inputs:
            if not hasattr(socket,'default_value'):continue
            value=socket.default_value
            if isinstance(value,(float,int,str,bool)) or value is None:pass
            else:
                try:value=list(value)
                except TypeError:value=str(type(value))
            inp.append((socket.name,value))
        state.append({'type':n.bl_idname,'name':n.name,'inputs':inp,'packed_image_sha256':hashlib.sha256(n.image.packed_file.data).hexdigest() if n.type=='TEX_IMAGE' and n.image and n.image.packed_file else None})
    links=sorted((l.from_node.name,l.from_socket.name,l.to_node.name,l.to_socket.name) for l in m.node_tree.links)
    return {'geometry_world_sha256':h.hexdigest(),'material_sha256':hashlib.sha256(json.dumps({'nodes':state,'links':links},sort_keys=True).encode()).hexdigest()}
result={}
for stem in ['23_petg_transparent_r5_hero','23_petg_transparent_r6_sideproof','23_petg_transparent_r7_backlit_control']:
    bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes'/(stem+'.blend')));result[stem]=signature()
assert len({json.dumps(v,sort_keys=True) for v in result.values()})==1,result
(ROOT/'tests/clear_petg_controls_identical_specimen.json').write_text(json.dumps({'unchanged_specimen_verified':True,'scenes':result,'scope':'Mesh topology/coordinates, world transform, shader inputs/links and packed image bytes'},indent=2)+'\n');print('CLEAR_CONTROLS_IDENTICAL_SPECIMEN_PASS',flush=True)
