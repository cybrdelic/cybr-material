# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
import bpy, json, hashlib, resource, time
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(str(_required_research_input("baseline-scene")))
NAMES = [f'Carpet / spun wool loop cores and fine filaments {i}' for i in range(4)]

def read(col, prop, width=1, dtype='f4'):
    a=np.empty(len(col)*width,dtype=dtype);col.foreach_get(prop,a)
    return a.reshape(-1,width) if width>1 else a

def mesh_sha(ob):
    h=hashlib.sha256();me=ob.data
    for col,prop,width,dtype in ((me.vertices,'co',3,'f4'),(me.loops,'vertex_index',1,'i4'),(me.polygons,'loop_start',1,'i4'),(me.polygons,'loop_total',1,'i4'),(me.polygons,'material_index',1,'i4'),(me.polygons,'use_smooth',1,'i4')):
        h.update(read(col,prop,width,dtype).tobytes())
    for uv in me.uv_layers:h.update(read(uv.data,'uv',2).tobytes())
    return h.hexdigest()

def material_sha(m):
    data=[]
    for n in m.node_tree.nodes:
        inputs=[]
        for s in n.inputs:
            if not hasattr(s,'default_value'):continue
            v=s.default_value
            if not isinstance(v,(int,float,str,bool)):
                try:v=list(v)
                except TypeError:v=str(v)
            inputs.append([s.name,v])
        data.append([n.name,n.bl_idname,getattr(n,'operation',None),inputs])
    links=[[l.from_node.name,l.from_socket.name,l.to_node.name,l.to_socket.name] for l in m.node_tree.links]
    return hashlib.sha256(json.dumps([data,links],sort_keys=True).encode()).hexdigest()

def run():
    start=time.time();bpy.ops.wm.open_mainfile(filepath=str(SOURCE))
    report={'objects':{},'materials':{},'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
    centers=[];radii=[];frames=[];colors=[];phases=[]
    for ob in bpy.context.scene.objects:
        if ob.type!='MESH' or not ob.name.startswith('Carpet /'):continue
        report['objects'][ob.name]={'vertices':len(ob.data.vertices),'polygons':len(ob.data.polygons),'geometry_sha256':mesh_sha(ob),'matrix_world':list(map(list,ob.matrix_world))}
        for m in ob.data.materials:report['materials'][m.name]=material_sha(m)
    for color,name in enumerate(NAMES):
        ob=bpy.data.objects[name];co=read(ob.data.vertices,'co',3)
        assert len(co)%1470==0 and len(ob.data.polygons)==len(co)//1470*1454
        lv=read(ob.data.loops,'vertex_index',dtype='i4')
        for i,block in enumerate(co.reshape(-1,1470,3)):
            rings=block[:294].reshape(49,6,3).astype('f8');center=rings.mean(1)
            radius=np.linalg.norm(rings-center[:,None],axis=-1).mean(1)
            a=(rings[:,0]-center)/radius[:,None]
            b=(rings[:,1]-center-.5*radius[:,None]*a)/(np.sqrt(3)*.5*radius[:,None])
            poly=ob.data.polygons[i*1454]
            assert tuple(lv[poly.loop_start:poly.loop_start+poly.loop_total])==(i*1470,i*1470+1,i*1470+7,i*1470+6)
            poly=ob.data.polygons[i*1454+290]
            assert tuple(lv[poly.loop_start:poly.loop_start+poly.loop_total])==(i*1470+294,i*1470+295,i*1470+299,i*1470+298)
            f=block[294:490].reshape(49,4,3).mean(1)-center
            phase=np.unwrap(np.arctan2(np.sum(f*b,axis=1),np.sum(f*a,axis=1)))
            centers.append(center);radii.append(radius);frames.append(np.stack((a,b),axis=1));colors.append(color);phases.append(phase)
    assert len(centers)==1764
    np.savez_compressed(ROOT/'receipts/baseline_loops.npz',center=centers,radius=radii,frame=frames,color=colors,phase=phases)
    report['verified_blocks']={'loops':1764,'vertices_per_loop':1470,'polygons_per_loop':1454,'core_vertices':294,'core_polygons':290,'retained_surface_vertices':1176,'retained_surface_polygons':1164,'core_rings':49,'core_sides':6,'surface_strands':6,'surface_sides':4}
    report['rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024;report['seconds']=time.time()-start
    (ROOT/'receipts/baseline_inspection.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'loops':len(centers),'rss_mib':report['rss_mib'],'seconds':report['seconds'],'hair_curves_api':hasattr(bpy.data,'hair_curves')}),flush=True)

if __name__=='__main__':run()
