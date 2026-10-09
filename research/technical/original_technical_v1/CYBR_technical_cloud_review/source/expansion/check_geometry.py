"""Blender CPU topology/determinism checks; no rendering or device acquisition."""
import argparse,collections,hashlib,json,math,sys
from pathlib import Path
import bpy,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from expansion.geometry import build
from expansion.catalog import SPECS
ROOT=Path(__file__).resolve().parents[2]

def stats(objects):
    report=[];digest=hashlib.sha256()
    for o in objects:
        if o.type!='MESH':continue
        verts=np.empty(len(o.data.vertices)*3,np.float32);o.data.vertices.foreach_get('co',verts)
        assert np.isfinite(verts).all(),o.name;digest.update(verts.tobytes())
        edges=collections.Counter()
        for p in o.data.polygons:
            ids=list(p.vertices);digest.update(np.asarray(ids,np.int32).tobytes())
            for a,b in zip(ids,ids[1:]+ids[:1]):edges[tuple(sorted((a,b)))]+=1
        boundary=sum(v==1 for v in edges.values());nonmanifold=sum(v>2 for v in edges.values())
        if 'closed interlocking rings' in o.name or 'fused bead wall' in o.name:
            assert boundary==0 and nonmanifold==0,(o.name,boundary,nonmanifold)
        report.append(dict(name=o.name,vertices=len(o.data.vertices),faces=len(o.data.polygons),
                           boundary_edges=boundary,nonmanifold_edges=nonmanifold,
                           dimensions_m=list(o.dimensions),components=o.get('procedural_components')))
    return digest.hexdigest(),report

def chain_clearance():
    a=np.linspace(0,math.tau,1024,endpoint=False)
    main=np.stack((.004*np.cos(a),.004*np.sin(a),np.zeros_like(a)),1)
    bridge=np.stack((.005+.0035*np.cos(a),np.zeros_like(a),.0035*np.sin(a)),1)
    minimum=float(np.sqrt(((main[:,None]-bridge[None,:])**2).sum(2)).min())-.0009
    assert minimum>0,minimum
    return dict(minimum_centerline_minus_wire_diameters_m=minimum,
                approximation='1024-sample centerline separation check; ring linkage defined by perpendicular bridge through main holes')

def main():
    p=argparse.ArgumentParser();p.add_argument('--only',nargs='*');p.add_argument('--repeat',action='store_true');p.add_argument('--extremes',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    path=ROOT/'expansion/geometry_checks.json'
    output=json.loads(path.read_text()) if path.exists() else {}
    for s in SPECS:
        if a.only and s['id'] not in a.only:continue
        bpy.ops.wm.read_factory_settings(use_empty=True);_,objects,_=build(s['id']);h,topology=stats(objects)
        if a.repeat:
            bpy.ops.wm.read_factory_settings(use_empty=True);_,objects,_=build(s['id']);other,_=stats(objects);assert h==other,s['id']
        output[s['id']]=dict(geometry_hash=h,repeatable=a.repeat,objects=topology)
        if s['family']=='Chainmail':output[s['id']]['clearance']=chain_clearance()
        if a.extremes and s['family'] in ('PETG','Chainmail'):
            sets=([{'layer_height_m':.0001,'wall_thickness_m':.0006,'bead_width_m':.00025},
                   {'layer_height_m':.0004,'wall_thickness_m':.0024,'bead_width_m':.0007}]
                  if s['family']=='PETG' else [{'pitch_m':.0095,'wire_radius_m':.0003},{'pitch_m':.012,'wire_radius_m':.0006}])
            variants=[]
            for params in sets:
                bpy.ops.wm.read_factory_settings(use_empty=True);_,objects,_=build(s['id'],{s['id']:params})
                vhash,vstats=stats(objects);variants.append(dict(parameters=params,geometry_hash=vhash,objects=vstats))
            output[s['id']]['parameter_extremes']=variants
        print('GEOMETRY_PASS',s['id'],sum(v['vertices'] for v in topology),flush=True)
        path=ROOT/'expansion/geometry_checks.json';path.write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':main()
