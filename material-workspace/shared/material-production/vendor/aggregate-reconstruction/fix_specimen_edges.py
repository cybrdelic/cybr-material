"""Replaces normal-offset Solidify skirts with flat-backed vertical-sided solids.
The complete accepted top mesh, UV coordinates, material and packed maps remain unchanged.
"""
import bpy, sys, math, hashlib, json
from pathlib import Path
args=sys.argv[sys.argv.index('--')+1:]
for source in args:
    path=Path(source);bpy.ops.wm.open_mainfile(filepath=str(path))
    changed=[]
    for ob in list(bpy.context.scene.objects):
        modifiers=[m for m in ob.modifiers if m.type=='SOLIDIFY' and m.name=='Solid specimen edge']
        if not modifiers:continue
        old=ob.data;n=round(math.sqrt(len(old.vertices)))-1
        assert (n+1)**2==len(old.vertices)
        top_count=len(old.polygons);v=[tuple(p.co) for p in old.vertices];f=[tuple(p.vertices) for p in old.polygons];uv=[(i/n,j/n) for j in range(n+1) for i in range(n+1)]
        boundary=list(range(n+1))+[j*(n+1)+n for j in range(1,n+1)]+[n*(n+1)+i for i in range(n-1,-1,-1)]+[j*(n+1) for j in range(n-1,0,-1)]
        bottom=[]
        for idx in boundary:
            bottom.append(len(v));x,y,z=v[idx];v.append((x,y,-.00095));uv.append(uv[idx])
        for i,a in enumerate(boundary):
            k=(i+1)%len(boundary);f.append((a,bottom[i],bottom[k],boundary[k]))
        f.append(tuple(reversed(bottom)))
        materials=list(old.materials);new=bpy.data.meshes.new(old.name+' / flat backing');new.from_pydata(v,[],f);new.update()
        for material in materials:new.materials.append(material)
        core=bpy.data.materials.new('Specimen / clean cut core');core.use_nodes=True;p=core.node_tree.nodes.get('Principled BSDF');asp='Asphalt' in ob.name
        p.inputs['Base Color'].default_value=(.025,.028,.027,1) if asp else (.19,.185,.165,1);p.inputs['Roughness'].default_value=.91
        new.materials.append(core);layer=new.uv_layers.new(name='UVMap')
        for poly in new.polygons:
            poly.use_smooth=poly.index<top_count
            if poly.index>=top_count:poly.material_index=len(materials)
            for li in poly.loop_indices:layer.data[li].uv=uv[new.loops[li].vertex_index]
        ob.data=new
        for mod in modifiers:ob.modifiers.remove(mod)
        ob['sidewall_construction']='Vertical perimeter walls to planar underside at -0.95 mm, no normal-extrusion skirts'
        ob['surface_top_unchanged']=True;changed.append(ob.name)
    assert changed,path
    target=path.with_name(path.stem+'_edge.blend');bpy.ops.wm.save_as_mainfile(filepath=str(target))
    receipt=dict(source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),result=str(target),sha256=hashlib.sha256(target.read_bytes()).hexdigest(),changed_objects=changed,top_mesh_unchanged=True,maps_unchanged=True)
    target.with_suffix('.edge.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
