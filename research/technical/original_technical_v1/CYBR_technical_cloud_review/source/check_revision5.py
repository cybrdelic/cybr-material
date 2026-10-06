import bpy,collections,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def closed(ob):
    c=collections.Counter()
    for f in ob.data.polygons:
        v=list(f.vertices)
        for a,b in zip(v,v[1:]+v[:1]):c[tuple(sorted((a,b)))]+=1
    assert all(n==2 for n in c.values()),ob.name
def geom():
    h=hashlib.sha256()
    for o in bpy.context.scene.objects:
        if o.type!='MESH':continue
        a=np.empty(len(o.data.vertices)*3,np.float32);o.data.vertices.foreach_get('co',a);h.update(a.tobytes())
        a=np.empty(len(o.data.loops),np.int32);o.data.loops.foreach_get('vertex_index',a);h.update(a.tobytes())
        h.update(np.asarray(o.matrix_world,np.float32).tobytes())
    return h.hexdigest()
results={}
for id in ('18_lcd','19_crt'):
    a=ROOT/'scenes'/f'{id}_r4_hero.blend';b=ROOT/'scenes'/f'{id}_r5_hero.blend'
    bpy.ops.wm.open_mainfile(filepath=str(a));before=geom();bpy.ops.wm.open_mainfile(filepath=str(b));after=geom();assert before==after
    cover=next(o for o in bpy.context.scene.objects if 'coverglass' in o.name or 'faceplate' in o.name);closed(cover)
    p=cover.data.materials[0].node_tree.nodes['Principled BSDF'];assert p.inputs['Transmission Weight'].default_value==1 and abs(p.inputs['IOR'].default_value-1.52)<1e-5
    results[id]={'r4_r5_geometry_identical':True,'geometry_hash':after,'ior':p.inputs['IOR'].default_value,'glass_thickness_m':cover['glass_thickness_m']}
for id in ('20_plastic_smooth','21_plastic_texture'):
    bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes'/f'{id}_r5_hero.blend'))
    objs=[o for o in bpy.context.scene.objects if o.name.startswith('Plastic /')]
    for ob in objs:closed(ob)
    upper=bpy.data.objects['Plastic / upper half rounded profile'];lower=bpy.data.objects['Plastic / lower half rounded profile'];seam=bpy.data.objects['Plastic / recessed seam connector']
    umax=max(v.co.x for v in upper.data.vertices);lmax=max(v.co.x for v in lower.data.vertices);smax=max(v.co.x for v in seam.data.vertices);assert abs(umax-lmax)<1e-7 and umax-smax>.00029
    zupper=min(v.co.z for v in upper.data.vertices);zlower=max(v.co.z for v in lower.data.vertices);assert abs(zupper-zlower-.00035)<1e-7
    nodes=upper.data.materials[0].node_tree.nodes;noise=next(n for n in nodes if n.type=='TEX_NOISE');bump=next(n for n in nodes if n.type=='BUMP')
    assert noise.inputs['Scale'].default_value==(8500 if id=='20_plastic_smooth' else 1850)
    assert abs(bump.inputs['Distance'].default_value-(.000008 if id=='20_plastic_smooth' else .000065))<1e-10
    results[id]={'all_three_parts_closed':True,'shared_half_width_m':umax,'seam_inset_m':umax-smax,'physical_gap_m':zupper-zlower,'shader_parameters_unchanged':True}
for id,rec in results.items():print('REVISION5_PASS',id,flush=True)
(ROOT/'tests/revision5_checks.json').write_text(json.dumps(results,indent=2)+'\n')
