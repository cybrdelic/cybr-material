"""Assemble a gated three-unit material comparison inside the exact selected r5.

This file does not construct fibres or render. Missing/unqualified input fails
closed. The control's frozen absolute rig is reused; no candidate-bound refit.
"""
# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
from pathlib import Path
import gc, hashlib, json, resource, sys, time
import bpy
import numpy as np
from mathutils import Matrix

ROOT=Path(__file__).resolve().parents[1]
R1=ROOT.parent/'carpet-r5-packed-core'
sys.path.insert(0,str(R1/'revisions/r1_frozen/source'))
import inspect_baseline as base
import build_packed_panel as rig
rig.ROOT=R1
BASE=Path(str(_required_research_input("baseline-scene")))
BASE_SHA='55df658fd44a0fb17bab14b46692c1c252b2d2658a967910e8c281f82e743879'
PROFILE=Path(str(_required_research_input("capture-profile")))
PROFILE_SHA='6339495108bf74c7a6a1e4c854cf61e656d9d188a7995362a0d54014a3c06a9a'
SELECTED={968,955,356}
START=None

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def guard(label):
    rss=int(next(line.split()[1] for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmRSS:')))/1024
    if rss>1536 or time.monotonic()-START>180:raise RuntimeError('Build resource limit exceeded at '+label)
    print(json.dumps({'stage':label,'rss_mib':rss}),flush=True)

def checked_inputs():
    manifest_path=ROOT/'receipts/qualified_inputs.json'
    manifest=json.loads(manifest_path.read_text())
    assert manifest['status']=='qualified_native_three_unit_input'
    assert {int(row['loop']) for row in manifest['loops']}==SELECTED
    assert len(manifest['loops'])==3
    assert manifest['scope']=='three central units only; no whole-panel or full manufacturing claim'
    loaded=[]
    for row in manifest['loops']:
        assert row['native_geometry_qualified'] is True
        path=Path(row['arrays']);receipt=Path(row['qualification_receipt'])
        assert sha(path)==row['arrays_sha256'] and sha(receipt)==row['qualification_receipt_sha256']
        state=np.load(path,allow_pickle=False)
        for key,count in [('body',187),('wraps',6)]:
            points=np.asarray(state[key]);assert points.ndim==3 and points.shape[0]==count and points.shape[2]==3
            assert points.shape[1]>=4 and np.isfinite(points).all()
        br=float(state['body_radius']);wr=float(state['wrap_radius'])
        assert 12e-6<br<16e-6 and abs(wr-16e-6)<1e-12
        loaded.append((row,state))
    return manifest,sha(manifest_path),loaded

def replace_old_units(loop_data):
    reports=[]
    for colour,name in enumerate(base.NAMES):
        indices=np.flatnonzero(loop_data['color']==colour)
        chosen=np.isin(indices,list(SELECTED))
        if not chosen.any():continue
        ob=bpy.data.objects[name];old=ob.data;n=len(indices)
        assert len(old.vertices)==n*1470 and len(old.polygons)==n*1454 and len(old.loops)==n*5820
        assert len(old.uv_layers)==1
        co=base.read(old.vertices,'co',3);old_vi=base.read(old.loops,'vertex_index',dtype='i4')
        old_uv=base.read(old.uv_layers.active.data,'uv',2)
        totals=base.read(old.polygons,'loop_total',dtype='i4')
        smooth=base.read(old.polygons,'use_smooth',dtype='i4');mi=base.read(old.polygons,'material_index',dtype='i4')
        kp=np.repeat(~chosen,1454);kl=np.repeat(~chosen,5820)
        vi=old_vi[kl];uv=old_uv[kl];count=totals[kp]
        assert len(vi)==int(count.sum())
        mesh=bpy.data.meshes.new(name+' / retained r5 faces')
        mesh.vertices.add(len(co));mesh.vertices.foreach_set('co',co.ravel())
        mesh.loops.add(len(vi));mesh.loops.foreach_set('vertex_index',vi)
        mesh.polygons.add(len(count));mesh.polygons.foreach_set('loop_start',np.r_[0,np.cumsum(count)[:-1]].astype('i4'))
        mesh.polygons.foreach_set('loop_total',count);mesh.polygons.foreach_set('use_smooth',smooth[kp]);mesh.polygons.foreach_set('material_index',mi[kp])
        layer=mesh.uv_layers.new(name=old.uv_layers.active.name);layer.data.foreach_set('uv',uv.ravel())
        for material in old.materials:mesh.materials.append(material)
        mesh.update(calc_edges=True)
        assert np.array_equal(co,base.read(mesh.vertices,'co',3))
        assert np.array_equal(vi,base.read(mesh.loops,'vertex_index',dtype='i4'))
        assert np.array_equal(uv,base.read(mesh.uv_layers.active.data,'uv',2))
        assert np.array_equal(count,base.read(mesh.polygons,'loop_total',dtype='i4'))
        assert np.array_equal(smooth[kp],base.read(mesh.polygons,'use_smooth',dtype='i4'))
        assert np.array_equal(mi[kp],base.read(mesh.polygons,'material_index',dtype='i4'))
        ob.data=mesh;bpy.data.meshes.remove(old)
        reports.append({'object':name,'removed_loop_indices':indices[chosen].tolist(),'removed_polygons':int(chosen.sum())*1454,'retained_vertices_faces_uv_material_indices_verified':True,'matrix_world':list(map(list,ob.matrix_world))})
        gc.collect();guard('retained colour '+str(colour))
    return reports

def append_curves(name,points,radius,material,transform):
    points=np.asarray(points,dtype='f4');count,keys,_=points.shape
    data=bpy.data.hair_curves.new(name);data.add_curves([keys]*count)
    data.position_data.foreach_set('vector',points.ravel())
    data.attributes.new('radius','FLOAT','POINT').data.foreach_set('value',np.full(count*keys,radius,dtype='f4'))
    data.attributes.new('curve_type','INT8','CURVE').data.foreach_set('value',np.zeros(count,dtype='i4'))
    data.attributes.new('resolution','INT','CURVE').data.foreach_set('value',np.ones(count,dtype='i4'))
    data.materials.append(material);ob=bpy.data.objects.new(name,data);bpy.context.scene.collection.objects.link(ob);ob.matrix_world=transform.copy()
    result=np.empty(points.size,dtype='f4');data.position_data.foreach_get('vector',result)
    assert np.array_equal(points.ravel(),result)
    read_radius=np.empty(count*keys,dtype='f4');data.attributes['radius'].data.foreach_get('value',read_radius)
    assert np.array_equal(read_radius,np.full(count*keys,radius,dtype='f4'))
    return {'object':name,'curves':count,'keys_per_curve':keys,'raw_positions_sha256':hashlib.sha256(result.tobytes()).hexdigest(),'radius_m':float(np.float32(radius)),'raw_radius_sha256':hashlib.sha256(read_radius.tobytes()).hexdigest(),'matrix_world':list(map(list,ob.matrix_world))}

def run():
    global START
    START=time.monotonic();manifest,input_hash,loaded=checked_inputs()
    assert sha(BASE)==BASE_SHA and sha(PROFILE)==PROFILE_SHA
    bpy.ops.wm.open_mainfile(filepath=str(BASE));guard('opened exact r5')
    before={o.name:base.mesh_sha(o) for o in bpy.context.scene.objects if o.type=='MESH' and o.name.startswith('Carpet /')}
    matrices={o.name:list(map(list,o.matrix_world)) for o in bpy.context.scene.objects if o.name.startswith('Carpet /')}
    materials={m.name:base.material_sha(m) for o in bpy.context.scene.objects if o.name.startswith('Carpet /') for m in o.data.materials}
    loop=np.load(R1/'receipts/baseline_loops.npz');changed=replace_old_units(loop);added=[]
    for row,state in loaded:
        index=int(row['loop']);colour=int(loop['color'][index]);old_ob=bpy.data.objects[base.NAMES[colour]]
        # Exact identity protects the already checked physical radii and coordinates.
        assert np.array_equal(np.asarray(old_ob.matrix_world),np.eye(4))
        material=bpy.data.materials[f'Carpet / fine wool yarn {colour}']
        for key,rkey in [('body','body_radius'),('wraps','wrap_radius')]:
            added.append(append_curves(f'Carpet / corrected loop {index} / {key}',state[key],float(state[rkey]),material,old_ob.matrix_world))
    changed_names={row['object'] for row in changed}
    retained={o.name:base.mesh_sha(o) for o in bpy.context.scene.objects if o.type=='MESH' and o.name.startswith('Carpet /') and o.name not in changed_names}
    assert all(before[name]==value for name,value in retained.items())
    assert all(list(map(list,bpy.data.objects[name].matrix_world))==matrix for name,matrix in matrices.items())
    assert all(base.material_sha(bpy.data.materials[name])==digest for name,digest in materials.items())
    rig_sha=rig.apply_rig();profile=json.loads(PROFILE.read_text());capture=profile['capture']
    bpy.context.scene.camera.matrix_world=Matrix(capture['camera_matrix_world'])
    bpy.context.scene.camera.data.type='ORTHO';bpy.context.scene.camera.data.ortho_scale=capture['camera_ortho_scale']
    bpy.context.scene.cycles_curves.shape='THICK';bpy.context.view_layer.update()
    bpy.context.scene['comparison_scope']='Only loops968,955,356 replaced. Original r5 surround/halo/backing retained. Construction comparison, not full-panel qualification.'
    guard('before save');bpy.context.preferences.filepaths.save_version=0
    target=ROOT/'scenes/24_carpet_three_corrected_units.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(target),compress=False);guard('saved')
    assert target.stat().st_size<=512*1024**2
    report={'source':str(target),'source_sha256':sha(target),'source_bytes':target.stat().st_size,'baseline_sha256':BASE_SHA,'qualified_input_manifest_sha256':input_hash,'inputs':manifest,'removed_blocks':changed,'added_curves':added,'retained_object_geometry_sha256':retained,'original_material_sha256':materials,'absolute_rig_sha256':rig_sha,'profile':str(PROFILE),'profile_sha256':PROFILE_SHA,'capture':capture,'scope':'Three central complete yarn units only, within the original r5 full specimen. Internal 8 mm material comparison.','limits':['No full-panel replacement or calibrated manufacturing simulation.','Native body and replacement wrap curves retain material graphs but do not reproduce the old mesh UV-driven bump field.','The two-point/control appearance comparison and any neighbour-contact limitations remain explicit in the input manifest.'],'wall_seconds':time.monotonic()-START,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'assembler_sha256':sha(Path(__file__))}
    (ROOT/'receipts/handoff.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:report[k] for k in ('source','source_sha256','source_bytes','wall_seconds','peak_rss_mib')}),flush=True)

if __name__=='__main__':run()
