"""Bounded carpet material comparison; exact r1 fibres, original r5 surround.

This is not a full-panel replacement or a manufacturing simulation.
The full original backing, flyaway halo, six surface strands and outer cores stay.
Only smooth cores whose loop centers fall inside the declared local region lose
their faces. Their vertices are retained but are unreferenced by render faces.
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
import sys, json, hashlib, gc, time, resource
import bpy
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R1=ROOT.parent/'carpet-r5-packed-core'
sys.path.insert(0,str(R1/'revisions/r1_frozen/source'))
import inspect_baseline as base
import build_packed_panel as original
original.ROOT=R1
BASE=Path(str(_required_research_input("baseline-scene")))
BASE_SHA='55df658fd44a0fb17bab14b46692c1c252b2d2658a967910e8c281f82e743879'
PROFILE=Path(str(_required_research_input("capture-profile")))
HALF_REGION=.012

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def memory(label):
    rss=int(next(s.split()[1] for s in Path('/proc/self/status').read_text().splitlines() if s.startswith('VmRSS:')))/1024
    print(label,round(rss,2),'MiB',flush=True)
    if rss>1536: raise RuntimeError('Build RSS admission exceeded')
    return rss

def run():
    started=time.time()
    assert sha(BASE)==BASE_SHA
    bpy.ops.wm.open_mainfile(filepath=str(BASE));memory('baseline')
    before={o.name:base.mesh_sha(o) for o in bpy.context.scene.objects if o.type=='MESH' and o.name.startswith('Carpet /')}
    materials={m.name:base.material_sha(m) for o in bpy.context.scene.objects if o.name.startswith('Carpet /') for m in o.data.materials}
    loop=np.load(R1/'receipts/baseline_loops.npz')
    c=loop['center'][:,24,:]
    selected=(np.abs(c[:,0])<=HALF_REGION)&(np.abs(c[:,1])<=HALF_REGION)
    records=[]
    profile=json.loads(PROFILE.read_text());camera=np.array(profile['capture']['camera_matrix_world'])
    camera_origin=camera[:3,3];camera_axes=camera[:3,:2];half_frame=profile['capture']['camera_ortho_scale']/2
    projected_guards=[]
    for color,name in enumerate(base.NAMES):
        ids=np.flatnonzero(loop['color']==color);local=selected[ids]
        ob=bpy.data.objects[name];old=ob.data;n=len(ids)
        assert len(old.vertices)==n*1470 and len(old.polygons)==n*1454
        co=base.read(old.vertices,'co',3)
        transform=np.array(ob.matrix_world)
        outer_core=co.reshape(n,1470,3)[~local,:294,:].astype('f8')
        outer_core=outer_core@transform[:3,:3].T+transform[:3,3]
        projected=(outer_core-camera_origin)@camera_axes
        box_min=projected.min(axis=1);box_max=projected.max(axis=1)
        axis_gap=np.maximum(np.maximum(box_min-half_frame,-half_frame-box_max),0)
        guard=np.linalg.norm(axis_gap,axis=1)
        assert np.all(guard>0), 'Unreplaced smooth core intersects the macro projection'
        projected_guards.append(float(guard.min()))
        del outer_core,projected,box_min,box_max,axis_gap,guard
        oldvi=base.read(old.loops,'vertex_index',dtype='i4')
        olduv=base.read(old.uv_layers.active.data,'uv',2)
        lt=base.read(old.polygons,'loop_total',dtype='i4')
        smooth=base.read(old.polygons,'use_smooth',dtype='i4')
        material=base.read(old.polygons,'material_index',dtype='i4')
        keep_poly=np.ones((n,1454),dtype=bool);keep_poly[local,:290]=False;keep_poly=keep_poly.ravel()
        keep_loop=np.ones((n,5820),dtype=bool);keep_loop[local,:1164]=False;keep_loop=keep_loop.ravel()
        vi=oldvi[keep_loop];uv=olduv[keep_loop];totals=lt[keep_poly]
        assert len(vi)==int(totals.sum())
        me=bpy.data.meshes.new('Carpet / exact r5 surround and retained strands '+str(color))
        me.vertices.add(len(co));me.vertices.foreach_set('co',co.ravel())
        me.loops.add(len(vi));me.loops.foreach_set('vertex_index',vi)
        me.polygons.add(len(totals));me.polygons.foreach_set('loop_start',np.r_[0,np.cumsum(totals)[:-1]].astype('i4'))
        me.polygons.foreach_set('loop_total',totals);me.polygons.foreach_set('use_smooth',smooth[keep_poly]);me.polygons.foreach_set('material_index',material[keep_poly])
        layer=me.uv_layers.new(name=old.uv_layers.active.name);layer.data.foreach_set('uv',uv.ravel())
        for m in old.materials:me.materials.append(m)
        me.update(calc_edges=True)
        assert np.array_equal(co,base.read(me.vertices,'co',3))
        assert np.array_equal(vi,base.read(me.loops,'vertex_index',dtype='i4'))
        assert np.array_equal(uv,base.read(me.uv_layers.active.data,'uv',2))
        ob.data=me;bpy.data.meshes.remove(old)
        ob['scope']='Exact r5 surround. Local smooth core faces only removed; retained strands and every kept polygon/UV are exact.'
        del co,oldvi,olduv,lt,smooth,material,keep_poly,keep_loop,vi,uv,totals;gc.collect()
        positions=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[local].copy()
        rad=np.load(R1/f'arrays/colour_{color}_radius_per_loop.npy')[local]
        count=len(rad)*187
        data=bpy.data.hair_curves.new('Carpet / exact local r1 body '+str(color));data.add_curves([144]*count)
        data.position_data.foreach_set('vector',positions.ravel())
        radius=data.attributes.new('radius','FLOAT','POINT');radius.data.foreach_set('value',np.repeat(rad,187*144))
        kind=data.attributes.new('curve_type','INT8','CURVE');kind.data.foreach_set('value',np.zeros(count,dtype='i4'))
        resolution=data.attributes.new('resolution','INT','CURVE');resolution.data.foreach_set('value',np.ones(count,dtype='i4'))
        # Compare raw authored keys/radii, not an interpolated viewport cache.
        check=np.empty(positions.size,dtype='f4');data.position_data.foreach_get('vector',check)
        assert np.array_equal(positions.ravel(),check)
        raw_sha=hashlib.sha256(check.tobytes()).hexdigest()
        data.materials.append(bpy.data.materials[f'Carpet / fine wool yarn {color}'])
        obj=bpy.data.objects.new('Carpet / exact local packed body '+str(color),data);bpy.context.scene.collection.objects.link(obj)
        obj['scope']='Exact frozen r1 constructed fibres. Local macro comparison only. Residual fold/overlap and native UV-bump limitations remain.'
        records.append({'color':color,'original_loops':n,'replaced_loop_indices':ids[local].tolist(),'replaced_loops':int(local.sum()),'fibres':count,'keys':count*144,'position_sha256':raw_sha,'source_position_file':str(R1/f'arrays/colour_{color}_positions.npy'),'raw_positions_exact':True,'radii_min_m':float(rad.min()),'radii_max_m':float(rad.max()),'kept_mesh_vertices_topology_uv_exact':True})
        del positions,rad,check;gc.collect();memory('group '+str(color))
    retained={o.name:base.mesh_sha(o) for o in bpy.context.scene.objects if o.type=='MESH' and o.name.startswith('Carpet /') and o.name not in base.NAMES}
    assert all(before[k]==v for k,v in retained.items())
    assert materials=={k:base.material_sha(bpy.data.materials[k]) for k in materials}
    rig_sha=original.apply_rig();bpy.context.scene.cycles_curves.shape='THICK'
    bpy.context.scene['local_material_scope']='24 mm exact constructed-core region in the full retained r5 panel, for an 8 mm macro only. Full-panel appearance unqualified.'
    bpy.context.view_layer.update();memory('pre-save')
    out=ROOT/'scenes/24_carpet_local_exact_r1_macro.blend'
    bpy.context.preferences.filepaths.save_version=0
    bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=False)
    memory('saved');assert out.stat().st_size<=512*1024**2
    report={'source':str(out),'source_sha256':sha(out),'source_bytes':out.stat().st_size,'baseline_sha256':BASE_SHA,'r1_scene_sha256':'8341139fb3427aedcf609d1eff2f3cdf9814819e9db048dddb2ef1bb4e2279ec','construction_changed':False,'region_half_width_m':HALF_REGION,'replaced_loops':int(selected.sum()),'groups':records,'retained_geometry_sha256':retained,'material_sha256':materials,'rig_sha256':rig_sha,'profile':str(PROFILE),'profile_sha256':sha(PROFILE),'capture':profile['capture'],'projection_guard':{'minimum_unreplaced_core_AABB_distance_from_macro_field_m':min(projected_guards),'by_color_m':projected_guards,'method':'Project every original outer core vertex into the exact frozen orthographic camera, then conservatively compare each projected core AABB to the 8 mm square field. All four groups have positive separation.','outside_material':'Original r5 cores remain, so this is a local construction comparison with baseline surrounding transport.'},'scope':'Actual material-scale local comparison only; every surrounding r5 object/core stays. Not a full-panel replacement or mechanics simulation.','limits':['Native body UV bump semantics differ from old mesh core.','Frozen r1 residual overlaps and tight crown folds remain; no new construction fix.','Only local macro is qualified by this source. Far surrounding cores retain r5 construction.'],'seconds':time.time()-started,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'source_code_sha256':sha(Path(__file__))}
    (ROOT/'receipts/handoff.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('source','source_sha256','source_bytes','replaced_loops','seconds','peak_rss_mib')}),flush=True)

if __name__=='__main__':run()
