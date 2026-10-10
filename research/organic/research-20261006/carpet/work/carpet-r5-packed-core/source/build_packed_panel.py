"""Core-only full-panel replacement. No renderer or independently fitted studio."""
# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
import bpy,numpy as np,json,hashlib,resource,time,gc,sys
from pathlib import Path
from mathutils import Matrix
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'source'))
from inspect_baseline import read,mesh_sha,material_sha,NAMES
LIMIT_MIB=3072
def rss():return int(next(x.split()[1] for x in Path('/proc/self/status').read_text().splitlines() if x.startswith('VmRSS:')))/1024
def guard(label):
    value=rss();print(label,'rss_mib',value,flush=True)
    if value>LIMIT_MIB:raise RuntimeError(f'3 GiB build admission exceeded at {label}: {value} MiB')
    return value
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def apply_rig():
    path=ROOT/'receipts/reconstructed_baseline_rig.json'
    assert sha(path)=='b8afacac3814045f28b1ed835b2222bdfb50658aaceb1b4e53748310771ceccc'
    rig=json.loads(path.read_text());scene=bpy.context.scene
    for ob in list(scene.objects):
        if ob.type in {'LIGHT','CAMERA'} or ob.name.startswith(('Slab /','Studio /','Capture /')):bpy.data.objects.remove(ob,do_unlink=True)
    for name,desc in rig['objects'].items():
        typ=desc['type']
        if typ=='CAMERA':data=bpy.data.cameras.new(name)
        elif typ=='LIGHT':data=bpy.data.lights.new(name,desc['light']['type'])
        else:
            data=bpy.data.meshes.new(name);data.from_pydata(desc['vertices'],[],desc['faces'])
            m=bpy.data.materials.new(desc['material']['name']);m.use_nodes=True;p=m.node_tree.nodes.get('Principled BSDF')
            for key in ['Base Color','Roughness']:p.inputs[key].default_value=desc['material'][key]
            data.materials.append(m)
        ob=bpy.data.objects.new(name,data);scene.collection.objects.link(ob);ob.matrix_world=Matrix(desc['matrix_world'])
        if typ in {'CAMERA','LIGHT'}:
            for key,value in desc[typ.lower()].items():setattr(data,key,value)
        if typ=='CAMERA':scene.camera=ob
    scene.world=bpy.data.worlds.new(rig['world']['name']);scene.world.use_nodes=True;bg=scene.world.node_tree.nodes['Background'];bg.inputs[0].default_value=rig['world']['color'];bg.inputs[1].default_value=rig['world']['strength']
    for target,section in [(scene.view_settings,'view'),(scene.render,'render'),(scene.cycles,'cycles'),(scene.unit_settings,'unit'),(scene.render.image_settings,'image')]:
        for key,value in rig[section].items():setattr(target,key,value)
    scene.cycles.samples=128
    scene.cycles_curves.shape='THICK'
    return sha(path)
def run():
    start=time.time();baseline=Path(str(_required_research_input("baseline-scene")))
    assert sha(baseline)=='55df658fd44a0fb17bab14b46692c1c252b2d2658a967910e8c281f82e743879'
    bpy.ops.wm.open_mainfile(filepath=str(baseline));guard('baseline loaded')
    before={ob.name:mesh_sha(ob) for ob in bpy.context.scene.objects if ob.type=='MESH' and ob.name.startswith('Carpet /')}
    mats={m.name:material_sha(m) for ob in bpy.context.scene.objects if ob.name.startswith('Carpet /') for m in ob.data.materials}
    thin=[];changed=[]
    for color,name in enumerate(NAMES):
        ob=bpy.data.objects[name];old=ob.data;n=len(old.vertices)//1470
        coords=read(old.vertices,'co',3).reshape(n,1470,3)[:,294:,:].copy().reshape(-1,3)
        oldvi=read(old.loops,'vertex_index',dtype='i4').reshape(n,5820)[:,1164:].copy().ravel()
        vi=(oldvi//1470)*1176+oldvi%1470-294
        uv=read(old.uv_layers.active.data,'uv',2).reshape(n,5820,2)[:,1164:,:].copy().reshape(-1,2)
        smooth=read(old.polygons,'use_smooth',dtype='i4').reshape(n,1454)[:,290:].copy().ravel()
        mi=read(old.polygons,'material_index',dtype='i4').reshape(n,1454)[:,290:].copy().ravel()
        me=bpy.data.meshes.new(f'Carpet / retained exact six surface strands {color}');me.vertices.add(len(coords));me.vertices.foreach_set('co',coords.ravel());me.loops.add(len(vi));me.loops.foreach_set('vertex_index',vi)
        me.polygons.add(n*1164);me.polygons.foreach_set('loop_start',np.arange(n*1164,dtype='i4')*4);me.polygons.foreach_set('loop_total',np.full(n*1164,4,dtype='i4'));me.polygons.foreach_set('use_smooth',smooth);me.polygons.foreach_set('material_index',mi)
        layer=me.uv_layers.new(name=old.uv_layers.active.name);layer.data.foreach_set('uv',uv.ravel())
        for m in old.materials:me.materials.append(m)
        me.update(calc_edges=True);ob.data=me;bpy.data.meshes.remove(old)
        assert np.array_equal(coords,read(me.vertices,'co',3));assert np.array_equal(vi,read(me.loops,'vertex_index',dtype='i4'));assert np.array_equal(uv,read(me.uv_layers.active.data,'uv',2))
        ob['construction']='Exact retained r5 six surface strands only; smooth six-sided core removed. Packed native body is a separate object.'
        thin.append({'object':name,'loops':n,'vertices':len(coords),'polygons':n*1164,'retained_geometry_sha256':mesh_sha(ob),'exact_vertex_topology_uv_arrays_verified':True})
        changed.append(name);del coords,oldvi,vi,uv,smooth,mi;gc.collect();guard(f'core removal {color}')
    construction=json.loads((ROOT/'receipts/full_construction.json').read_text());curve_objects=[]
    for color in range(4):
        points=np.load(ROOT/f'arrays/colour_{color}_positions.npy',mmap_mode='r');n,f,k,_=points.shape;count=n*f
        radii=np.load(ROOT/f'arrays/colour_{color}_radius_per_loop.npy');name=f'Carpet / packed crimped body fibres {color}'
        data=bpy.data.hair_curves.new(name);data.add_curves([k]*count);data.position_data.foreach_set('vector',points.ravel());del points
        ar=data.attributes.new('radius','FLOAT','POINT');ar.data.foreach_set('value',np.repeat(radii,f*k));del radii
        typ=data.attributes.new('curve_type','INT8','CURVE');typ.data.foreach_set('value',np.zeros(count,dtype='i4'));res=data.attributes.new('resolution','INT','CURVE');res.data.foreach_set('value',np.full(count,3,dtype='i4'))
        ob=bpy.data.objects.new(name,data);bpy.context.scene.collection.objects.link(ob);data.materials.append(bpy.data.materials[f'Carpet / fine wool yarn {color}'])
        ob['construction']='Packed continuous constructed yarn; corrected inner crown and intrinsic crimp; not a solved mechanical state; measured residual local overlap remains.'
        ob['body_fibres_per_loop']=f;ob['points_per_fibre']=k;ob['fiber_uv_note']='Original Principled nodes unchanged. Native body curves have no original circumferential mesh UV; original UV bump semantics are not preserved on new bodies.'
        curve_objects.append({'object':name,'curves':count,'points':count*k,'geometry_role':'native Catmull-Rom fibre volume'});gc.collect();guard(f'curves {color}')
    retained={ob.name:mesh_sha(ob) for ob in bpy.context.scene.objects if ob.type=='MESH' and ob.name.startswith('Carpet /') and ob.name not in NAMES}
    assert all(before[name]==digest for name,digest in retained.items())
    after_mats={name:material_sha(bpy.data.materials[name]) for name in mats};assert mats==after_mats
    rig_sha=apply_rig();guard('rig applied')
    bpy.context.scene['constructed_yarn_scope']='Full 6 cm r5 panel; only body cores replaced. All body endpoints buried. Local residual packing violations remain. No full manufacturing simulation.'
    out=ROOT/'scenes/24_carpet_r5_packed_core_constructed_r1.blend';bpy.context.preferences.filepaths.save_version=0
    bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=False);guard('saved')
    report={'source':str(out),'source_sha256':sha(out),'baseline_sha256':sha(baseline),'rig_sha256':rig_sha,'retained_objects_geometry_sha256':retained,'retained_materials_sha256':after_mats,'retained_surface_strands':thin,'changed_original_objects':changed,'added_objects':curve_objects,'removed_core_vertices':1764*294,'removed_core_polygons':1764*290,'construction':construction,'packing_preflight':json.loads((ROOT/'receipts/packing_preflight.json').read_text()),'sampling_convergence':json.loads((ROOT/'receipts/sampling_convergence.json').read_text()),'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'last_rss_mib':rss(),'source_hashes':{p.name:sha(p) for p in (ROOT/'source').glob('*.py')},'cycles_curve_shape':bpy.context.scene.cycles_curves.shape,'cycles_curve_subdivisions':bpy.context.scene.cycles_curves.subdivisions,'material_limit':'Original Principled response and all material node graphs retained exactly. New native bodies do not carry original mesh circumferential UV semantics. No Hair BSDF.','output_stem':'24_carpet_r5_packed_core_constructed_r1'}
    report['all_endpoints_buried']=json.loads((ROOT/'receipts/all_endpoints_buried.json').read_text())
    report['uv_shader_dependency']={'chain':'Texture Coordinate.UV -> Vector Math MULTIPLY(32,2,1) -> Noise Texture(scale1.5,detail3) -> Bump(strength0.15,distance8micrometres) -> Principled.Normal','affected_input':'Principled.Normal only','unaffected_constant_inputs':'Original four Base Color values, Roughness0.96, Sheen Weight0.32 and all other original constant Principled inputs remain unchanged.','new_curve_uv':'No original circumferential mesh UV attribute on native body fibres; exact spatial bump evaluation is therefore not equivalent. Retained six surface strands keep exact UV arrays.'}
    (ROOT/'receipts/packed_core_handoff.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:report[k] for k in ['source','source_sha256','seconds','peak_rss_mib']}),flush=True)
if __name__=='__main__':run()
