"""Non-render API probe: resolution cache versus unchanged native curve keys."""
import bpy,numpy as np,hashlib,json,sys,time,resource,gc
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
resolution=int(sys.argv[sys.argv.index('--')+1]);start=time.time()
def rss():return int(next(x.split()[1] for x in Path('/proc/self/status').read_text().splitlines() if x.startswith('VmRSS:')))/1024
def arrays(data):
    p=np.empty(len(data.points)*3,dtype='f4');data.position_data.foreach_get('vector',p)
    r=np.empty(len(data.points),dtype='f4');data.attributes['radius'].data.foreach_get('value',r)
    return p,r
def digest(p,r):return hashlib.sha256(p.tobytes()+r.tobytes()).hexdigest()
for ob in list(bpy.data.objects):bpy.data.objects.remove(ob,do_unlink=True)
baseline_rss=rss();data0=np.load(ROOT/'arrays/colour_0_positions.npy',mmap_mode='r')
indices=np.array(list(range(30))+[230,672]);positions=data0[indices].copy();radii=np.load(ROOT/'arrays/colour_0_radius_per_loop.npy')[indices]
del data0
n,f,k,_=positions.shape;count=n*f;expected_p=positions.ravel();expected_r=np.repeat(radii,f*k);expected_hash=digest(expected_p,expected_r)
data=bpy.data.hair_curves.new('Native cache probe');data.add_curves([k]*count);data.position_data.foreach_set('vector',expected_p)
data.attributes.new('radius','FLOAT','POINT').data.foreach_set('value',expected_r)
data.attributes.new('curve_type','INT8','CURVE').data.foreach_set('value',np.zeros(count,dtype='i4'))
data.attributes.new('resolution','INT','CURVE').data.foreach_set('value',np.full(count,resolution,dtype='i4'))
ob=bpy.data.objects.new('Native cache probe',data);bpy.context.scene.collection.objects.link(ob)
native_rss=rss();bpy.context.view_layer.update();dg=bpy.context.evaluated_depsgraph_get();ev=ob.evaluated_get(dg)
bound=np.array(ev.bound_box,dtype='f8');bbox_rss=rss();ep,er=arrays(ev.data);assert digest(ep,er)==expected_hash
del ep,er
# Evaluate the normal field on the point domain, the same field used by Cycles'
# optional ATTR_STD_VERTEX_NORMAL export. This is not a material or scene change.
group=bpy.data.node_groups.new('Probe point normals','GeometryNodeTree')
group.interface.new_socket(name='Geometry',in_out='INPUT',socket_type='NodeSocketGeometry')
group.interface.new_socket(name='Geometry',in_out='OUTPUT',socket_type='NodeSocketGeometry')
inn=group.nodes.new('NodeGroupInput');out=group.nodes.new('NodeGroupOutput')
store=group.nodes.new('GeometryNodeStoreNamedAttribute');store.domain='POINT';store.data_type='FLOAT_VECTOR';store.inputs['Name'].default_value='probe_normal'
normal=group.nodes.new('GeometryNodeInputNormal');group.links.new(inn.outputs['Geometry'],store.inputs['Geometry']);group.links.new(normal.outputs['Normal'],store.inputs['Value']);group.links.new(store.outputs['Geometry'],out.inputs['Geometry'])
modifier=ob.modifiers.new('Probe only','NODES');modifier.node_group=group
bpy.context.view_layer.update();dg=bpy.context.evaluated_depsgraph_get();ev=ob.evaluated_get(dg)
normal_values=np.empty(len(ev.data.points)*3,dtype='f4');ev.data.attributes['probe_normal'].data.foreach_get('vector',normal_values)
np.save(HERE/f'point_normals_resolution_{resolution}.npy',normal_values.reshape(-1,3))
ep,er=arrays(ev.data);assert digest(ep,er)==expected_hash
report={'resolution':resolution,'loops':n,'fibres':count,'authored_native_points':len(data.points),'post_evaluation_native_points':len(ev.data.points),'raw_positions_and_radii_sha256':expected_hash,'evaluated_positions_and_radii_sha256':digest(ep,er),'native_export_inputs_bitwise_invariant':True,'cycles_export_evidence':'Blender4.3.2 export_hair_curves reads b_curves.positions() and radius attribute directly; no resolution-dependent evaluated-position substitution. This probe checks the exact evaluated Curves data consumed by that source path, not private Cycles device buffers.','bounds_min_m':bound.min(0).tolist(),'bounds_max_m':bound.max(0).tolist(),'rss_start_mib':baseline_rss,'rss_native_mib':native_rss,'rss_after_bounds_mib':bbox_rss,'bounds_added_mib':bbox_rss-native_rss,'rss_after_normal_field_mib':rss(),'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'seconds':time.time()-start,'no_render':True}
assert report['peak_rss_mib']<1024 and report['seconds']<120
(HERE/f'resolution_cache_{resolution}.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
