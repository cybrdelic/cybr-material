import bpy,numpy as np,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def rss():return int(next(x.split()[1] for x in Path('/proc/self/status').read_text().splitlines() if x.startswith('VmRSS:')))/1024
for ob in list(bpy.data.objects):bpy.data.objects.remove(ob,do_unlink=True)
p=np.load(ROOT/'arrays/colour_0_positions.npy',mmap_mode='r')[:32].copy();n,f,k,_=p.shape
d=bpy.data.hair_curves.new('In-place cache probe');d.add_curves([k]*(n*f));d.position_data.foreach_set('vector',p.ravel())
rr=np.repeat(np.load(ROOT/'arrays/colour_0_radius_per_loop.npy')[:32],f*k)
d.attributes.new('radius','FLOAT','POINT').data.foreach_set('value',rr)
d.attributes.new('curve_type','INT8','CURVE').data.foreach_set('value',np.zeros(n*f,dtype='i4'))
at=d.attributes.new('resolution','INT','CURVE');at.data.foreach_set('value',np.full(n*f,3,dtype='i4'))
ob=bpy.data.objects.new('In-place cache probe',d);bpy.context.scene.collection.objects.link(ob);bpy.context.view_layer.update()
bb0=np.array(ob.evaluated_get(bpy.context.evaluated_depsgraph_get()).bound_box);rss0=rss()
d.attributes.remove(at);at=d.attributes.new('resolution','INT','CURVE');at.data.foreach_set('value',np.ones(n*f,dtype='i4'));d.update_tag();ob.update_tag();bpy.context.view_layer.update()
bb1=np.array(ob.evaluated_get(bpy.context.evaluated_depsgraph_get()).bound_box);rss1=rss()
q=np.empty(p.size,dtype='f4');d.position_data.foreach_get('vector',q);assert np.array_equal(q,p.ravel())
rnew=np.empty(len(d.points),dtype='f4');d.attributes['radius'].data.foreach_get('value',rnew);assert np.array_equal(rr,rnew)
types=np.empty(n*f,dtype='i4');d.attributes['curve_type'].data.foreach_get('value',types);assert np.all(types==0)
expected_hash=hashlib.sha256(p.tobytes()+rr.tobytes()+np.zeros(n*f,dtype='i4').tobytes()).hexdigest();actual_hash=hashlib.sha256(q.tobytes()+rnew.tobytes()+types.tobytes()).hexdigest();assert expected_hash==actual_hash
report={'change':'remove/recreate resolution attribute from3 to1 after evaluating bounds, then data.update_tag and object.update_tag','points_unchanged':True,'radii_unchanged':True,'curve_types_unchanged':True,'before_sha256':expected_hash,'after_sha256':actual_hash,'bbox_difference_um':float(np.abs(bb1-bb0).max()*1e6),'rss_before_mib':rss0,'rss_after_mib':rss1,'rss_delta_mib':rss1-rss0,'bounds_changed':not np.array_equal(bb0,bb1)}
(HERE/'inplace_resolution_recreated_attribute.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
