"""Closed flat-interior coupons isolating geometry sampling, with fixed state/optics.
Reconstructed new studies; not exact recovery of a historical scene or a full molded part.
"""
import bpy,sys,json,hashlib,time,resource
import numpy as np
from pathlib import Path
from mathutils import Vector
R=Path(__file__).resolve().parents[1];ROOT=R.parents[2];sys.path.insert(0,str(R/'source'));from studio import configure
CASE=sys.argv[sys.argv.index('--')+1];assert CASE in {'coarse','native'}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();rec=json.loads((R/'receipts/native_recovery.json').read_text());assert rec['tool_exact_archive_bytes'] and rec['production_state_exact_archive_bytes'] and rec['both_native_maps_exact_archive_bytes']
t=time.monotonic();data=np.load(R/'state/interior_geometry_fields.npz');key='coarse' if CASE=='coarse' else 'fine';a=data[key+'_axis_m'];h=data[key+'_height_m'];n=len(a);assert h.shape==(n,n);P=4*(n-1);nt=n*n
bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene
src=R/'inputs/21_plastic_texture_r5_hero.blend'
with bpy.data.libraries.load(str(src),link=False) as (old,new):new.materials=['Plastic / Molded stipple']
mat=new.materials[0];assert mat;mat.name='Plastic / fixed production state / unchanged optical closure';nodes=mat.node_tree.nodes;links=mat.node_tree.links;p=nodes['Principled BSDF']
scalar={i.name:list(i.default_value) if hasattr(i.default_value,'__len__') else i.default_value for i in p.inputs if hasattr(i,'default_value')}
assert scalar['Base Color']==[.10999999940395355,.14000000059604645,.1599999964237213,1.]
for node in list(nodes):
 if node.type not in {'OUTPUT_MATERIAL','BSDF_PRINCIPLED'}:nodes.remove(node)
imgpath=R/'maps/production_state_roughness_native4096.png';assert sha(imgpath)==rec['native_roughness_sha256'];tex=nodes.new('ShaderNodeTexImage');tex.name='Fixed native4096 source roughness';tex.image=bpy.data.images.load(str(imgpath));tex.image.colorspace_settings.name='Non-Color';tex.interpolation='Linear';tex.extension='EXTEND';assert list(tex.image.size)==[4096,4096]
geo=nodes.new('ShaderNodeNewGeometry');scale=nodes.new('ShaderNodeVectorMath');scale.operation='SCALE';scale.inputs[3].default_value=12.5;offset=nodes.new('ShaderNodeVectorMath');offset.operation='ADD';offset.inputs[1].default_value=(.5,.5,0);links.new(geo.outputs['Position'],scale.inputs[0]);links.new(scale.outputs[0],offset.inputs[0]);links.new(offset.outputs[0],tex.inputs['Vector']);links.new(tex.outputs['Color'],p.inputs['Roughness'])
# One connected watertight mesh. Structured top quads; vertical cut sides; bottom fan.
verts=np.empty((nt+P+1,3),'f4');verts[:nt,0]=np.tile(a,n);verts[:nt,1]=np.repeat(a,n);verts[:nt,2]=(.008+h).ravel()
per=np.concatenate([np.arange(n),np.arange(2*n-1,nt,n),np.arange(nt-2,nt-n-1,-1),np.arange(nt-2*n,0,-n)]).astype('i4');assert len(per)==P and len(np.unique(per))==P
verts[nt:nt+P]=verts[per];verts[nt:nt+P,2]=.005;verts[-1]=(0,0,.005)
i=np.arange(n-1,dtype='i4')[None,:]+np.arange(n-1,dtype='i4')[:,None]*n;i=i.ravel();q=np.stack([i,i+1,i+n+1,i+n],axis=1);bot=np.arange(nt,nt+P,dtype='i4');side=np.stack([per,bot,np.roll(bot,-1),np.roll(per,-1)],axis=1);bottom=np.stack([np.full(P,nt+P,'i4'),np.roll(bot,-1),bot],axis=1)
indices=np.concatenate([q.ravel(),side.ravel(),bottom.ravel()]);counts=np.concatenate([np.full(len(q)+P,4,'i4'),np.full(P,3,'i4')]);starts=np.cumsum(np.r_[0,counts[:-1]],dtype='i4')
mesh=bpy.data.meshes.new('Fixed physical relief / '+CASE);mesh.vertices.add(len(verts));mesh.vertices.foreach_set('co',verts.ravel());mesh.loops.add(len(indices));mesh.loops.foreach_set('vertex_index',indices);mesh.polygons.add(len(counts));mesh.polygons.foreach_set('loop_start',starts);mesh.polygons.foreach_set('loop_total',counts);smooth=np.zeros(len(counts),dtype=bool);smooth[:len(q)]=True;mesh.polygons.foreach_set('use_smooth',smooth);mesh.update(calc_edges=True);assert not mesh.validate(verbose=False)
ob=bpy.data.objects.new('21 / flat production interior / '+CASE,mesh);s.collection.objects.link(ob);mesh.materials.append(mat)
# Half-edge counts prove the actual mesh closes. Opposite directed incidences prove orientation.
edges=np.empty((len(indices),2),'i4');edges[:,0]=indices;ends=starts+counts-1;edges[:,1]=np.roll(indices,-1);edges[ends,1]=indices[starts];lo=np.minimum(edges[:,0],edges[:,1]).astype('i8');hi=np.maximum(edges[:,0],edges[:,1]).astype('i8');codes=lo*len(verts)+hi;order=np.argsort(codes);sortedcodes=codes[order];assert np.all(sortedcodes[::2]==sortedcodes[1::2]);assert np.all(sortedcodes[1:-1:2]!=sortedcodes[2::2]);assert np.all(edges[order[::2],0]==edges[order[1::2],1]);del edges,lo,hi,codes,order,sortedcodes
norm=np.empty((len(mesh.polygons),3),'f4');mesh.polygons.foreach_get('normal',norm.ravel());assert np.isfinite(norm).all() and np.min(np.linalg.norm(norm,axis=1))>.999;assert norm[:len(q),2].min()>0 and np.all(norm[-P:,2]<-.999)
# Use fixed state bounds for identical studios even if the coarse mesh misses a peak.
loz=.005;hiz=.008+float(data['fine_height_m'].max());width=float(a[-1]-a[0]);bpy.ops.mesh.primitive_cube_add(size=1,location=(0,0,(loz+hiz)/2));bound=bpy.context.object;bound.name='Temporary fixed state bounds';bound.dimensions=(width,width,hiz-loz);bpy.ops.object.transform_apply(location=False,rotation=False,scale=True);report=configure('21_plastic_texture',[bound]);bpy.data.objects.remove(bound,do_unlink=True)
s.camera.name='Plastic / matched physical interior';s.camera.data.ortho_scale=.028
for obj in s.objects:
 if obj.type=='LIGHT' and obj.data.type=='AREA':obj.rotation_euler=(0,0,0)
s.render.dither_intensity=0;s.render.image_settings.color_depth='16';s.cycles.samples=128;s.cycles.use_adaptive_sampling=True;s.cycles.adaptive_threshold=.01
ob['catalog_id']='21_plastic_texture';ob['geometry_pitch_m']=float(a[1]-a[0]);ob['native_field_domain_m']=.08;ob['physical_state_sha256']=rec['production_state_sha256'];ob['scope']='Flat interior sampling study; cut edges are artificial coupon boundaries';mat['known_limitation']=rec['roughness_limitation'];s['study_case']=CASE;s['selected']=False;s['physical_qualification']=False
rig={o.name:{'type':o.type,'matrix':[v for row in o.matrix_world for v in row],'data':{'energy':o.data.energy,'size':o.data.size,'color':list(o.data.color),'shape':o.data.shape} if o.type=='LIGHT' else {'ortho_scale':o.data.ortho_scale,'lens':o.data.lens} if o.type=='CAMERA' else {}} for o in s.objects if o!=ob}
# Stored float32 heights and source sample precision checked independently from saved coordinates.
maxerr=float(abs(verts[:nt,2].astype('f8')-(.008+h.ravel())).max());assert maxerr<1e-9
report.update(case=CASE,source_material_scene=str(src),source_material_scene_sha256=sha(src),source_material_principled_defaults=scalar,roughness_path=str(imgpath),roughness_sha256=sha(imgpath),roughness_limitation=rec['roughness_limitation'],geometry_grid=[n,n],geometry_pitch_m=float(a[1]-a[0]),coupon_xy_extent_m=width,nominal_cut_thickness_m=.003,top_height_origin_m=.008,vertex_count=len(verts),face_count=len(counts),closed_oriented_manifold=True,height_float32_max_error_m=maxerr,physical_state_sha256=rec['production_state_sha256'],rig_signature=rig,source_recovery='New reconstructed flat-interior study from byte-exact recovered original native state; not historical blend recovery',quality_acceptance=False)
(R/'scenes').mkdir(exist_ok=True);path=R/'scenes'/('21_plastic_fixed_state_'+CASE+'.blend');bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True);report.update(scene=str(path),scene_sha256=sha(path),wall_seconds=time.monotonic()-t,peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024);(R/'receipts'/('scene_'+CASE+'.json')).write_text(json.dumps(report,indent=2));print('PLASTIC_READY',CASE,report['scene_sha256'],report['wall_seconds'],flush=True)
