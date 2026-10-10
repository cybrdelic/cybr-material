"""Save-only bounded actual-fibre material scene. No scene renderer is called."""
# Explicit external inputs replace machine-specific historical file locations.
from pathlib import Path as _ResearchInputPath
import sys as _research_input_sys
for _research_parent in _ResearchInputPath(__file__).resolve().parents:
    if (_research_parent / "research_inputs.py").is_file():
        _research_input_sys.path.insert(0, str(_research_parent))
        break
from research_inputs import required_input as _required_research_input
from pathlib import Path
import bpy,numpy as np,math,json,hashlib,time,resource
from mathutils import Vector
P=Path(__file__).resolve().parents[1];start=time.monotonic();meta=json.loads((P/'receipts/qualified_seven_geometry.json').read_text());z=np.load(meta['state']);points=z['points'];r=float(z['radius_m']);N=len(points);assert N==7 and meta['normal_contact']['max_penetration_m']<1e-6;bpy.ops.wm.read_factory_settings(use_empty=True);bpy.context.preferences.filepaths.save_version=0
selected=Path(str(_required_research_input("baseline-scene")))
with bpy.data.libraries.load(str(selected),link=False) as (a,b):b.materials=[n for n in a.materials if n.startswith('Carpet / fine wool yarn ')]
mats=[]
for i in range(4):
 m=bpy.data.materials.get('Carpet / fine wool yarn '+str(i));assert m;bs=m.node_tree.nodes.get('Principled BSDF')
 for link in list(bs.inputs['Normal'].links):m.node_tree.links.remove(link)
 m['provenance']='Selected r5 colour, roughness and sheen. Explicit fibre geometry replaces its yarn-scale normal overlay.';mats.append(m)
verts=[];faces=[];mi=[];sides=24;phi=np.arange(sides)*2*np.pi/sides
# The visible solid is the union representation of solved capsules. Cylinders
# and joint spheres overlap internally; no opaque bundle core exists.
for i,p in enumerate(points):
 colour=(i*17+i//9)%4
 for a,b in zip(p[:-1],p[1:]):
  t=b-a;t/=np.linalg.norm(t);axis=np.eye(3)[np.argmin(abs(t))];u=np.cross(t,axis);u/=np.linalg.norm(u);v=np.cross(t,u);ring=r*(np.cos(phi)[:,None]*u+np.sin(phi)[:,None]*v);base=len(verts);verts.extend((a+ring).tolist());verts.extend((b+ring).tolist())
  for k in range(sides):faces.append((base+k,base+(k+1)%sides,base+sides+(k+1)%sides,base+sides+k));mi.append(colour)
 for a in p:
  base=len(verts);verts.append((a+np.array([0,0,r])).tolist())
  for theta in np.arange(1,12)*np.pi/12:verts.extend((a+r*np.c_[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.full(sides,np.cos(theta))]).tolist())
  bottom=len(verts);verts.append((a-np.array([0,0,r])).tolist())
  for k in range(sides):faces.append((base,base+1+k,base+1+(k+1)%sides));mi.append(colour)
  for j in range(10):
   first=base+1+j*sides;other=first+sides
   for k in range(sides):faces.append((first+k,other+k,other+(k+1)%sides,first+(k+1)%sides));mi.append(colour)
  first=base+1+10*sides
  for k in range(sides):faces.append((first+k,bottom,first+(k+1)%sides));mi.append(colour)
mesh=bpy.data.meshes.new('Seven exact solved capsule centreline sets');mesh.from_pydata(verts,[],faces);mesh.update();mesh.polygons.foreach_set('use_smooth',np.ones(len(faces),dtype=bool));mesh.polygons.foreach_set('material_index',np.array(mi,dtype=np.int32))
for m in mats:mesh.materials.append(m)
o=bpy.data.objects.new('Seven actual crimped fibres; one free right end',mesh);bpy.context.collection.objects.link(o);o['radius_m']=r;o['state_sha256']=meta['state_sha256'];o['capsule_tessellation_radial_bound_m']=r*(1-math.cos(math.pi/24));o['opaque_bundle_core']=False

def mat(name,c,rough):
 m=bpy.data.materials.new(name);m.use_nodes=True;bs=m.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(*c,1);bs.inputs['Roughness'].default_value=rough;return m
floor=mat('Neutral fixture floor',(.16,.151,.136),.86);clamp=mat('Modelled fixed-root supports',(.085,.09,.081),.66)
def box(name,loc,dim,m):
 bpy.ops.mesh.primitive_cube_add(size=1,location=loc);o=bpy.context.object;o.name=name;o.dimensions=dim;bpy.ops.object.transform_apply(location=False,rotation=False,scale=True);o.data.materials.append(m);return o
box('Backing contact surface at solved z = 0',(.0012,0,-.00011),(.012,.012,.00022),floor)
for end,label in [(0,'left'),(-1,'right')]:
 ids=np.arange(N) if end==0 else np.flatnonzero(z['endpoint_fixed']);p=points[ids,end];lo=p.min(0)-r;hi=p.max(0)+r;depth=.00011;center=(lo+hi)/2;center[0]=float(np.mean(p[:,0]))+(-depth/2+.000008 if end==0 else depth/2-.000008);dims=hi-lo;dims[0]=depth;dims[1:]+=.000025;box('Explicit '+label+' root clamp support',center,dims,clamp)
scene=bpy.context.scene;scene.name='Qualified short wool-fibre yarn diagnostic';scene.render.engine='CYCLES';scene.cycles.samples=128;scene.cycles.use_denoising=False;scene.cycles.use_light_tree=False;scene.render.resolution_x=960;scene.render.resolution_y=768;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='OPEN_EXR_MULTILAYER';scene.render.image_settings.color_depth='32';scene.render.dither_intensity=0;scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.unit_settings.system='METRIC';scene.unit_settings.scale_length=1;scene.view_layers[0].cycles.denoising_store_passes=True
world=bpy.data.worlds.new('Neutral fibre studio');world.use_nodes=True;world.node_tree.nodes['Background'].inputs[0].default_value=(.18,.21,.25,1);world.node_tree.nodes['Background'].inputs[1].default_value=.22;scene.world=world
def aim(o,target):o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
target=(.0012,-.00018,.00007);scale=.003
for name,loc,size,power,color in [('Key',(-.0009,-.002,.0042),.00315,36*scale*scale,(1,.92,.82)),('Rim',(.00285,.00075,.0018),.00225,22*scale*scale,(.8,.9,1)),('Fill',(.0036,-.00165,.00165),.0036,5*scale*scale,(1,1,1))]:
 d=bpy.data.lights.new('Capture / '+name,'AREA');d.energy=power;d.size=size;d.shape='DISK';d.color=color;o=bpy.data.objects.new(d.name,d);bpy.context.collection.objects.link(o);o.location=loc;aim(o,target)
d=bpy.data.cameras.new('Capture / actual short yarn');d.type='ORTHO';d.ortho_scale=scale;d.clip_start=.00001;d.clip_end=1.;cam=bpy.data.objects.new(d.name,d);bpy.context.collection.objects.link(cam);cam.location=(.0012,-.00155,.0055);aim(cam,target);scene.camera=cam;scene['scope']=meta['scope'];scene['field_m']=scale;scene['fibre_diameter_m']=2*r;scene['preserved_r5_source_sha256']=hashlib.sha256(selected.read_bytes()).hexdigest();source=P/'assets/seven_fibre_short_yarn_r01.blend';scene.render.filepath=str(P/'assets/SHORT_YARN_guides.exr');text=bpy.data.texts.new('Model provenance.json');text.write(json.dumps(meta,indent=2));bpy.ops.wm.save_as_mainfile(filepath=str(source),compress=True);sha=hashlib.sha256(source.read_bytes()).hexdigest();job=dict(name='CARPET_SHORT_YARN_SEVEN',source=str(source),source_sha256=sha,scene=scene.name,camera=cam.name,engine='CYCLES',width=960,height=768,samples=128,save_guides=True,light_tree=False,output=scene.render.filepath,orthographic_scale_m=scale,intervention='Solved intrinsic crimp, director-aware local contact friction and actual connected free endpoint; selected r5 dye scalars on actual fibre surfaces.',shared_controls='One source and one physical state. No exact A/B claim.');out=dict(status='ready_for_capture',source=str(source),source_sha256=sha,source_bytes=source.stat().st_size,geometry_state_sha256=meta['state_sha256'],fibres=N,radius_m=r,stock_length_m=.0024,field_m=scale,geometry=meta,source_frozen=True,rendered=False,scope=meta['scope'],jobs=[job],build_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,mesh_vertices=len(mesh.vertices),mesh_faces=len(mesh.polygons));(P/'receipts/seven_short_yarn_ready.json').write_text(json.dumps(out,indent=2));print('READY',json.dumps(out),flush=True)
