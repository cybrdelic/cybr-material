"""SI-unit neutral specimen studio. Never scales source geometry or material coordinates."""
import bpy,hashlib
from pathlib import Path
from mathutils import Vector

def aim(o,t):o.rotation_euler=(Vector(t)-o.location).to_track_quat('-Z','Y').to_euler()
def configure(specimen_id,objects=None):
 s=bpy.context.scene
 for o in list(s.objects):
  if o.type in {'LIGHT','CAMERA'} or o.name.startswith(('Studio /','Capture /')):bpy.data.objects.remove(o,do_unlink=True)
 obs=objects or [o for o in s.objects if o.type in {'MESH','CURVE','CURVES','SURFACE'} and not o.hide_render]
 bpy.context.view_layer.update()
 pts=[o.matrix_world@Vector(p) for o in obs for p in o.bound_box]
 assert pts,'No specimen bounds'
 lo=Vector(tuple(min(p[i] for p in pts) for i in range(3)));hi=Vector(tuple(max(p[i] for p in pts) for i in range(3)));center=(lo+hi)/2;w=max(hi.x-lo.x,hi.y-lo.y,hi.z-lo.z)
 assert .001<w<10,(specimen_id,w)
 # Same relative camera angle and frame, preserving each physical specimen.
 bpy.ops.object.camera_add(location=center+Vector((.82,-1.10,1.02))*w);cam=bpy.context.object;cam.name='Slab / shared camera';cam.data.type='ORTHO';cam.data.ortho_scale=w*1.62;cam.data.clip_start=.00001;cam.data.clip_end=100;aim(cam,center);s.camera=cam
 bpy.ops.mesh.primitive_plane_add(size=w*200,location=(center.x,center.y,lo.z-0.000001));floor=bpy.context.object;floor.name='Slab / neutral floor'
 mat=bpy.data.materials.new('Slab / neutral floor');mat.use_nodes=True;p=mat.node_tree.nodes.get('Principled BSDF');p.inputs['Base Color'].default_value=(.18,.18,.18,1);p.inputs['Roughness'].default_value=.85;floor.data.materials.append(mat)
 for name,loc,size,power in [('key',(-.70,-.60,1.30),.95,28),('rim',(.65,.35,.80),.65,14),('fill',(.6,-.65,.50),1.2,5)]:
  d=bpy.data.lights.new('Slab / '+name,'AREA');d.energy=power*w*w;d.shape='DISK';d.size=size*w;d.color=(1,1,1);o=bpy.data.objects.new(d.name,d);s.collection.objects.link(o);o.location=center+Vector(loc)*w;aim(o,center)
 s.world=bpy.data.worlds.new('Slab / neutral world');s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs[0].default_value=(.18,.18,.18,1);s.world.node_tree.nodes['Background'].inputs[1].default_value=.22
 s.view_settings.view_transform='AgX';s.view_settings.look='AgX - Medium High Contrast';s.view_settings.exposure=0;s.view_settings.gamma=1
 s.render.engine='CYCLES';s.cycles.device='CPU';s.cycles.use_denoising=False;s.cycles.use_light_tree=False;s.cycles.max_bounces=12;s.cycles.transmission_bounces=12
 s.render.resolution_x=960;s.render.resolution_y=960;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.render.image_settings.color_mode='RGB';s.render.image_settings.color_depth='16';s.render.film_transparent=False;s.unit_settings.system='METRIC';s.unit_settings.scale_length=1
 s['slab_studio_revision']='r2: ground at specimen minimum minus1micrometre';s['slab_catalog_id']=specimen_id;s['slab_bounds_m']=[*lo,*hi];s['slab_scope']='Actual SI-scale selected material specimen, unified neutral studio; no physical microstructure scaling'
 return {'studio_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'ground_clearance_m':0.000001,'id':specimen_id,'bounds_min_m':list(lo),'bounds_max_m':list(hi),'dimensions_m':list(hi-lo),'frame_width_m':cam.data.ortho_scale,'specimen_objects':[o.name for o in obs]}
