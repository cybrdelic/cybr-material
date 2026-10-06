"""Non-mutating EEVEE/Cycles review renderer. -- JSON job or list. Source geometry/material preserved unless job records an explicit adaptation."""
import bpy,sys,json,time,hashlib,pathlib,math
from mathutils import Vector
config=json.load(open(sys.argv[sys.argv.index('--')+1]));jobs=config if isinstance(config,list) else [config]
for job in jobs:
 source=pathlib.Path(job['source']);before=hashlib.sha256(source.read_bytes()).hexdigest();bpy.ops.wm.open_mainfile(filepath=str(source))
 if job.get('scene'):bpy.context.window.scene=bpy.data.scenes[job['scene']]
 s=bpy.context.scene
 if job.get('camera'):s.camera=bpy.data.objects[job['camera']]
 assert s.camera,'No active camera'
 missing=[im.filepath for im in bpy.data.images if im.source=='FILE' and im.filepath and not im.packed_file and not pathlib.Path(bpy.path.abspath(im.filepath)).is_file()]
 if missing:raise RuntimeError('Missing images: '+str(missing))
 adaptations=['Engine/resolution/samples changed in memory only']
 if job.get('camera_focus_prefix'):
  bpy.context.view_layer.update();obs=[o for o in s.objects if o.name.startswith(job['camera_focus_prefix'])];assert obs,'No focus objects'
  pts=[o.matrix_world@Vector(p) for o in obs for p in o.bound_box];center=sum(pts,Vector())/len(pts);s.camera.location=center+Vector(job.get('focus_offset',[.3,-1.05,.3]));s.camera.rotation_euler=(center-s.camera.location).to_track_quat('-Z','Y').to_euler();adaptations.append('Object-bounds detail camera framing')
 if 'ortho_scale' in job:s.camera.data.ortho_scale=job['ortho_scale'];adaptations.append('Orthographic detail framing')
 if 'camera_location' in job:s.camera.location=job['camera_location'];s.camera.rotation_euler=(Vector(job['camera_target'])-s.camera.location).to_track_quat('-Z','Y').to_euler();adaptations.append('Preview camera pose changed')
 if 'lens' in job:s.camera.data.lens=job['lens']
 if 'key_location' in job:
  for o in s.objects:
   if o.type=='LIGHT' and 'key' in o.name.lower():o.location=job['key_location'];o.rotation_euler=(Vector(job.get('key_target',[0,0,0]))-o.location).to_track_quat('-Z','Y').to_euler();adaptations.append('Key moved for raking-light proof')
 if 'light_size' in job:
  for o in s.objects:
   if o.type=='LIGHT' and 'key' in o.name.lower():o.data.size=job['light_size'];adaptations.append('Key size changed for hard-light proof')
 s.render.engine=job.get('engine','BLENDER_EEVEE_NEXT')
 if s.render.engine=='CYCLES':
  s.cycles.device='CPU';s.cycles.samples=job.get('samples',128);s.cycles.use_denoising=False;s.cycles.use_adaptive_sampling=True;s.cycles.adaptive_threshold=.02;s.cycles.adaptive_min_samples=64;adaptations.append('Cycles CPU optical proof, denoising disabled because unavailable')
 else:s.eevee.taa_render_samples=job.get('samples',64)
 if hasattr(s.eevee,'use_raytracing'):s.eevee.use_raytracing=True
 for m in bpy.data.materials:
  if m.use_nodes and any(n.type=='BSDF_PRINCIPLED' and n.inputs.get('Transmission Weight') and n.inputs['Transmission Weight'].default_value>0 for n in m.node_tree.nodes):
   if hasattr(m,'use_screen_refraction'):m.use_screen_refraction=True
 width=job.get('width',960);height=job.get('height',round(width*s.render.resolution_y/s.render.resolution_x));s.render.resolution_x=width;s.render.resolution_y=height;s.render.resolution_percentage=100
 out=pathlib.Path(job['output']);out.parent.mkdir(parents=True,exist_ok=True);s.render.filepath=str(out);s.render.image_settings.file_format='PNG'
 start=time.monotonic();bpy.ops.render.render(write_still=True);elapsed=time.monotonic()-start
 after=hashlib.sha256(source.read_bytes()).hexdigest();assert after==before
 report={'job':job,'source_sha256':before,'render_seconds':elapsed,'engine':s.render.engine,'blender_version':bpy.app.version_string,'adaptations':adaptations,'camera':s.camera.name,'camera_location':list(s.camera.location),'camera_rotation':list(s.camera.rotation_euler),'missing_images':missing,'objects':len(s.objects),'resolution':[width,height]}
 out.with_suffix(out.suffix+'.json').write_text(json.dumps(report,indent=2));print('VIEW_COMPLETE '+str(out)+' '+str(round(elapsed,2))+' s',flush=True)
