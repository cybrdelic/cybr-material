import bpy,json,hashlib
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1];src=ROOT/'scenes/25_tile_ceramic_r4_hero.blend'
bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
cam=s.camera;cam.location=(.16,-.245,.132);target=(0,-.025,.009);cam.rotation_euler=(Vector(target)-cam.location).to_track_quat('-Z','Y').to_euler();cam.data.ortho_scale=.25
light=bpy.data.lights.new('Capture / narrow glaze strip','AREA');light.shape='RECTANGLE';light.size=.19;light.size_y=.013;light.energy=.65
ob=bpy.data.objects.new(light.name,light);s.collection.objects.link(ob);ob.location=(-.14,.24,.13);ob.rotation_euler=(Vector(target)-ob.location).to_track_quat('-Z','Y').to_euler()
s['revision']='r5 glaze-proof camera and narrow strip light; unchanged materials/geometry';s.render.resolution_x=960;s.render.resolution_y=800
out=ROOT/'scenes/25_tile_ceramic_r5_grazing.blend';s.render.filepath=str(ROOT/'evidence/25_tile_ceramic_r5_grazing.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
(ROOT/'tests/ceramic_glaze_r5_patch.json').write_text(json.dumps({'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'shader_changed':False,'geometry_changed':False,'camera_changed':True,'light_added':'0.65W narrow strip'},indent=2)+'\n')
