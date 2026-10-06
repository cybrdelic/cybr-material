import bpy,math
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1];bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes/22_petg_r5_hero.blend'));s=bpy.context.scene
a=-1.72;n=Vector((math.cos(a),math.sin(a),0));t=Vector((-math.sin(a),math.cos(a),0));target=n*.014+Vector((0,0,.013))
cam=s.camera;cam.location=n*.069+Vector((0,0,.016));cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler();cam.data.ortho_scale=.012
key=s.objects['Capture / key'];key.location=target+t*.045+n*.013+Vector((0,0,.025));key.rotation_euler=(target-key.location).to_track_quat('-Z','Y').to_euler();key.data.size=.012;key.data.size_y=.026;key.data.energy=.075
s.objects['Capture / fill'].data.energy*=.25;s.render.resolution_x=800;s.render.resolution_y=800;s['revision']='r6: seam-centered12mm macro with raking physical light; unchanged55um seam geometry and shader'
out=ROOT/'scenes/22_petg_r6_seamproof.blend';s.render.filepath=str(ROOT/'evidence/22_petg_r6_seamproof.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
