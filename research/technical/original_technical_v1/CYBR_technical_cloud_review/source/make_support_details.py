import bpy
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1]
for id,w,h in [('22_petg',.045,.012),('14_chainmail',.109,.010)]:
    bpy.ops.wm.open_mainfile(filepath=str(ROOT/'scenes'/f'{id}_r5_hero.blend'));s=bpy.context.scene;cam=s.camera
    cam.location=(w*.27,-w*.37,w*.45+h);cam.data.ortho_scale=w*.46;target=(0,-w*.13,h);cam.rotation_euler=(Vector(target)-cam.location).to_track_quat('-Z','Y').to_euler()
    out=ROOT/'scenes'/f'{id}_r5_detail.blend';s.render.filepath=str(ROOT/'evidence'/f'{id}_r5_detail.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
