"""Preserve r4 optical models; move the broad reflected fill off the chart."""
import bpy,json,hashlib
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1]
for id in ('18_lcd','19_crt'):
    path=ROOT/'scenes'/f'{id}_r4_hero.blend';bpy.ops.wm.open_mainfile(filepath=str(path));s=bpy.context.scene
    fill=bpy.data.objects['Capture / fill'];old=list(fill.location);w=s.camera.data.ortho_scale/1.36;target=(0,0,.004)
    fill.location=(w*1.1,w*.5,w*.9+.004);fill.rotation_euler=(Vector(target)-fill.location).to_track_quat('-Z','Y').to_euler();fill.data.energy*=.5
    s.cycles.use_denoising=False;s.cycles.samples=256;s.cycles.adaptive_min_samples=64;s.cycles.adaptive_threshold=.012
    s.render.resolution_x=1536;s.render.resolution_y=1280
    out=ROOT/'scenes'/f'{id}_r5_hero.blend';s.render.filepath=str(ROOT/'evidence'/f'{id}_r5_hero.png')
    s['revision']='r5: fill-light position/power and render sampling only; geometry and materials unchanged'
    bpy.ops.wm.save_as_mainfile(filepath=str(out))
    (ROOT/'tests'/f'{id}_r5_patch.json').write_text(json.dumps({'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'output':str(out),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'fill_before':old,'fill_after':list(fill.location),'optical_materials_changed':False,'geometry_changed':False},indent=2)+'\n')
    print('READY_SCENE',out,flush=True)
