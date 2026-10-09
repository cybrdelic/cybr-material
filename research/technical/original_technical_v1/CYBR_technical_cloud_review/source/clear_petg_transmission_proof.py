"""A side-on external chart tests transmission through both real printed walls."""
import bpy,sys,json,hashlib
from pathlib import Path
from mathutils import Vector
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import plain,box
ROOT=Path(__file__).resolve().parents[1];src=ROOT/'scenes/23_petg_transparent_r5_hero.blend';bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
for o in list(s.objects):
    if o.name.startswith('Optics / external background'):bpy.data.objects.remove(o,do_unlink=True)
for i,c in enumerate([(.009,.012,.016),(.78,.75,.68),(.42,.025,.012)]):
    box('Optics / side-on contrast chart '+str(i),((i-1)*.018,.026,.021),(.018,.001,.042),plain('Optics / matte chart ink '+str(i),c,.86))
ink=plain('Optics / fine dark chart bars',(.007,.01,.012),.92)
for i in range(7):box('Optics / chart horizontal bar '+str(i),(0,.0254,.004+i*.005),(.054,.0001,.0008),ink)
s.camera.location=(.003,-.09,.023);target=(0,0,.014);s.camera.rotation_euler=(Vector(target)-s.camera.location).to_track_quat('-Z','Y').to_euler();s.camera.data.ortho_scale=.06
s.render.resolution_x=960;s.render.resolution_y=720;s.cycles.samples=512;s.cycles.adaptive_min_samples=128;s.cycles.adaptive_threshold=.008;s.cycles.use_denoising=False
s['revision']='r6: side-on external reference-chart geometry/camera only, unchanged specimen and optical shader'
s['optical_expectation']='FDM layered PETG is rough and refracting. Test external color/bar transmission, not sharp optical-glass imaging.'
out=ROOT/'scenes/23_petg_transparent_r6_sideproof.blend';s.render.filepath=str(ROOT/'evidence/23_petg_transparent_r6_sideproof.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
(ROOT/'tests/clear_petg_r6_patch.json').write_text(json.dumps({'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output':str(out),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'specimen_geometry_changed':False,'material_changed':False,'roughness_map_range':[.160784,.317647],'roughness_map_mean':.18898,'changed':'Side-on camera and wide upright reference chart only'},indent=2)+'\n')
