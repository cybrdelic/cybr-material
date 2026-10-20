"""Add a supplementary native Cycles view from an existing editable scene.

Never modifies the input .blend, accepted PNGs, map set or selection registry.
Example (Blender 4.5.3 LTS):
  blender -b blender/CYBR_Cycles_Details.blend --python tools/gallery/render_variant.py -- \
    --material 05_champagne_brass --size 4096 --samples 128 --minimum 32
"""
import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import bpy
from mathutils import Matrix

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--material', default='05_champagne_brass')
parser.add_argument('--size', type=int, default=4096)
parser.add_argument('--samples', type=int, default=128)
parser.add_argument('--minimum', type=int, default=32)
args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
assert args.size >= 1024 and args.samples >= args.minimum >= 1
source = Path(bpy.data.filepath)
assert source.resolve() == (ROOT/'blender/CYBR_Cycles_Details.blend').resolve()
spec = next(s for s in json.loads((ROOT/'materials/manifest.json').read_text())['materials'] if s['id']==args.material)
scene = next(s for s in bpy.data.scenes if s.name=='CYCLES / Detail / '+spec['name'])
assert scene.render.engine=='CYCLES'
bpy.context.window.scene=scene
rotation=Matrix.Rotation(math.radians(90),4,'Z')
lights=[]
for ob in scene.objects:
    if ob.type=='LIGHT':
        lights.append(dict(name=ob.name,before=[list(row) for row in ob.matrix_world]))
        ob.matrix_world=rotation@ob.matrix_world
        lights[-1]['after']=[list(row) for row in ob.matrix_world]
assert lights, 'No existing scene lights; do not substitute a new rig'
dependencies=[]
for im in bpy.data.images:
    if im.source=='FILE' and im.filepath:
        p=Path(bpy.path.abspath(im.filepath))
        assert p.is_file(), str(p)
        assert not im.packed_file, 'Expected portable relative map dependencies'
        dependencies.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
scene.cycles.device='CPU'
scene.cycles.samples=args.samples
scene.cycles.use_adaptive_sampling=True
scene.cycles.adaptive_min_samples=args.minimum
scene.cycles.adaptive_threshold=.01
scene.cycles.use_denoising=True
scene.cycles.denoiser='OPENIMAGEDENOISE'
scene.cycles.denoising_input_passes='RGB_ALBEDO_NORMAL'
scene.cycles.denoising_prefilter='ACCURATE'
scene.render.threads_mode='AUTO'
scene.render.resolution_x=args.size
scene.render.resolution_y=args.size
scene.render.resolution_percentage=100
scene.render.image_settings.file_format='PNG'
scene.render.image_settings.color_depth='16'
scene.render.image_settings.color_mode='RGB'
out=ROOT/'gallery/renders'/f'{args.material}-rotated-light-{args.size}'
out.parent.mkdir(parents=True,exist_ok=True)
assert not out.with_suffix('.png').exists(), 'Refusing to overwrite a prior supplementary image'
scene.render.filepath=str(out.with_suffix('.png'))
start=time.monotonic()
bpy.ops.render.render(write_still=True)
p=out.with_suffix('.png')
receipt=dict(schema='cybr.material.gallery.supplementary.v1',material_id=args.material,
    selection='Supplementary view; does not replace or promote an accepted baseline',
    source_scene=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    blender=bpy.app.version_string,engine='CYCLES',device='CPU',resolution=[args.size,args.size],
    maximum_samples=args.samples,minimum_samples=args.minimum,adaptive_threshold=.01,
    denoising='OPENIMAGEDENOISE / RGB_ALBEDO_NORMAL / ACCURATE',
    postprocessing='Native Cycles output and AgX; no resizing or sharpening',
    view_transform=scene.view_settings.view_transform,look=scene.view_settings.look,
    exposure=scene.view_settings.exposure,light_rotation_degrees=90,lights=lights,
    camera_matrix=[list(row) for row in scene.camera.matrix_world],
    lens_mm=scene.camera.data.lens,dependencies=dependencies,elapsed_seconds=time.monotonic()-start,
    file=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size)
out.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('SUPPLEMENTARY_RENDER_COMPLETE',p)
