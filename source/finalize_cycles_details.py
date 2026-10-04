"""Normalize portable outputs in an existing Cycles details file."""
from pathlib import Path
import bpy

for scene in list(bpy.data.scenes):
    if not scene.objects:
        bpy.data.scenes.remove(scene)
    else:
        assert scene.render.engine == 'CYCLES'
        scene.render.filepath = '//../path_traced/renders/' + Path(scene.render.filepath).name

for image in bpy.data.images:
    if image.source == 'FILE':
        image.filepath = '//../materials/' + Path(image.filepath).parent.name + '/' + Path(image.filepath).name

bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath, compress=True)
print('CYCLES_DETAILS_PORTABLE', flush=True)
