import bpy
from pathlib import Path
root=Path(__file__).resolve().parents[1]
for id in ('12_lattice_wire','14_chainmail','22_petg','23_petg_transparent'):
    path=root/'scenes'/f'{id}_r4_hero.blend';bpy.ops.wm.open_mainfile(filepath=str(path))
    for image in bpy.data.images:
        if image.source=='FILE':image.pack()
    bpy.ops.wm.save_as_mainfile(filepath=str(path))
    print('PACKED',id,flush=True)
