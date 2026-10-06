"""Explicit backlit control: illuminate through the unchanged rough printed cup."""
import bpy,json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];src=ROOT/'scenes/23_petg_transparent_r6_sideproof.blend';bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
sys.path.insert(0,str(Path(__file__).resolve().parent))
from expansion.geometry import box,plain
# Put all three primary contrast fields behind the central wall, rather than only at its edges.
for o in list(s.objects):
    if o.name.startswith('Optics / side-on contrast chart'):bpy.data.objects.remove(o,do_unlink=True)
for i,c in enumerate([(.78,.75,.68),(.009,.012,.016),(.78,.75,.68),(.42,.025,.012),(.78,.75,.68)]):
    box('Optics / side-on contrast chart '+str(i),((i-2)*.009,.026,.021),(.009,.001,.042),plain('Optics / narrow backlit field '+str(i),c,.86))
for o in s.objects:
    if o.name.startswith('Optics / side-on contrast chart'):
        m=o.data.materials[0];p=m.node_tree.nodes['Principled BSDF'];p.inputs['Emission Color'].default_value=p.inputs['Base Color'].default_value;p.inputs['Emission Strength'].default_value=3
        o['role']='Emissive physical transmission-control panel, not ordinary unlit printed card'
    if o.type=='LIGHT':o.data.energy*=.25
    if o.name.startswith('Optics / chart horizontal bar'):o.dimensions.x=.045
s.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.08
s.cycles.samples=128;s.cycles.adaptive_min_samples=32;s.cycles.adaptive_threshold=.02
s['revision']='r7: explicit backlit control only, unchanged cup mesh and optical material'
s['optical_expectation']='Diagnostic evidence of transmission and blur through a rough FDM print; not a clear optical-window or finished hero claim.'
out=ROOT/'scenes/23_petg_transparent_r7_backlit_control.blend';s.render.filepath=str(ROOT/'evidence/23_petg_transparent_r7_backlit_control.png');bpy.ops.wm.save_as_mainfile(filepath=str(out));print('READY_SCENE',out,flush=True)
(ROOT/'tests/clear_petg_r7_patch.json').write_text(json.dumps({'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'cup_geometry_changed':False,'cup_material_changed':False,'chart_emission_strength':3,'front_area_light_multiplier':.25,'world_strength':.08,'claim_limit':'Explicit backlit transmission diagnostic, not ordinary front-lit transparency acceptance'},indent=2)+'\n')
