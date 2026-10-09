"""Validate native clear-PETG scene shader and physical shell without rendering."""
import json,sys
from pathlib import Path
import bpy
root=Path(__file__).resolve().parents[2]
path=root/'expansion/scenes/23_petg_transparent_hero_soft.blend'
bpy.ops.wm.open_mainfile(filepath=str(path))
wall=next(o for o in bpy.context.scene.objects if o.name.startswith('PETG / fused'))
mat=wall.data.materials[0];nodes=mat.node_tree.nodes;p=nodes.get('Principled BSDF')
assert p.inputs['Transmission Weight'].default_value>.9
assert p.inputs['Alpha'].default_value==1
assert p.inputs['IOR'].default_value>1.5
assert abs(wall['wall_thickness_m']-.0012)<1e-7
assert any(n.type=='VOLUME_ABSORPTION' for n in nodes)
assert any('Roughness'==n.label for n in nodes if n.type=='TEX_IMAGE')
report=dict(scene=str(path),actual_transmission=float(p.inputs['Transmission Weight'].default_value),
            alpha=float(p.inputs['Alpha'].default_value),ior=float(p.inputs['IOR'].default_value),
            wall_thickness_m=wall['wall_thickness_m'],layer_pitch_m=wall['layer_height_m'],
            geometry='closed hollow solid wall; printed base separate fused overlap',
            layer_boundary='roughness map follows UV layer pitch; actual shoulders in geometry',
            volume_absorption_density=next(n.inputs['Density'].default_value for n in nodes if n.type=='VOLUME_ABSORPTION'),
            limitation='Node/geometry checks are not visual transparency evidence; real Cycles review remains required')
(root/'expansion/transparent-optics-check.json').write_text(json.dumps(report,indent=2)+'\n')
print('CLEAR_PETG_OPTICS_CHECK',report)
