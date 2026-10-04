"""Validate render-ready Cycles scenes and their portable physical materials."""
import json
from pathlib import Path
import bpy

ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'materials/manifest.json').read_text())
specs={s['id']:s for s in manifest['materials']}
scenes=[s for s in bpy.data.scenes if s.objects]
is_detail='Details' in Path(bpy.data.filepath).name
assert len(scenes)==(10 if is_detail else 3),[s.name for s in scenes]
seen=set()
displaced=0
for scene in scenes:
    assert scene.render.engine=='CYCLES' and scene.cycles.device=='CPU'
    assert scene.cycles.samples==768
    assert scene.cycles.use_denoising is False
    assert scene.cycles.use_adaptive_sampling and scene.cycles.adaptive_min_samples==128
    assert scene.render.image_settings.file_format=='PNG'
    assert scene.render.image_settings.color_depth=='16'
    assert scene.render.filepath.startswith('//../path_traced/renders/')
    assert scene.unit_settings.system=='METRIC'
    assert scene.camera and scene.camera.data.clip_start>0
    assert scene.view_settings.view_transform=='AgX'
    assert not scene.use_nodes
    assert len([o for o in scene.objects if o.type=='LIGHT'])==3
    for obj in scene.objects:
        if obj.type!='MESH':continue
        for mat in obj.data.materials:
            if 'aurel_id' not in mat:continue
            id=mat['aurel_id'];seen.add(id)
            spec=specs[id]
            group=next(n for n in mat.node_tree.nodes if n.type=='GROUP')
            assert group.inputs['Normal Strength'].default_value==1
            assert all(abs(v-1)<1e-6 for v in group.inputs['UV Scale'].default_value)
            displace=[m for m in obj.modifiers if m.type=='DISPLACE']
            if displace:
                displaced+=1
                assert group.inputs['Displacement Mode'].default_value
                for mod in displace:
                    assert abs(mod.strength-spec['height_scale_m'])<1e-7
                    assert mod.texture_coords=='UV'
                    assert mod.texture.image.filepath.endswith('Height_Macro.png')
            elif group.inputs['Displacement Mode'].default_value:
                raise AssertionError('Residual normal without geometric height: '+obj.name)
            for node in group.node_tree.nodes:
                if node.type=='TEX_IMAGE':
                    assert node.image.filepath.startswith('//../materials/')
                    assert Path(bpy.path.abspath(node.image.filepath)).is_file()
                    assert node.image.colorspace_settings.name==('sRGB' if node.image.filepath.endswith('BaseColor.png') else 'Non-Color')
assert seen==set(specs),seen
assert displaced>=(10 if is_detail else 5),displaced
result=dict(file=Path(bpy.data.filepath).name,scenes=len(scenes),materials=len(seen),
            displaced_objects=displaced,engine='CYCLES',portable_texture_links=True,
            denoising=False,passed=True)
folder=ROOT/'path_traced/verification';folder.mkdir(parents=True,exist_ok=True)
(folder/(Path(bpy.data.filepath).stem+'.json')).write_text(json.dumps(result,indent=2)+'\n')
print('PATH_TRACE_SCENES_VERIFIED '+json.dumps(result),flush=True)
