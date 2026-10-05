"""Validate render-ready Cycles scenes and their portable physical materials."""
import json
import sys
from pathlib import Path
import bpy
from mathutils import Vector, Matrix
from mathutils.bvhtree import BVHTree
import random

ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'materials/manifest.json').read_text())
specs={s['id']:s for s in manifest['materials']}
geometry_only='--geometry-only' in sys.argv
scenes=[s for s in bpy.data.scenes if s.objects]
is_detail='Details' in Path(bpy.data.filepath).name
is_architecture='Architecture' in Path(bpy.data.filepath).name
assert len(scenes)==(4 if is_architecture else 10 if is_detail else 3),[s.name for s in scenes]
seen=set()
displaced=0
for scene in scenes:
    assert scene.render.engine=='CYCLES' and scene.cycles.device=='CPU'
    guided=is_architecture and 'Joints and Surface Wear' not in scene.name
    assert scene.cycles.samples==(512 if guided else 768)
    if not geometry_only:assert scene.cycles.use_denoising is guided
    if guided and not geometry_only:
        assert scene.cycles.denoiser=='OPENIMAGEDENOISE'
        assert scene.cycles.denoising_input_passes=='RGB_ALBEDO_NORMAL'
        assert scene.cycles.denoising_prefilter=='ACCURATE'
    assert scene.cycles.use_adaptive_sampling and scene.cycles.adaptive_min_samples==128
    assert scene.render.image_settings.file_format=='PNG'
    assert scene.render.image_settings.color_depth=='16'
    assert scene.render.filepath.startswith('//../path_traced/renders/')
    assert scene.unit_settings.system=='METRIC'
    assert scene.camera and scene.camera.data.clip_start>0
    assert scene.view_settings.view_transform=='AgX'
    assert not scene.use_nodes
    assert len([o for o in scene.objects if o.type=='LIGHT'])==3
    if is_architecture:
        assert scene['source_image_inputs']==0
        assert scene['architecture_object_count']>50
        if 'floor_material' in scene and 'Joints' not in scene.name:
            assert scene['floor_board_count']>300
        if 'room_dimensions_m' in scene:
            assert scene['room_dimensions_m'][0]>6 and scene['room_dimensions_m'][1]>5
        if 'Oak Reading Room' in scene.name:
            def z_bounds(ob):
                transform=Matrix.LocRotScale(ob.location,ob.rotation_euler.to_quaternion(),ob.scale)
                values=[(transform@v.co).z for v in ob.data.vertices]
                return min(values),max(values)
            shelf_tops=[z_bounds(ob)[1] for ob in scene.objects if ob.name.startswith('Library / steel shelf')]
            spines=[ob for ob in scene.objects if ob.name.startswith('Library / bound spine')]
            assert len(spines)==240 and len(shelf_tops)==5
            for ob in spines:
                bottom=z_bounds(ob)[0]
                gap=min((bottom-top for top in shelf_tops if top<=bottom+.0001),default=1)
                assert -.0001<=gap<=.001,(ob.name,'Book/shelf contact gap',gap)
        if 'floor_material' in scene:
            # Check the actual mesh layout, independently of its tiling recipe.
            # This catches uncovered patches and overlapping basket modules.
            verts=[];faces=[]
            for ob in scene.objects:
                if not ob.name.startswith('Parquet / individually'):continue
                # Inactive scenes can retain unevaluated matrix_world values
                # in a build-only file. Boards are unparented and unconstrained;
                # compose their saved transforms without a render dependency.
                assert ob.parent is None and not ob.constraints
                transform=Matrix.LocRotScale(ob.location,ob.rotation_euler.to_quaternion(),ob.scale)
                start=len(verts);verts.extend(transform@v.co for v in ob.data.vertices)
                faces.extend(tuple(start+i for i in p.vertices) for p in ob.data.polygons)
            tree=BVHTree.FromPolygons(verts,faces)
            dims=scene.get('room_dimensions_m',(1.85,1.65));rng=random.Random(140551)
            hits=0;count=12000
            for i in range(count):
                x=rng.uniform(-dims[0]/2+.03,dims[0]/2-.03);y=rng.uniform(-dims[1]/2+.03,dims[1]/2-.03)
                if tree.ray_cast(Vector((x,y,.15)),Vector((0,0,-1)),.30)[0] is not None:hits+=1
            assert hits/count>.980,(scene.name,'Parquet coverage',hits/count)
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
assert displaced>=(1 if is_architecture else 10 if is_detail else 5),displaced
result=dict(file=Path(bpy.data.filepath).name,scenes=len(scenes),materials=len(seen),
            displaced_objects=displaced,engine='CYCLES',portable_texture_links=True,
            guided_interior_denoising=any(s.cycles.use_denoising for s in scenes),
            verification='geometry_only' if geometry_only else 'complete',passed=True)
folder=ROOT/'path_traced/verification';folder.mkdir(parents=True,exist_ok=True)
(folder/(Path(bpy.data.filepath).stem+('_geometry' if geometry_only else '')+'.json')).write_text(json.dumps(result,indent=2)+'\n')
print('PATH_TRACE_SCENES_VERIFIED '+json.dumps(result),flush=True)
