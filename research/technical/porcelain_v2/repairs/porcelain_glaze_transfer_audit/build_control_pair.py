"""Build admitted full-material control/candidate pairs. No rendering.

Prerequisite in both: previously verified continuous visible UV atlas.
Candidate-only intervention: explicit shell macro normals and normal-region split.
"""
import bpy, hashlib, json, math, sys
from pathlib import Path
from mathutils import Matrix, Vector

argv = sys.argv[sys.argv.index('--')+1:]
assert len(argv) == 3 and argv[2] == '--admitted', 'Explicit resource admission flag required'
SOURCE, OUT = Path(argv[0]).resolve(), Path(argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
EXPECTED_SOURCE = '5837b86217a72b5d380321a7c2583feea186a181921801ab3538fdf2ba5bd0c1'
NAME = '07 new / finite60micron mean glaze shell'
IX, ZC, PERIOD = .0474, .0012, .0192

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while block := f.read(1024*1024):
            h.update(block)
    return h.hexdigest()

def scalar_rna(value, exclude=()):
    result = {}
    for prop in value.bl_rna.properties:
        key = prop.identifier
        if key in exclude or key in {'rna_type', 'name', 'name_full', 'is_updated', 'is_updated_data',
                                    'is_evaluated', 'is_runtime_data', 'session_uid', 'users', 'tag'}:
            continue
        if prop.type not in {'BOOLEAN','INT','FLOAT','STRING','ENUM'}:
            continue
        try:
            v = getattr(value, key)
            if hasattr(v, 'to_list'):
                v = v.to_list()
            elif not isinstance(v, (str, bool, int, float)):
                v = list(v)
            json.dumps(v, allow_nan=False)
            result[key] = v
        except (TypeError, ValueError, AttributeError):
            pass
    return result

def tree(nt):
    return {'nodes': {n.name: {
        'type': n.bl_idname, 'properties': scalar_rna(n),
        'inputs': [{**scalar_rna(i), **({'value': list(i.default_value) if hasattr(i.default_value, '__len__') and not isinstance(i.default_value, str) else i.default_value} if hasattr(i, 'default_value') else {})} for i in n.inputs],
        'image': ({'sha256': sha(bpy.path.abspath(n.image.filepath)),
                   'colorspace': n.image.colorspace_settings.name,
                   'alpha_mode': n.image.alpha_mode, 'source': n.image.source}
                  if n.type == 'TEX_IMAGE' and n.image else None)} for n in nt.nodes},
        'links': sorted([l.from_node.name, l.from_socket.identifier, l.to_node.name, l.to_socket.identifier] for l in nt.links)}

def snapshot():
    s = bpy.context.scene
    for layer in s.view_layers:
        layer.update()
    result = {'objects': {}, 'materials': {m.name: {'properties': scalar_rna(m), 'tree': tree(m.node_tree) if m.use_nodes else None} for m in bpy.data.materials},
              'world': {'properties': scalar_rna(s.world), 'tree': tree(s.world.node_tree)},
              'view_settings': scalar_rna(s.view_settings), 'display_settings': scalar_rna(s.display_settings),
              'unit_settings': scalar_rna(s.unit_settings), 'render': scalar_rna(s.render),
              'cycles': scalar_rna(s.cycles), 'active_camera': s.camera.name}
    for ob in s.objects:
        d = {'type': ob.type, 'matrix': [list(v) for v in ob.matrix_world],
             'visibility': [ob.hide_render, ob.hide_viewport], 'modifiers': [m.type for m in ob.modifiers],
             'constraints': [c.type for c in ob.constraints]}
        if ob.type == 'MESH':
            me = ob.data
            d.update(vertices=[list(v.co) for v in me.vertices], edges=[list(e.vertices) for e in me.edges],
                     polygons=[{'vertices': list(p.vertices), 'material_index': p.material_index} for p in me.polygons],
                     loops=[l.vertex_index for l in me.loops], materials=[m.name for m in me.materials],
                     uv_layers={uv.name: {'active_render': uv.active_render, 'uv': [list(x.uv) for x in uv.data]} for uv in me.uv_layers})
        elif ob.type in {'LIGHT','CAMERA'}:
            d['data'] = scalar_rna(ob.data)
            if ob.type == 'CAMERA':
                d['dof'] = scalar_rna(ob.data.dof)
        result['objects'][ob.name] = d
    return result

def normal_state():
    return {o.name: {'smooth': [p.use_smooth for p in o.data.polygons],
                     'sharp': [e.use_edge_sharp for e in o.data.edges],
                     'normals': [list(v.vector) for v in o.data.corner_normals],
                     'has_custom_normals': o.data.has_custom_normals}
            for o in bpy.context.scene.objects if o.type == 'MESH'}

def atlas(me):
    uv = me.uv_layers.active
    flat = {k: list(uv.data[k].uv) for p in me.polygons if p.normal.z > .999 for k in p.loop_indices}
    cache, changed = {}, 0
    for p in me.polygons:
        if p.normal.z <= .02:
            continue
        for k in p.loop_indices:
            vid = me.loops[k].vertex_index
            if vid not in cache:
                x,y,z = me.vertices[vid].co
                cx,cy = min(max(x,-IX),IX), min(max(y,-IX),IX)
                dx,dy = x-cx,y-cy
                rho = math.hypot(dx,dy)
                if rho < 1e-9:
                    u,v = x,y
                else:
                    rise=z-ZC
                    scale=math.hypot(rho,rise)*math.atan2(rho,rise)/rho
                    u,v=cx+dx*scale,cy+dy*scale
                cache[vid] = (u/PERIOD+.5,v/PERIOD+.5)
            before = tuple(uv.data[k].uv)
            uv.data[k].uv = cache[vid]
            changed += tuple(uv.data[k].uv) != before
    assert all(list(uv.data[k].uv) == v for k,v in flat.items())
    # Exact per-vertex assignment guarantees continuous outer chart.
    for p in me.polygons:
        if p.normal.z > .02:
            for k in p.loop_indices:
                assert tuple(uv.data[k].uv) == tuple(Vector(cache[me.loops[k].vertex_index]))
    return {'changed_loop_count': changed, 'flat_top_uv_bit_exact': True,
            'continuous_outer_chart': True, 'method': 'Existing radial continuous-atlas prerequisite; corners remain non-isometric'}

def nominal(inner):
    v = Vector((inner.x-min(max(inner.x,-IX),IX), inner.y-min(max(inner.y,-IX),IX), inner.z-ZC))
    for i in range(3):
        if abs(v[i]) < 1e-8:
            v[i] = 0
    v.normalize()
    return v

def correct_normals(me):
    half = len(me.vertices)//2
    assert len(me.vertices) == 488
    target, regions, edge_regions = [None]*len(me.loops), {}, {}
    for p in me.polygons:
        vertices = [me.loops[k].vertex_index for k in p.loop_indices]
        if p.normal.z > .02:
            assert all(v < half for v in vertices)
            kind = 'outer'
        elif all(v >= half for v in vertices):
            kind = 'inner'
        else:
            kind = 'termination'
        regions[p.index] = kind
        # Keep planar faces flat so their actual evaluated normals remain exact;
        # custom-normal packing otherwise introduces a tiny nonzero XY tilt.
        p.use_smooth = kind != 'termination' and abs(p.normal.z) < .999
        for k in p.loop_indices:
            vid = me.loops[k].vertex_index
            if kind == 'outer':
                assert (me.vertices[vid+half].co-me.vertices[vid].co).length < .00007
                normal = nominal(me.vertices[vid+half].co)
            elif kind == 'inner':
                normal = -nominal(me.vertices[vid].co)
            else:
                normal = p.normal.copy()
            target[k] = tuple(normal)
            edge_regions.setdefault(me.loops[k].edge_index,set()).add(kind)
    for e in me.edges:
        e.use_edge_sharp = 'termination' in edge_regions[e.index] and len(edge_regions[e.index]) > 1
    me.normals_split_custom_set(target)
    me.update()
    return target, regions

def verify_normals(me, target, regions):
    actual = [v.vector.copy() for v in me.corner_normals]
    def precise_angle(a,b):
        cross=(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
        return math.degrees(math.atan2(math.sqrt(sum(x*x for x in cross)),sum(x*y for x,y in zip(a,b))))
    errors = [precise_angle(a,b) for a,b in zip(actual,target)]
    # Blender's packed custom normal representation is quantized.
    assert max(errors) < .02, max(errors)
    for p in me.polygons:
        if p.normal.z > .999:
            assert all((actual[k]-Vector((0,0,1))).length < 1e-6 for k in p.loop_indices)
        if regions[p.index] == 'outer':
            assert all(actual[k].z >= -2e-4 for k in p.loop_indices), 'Outer normal exceeds packed-normal quantization tolerance'
    assert me.has_custom_normals
    outer_z = [actual[k].z for p in me.polygons if regions[p.index] == 'outer' for k in p.loop_indices]
    return {'max_saved_custom_normal_error_degrees': max(errors), 'exact_flat_top_normals': True,
            'minimum_outer_normal_z': min(outer_z), 'outer_z_quantization_tolerance': 2e-4,
            'split_termination_edges': sum(e.use_edge_sharp for e in me.edges)}

def rotate_lights():
    turn = Matrix(((0,-1,0,0),(1,0,0,0),(0,0,1,0),(0,0,0,1)))
    names = []
    for ob in bpy.context.scene.objects:
        if ob.type == 'LIGHT' and ob.data.type == 'AREA' and ob.name.startswith('Slab /'):
            ob.matrix_world = turn @ ob.matrix_world
            names.append(ob.name)
    bpy.context.view_layer.update()
    assert names
    return sorted(names)

assert sha(SOURCE) == EXPECTED_SOURCE
receipts, pair_snapshots, source_state = [], {}, None
for rig in ['main','second_light']:
    for treatment in ['control','candidate']:
        bpy.ops.wm.open_mainfile(filepath=str(SOURCE),load_ui=False)
        bpy.context.preferences.filepaths.save_version=0
        initial = snapshot()
        initial_normals = normal_state()
        if source_state is None:
            source_state = initial
        else:
            assert initial == source_state
        s = bpy.context.scene
        assert s['deposition_layout_control'] == 'control'
        me = s.objects[NAME].data
        prerequisite = atlas(me)
        after_atlas = snapshot()
        # Only glaze UVs may change in the shared prerequisite.
        audit = json.loads(json.dumps(after_atlas))
        audit['objects'][NAME]['uv_layers'] = initial['objects'][NAME]['uv_layers']
        assert audit == initial
        target = regions = None
        if treatment == 'candidate':
            target,regions = correct_normals(me)
        changed_lights = rotate_lights() if rig == 'second_light' else []
        before_save = snapshot()
        expected_normals = normal_state()
        if rig == 'main':
            assert before_save == after_atlas
        else:
            check = json.loads(json.dumps(before_save))
            for name in changed_lights:
                check['objects'][name]['matrix'] = after_atlas['objects'][name]['matrix']
            assert check == after_atlas
        if treatment == 'control':
            assert expected_normals == initial_normals
        else:
            assert {k:v for k,v in expected_normals.items() if k != NAME} == {k:v for k,v in initial_normals.items() if k != NAME}
        path = OUT/f'07_porcelain_control_atlas_{treatment}_{rig}.blend'
        assert not path.exists(), f'Refusing overwrite: {path}'
        s['macro_normal_treatment'] = treatment
        s['light_rig_variant'] = rig
        bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True)
        # A fresh reopen verifies persisted data, including evaluated loop normals.
        bpy.ops.wm.open_mainfile(filepath=str(path),load_ui=False)
        reloaded = snapshot()
        assert reloaded == before_save, 'Fresh reopen changed a frozen datum'
        assert normal_state() == expected_normals, 'Saved normal/smooth/sharp data differ'
        verification = verify_normals(bpy.context.scene.objects[NAME].data,target,regions) if target is not None else {'original_normals_smooth_sharp_preserved': True}
        if rig in pair_snapshots:
            assert reloaded == pair_snapshots[rig], 'Control/candidate differ beyond declared normal attributes'
        else:
            pair_snapshots[rig] = reloaded
        receipts.append({'treatment':treatment,'rig':rig,'scene':str(path),'sha256':sha(path),
                         'bytes':path.stat().st_size,'prerequisite':prerequisite,'normal_verification':verification,
                         'cameras':[o.name for o in bpy.context.scene.objects if o.type=='CAMERA'],
                         'second_light_names':changed_lights,'fresh_reopen_equal':True})
        print('FROZEN_SCENE', json.dumps(receipts[-1]),flush=True)
assert sha(SOURCE) == EXPECTED_SOURCE
assert sum(p.stat().st_size for p in OUT.glob('*.blend')) < 8*1024*1024
report = {'source':str(SOURCE),'source_sha256':EXPECTED_SOURCE,
          'source_unchanged':True,'arrival_case':'control','mean_film_m':6.362943104938368e-5,
          'candidate_only_changes':['shell explicit custom loop normals','shell smooth flags','termination sharp edges'],
          'matched_pair_geometry_uv_materials_lights_cameras_world_render_settings_equal':True,
          'second_light_protocol':'Every existing Slab area light rotated +90 degrees world Z about x=y=0; height and light data unchanged',
          'historical_reference_is_not_claimed_as_same_binary':True,
          'known_limits':['Opaque pigmented surface with coat lobe; no glaze/body layer transport','Native retained field slopes remain','Corner radial chart is not isometric'],
          'rendered':False,'visual_acceptance':False,'scenes':receipts}
(OUT/'build_receipt.json').write_text(json.dumps(report,indent=2))
print('MATCHED_PAIR_FROZEN',str(OUT/'build_receipt.json'),flush=True)
