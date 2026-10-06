"""Admitted source-only effective-layer experiment. Never calls a renderer.

The coating remains an opaque effective closure. This is an optical-model
change, not coefficient-identical repair or finite-volume glass transport.
"""
import ast
import copy
import hashlib
import json
import math
import sys
from pathlib import Path

import bpy

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = ROOT.parents[1]
AUDIT = PRODUCTION / 'experiments/porcelain_glaze_transfer_audit'
OUT = ROOT / 'scenes'
SHELL = '07 new / finite60micron mean glaze shell'
BODY = '07 new / unglazed porcelain body'
GLAZE = '07 new / mean-calibrated pigmented clearcoat surface'
ADDED = ['Layer substrate height', 'Layer substrate bump']
SOURCES = {
    'main': ('07_porcelain_control_atlas_candidate_main.blend', '36e6117e953550b50636ef8fcbb0ab5494edaafc69db8ab0c656bdc41601433a'),
    'second_light': ('07_porcelain_control_atlas_candidate_second_light.blend', 'ca652a0e4db67a9e34b0818f1ac559188cbe4b1d53269f02a588cd4606bd771e'),
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


# Reuse only audited observation helpers, never the old builder's executable body.
AUDIT_HELPERS_SHA = '0ea5f1500b25e28ab7c5b3b3c505054711273f4c4df40ce822e709a68fd44087'
helper_path = AUDIT / 'build_control_pair.py'
assert sha(helper_path) == AUDIT_HELPERS_SHA
parsed = ast.parse(helper_path.read_text())
definitions = [n for n in parsed.body if isinstance(n, ast.FunctionDef) and
               n.name in {'scalar_rna', 'tree', 'snapshot', 'normal_state'}]
assert len(definitions) == 4
exec(compile(ast.Module(body=definitions, type_ignores=[]), str(helper_path), 'exec'))


def id_properties(block):
    def plain(v):
        if hasattr(v, 'to_dict'):
            return {k: plain(x) for k, x in v.to_dict().items()}
        if hasattr(v, 'to_list'):
            return [plain(x) for x in v.to_list()]
        return v
    return {k: plain(v) for k, v in block.items()}


def frozen_state():
    s = bpy.context.scene
    result = snapshot()
    result['normal_state'] = normal_state()
    result['scene_properties'] = id_properties(s)
    result['image_settings'] = scalar_rna(s.render.image_settings)
    result['frame'] = s.frame_current
    result['scene_scalar'] = scalar_rna(s)
    result['compositor'] = tree(s.node_tree) if s.use_nodes else None
    result['view_layers'] = {v.name: {'properties': scalar_rna(v), 'cycles': scalar_rna(v.cycles)} for v in s.view_layers}
    result['material_properties'] = {m.name: id_properties(m) for m in bpy.data.materials}
    result['object_properties'] = {o.name: id_properties(o) for o in s.objects}
    result['world_properties'] = id_properties(s.world)
    result['images'] = {i.name: {'file': str(Path(bpy.path.abspath(i.filepath)).resolve()),
                               'sha256': sha(bpy.path.abspath(i.filepath)),
                               'colorspace': i.colorspace_settings.name,
                               'source': i.source, 'alpha_mode': i.alpha_mode,
                               'packed': i.packed_file is not None}
                        for i in bpy.data.images if i.source == 'FILE'}
    for ob in s.objects:
        entry = result['objects'][ob.name]
        entry['properties'] = scalar_rna(ob)
        entry['modifiers'] = [scalar_rna(x) for x in ob.modifiers]
        entry['constraints'] = [scalar_rna(x) for x in ob.constraints]
        if ob.type == 'MESH':
            entry['mesh_properties'] = scalar_rna(ob.data)
            entry['mesh_custom_properties'] = id_properties(ob.data)
            entry['active_uv'] = ob.data.uv_layers.active_index
    return result


def close(a, b):
    assert math.isclose(a, b, rel_tol=0, abs_tol=1e-7), (a, b)


def incoming(socket):
    return [(l.from_node.name, l.from_socket.identifier) for l in socket.links]


def copy_node(source, tree, name):
    """Copy a single supported body node without touching its source material."""
    target = tree.nodes.new(source.bl_idname)
    for prop in source.bl_rna.properties:
        if prop.is_readonly or prop.identifier in {'rna_type', 'name', 'parent'}:
            continue
        if prop.type in {'BOOLEAN', 'INT', 'FLOAT', 'STRING', 'ENUM'}:
            try:
                setattr(target, prop.identifier, getattr(source, prop.identifier))
            except (AttributeError, TypeError):
                pass
    target.name = name
    if source.type == 'TEX_IMAGE':
        target.image = source.image
        target.image_user.frame_duration = source.image_user.frame_duration
        target.image_user.frame_start = source.image_user.frame_start
        target.image_user.frame_offset = source.image_user.frame_offset
        target.image_user.use_auto_refresh = source.image_user.use_auto_refresh
        target.image_user.use_cyclic = source.image_user.use_cyclic
    for a, b in zip(source.inputs, target.inputs):
        assert a.identifier == b.identifier
        if hasattr(a, 'default_value'):
            b.default_value = a.default_value
    # Node presentation is copied as well, and excluded from no equality tests.
    return target


def assert_node_copy(a, b):
    old, new = tree(a.id_data)['nodes'][a.name], tree(b.id_data)['nodes'][b.name]
    # Names are keyed externally; connectedness is checked by the explicit graph below.
    for entry in (old, new):
        for socket in entry['inputs']:
            socket.pop('is_linked', None)
    assert old == new, ('Body node copy differs', a.name, b.name)


def assign_layers():
    glaze, body = bpy.data.materials[GLAZE], bpy.data.materials[BODY]
    gt, bt = glaze.node_tree, body.node_tree
    p, b = gt.nodes['Principled BSDF'], bt.nodes['Principled BSDF']
    assert list(bpy.context.scene.objects[SHELL].data.materials) == [glaze]
    assert incoming(b.inputs['Normal']) == [('Bump', 'Normal')]
    source_bump = bt.nodes['Bump']
    assert incoming(source_bump.inputs['Height']) == [('Image Texture', 'Color')]
    source_height = bt.nodes['Image Texture']
    assert Path(bpy.path.abspath(source_height.image.filepath)).name == 'Substrate_Height.png'
    close(source_bump.inputs['Distance'].default_value, 20e-6)
    close(source_bump.inputs['Strength'].default_value, 1)
    assert not source_bump.inputs['Normal'].is_linked
    assert not b.inputs['Base Color'].is_linked and not b.inputs['Roughness'].is_linked
    close(b.inputs['Roughness'].default_value, .67)
    close(b.inputs['IOR'].default_value, 1.49)
    close(p.inputs['IOR'].default_value, 1.49)
    close(p.inputs['Coat Weight'].default_value, .32)
    close(p.inputs['Coat Roughness'].default_value, .12)
    for key in ['Transmission Weight', 'Metallic', 'Subsurface Weight', 'Sheen Weight', 'Thin Film Thickness']:
        close(p.inputs[key].default_value, 0)
        assert not p.inputs[key].is_linked
    assert incoming(p.inputs['Normal']) == [('Normal Map', 'Normal')]
    assert incoming(p.inputs['Coat Normal']) == [('Normal Map', 'Normal')]
    assert incoming(p.inputs['Roughness']) == [('Image Texture.001', 'Color')]
    assert not p.inputs['Coat Roughness'].is_linked
    assert Path(bpy.path.abspath(gt.nodes['Image Texture.001'].image.filepath)).name == 'Glaze_Roughness.png'
    assert Path(bpy.path.abspath(gt.nodes['Image Texture'].image.filepath)).name == 'Glaze_Normal_OpenGL_RGB16.png'
    close(gt.nodes['Normal Map'].inputs['Strength'].default_value, 1)
    assert gt.nodes['Normal Map'].space == 'TANGENT' and gt.nodes['Normal Map'].uv_map == 'UVMap'
    for node in [source_height, gt.nodes['Image Texture'], gt.nodes['Image Texture.001']]:
        assert tuple(node.image.size) == (4096, 4096)
        assert node.image.colorspace_settings.name == 'Non-Color'
    assert tuple(p.inputs['Coat Tint'].default_value) == (1, 1, 1, 1)
    old_glaze_index = p.inputs['IOR'].default_value
    body_index = b.inputs['IOR'].default_value
    old_coat_index = p.inputs['Coat IOR'].default_value
    old_color = list(p.inputs['Base Color'].default_value)
    old_rough = p.inputs['Roughness'].default_value
    height = copy_node(source_height, gt, ADDED[0])
    bump = copy_node(source_bump, gt, ADDED[1])
    p.inputs['Base Color'].default_value = b.inputs['Base Color'].default_value
    p.inputs['Roughness'].default_value = b.inputs['Roughness'].default_value
    p.inputs['IOR'].default_value = body_index / old_glaze_index
    p.inputs['Coat Weight'].default_value = 1
    p.inputs['Coat IOR'].default_value = old_glaze_index
    for name in ['Normal', 'Roughness']:
        for link in list(p.inputs[name].links):
            gt.links.remove(link)
    gt.links.new(height.outputs['Color'], bump.inputs['Height'])
    gt.links.new(bump.outputs['Normal'], p.inputs['Normal'])
    gt.links.new(gt.nodes['Image Texture.001'].outputs['Color'], p.inputs['Coat Roughness'])
    assert_node_copy(source_height, height)
    assert_node_copy(source_bump, bump)
    assert incoming(p.inputs['Normal']) == [(ADDED[1], 'Normal')]
    assert incoming(p.inputs['Coat Normal']) == [('Normal Map', 'Normal')]
    assert incoming(p.inputs['Coat Roughness']) == [('Image Texture.001', 'Color')]
    assert not p.inputs['Roughness'].is_linked
    assert p.inputs['IOR'].default_value == 1.0
    return {'stored_glaze_ior': old_glaze_index, 'stored_body_ior': body_index,
            'old_coat_ior': old_coat_index, 'buried_relative_ior': p.inputs['IOR'].default_value,
            'old_base_color': old_color, 'new_base_color_from_body': list(p.inputs['Base Color'].default_value),
            'old_unlinked_base_roughness_default': old_rough,
            'new_base_roughness_from_body': p.inputs['Roughness'].default_value,
            'coat_roughness': 'unchanged native Glaze_Roughness; old 0.12 default is now inactive',
            'coat_weight': 1, 'transmission_weight': 0,
            'substrate_normal': 'exact copied body Substrate_Height -> Bump; shell UV/macro normal unchanged'}


def verify_only_declared_changes(before, after):
    normalized = copy.deepcopy(after)
    old = before['materials'][GLAZE]['tree']
    new = normalized['materials'][GLAZE]['tree']
    assert set(new['nodes']) - set(old['nodes']) == set(ADDED)
    for name in ADDED:
        del new['nodes'][name]
    allowed_inputs = {'Base Color', 'Roughness', 'IOR', 'Coat Weight', 'Coat IOR', 'Normal', 'Coat Roughness'}
    old_inputs = old['nodes']['Principled BSDF']['inputs']
    new_inputs = new['nodes']['Principled BSDF']['inputs']
    changes = []
    for i, (a, b) in enumerate(zip(old_inputs, new_inputs)):
        if a != b:
            assert a['identifier'] in allowed_inputs, a['identifier']
            changes.append({'input': a['identifier'], 'before': a, 'after': b})
            new_inputs[i] = copy.deepcopy(a)
    removed = {('Normal Map', 'Normal', 'Principled BSDF', 'Normal'),
               ('Image Texture.001', 'Color', 'Principled BSDF', 'Roughness')}
    added = {(ADDED[0], 'Color', ADDED[1], 'Height'),
             (ADDED[1], 'Normal', 'Principled BSDF', 'Normal'),
             ('Image Texture.001', 'Color', 'Principled BSDF', 'Coat Roughness')}
    old_links, new_links = set(map(tuple, old['links'])), set(map(tuple, new['links']))
    assert old_links - new_links == removed
    assert new_links - old_links == added
    new['links'] = copy.deepcopy(old['links'])
    assert normalized == before, 'Undeclared scene, field, material, node, slot or capture change'
    return {'socket_changes': changes, 'removed_links': sorted(removed), 'added_links': sorted(added),
            'added_nodes': ADDED, 'material_slot_changes': [], 'all_other_observed_data_equal': True}


def main():
    args = sys.argv[sys.argv.index('--') + 1:]
    assert args == ['--admitted'], 'Explicit 1 GiB / 120 s / 8 MiB source-build admission required'
    assert bpy.app.version[:3] == (4, 3, 2)
    OUT.mkdir(exist_ok=True)
    coverage = json.loads((AUDIT / 'positive_film_coverage.json').read_text())
    assert coverage['all_positive'] and coverage['positive_coverage_fraction'] == 1.0
    assert coverage['film_m'][0] > 0
    ref = json.loads((AUDIT / 'inputs/primary_layer_model_reference.json').read_text())
    excerpt = AUDIT / 'inputs/blender_4_3_2_closure_excerpt.h'
    assert sha(excerpt) == ref['local_sha256']
    code = excerpt.read_text()
    assert 'bsdf->ior = coat_ior;' in code and 'float eta = ior;' in code
    assert 'eta != 1.0f || thinfilm_thickness > 0.1f' in code
    maproot = PRODUCTION / 'experiments/porcelain_deposition_arrivals/maps'
    mapbytes = {str(p.relative_to(PRODUCTION)): sha(p) for p in sorted(maproot.rglob('*')) if p.is_file()}
    receipts = []
    for rig, (filename, expected) in SOURCES.items():
        source = PRODUCTION / 'experiments/porcelain_macro_normals_control/scenes' / filename
        assert sha(source) == expected
        destination = OUT / f'07_porcelain_effective_layer_{rig}.blend'
        assert not destination.exists(), f'Refusing overwrite: {destination}'
        bpy.ops.wm.open_mainfile(filepath=str(source), load_ui=False)
        bpy.context.preferences.filepaths.save_version = 0
        before = frozen_state()
        roles = assign_layers()
        after = frozen_state()
        delta = verify_only_declared_changes(before, after)
        bpy.ops.wm.save_as_mainfile(filepath=str(destination), compress=True)
        bpy.ops.wm.open_mainfile(filepath=str(destination), load_ui=False)
        reopened = frozen_state()
        assert reopened == after, 'Fresh reopen differs from pre-save candidate'
        assert verify_only_declared_changes(before, reopened) == delta
        assert sha(source) == expected
        receipts.append({'rig': rig, 'control': str(source), 'control_sha256': expected,
                         'scene': str(destination), 'sha256': sha(destination), 'bytes': destination.stat().st_size,
                         'cameras': [o.name for o in bpy.context.scene.objects if o.type == 'CAMERA'],
                         'active_camera': bpy.context.scene.camera.name, 'field_roles': roles,
                         'delta': delta, 'before_state_sha256': digest(before), 'after_state_sha256': digest(after),
                         'fresh_reopen_equal': True, 'source_unchanged': True})
        print('FROZEN_LAYER_SOURCE', json.dumps({k: receipts[-1][k] for k in ['rig','scene','sha256','bytes']}), flush=True)
    assert {str(p.relative_to(PRODUCTION)): sha(p) for p in sorted(maproot.rglob('*')) if p.is_file()} == mapbytes
    analytic = {'air_glaze_normal_incidence_fresnel': ((1.49 - 1) / (1.49 + 1)) ** 2,
                'matched_buried_normal_incidence_fresnel': ((1.49 - 1.49) / (1.49 + 1.49)) ** 2,
                'matched_buried_relative_index': 1.49 / 1.49,
                'zero_film_limit': 'Not exercised: this all-positive-film candidate is not a varying-coverage model. Removing the coat requires restoring the substrate-air IOR.',
                'coat_tint_white': 'No additional absorption is invented.',
                'finite_volume_transport': False}
    assert analytic['matched_buried_normal_incidence_fresnel'] == 0
    assert 0 < analytic['air_glaze_normal_incidence_fresnel'] < 1
    report = {'experiment': 'porcelain_glaze_layer_assignment', 'optical_response_changed': True,
              'coefficient_identical': False, 'effective_layered_closure': True,
              'coverage_receipt': str(AUDIT / 'positive_film_coverage.json'), 'coverage': coverage,
              'primary_reference': ref, 'analytic_limits': analytic,
              'original_map_files_sha256': mapbytes, 'all_map_bytes_unchanged': True,
              'geometry_custom_normals_uv_lights_cameras_world_exposure_unchanged': True,
              'original_body_material_unchanged': True, 'material_slot_changes': [],
              'source_builder_sha256': sha(__file__), 'frozen_observation_helpers_sha256': AUDIT_HELPERS_SHA,
              'scenes': receipts, 'rendered': False, 'visual_acceptance': False, 'selected': False,
              'limitations': ['Opaque substrate response evaluated at the shell; no refracted finite-volume paths.',
                             'Positive film supports this explicit full-coat approximation, not a uniquely calibrated optical model.',
                             'Existing pigment, roughness and stored 1.49 indices are inherited approximations, not measured optical constants.',
                             'No spectral absorption, glaze-body scattering or firing chemistry is introduced.',
                             'Retained map slopes and non-isometric corner atlas remain.']}
    (OUT / 'build_receipt.json').write_text(json.dumps(report, indent=2))
    assert sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()) < 8 * 1024 * 1024
    print('LAYER_SOURCE_ONLY_COMPLETE', str(OUT / 'build_receipt.json'), flush=True)


if __name__ == '__main__':
    main()
