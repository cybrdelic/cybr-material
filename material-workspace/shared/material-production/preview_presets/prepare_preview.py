#!/usr/bin/env python3
"""Prepare a fresh derivative .blend + guided job; deliberately cannot render."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from grazing_preview import (apply_preset, assert_source_unchanged, build_job,
                             load_preset, sha256, source_guard)


def arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Immutable source .blend')
    parser.add_argument('--out', required=True, help='Brand-new directory; parent must exist')
    parser.add_argument('--specimen-width-m', required=True, type=float,
                        help='Actual specimen width in metres; never rescales geometry')
    parser.add_argument('--specimen-pivot-m', required=True, nargs=3, type=float,
                        metavar=('X', 'Y', 'Z'), help='Surface-centre datum in world metres')
    parser.add_argument('--composition', choices=('reference', 'full', 'detail'))
    parser.add_argument('--specimen-object', action='append', default=[],
                        help='Exact geometry object name; repeat for full-specimen fit')
    parser.add_argument('--specimen-bounds-m', nargs=6, type=float,
                        metavar=('MIN_X', 'MIN_Y', 'MIN_Z', 'MAX_X', 'MAX_Y', 'MAX_Z'))
    parser.add_argument('--crop-width-m', type=float, help='Detail horizontal field of view')
    parser.add_argument('--crop-center-m', nargs=3, type=float)
    parser.add_argument('--padding', type=float, help='Full-specimen framing margin, >=1')
    parser.add_argument('--width', type=int)
    parser.add_argument('--height', type=int)
    parser.add_argument('--blender', default='blender', help='Blender executable for launcher')
    parser.add_argument('--expect-source-sha256', help='Optional known immutable source digest')
    parser.add_argument('--protect-file', action='append', default=[],
                        help='Optional registry or additional immutable file to hash-check')
    return parser.parse_args(argv)


def inspect_scene(scene):
    camera = scene.camera
    key = next(o for o in scene.objects if o.get('cybr_preview_role') == 'key')
    bg = scene.world.node_tree.nodes.get('Background')
    return {'camera_location_m': list(camera.location),
            'camera_rotation_euler': list(camera.rotation_euler),
            'camera_matrix_world': [list(row) for row in camera.matrix_world],
            'camera_ortho_scale_m': camera.data.ortho_scale,
            'key_location_m': list(key.location), 'key_energy_w': key.data.energy,
            'key_diameter_m': key.data.size, 'key_color': list(key.data.color),
            'other_light_energies': {o.name: o.data.energy for o in scene.objects if o.type == 'LIGHT' and o != key},
            'world_rgba': list(bg.inputs['Color'].default_value),
            'world_strength': bg.inputs['Strength'].default_value,
            'cycles': {name: getattr(scene.cycles, name) for name in load_preset()['cycles']},
            'display': {name: getattr(scene.view_settings, name) for name in ('view_transform', 'look', 'exposure', 'gamma')},
            'png': {'depth': scene.render.image_settings.color_depth,
                    'mode': scene.render.image_settings.color_mode,
                    'dither': scene.render.dither_intensity},
            'rendered': False}


def run_inside_blender(args):
    import bpy
    source, out, digest = source_guard(args.source, args.out)
    if args.expect_source_sha256 and digest != args.expect_source_sha256:
        raise ValueError('Source does not match expected SHA256')
    protected = {str(Path(p).resolve(strict=True)): sha256(p) for p in args.protect_file}
    out.mkdir(exist_ok=False)
    bpy.ops.wm.open_mainfile(filepath=str(source), load_ui=False, use_scripts=False)
    scene = bpy.context.scene
    # A small preservation receipt checks the adapter did not alter specimen placement.
    original_geometry = {o.name: {'type': o.type, 'matrix': [list(r) for r in o.matrix_world],
                         'data_name': o.data.name if o.data else None,
                         'materials': [m.name if m else None for m in getattr(o.data, 'materials', [])]}
                         for o in scene.objects if o.type not in ('LIGHT', 'CAMERA')}
    bounds = args.specimen_bounds_m
    plan = apply_preset(scene, specimen_width_m=args.specimen_width_m,
                        specimen_pivot_m=args.specimen_pivot_m,
                        composition=args.composition, specimen_objects=args.specimen_object,
                        specimen_bounds_m=[bounds[:3], bounds[3:]] if bounds else None,
                        crop_width_m=args.crop_width_m, crop_center_m=args.crop_center_m,
                        padding=args.padding, width=args.width, height=args.height)
    for name, before in original_geometry.items():
        ob = scene.objects[name]
        after = {'type': ob.type, 'matrix': [list(r) for r in ob.matrix_world],
                 'data_name': ob.data.name if ob.data else None,
                 'materials': [m.name if m else None for m in getattr(ob.data, 'materials', [])]}
        if before != after:
            raise RuntimeError('Source geometry placement or material slots changed: ' + name)
    derivative = out / 'preview.blend'
    if derivative.exists():
        raise FileExistsError('Refusing derivative overwrite')
    bpy.ops.wm.save_as_mainfile(filepath=str(derivative), compress=True,
                               relative_remap=True, copy=True, check_existing=True)
    # Verify settings from the saved derivative, not just pre-save memory.
    bpy.ops.wm.open_mainfile(filepath=str(derivative), load_ui=False, use_scripts=False)
    scene = bpy.context.scene
    assert_source_unchanged(source, digest)
    for path, expected in protected.items():
        assert_source_unchanged(path, expected)
    job = build_job(derivative, out, plan)
    receipt = {'preset_id': plan['preset_id'], 'preset_sha256': plan['preset_sha256'],
               'immutable_source': str(source), 'source_sha256_before': digest,
               'source_sha256_after': sha256(source), 'source_unchanged': True,
               'derivative': str(derivative), 'derivative_sha256': sha256(derivative),
               'protected_files_sha256': protected,
               'geometry_placement_and_material_slots_unchanged': True,
               'plan': plan, 'inspection': inspect_scene(scene),
               'blender_version': bpy.app.version_string,
               'rendered': False, 'selected_registry_modified': False,
               'limitations': ['Presentation changes lighting; only call comparisons matched when both use the same preset and framing.',
                              'Geometry and physical feature dimensions are not scaled.',
                              'Emissive source geometry remains unchanged; the preset only disables other explicit lights.',
                              'No render or visual equivalence validation performed.']}
    (out / 'job.json').write_text(json.dumps(job, indent=2) + '\n')
    (out / 'preparation_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'prepared': str(out), 'job': str(out / 'job.json'), 'rendered': False}))


def main():
    inside = '--' in sys.argv
    argv = sys.argv[sys.argv.index('--') + 1:] if inside else sys.argv[1:]
    args = arguments(argv)
    if inside:
        run_inside_blender(args)
        return
    # Validate first; no source load or output write on known collisions.
    source, out, digest = source_guard(args.source, args.out)
    if args.expect_source_sha256 and digest != args.expect_source_sha256:
        raise ValueError('Source does not match expected SHA256')
    cmd = [args.blender, '--background', '--factory-startup', '--disable-autoexec',
           '-t', '1', '--python-exit-code', '1', '--python', str(Path(__file__).resolve()),
           '--', *argv]
    subprocess.run(cmd, check=True)
    assert_source_unchanged(source, digest)


if __name__ == '__main__':
    main()
