"""Local, no-render grazing preview adapter. Pure numeric functions import without bpy."""
from copy import deepcopy
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def defaults():
    return json.loads((ROOT / 'preview_defaults.json').read_text())

def load_preset():
    config = defaults()
    path = (ROOT / config['preset_file']).resolve()
    if path.parent != ROOT:
        raise ValueError('Preset must be local to this adapter')
    preset = json.loads(path.read_text())
    if preset['id'] != config['default_preview_preset']:
        raise ValueError('Default and preset identity disagree')
    return preset

def vector3(value, name):
    if len(value) != 3 or not all(math.isfinite(float(v)) for v in value):
        raise ValueError(name + ' must contain three finite numbers')
    return [float(v) for v in value]

def positive(value, name):
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(name + ' must be positive and finite') from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + ' must be positive and finite')
    return value

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def subtract(a, b):
    return [x - y for x, y in zip(a, b)]

def plan_preset(*, specimen_width_m=.25, specimen_pivot_m=(0., 0., .008),
                composition=None, crop_width_m=None, crop_center_m=None,
                specimen_bounds_m=None, padding=None, width=None, height=None):
    """Distances are physical metres. Width scales the rig, NEVER the specimen.

    Pivot is the specimen surface-centre datum corresponding to (0, 0, .008)
    in the 0.25 m asphalt reference. Reference mode preserves camera offset.
    full: project explicit world-space bounds; detail: explicit view width.
    """
    p = load_preset()
    config = defaults()
    mode = composition or config['default_composition']
    if mode not in ('reference', 'full', 'detail'):
        raise ValueError('Composition must be reference, full, or detail')
    extent = positive(specimen_width_m, 'specimen_width_m')
    pivot = vector3(specimen_pivot_m, 'specimen_pivot_m')
    factor = extent / p['reference_specimen_width_m']
    width = config['width'] if width is None else width
    height = config['height'] if height is None else height
    if any(type(v) is not int or v <= 0 for v in (width, height)):
        raise ValueError('Output dimensions must be positive integers')
    origin = p['reference_specimen_pivot_m']
    def transform(position):
        return [pivot[i] + (position[i] - origin[i]) * factor for i in range(3)]
    camera = deepcopy(p['camera'])
    key = deepcopy(p['key'])
    camera['location_m'] = transform(camera['location_m'])
    camera['ortho_scale'] *= factor
    camera['clip_start'] *= factor
    camera['clip_end'] *= factor
    key['location_m'] = transform([key['matrix_world'][i][3] for i in range(3)])
    key['target_m'] = pivot
    key['energy'] *= factor ** 2
    key['size'] *= factor
    key['size_y'] *= factor
    key['shadow_soft_size'] *= factor
    matrix = camera['matrix_world']
    axes = [[matrix[i][j] for i in range(3)] for j in (0, 1)]
    axes = [[v / math.sqrt(dot(a, a)) for v in a] for a in axes]
    crop = {'mode': mode, 'crop_width_m': crop_width_m, 'crop_center_m': crop_center_m}
    if mode == 'reference':
        if crop_width_m is not None or crop_center_m is not None:
            raise ValueError('Reference mode cannot change crop; use detail or full')
    else:
        if mode == 'full':
            if crop_width_m is not None or crop_center_m is not None:
                raise ValueError('Full mode derives its crop from specimen bounds')
            if specimen_bounds_m is None or len(specimen_bounds_m) != 2:
                raise ValueError('Full mode requires explicit specimen bounds or objects')
            low, high = [vector3(v, 'specimen_bounds_m') for v in specimen_bounds_m]
            if any(a > b for a, b in zip(low, high)):
                raise ValueError('Specimen bounds are reversed')
            corners = list(itertools.product(*zip(low, high)))
            projected = [[dot(c, axis) for c in corners] for axis in axes]
            spans = [max(a) - min(a) for a in projected]
            pad = positive(config['full_specimen_padding'] if padding is None else padding, 'padding')
            if pad < 1:
                raise ValueError('Full specimen padding cannot be below 1')
            # Blender ortho_scale is the larger image dimension.
            view_width = max(spans[0], spans[1] * width / height) * pad
            positive(view_width, 'projected specimen width')
            center = [(a + b) / 2 for a, b in zip(low, high)]
            crop.update({'specimen_bounds_m': [low, high], 'padding': pad})
        else:
            view_width = positive(crop_width_m, 'crop_width_m')
            center = pivot if crop_center_m is None else vector3(crop_center_m, 'crop_center_m')
        # Shift in the image plane only. Preserve camera orientation and depth.
        delta = subtract(center, camera['location_m'])
        offsets = [dot(delta, axis) for axis in axes]
        camera['location_m'] = [camera['location_m'][i] + sum(offsets[j] * axes[j][i] for j in range(2)) for i in range(3)]
        camera['ortho_scale'] = view_width * max(1., height / width)
        crop.update({'crop_width_m': view_width, 'crop_center_m': center})
    for block in (camera, key):
        for i in range(3):
            block['matrix_world'][i][3] = block['location_m'][i]
    return {'preset_id': p['id'], 'preset_sha256': sha256(ROOT / config['preset_file']),
            'specimen_width_m': extent, 'specimen_pivot_m': pivot, 'scale': factor,
            'composition': crop, 'camera': camera, 'key': key, 'world': p['world'],
            'cycles': p['cycles'], 'display': p['display'], 'render': p['render'],
            'width': width, 'height': height, 'source_geometry_rescaled': False}

def source_guard(source, output_dir):
    """Validate immutable source + brand-new output directory; no writes here."""
    source = Path(source).expanduser().resolve(strict=True)
    raw_output = Path(output_dir).expanduser()
    if raw_output.exists() or raw_output.is_symlink():
        raise FileExistsError('Output directory must not exist: ' + str(raw_output))
    output = raw_output.resolve()
    if source.suffix.lower() != '.blend' or not source.is_file():
        raise ValueError('Source must be an existing .blend file')
    if source == output or output in source.parents:
        raise ValueError('Output cannot contain or overwrite source')
    if not output.parent.is_dir():
        raise ValueError('Output parent must already exist')
    return source, output, sha256(source)

def assert_source_unchanged(source, digest):
    if sha256(source) != digest:
        raise RuntimeError('Immutable source hash changed; reject derivative')

def apply_preset(scene, *, specimen_objects=None, **options):
    """Apply to loaded scene IN MEMORY ONLY. Never saves or renders a file."""
    import bpy
    from mathutils import Vector
    if not math.isclose(scene.unit_settings.scale_length, 1., rel_tol=0., abs_tol=1e-9):
        raise ValueError('Expected metre-based scene (unit scale 1); no automatic geometry/unit conversion')
    if specimen_objects and options.get('specimen_bounds_m') is not None:
        raise ValueError('Choose specimen objects OR explicit bounds, not both')
    if specimen_objects:
        for layer in scene.view_layers:
            layer.update()
        objects = [scene.objects.get(name) for name in specimen_objects]
        if any(o is None or o.type not in {'MESH', 'CURVE', 'CURVES', 'SURFACE', 'FONT'} for o in objects):
            raise ValueError('Every specimen object must be an existing geometry object')
        points = [o.matrix_world @ Vector(c) for o in objects for c in o.bound_box]
        options['specimen_bounds_m'] = [[min(p[i] for p in points) for i in range(3)], [max(p[i] for p in points) for i in range(3)]]
    plan = plan_preset(**options)
    scene.frame_set(plan['render']['frame'])
    # Reuse only this adapter's tagged rig. Existing source cameras are untouched.
    owned = {o.get('cybr_preview_role'): o for o in scene.objects if o.get('cybr_preview_role')}
    key = owned.get('key')
    if key is None:
        data = bpy.data.lights.new('Preview / warm grazing key', 'AREA')
        key = bpy.data.objects.new('Preview / warm grazing key', data)
        scene.collection.objects.link(key)
        key['cybr_preview_role'] = 'key'
    for ob in scene.objects:
        if ob.type == 'LIGHT' and ob != key:
            ob.data = ob.data.copy()
            ob.data.animation_data_clear()
            ob.data.energy = 0.
            ob.data.use_nodes = False
    key.parent = None
    key.animation_data_clear()
    key.constraints.clear()
    key.data.animation_data_clear()
    key.location = plan['key']['location_m']
    key.rotation_mode = 'XYZ'
    key.rotation_euler = plan['key']['rotation_euler']
    key.scale = (1., 1., 1.)
    key.hide_render = False
    for name in ('type', 'energy', 'color', 'shape', 'size', 'size_y', 'shadow_soft_size', 'spread', 'use_shadow', 'use_nodes'):
        setattr(key.data, name, plan['key'][name])
    camera = owned.get('camera')
    if camera is None:
        data = bpy.data.cameras.new('Preview / grazing reference camera')
        camera = bpy.data.objects.new('Preview / grazing reference camera', data)
        scene.collection.objects.link(camera)
        camera['cybr_preview_role'] = 'camera'
    camera.parent = None
    camera.constraints.clear()
    camera.animation_data_clear()
    camera.data.animation_data_clear()
    camera.location = plan['camera']['location_m']
    camera.rotation_mode = plan['camera']['rotation_mode']
    camera.rotation_euler = plan['camera']['rotation_euler']
    camera.scale = (1., 1., 1.)
    for name in ('type', 'ortho_scale', 'lens', 'clip_start', 'clip_end', 'sensor_width', 'sensor_height', 'shift_x', 'shift_y'):
        setattr(camera.data, name, plan['camera'][name])
    camera.data.dof.use_dof = False
    scene.camera = camera
    # Dedicated world prevents source HDRIs, animation, or linked worlds leaking in.
    world = bpy.data.worlds.new('Preview / grazing constant world')
    world.use_nodes = True
    world.node_tree.nodes.clear()
    background = world.node_tree.nodes.new('ShaderNodeBackground')
    background.inputs['Color'].default_value = plan['world']['color_rgba']
    background.inputs['Strength'].default_value = plan['world']['strength']
    output = world.node_tree.nodes.new('ShaderNodeOutputWorld')
    world.node_tree.links.new(background.outputs['Background'], output.inputs['Surface'])
    world.cycles.sampling_method = plan['world']['sampling_method']
    scene.world = world
    scene.render.engine = 'CYCLES'
    for name, value in plan['cycles'].items():
        if not hasattr(scene.cycles, name):
            raise RuntimeError('Blender lacks required Cycles setting: ' + name)
        setattr(scene.cycles, name, value)
    for name in ('view_transform', 'look', 'exposure', 'gamma', 'use_curve_mapping'):
        setattr(scene.view_settings, name, plan['display'][name])
    scene.display_settings.display_device = plan['display']['display_device']
    for name in ('film_transparent', 'pixel_aspect_x', 'pixel_aspect_y', 'dither_intensity'):
        setattr(scene.render, name, plan['render'][name])
    scene.render.resolution_x, scene.render.resolution_y = plan['width'], plan['height']
    scene.render.resolution_percentage = 100
    scene.render.use_border = False
    scene.render.use_crop_to_border = False
    scene.render.use_sequencer = False
    scene.use_nodes = True
    scene.node_tree.nodes.clear()
    scene.use_nodes = False
    scene.render.use_compositing = False
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGB'
    scene.render.image_settings.color_depth = '16'
    scene['preview_preset_id'] = plan['preset_id']
    scene['preview_preset_sha256'] = plan['preset_sha256']
    scene['preview_rig_scale'] = plan['scale']
    scene['preview_composition'] = plan['composition']['mode']
    for layer in scene.view_layers:
        layer.update()
    return plan

def build_job(derivative_source, output_dir, plan):
    p = load_preset()
    job = deepcopy(p['job_defaults'])
    job.update({'source': str(Path(derivative_source).resolve()),
                'source_sha256': sha256(derivative_source),
                'width': plan['width'], 'height': plan['height'],
                'output': str(Path(output_dir).resolve() / 'INTERNAL_PANEL.png'),
                'preview_preset': {'id': plan['preset_id'], 'sha256': plan['preset_sha256'],
                                   'scale': plan['scale'], 'composition': plan['composition']}})
    return job
