"""Build original native-4096 panel scenes, without rendering or image resampling.

Private panel builder based on build_original_native4k.py. The recovered images
stay external and byte-identical. SI-unit edge easing and noncollapsed side UVs
are added to the closed macro-height specimen mesh.
"""
import bpy
import bmesh
import sys
import json
import hashlib
import math
import resource
from pathlib import Path
import numpy as np
from PIL import Image

RECOVERED = Path('/workspace/shared/material-native4k/original-verified')
OUT = Path('/workspace/shared/material-slabs')
sys.path.insert(0, str(OUT / 'source'))
from studio import configure

CID = sys.argv[sys.argv.index('--') + 1]
assert CID[:2] in {'02', '03', '05', '06', '07', '08', '09'}, CID
P = RECOVERED / 'materials' / CID
COMMIT = 'de4146ac759164e29323b152d37f64417c00c3fb'
metadata = json.loads((P / 'material.json').read_text())
assert metadata['resolution'] == 4096
recovery = json.loads((RECOVERED / 'recovery_receipt.json').read_text())
assert recovery['commit'] == COMMIT
verified = {entry['path']: entry for entry in recovery['verified_files']}
tree = json.loads((RECOVERED / 'tree.json').read_text())
assert tree['commit'] == COMMIT
expected = {entry['path']: entry for entry in tree['files']}

sources = {}
for file in [P / name for name in metadata['maps']] + [P / 'material.json']:
    relative = str(file.relative_to(RECOVERED))
    recorded = verified[relative]
    data = file.read_bytes()
    sha256 = hashlib.sha256(data).hexdigest()
    git_blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    assert sha256 == recorded['sha256'], relative
    assert git_blob == recorded['git_blob_sha'] == expected[relative]['sha'], relative
    dimensions = None
    if file.suffix == '.png':
        with Image.open(file) as image:
            assert image.size == (4096, 4096), (relative, image.size)
            image.verify()
        dimensions = [4096, 4096]
    sources[file.name] = {'sha256': sha256, 'git_blob_sha': git_blob,
                          'bytes': len(data), 'dimensions': dimensions,
                          'source_url': recorded['url']}
    del data
print('VERIFIED_NATIVE_MAPS', CID, len(sources), flush=True)

bpy.ops.wm.read_factory_settings(use_empty=True)
s = bpy.context.scene
material = bpy.data.materials.new('Native4K / ' + metadata['name'])
material.use_nodes = True
nodes = material.node_tree.nodes
links = material.node_tree.links
principled = nodes.get('Principled BSDF')
for socket, key, default in [('IOR', 'ior', 1.5), ('Coat Weight', 'coat', 0),
                              ('Coat Roughness', 'coat_roughness', .3),
                              ('Anisotropic', 'anisotropy', 0),
                              ('Sheen Weight', 'sheen', 0)]:
    principled.inputs[socket].default_value = metadata.get(key, default)
for name, socket in [('BaseColor', 'Base Color'), ('Roughness', 'Roughness'),
                     ('Metallic', 'Metallic'), ('Opacity', 'Alpha')]:
    image = bpy.data.images.load(str(P / (name + '.png')))
    assert list(image.size) == [4096, 4096]
    image.colorspace_settings.name = 'sRGB' if name == 'BaseColor' else 'Non-Color'
    node = nodes.new('ShaderNodeTexImage')
    node.image = image
    node.label = name + ' / original 4096'
    links.new(node.outputs['Color'], principled.inputs[socket])
image = bpy.data.images.load(str(P / 'Normal_Micro_OpenGL.png'))
assert list(image.size) == [4096, 4096]
image.colorspace_settings.name = 'Non-Color'
node = nodes.new('ShaderNodeTexImage')
node.image = image
node.label = 'Residual native-4096 micro normal paired with Height_Macro geometry'
normal = nodes.new('ShaderNodeNormalMap')
normal.uv_map = 'UVMap'
links.new(node.outputs['Color'], normal.inputs['Color'])
links.new(normal.outputs['Normal'], principled.inputs['Normal'])

width = metadata['tile_m']
depth = metadata.get('tile_y_m', width)
thickness = {'02': .020, '03': .018, '05': .003, '06': .003,
             '07': .008, '08': .003, '09': .015}[CID[:2]]
edge_radius = {'02': .0005, '03': .00035, '05': .00015, '06': .00015,
               '07': .0003, '08': .0003, '09': .0005}[CID[:2]]
with Image.open(P / 'Height_Macro.png') as image:
    assert image.size == (4096, 4096)
    a = np.array(image, dtype='f4') / 65535
N = 256
indices = np.rint(np.linspace(0, 4095, N + 1)).astype(int)
height = (a[np.ix_(indices, indices)] - .5) * metadata['height_scale_m']
del a
vertices = []
uv = []
faces = []
for j in range(N + 1):
    for i in range(N + 1):
        vertices.append(((i / N - .5) * width, (.5 - j / N) * depth,
                         thickness + float(height[j, i])))
        uv.append((i / N, 1 - j / N))
for j in range(N):
    for i in range(N):
        k = j * (N + 1) + i
        faces.append((k, k + N + 1, k + N + 2, k + 1))
boundary = (list(range(N + 1)) + [j * (N + 1) + N for j in range(1, N + 1)]
            + [N * (N + 1) + i for i in range(N - 1, -1, -1)]
            + [j * (N + 1) for j in range(N - 1, 0, -1)])
bottom = []
for k in boundary:
    bottom.append(len(vertices))
    x, y, z = vertices[k]
    vertices.append((x, y, 0))
    uv.append(uv[k])
for i, k in enumerate(boundary):
    q = (i + 1) % len(boundary)
    faces.append((k, boundary[q], bottom[q], bottom[i]))
faces.append(tuple(bottom))
mesh = bpy.data.meshes.new(CID + ' / closed true-thickness slab')
mesh.from_pydata(vertices, [], faces)
mesh.update()
ob = bpy.data.objects.new(CID + ' / full native tile slab', mesh)
s.collection.objects.link(ob)
mesh.materials.append(material)
layer = mesh.uv_layers.new(name='UVMap')
layer.active_render = True
for face in mesh.polygons:
    for loop_index in face.loop_indices:
        layer.data[loop_index].uv = uv[mesh.loops[loop_index].vertex_index]
    face.use_smooth = face.index < N * N
    if N * N <= face.index < N * N + len(boundary):
        # Continue the physical tile scale down every thin side, rather than
        # collapsing its UVs to one border line as on a simple extrusion.
        va, vb, bbottom, abottom = [mesh.vertices[k].co for k in face.vertices]
        along_x = abs(vb.x - va.x) > abs(vb.y - va.y)
        coordinate = 1 if along_x else 0
        extent = depth if along_x else width
        side_sign = 1 if (va.y if along_x else va.x) > 0 else -1
        for loop_offset, top_point in [(2, vb), (3, va)]:
            loop_index = face.loop_indices[loop_offset]
            layer.data[loop_index].uv[coordinate] += side_sign * top_point.z / extent

# Check base geometry and UVs before adding the physical edge modifier.
bm = bmesh.new()
bm.from_mesh(mesh)
bm.normal_update()
assert all(len(edge.link_faces) == 2 for edge in bm.edges), 'Open mesh edge'
assert all(edge.is_contiguous for edge in bm.edges), 'Inconsistent winding'
signed_volume = bm.calc_volume(signed=True)
assert signed_volume > 0, signed_volume
bm.free()
assert all(p.normal.z > 0 for p in list(mesh.polygons)[:N * N])
assert mesh.polygons[-1].normal.z < -.999
min_uv_area = math.inf
for face in mesh.polygons:
    coordinates = [layer.data[i].uv for i in face.loop_indices]
    area = abs(sum(v.x * coordinates[(i + 1) % len(coordinates)].y
                   - coordinates[(i + 1) % len(coordinates)].x * v.y
                   for i, v in enumerate(coordinates))) / 2
    assert area > 1e-12, (face.index, area)
    min_uv_area = min(min_uv_area, area)

bevel = ob.modifiers.new('Physical %.2f mm eased edges' % (edge_radius * 1000), 'BEVEL')
bevel.width = edge_radius
bevel.segments = 3
bevel.limit_method = 'ANGLE'
bevel.angle_limit = math.radians(30)
bevel.use_clamp_overlap = True
bevel.affect = 'EDGES'

# Verify the actual evaluated render geometry is still closed and outward.
bpy.context.view_layer.update()
evaluated = ob.evaluated_get(bpy.context.evaluated_depsgraph_get())
evaluated_mesh = evaluated.to_mesh()
bm = bmesh.new()
bm.from_mesh(evaluated_mesh)
bm.normal_update()
assert all(len(edge.link_faces) == 2 for edge in bm.edges), 'Beveled mesh is open'
assert all(edge.is_contiguous for edge in bm.edges), 'Beveled winding inconsistent'
evaluated_volume = bm.calc_volume(signed=True)
assert evaluated_volume > 0
beveled_counts = {'vertices': len(bm.verts), 'edges': len(bm.edges), 'faces': len(bm.faces)}
evaluated_min_z = min(v.co.z for v in evaluated_mesh.vertices)
bm.free()
evaluated.to_mesh_clear()

ob['physical_tile_m'] = [width, depth]
ob['catalog_id'] = CID
ob['macro_geometry_grid'] = N
ob['maps_native_resolution'] = 4096
ob['no_double_relief'] = 'Height_Macro mesh plus residual Normal_Micro_OpenGL only'
ob['edge_radius_m'] = edge_radius
report = configure(CID, [ob])
floor = bpy.data.objects['Slab / neutral floor']
assert abs((evaluated_min_z - floor.location.z) - .000001) < 1e-8
variant = ('Exact recovered original native-4096 map material; distinct from later '
           'object-space procedural volume shader') if CID[:2] == '03' else 'Exact recovered original native-4096 material'
report.update(
    native_map_source=str(P), repository_commit=COMMIT,
    provenance='Byte-verified original native-4096 maps; no upsampling or regeneration',
    physical_height_scale_m=metadata['height_scale_m'], height_midlevel=.5,
    physical_thickness_m=thickness, physical_edge_bevel_m=edge_radius,
    bevel_segments=3, bevel_limit_degrees=30, bevel_overlap_clamp=True,
    shader_variant=variant, source_files=sources,
    source_hashes={name: entry['sha256'] for name, entry in sources.items()},
    builder_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    dependencies_external_not_packed=True,
    loaded_native_images={image.name: list(image.size) for image in bpy.data.images},
    manifold_two_faces_per_edge=True, normals_consistent_outward=True,
    signed_volume_m3=signed_volume, beveled_signed_volume_m3=evaluated_volume,
    evaluated_geometry_counts=beveled_counts, actual_geometry_min_z_m=evaluated_min_z,
    floor_gap_m=evaluated_min_z-floor.location.z,
    top_uv_domain=[0, 0, 1, 1], side_uvs='Physical-scale continuation from border; nonzero area',
    min_base_face_uv_area=min_uv_area, map_pixel_pitch_m=[width / 4096, depth / 4096],
    macro_geometry_grid=[N, N], macro_geometry_sample_pitch_m=[width / N, depth / N],
    relief_note='Macro geometry is sampled from the native 4096 map on a 256x256 grid; residual detail is shaded by native 4096 normal map.',
    estimated_build_memory_gib=2.0,
)
s['slab_source_version'] = variant
s.name = 'Slab / ' + CID + ' / original native4K'
target = OUT / 'scenes' / (CID + '_original_native4k_slab.blend')
target.parent.mkdir(exist_ok=True)
bpy.ops.wm.save_as_mainfile(filepath=str(target), compress=True)
report['scene'] = str(target)
report['scene_sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
report['build_peak_rss_mib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
(OUT / 'receipts').mkdir(exist_ok=True)
receipt = OUT / 'receipts' / (CID + '_original_native4k_slab.json')
receipt.write_text(json.dumps(report, indent=2))
print('NATIVE4K_SLAB_READY', json.dumps({'scene': str(target), 'receipt': str(receipt),
      'sha256': report['scene_sha256'], 'build_peak_rss_mib': report['build_peak_rss_mib']}), flush=True)
