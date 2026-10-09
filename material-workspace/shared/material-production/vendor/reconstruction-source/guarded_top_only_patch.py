"""UNEXECUTED Blender 4.3.2 adapter. Resource admission is mandatory.

Makes a NEW candidate .blend; never renders, rebakes, edits source maps, changes
physical parameters, or overwrites input scenes. Corrects only source top faces.
Generated bevel / source side / source bottom retain the original material.

Usage (ONLY after the rendering owner provides a valid admission file):
blender -b --factory-startup --python guarded_top_only_patch.py -- \
  --material 02_roman_travertine --admission-file /approved/admission.json

A hash in an admission file is a resource coordination guard, not new user
permission. This script must still only be run within authorized task scope.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import struct
import sys
import numpy as np

ROOT = Path('/workspace/shared/material-causal-rebuild/curved-normal-transfer')
SOURCE_ROOT = Path('/workspace/shared/material-slabs/scenes')
MAP_ROOT = Path('/workspace/shared/material-native4k/original-verified/materials')
SOURCES = {
    '01_calacatta_oro': '01_calacatta_oro_original_native4k_slab_r2.blend',
    '02_roman_travertine': '02_roman_travertine_original_native4k_slab.blend',
    '03_american_walnut': '03_american_walnut_original_native4k_slab.blend',
    '05_champagne_brass': '05_champagne_brass_original_native4k_slab_r2_uv_tangent.blend',
    '06_blackened_steel': '06_blackened_steel_original_native4k_slab.blend',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def bound_image_sha(image, bpy):
    if image.packed_file:
        return hashlib.sha256(image.packed_file.data).hexdigest()
    return sha(Path(bpy.path.abspath(image.filepath)))


def require(ok, message):
    if not ok: raise RuntimeError(message)


def resource_guard(path, src, cid):
    a = json.loads(Path(path).read_text())
    require(a.get('authorized') is True, 'Resource admission was not granted')
    require(a.get('owner') == '/root/finish_clean_material_panels', 'Wrong resource-admission owner')
    require(a.get('action') == 'blender_metadata_and_candidate_save_no_render', 'Wrong admitted action')
    require(a.get('material_id') == cid and a.get('source_scene_sha256') == sha(src), 'Admission does not bind this material/scene')
    require(a.get('maximum_output_bytes', 0) >= 64*1024*1024, 'Require a 64 MiB output envelope')
    expiry = dt.datetime.fromisoformat(a['expires_utc'].replace('Z', '+00:00'))
    require(expiry.tzinfo is not None and dt.datetime.now(dt.timezone.utc) < expiry, 'Admission expired or lacks timezone')
    require(shutil.disk_usage(ROOT).free >= 256*1024*1024, 'Under 256 MiB free; do not start')
    return a


def array(collection, name, width, dtype='f4'):
    a = np.empty((len(collection), width), dtype=dtype)
    collection.foreach_get(name, a.ravel())
    return a


def mesh_arrays(mesh):
    mesh.calc_loop_triangles()
    data = {
        'vertices': array(mesh.vertices, 'co', 3),
        'edges': array(mesh.edges, 'vertices', 2, 'i4'),
        'loop_vertices': array(mesh.loops, 'vertex_index', 1, 'i4'),
        'face_loop_start': array(mesh.polygons, 'loop_start', 1, 'i4'),
        'face_loop_total': array(mesh.polygons, 'loop_total', 1, 'i4'),
        'face_smooth': array(mesh.polygons, 'use_smooth', 1, '?'),
        'corner_normals': array(mesh.corner_normals, 'vector', 3),
        'triangles': array(mesh.loop_triangles, 'loops', 3, 'i4'),
        'triangle_faces': array(mesh.loop_triangles, 'polygon_index', 1, 'i4'),
    }
    for uv in mesh.uv_layers:
        data['uv:' + uv.name] = array(uv.data, 'uv', 2)
    return data


def array_digest(data):
    h = hashlib.sha256()
    for name, a in sorted(data.items()):
        h.update(name.encode()); h.update(str(a.shape).encode()); h.update(a.dtype.str.encode()); h.update(a.tobytes())
    return h.hexdigest()


def eval_snapshot(bpy, ob):
    bpy.context.view_layer.update()
    deps = bpy.context.evaluated_depsgraph_get()
    obj = ob.evaluated_get(deps)
    mesh = obj.to_mesh(preserve_all_data_layers=True, depsgraph=deps)
    try:
        return mesh_arrays(mesh), array(mesh.polygons, 'material_index', 1, 'i4')[:, 0]
    finally: obj.to_mesh_clear()


def primitive(value):
    if isinstance(value, (str, int, float, bool)) or value is None: return value
    try: return list(value)
    except TypeError: return str(value)


def scalar_rna(value, exclude=()):
    result = {}
    for prop in value.bl_rna.properties:
        if prop.identifier in {'rna_type', *exclude}: continue
        if prop.type in {'BOOLEAN', 'INT', 'FLOAT', 'STRING', 'ENUM'}:
            result[prop.identifier] = primitive(getattr(value, prop.identifier))
    return result


def shader_signature(mat, image_override=None, normal_override=None):
    rows = []
    for node in mat.node_tree.nodes:
        row = {'name':node.name, 'type':node.bl_idname,
               'settings':scalar_rna(node, ('select', 'dimensions', 'location', 'width', 'height')),
               'inputs':[(s.name, primitive(s.default_value)) for s in node.inputs if hasattr(s, 'default_value')]}
        if node.type == 'TEX_IMAGE':
            row['image'] = {'name':node.image.name, 'path':node.image.filepath,
                            'colorspace':node.image.colorspace_settings.name} if node.image else None
            if image_override and node.name == image_override[0]: row['image'] = image_override[1]
        if normal_override and node.name == normal_override[0]: row['settings']['space'] = normal_override[1]
        rows.append(row)
    return {'nodes':rows, 'links':sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier)
                                       for l in mat.node_tree.links),
            'material_settings':scalar_rna(mat, ('name', 'name_full', 'is_updated', 'is_updated_data', 'session_uid'))}


def scene_controls(bpy):
    s = bpy.context.scene
    return {'objects':{o.name:{'matrix': [list(r) for r in o.matrix_world],
                                      'visibility': [o.hide_render, o.hide_viewport]} for o in s.objects},
            'lights':{x.name:scalar_rna(x) for x in bpy.data.lights},
            'cameras':{x.name:scalar_rna(x) for x in bpy.data.cameras},
            'camera':s.camera.name if s.camera else None,
            'render':scalar_rna(s.render), 'cycles':scalar_rna(s.cycles),
            'view':scalar_rna(s.view_settings), 'display':scalar_rna(s.display_settings),
            'world':s.world.name if s.world else None,
            'world_nodes':[(n.name, scalar_rna(n), [(v.name,primitive(v.default_value)) for v in n.inputs if hasattr(v,'default_value')])
                           for n in s.world.node_tree.nodes] if s.world and s.world.use_nodes else None}


def verify_top_chart(ob, metadata):
    g = int(ob.get('macro_geometry_grid', -1))
    require(g == 256, 'Only inspected 256-square original slabs are supported')
    uv = ob.data.uv_layers.get('UVMap')
    require(uv is not None and uv.active_render, 'Expected active render UVMap')
    require(len(ob.data.materials) == 1 and all(p.material_index == 0 for p in ob.data.polygons), 'Expected one unchanged source material')
    width, depth = metadata['tile_m'], metadata.get('tile_y_m', metadata['tile_m'])
    for j in range(g):
        for i in range(g):
            k = j*(g+1)+i
            p = ob.data.polygons[j*g+i]
            require(tuple(p.vertices) == (k,k+g+1,k+g+2,k+1), 'Unexpected source top topology/winding')
            for li in p.loop_indices:
                co = ob.data.vertices[ob.data.loops[li].vertex_index].co
                u,v = uv.data[li].uv
                require(abs(co.x - (u-.5)*width) < 1e-6 and abs(co.y - (v-.5)*depth) < 1e-6,
                        'Full-normal object chart axes no longer match geometry/UV')
    return g*g


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--material', choices=SOURCES, required=True)
    ap.add_argument('--admission-file', required=True)
    args = ap.parse_args(argv)
    cid = args.material; src = SOURCE_ROOT/SOURCES[cid]; maps = MAP_ROOT/cid
    admission = resource_guard(args.admission_file, src, cid)
    # This import deliberately follows resource validation. Help/compilation
    # and analytic unit tests work without starting Blender.
    import bpy
    require(tuple(bpy.app.version) == (4,3,2), 'Only pinned Blender 4.3.2 is source-audited')
    outdir = ROOT/'candidates'/cid
    require(not outdir.exists(), 'Candidate directory exists; no overwrite/retry by default')
    source_sha = sha(src)
    map_hashes = {p.name:sha(p) for p in maps.glob('*.png')}
    fullpath = maps/'Normal_OpenGL.png'
    with fullpath.open('rb') as f: header=f.read(29)
    require(header[:8] == b'\x89PNG\r\n\x1a\n' and struct.unpack('>II',header[16:24]) == (4096,4096)
            and header[24:26] == bytes((8,2)), 'Expected original native 4096 RGB8 full normal')
    metadata_path = maps/'material.json'
    if not metadata_path.exists():
        metadata_path = Path('/workspace/shared/material-slabs/local-recovery-candidates')/cid/'material.json'
    metadata = json.loads(metadata_path.read_text())
    bpy.ops.wm.open_mainfile(filepath=str(src))
    ob = bpy.data.objects[cid+' / full native tile slab']
    require(not bpy.context.scene.render.use_motion_blur, 'Motion blur requires a separate time-dependent contract')
    require(not ob.animation_data and not ob.data.animation_data and not ob.data.shape_keys, 'Deformation/animation is outside scope')
    active = [m for m in ob.modifiers if m.show_render or m.show_viewport]
    require(len(active) == 1 and active[0].type == 'BEVEL', 'Unexpected modifier stack')
    bevel = active[0]
    require(bevel.show_viewport == bevel.show_render and bevel.show_render, 'Viewport/render modifier evaluation differs')
    require(bevel.limit_method == 'ANGLE' and bevel.affect == 'EDGES', 'Unexpected bevel rule')
    require(bevel.material in (-1,0), 'Existing bevel uses a different material')
    n_top = verify_top_chart(ob, metadata)
    mat = ob.data.materials[0]
    require(mat and mat.use_nodes and not mat.animation_data and not mat.node_tree.animation_data, 'Unsupported material')
    normal_nodes = [n for n in mat.node_tree.nodes if n.type == 'NORMAL_MAP']
    require(len(normal_nodes) == 1, 'Expected one original Normal Map node')
    normal = normal_nodes[0]
    require(normal.space == 'TANGENT' and not normal.inputs['Strength'].is_linked and normal.inputs['Strength'].default_value == 1,
            'Only unchanged tangent Normal Map strength 1 is supported')
    require(len(normal.inputs['Color'].links) == 1, 'Unexpected normal-color binding')
    image_node = normal.inputs['Color'].links[0].from_node
    require(image_node.type == 'TEX_IMAGE' and image_node.image and 'Normal_Micro_OpenGL' in image_node.image.name,
            'Normal Map input is not the original micro normal image')
    require(image_node.image.colorspace_settings.name == 'Non-Color' and image_node.projection == 'FLAT', 'Unexpected map interpretation')
    require(not image_node.inputs['Vector'].is_linked, 'Custom UV mapping requires explicit chart qualification')
    require(all(link.to_node.type == 'BSDF_PRINCIPLED' and link.to_socket.name == 'Normal' for link in normal.outputs['Normal'].links),
            'Normal output has unexpected consumers')
    require(len(normal.outputs['Normal'].links) > 0, 'Normal output is unused')
    for node in mat.node_tree.nodes:
        if node.type == 'TEX_IMAGE' and node.image:
            basename = Path(bpy.path.abspath(node.image.filepath)).name
            require(basename in map_hashes, 'Source image is outside the verified original map set')
            require(bound_image_sha(node.image, bpy) == map_hashes[basename], 'Bound source image bytes differ from original map')
    base_before = mesh_arrays(ob.data); evaluated_before, assignments_before = eval_snapshot(bpy, ob)
    require(np.all(assignments_before == 0), 'Existing evaluated material partition is unexpected')
    controls = scene_controls(bpy); original_signature = shader_signature(mat)
    modifier_before = scalar_rna(bevel)
    corrected = mat.copy(); corrected.name = mat.name+' / top-only exact full-normal candidate'
    nfix = corrected.node_tree.nodes[normal.name]
    ifix = corrected.node_tree.nodes[image_node.name]
    old_image = {'name':ifix.image.name, 'path':ifix.image.filepath, 'colorspace':ifix.image.colorspace_settings.name}
    # Existing data maps are untouched. The new shader loads original bytes.
    newimage = bpy.data.images.load(str(fullpath), check_existing=False)
    newimage.colorspace_settings.name = 'Non-Color'
    require(list(newimage.size) == [4096,4096], 'Full-normal image loading failed')
    ifix.image = newimage; nfix.space = 'OBJECT'
    ob.data.materials.append(corrected)
    for p in list(ob.data.polygons)[:n_top]: p.material_index = 1
    # Critical: otherwise bevel can inherit the top material and be flattened.
    bevel.material = 0
    base_after = mesh_arrays(ob.data); evaluated_after, assignments_after = eval_snapshot(bpy, ob)
    require(array_digest(base_before) == array_digest(base_after), 'Base geometry/normals/UVs changed')
    require(array_digest(evaluated_before) == array_digest(evaluated_after), 'Evaluated geometry/normals/UVs changed')
    require(set(assignments_after.tolist()) == {0,1}, 'Missing controlled surface partition')
    require(np.count_nonzero(assignments_after == 1) == n_top, 'Not all source top faces survived as corrected faces')
    require(shader_signature(mat) == original_signature, 'Original side/bevel material changed')
    require(shader_signature(corrected, (image_node.name,old_image), (normal.name,'TANGENT')) == original_signature,
            'A shader property changed outside the full-normal binding/space')
    modifier_after = scalar_rna(bevel); modifier_after['material'] = modifier_before['material']
    require(modifier_after == modifier_before, 'A bevel setting changed outside material routing')
    require(scene_controls(bpy) == controls, 'Camera/light/environment/render controls changed')
    require(sha(src) == source_sha and {p.name:sha(p) for p in maps.glob('*.png')} == map_hashes, 'Source scene/maps changed')
    # Actual evaluated triangle topology + UV barycentrics inputs. No guessed
    # grid diagonal; no unqualified tangent attributes advertised as Cycles.
    tri = evaluated_after['triangles']; uv = evaluated_after['uv:UVMap']
    face_ids = evaluated_after['triangle_faces'][:,0]
    export = {'triangle_loop_indices':tri, 'triangle_uv':uv[tri],
              'triangle_positions':evaluated_after['vertices'][evaluated_after['loop_vertices'][tri,0]],
              'triangle_corner_normals':evaluated_after['corner_normals'][tri],
              'triangle_material_indices':assignments_after[face_ids],
              'triangle_polygon_indices':face_ids,
              'object_to_world':np.array(ob.matrix_world)}
    outdir.mkdir(parents=True)
    np.savez_compressed(outdir/'evaluated_triangle_contract.npz', **export)
    dest = outdir/(cid+'_top_only_normal_candidate.blend')
    bpy.ops.wm.save_as_mainfile(filepath=str(dest), compress=True)
    require(sha(src) == source_sha, 'Pristine source scene was modified')
    require(sum(p.stat().st_size for p in outdir.iterdir()) < admission['maximum_output_bytes'], 'Output envelope exceeded')
    report = {'status':'candidate_built_not_render_verified', 'scope':'source top only; bevel/side/bottom excluded unchanged control',
              'material_id':cid,'source_scene':str(src),'source_sha256':source_sha,
              'candidate':str(dest),'candidate_sha256':sha(dest),'map_sha256':map_hashes,
              'base_geometry_uv_normal_sha256':array_digest(base_after),
              'evaluated_geometry_uv_normal_sha256':array_digest(evaluated_after),
              'corrected_evaluated_faces':int(np.count_nonzero(assignments_after==1)),
              'unchanged_control_evaluated_faces':int(np.count_nonzero(assignments_after==0)),
              'target_contract':'OBJECT Normal Map of original linear full-normal RGB; strength 1; original texture filtering',
              'frame_provenance':'no tangent frame used by top-only route; exported normals diagnostic only',
              'physical_shader_coefficients_changed':False, 'source_pixel_maps_changed':False,
              'camera_light_world_geometry_changed':False, 'render_run':False,
              'limitations':['Not a whole-slab correction','Shading normal target does not restore unresolved geometry/occlusion',
                             'Renderer closure restrictions remain; raw normal shader equality needs an admitted renderer fixture',
                             'Original RGB8 normal quantization remains','No optical/formation-model validation']}
    (outdir/'candidate_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    argv = sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else sys.argv[1:]
    main(argv)
