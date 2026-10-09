"""Read-only Blender data extraction. Does not render, load image pixels or save scenes."""
import bpy, hashlib, json, sys
from pathlib import Path

src = Path(sys.argv[sys.argv.index('--') + 1]).resolve()
out = Path(sys.argv[sys.argv.index('--') + 2]).resolve()
before = hashlib.sha256(src.read_bytes()).hexdigest()
bpy.ops.wm.open_mainfile(filepath=str(src), load_ui=False)
scene = bpy.context.scene
result = {'scene': str(src), 'sha256': before, 'blender': bpy.app.version_string,
          'objects': {}, 'materials': {}}
for ob in scene.objects:
    if ob.type != 'MESH' or not ob.name.startswith('07 new /'):
        continue
    me = ob.data
    me.calc_loop_triangles()
    me.calc_tangents(uvmap='UVMap')
    result['objects'][ob.name] = {
        'matrix': [list(v) for v in ob.matrix_world],
        'has_custom_normals': me.has_custom_normals,
        'vertices': [list(v.co) for v in me.vertices],
        'vertex_normals': [list(v.normal) for v in me.vertices],
        'polygons': [{'index': p.index, 'normal': list(p.normal),
                      'smooth': p.use_smooth, 'loops': list(p.loop_indices)}
                     for p in me.polygons],
        'loops': [{'vertex': l.vertex_index, 'normal': list(me.corner_normals[l.index].vector),
                   'tangent': list(l.tangent), 'bitangent_sign': l.bitangent_sign,
                   'uv': list(me.uv_layers.active.data[l.index].uv)} for l in me.loops],
        'triangles': [{'loops': list(t.loops), 'polygon': t.polygon_index}
                      for t in me.loop_triangles],
        'sharp_edges': sum(e.use_edge_sharp for e in me.edges)}
for mat in bpy.data.materials:
    if not mat.use_nodes:
        continue
    result['materials'][mat.name] = {
        'nodes': [{'name': n.name, 'type': n.bl_idname,
                   **({'space': n.space, 'uv_map': n.uv_map,
                       'strength': n.inputs['Strength'].default_value}
                      if n.type == 'NORMAL_MAP' else {}),
                   **({'filepath': n.image.filepath, 'colorspace': n.image.colorspace_settings.name,
                       'interpolation': n.interpolation, 'extension': n.extension}
                      if n.type == 'TEX_IMAGE' and n.image else {})}
                  for n in mat.node_tree.nodes],
        'links': [[l.from_node.name, l.from_socket.name, l.to_node.name, l.to_socket.name]
                  for l in mat.node_tree.links]}
assert hashlib.sha256(src.read_bytes()).hexdigest() == before
out.write_text(json.dumps(result, separators=(',', ':')))
assert out.stat().st_size < 5 * 1024 * 1024
print('READ_ONLY_SCENE_EXTRACTED', out, out.stat().st_size)
