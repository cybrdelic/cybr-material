"""Render the authored 4K surfaces with Blender Cycles CPU path tracing.

Example: blender -b --python source/render_path_traced.py -- --kind macros
No denoising, image blur, sharpening, artificial grain, or image upscaling.
"""
import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import bpy
from mathutils import Vector

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'path_traced'
argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
parser = argparse.ArgumentParser()
parser.add_argument('--kind', choices=('macros', 'scenes', 'all'), default='all')
parser.add_argument('--only', nargs='*', default=[])
parser.add_argument('--size', type=int, default=1280)
parser.add_argument('--samples', type=int, default=1024)
parser.add_argument('--threshold', type=float, default=.003)
parser.add_argument('--min-samples', type=int, default=128)
parser.add_argument('--draft', action='store_true')
parser.add_argument('--use-final-maps', action='store_true', help='Draft output with native material maps')
parser.add_argument('--build-only', action='store_true')
args = parser.parse_args(argv)
sys.argv = ['blender', '--'] + (['--draft'] if args.draft and not args.use_final_maps else [])
module = importlib.util.spec_from_file_location('macro_setup', ROOT / 'source/render_macros.py')
macro_setup = importlib.util.module_from_spec(module)
module.loader.exec_module(macro_setup)
atelier = macro_setup.atelier
macro_setup.args.size = args.size
macro_setup.args.samples = args.samples
atelier.args.samples = args.samples


def configure(scene, name):
    scene.name = name
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = args.samples
    scene.cycles.use_denoising = False
    scene.cycles.use_adaptive_sampling = True
    scene.cycles.adaptive_threshold = args.threshold
    scene.cycles.adaptive_min_samples = min(args.min_samples, args.samples)
    scene.cycles.max_bounces = 8
    scene.cycles.diffuse_bounces = 4
    scene.cycles.glossy_bounces = 4
    scene.cycles.transmission_bounces = 6
    scene.cycles.transparent_max_bounces = 8
    scene.cycles.use_light_tree = True
    scene.cycles.use_fast_gi = False
    scene.cycles.caustics_reflective = False
    scene.cycles.caustics_refractive = False
    scene.cycles.sample_clamp_indirect = 4
    scene.cycles.seed = 20261004
    scene.render.threads_mode = 'FIXED'
    scene.render.threads = 4
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGB'
    scene.render.image_settings.color_depth = '16'
    scene.render.image_settings.compression = 15
    scene.render.film_transparent = False
    scene.use_nodes = False
    scene['renderer'] = 'Blender Cycles / CPU / path tracing'
    scene['postprocessing'] = 'AgX display transform only; no denoising or image filtering'
    return scene


def make_macro(spec):
    scene = configure(macro_setup.macro(spec), 'CYCLES / Detail / ' + spec['name'])
    return scene


def surface(id, displaced=False):
    spec = next(s for s in atelier.SPECS if s['id'] == id)
    mat = atelier.material(spec)
    mat.asset_clear()
    node = next(n for n in mat.node_tree.nodes if n.type == 'GROUP')
    node.inputs['Displacement Mode'].default_value = displaced
    return mat, spec


def assign_mesh(name, verts, faces, mat, uv_fn):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    ob = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(ob)
    mesh.materials.append(mat)
    uv = mesh.uv_layers.new(name='UVMap')
    for poly in mesh.polygons:
        poly.use_smooth = True
        for loop_idx in poly.loop_indices:
            vertex_idx = mesh.loops[loop_idx].vertex_index
            uv.data[loop_idx].uv = uv_fn(vertex_idx, poly, loop_idx)
    return ob


def relief(ob, spec, direction='NORMAL'):
    tex = bpy.data.textures.new('Measured relief / ' + ob.name, type='IMAGE')
    tex.image = bpy.data.images.load(str(atelier.MAP_ROOT / spec['id'] / 'Height_Macro.png'), check_existing=True)
    tex.image.colorspace_settings.name = 'Non-Color'
    tex.use_interpolation = True
    tex.extension = 'EXTEND' if spec.get('tiling')=='finite-cut' else 'REPEAT'
    mod = ob.modifiers.new('Height_Macro / true geometric relief', 'DISPLACE')
    mod.texture = tex
    mod.texture_coords = 'UV'
    mod.uv_layer = ob.data.uv_layers.active.name
    mod.direction = direction
    mod.mid_level = .5
    mod.strength = spec['height_scale_m']


def patch(name, mat, spec, size, location, folds=False, center=(.5, .5), n=241):
    width, depth = size
    verts = []
    for j in range(n):
        y = (j / (n - 1) - .5) * depth
        for i in range(n):
            x = (i / (n - 1) - .5) * width
            z = 0.
            if folds:
                # Long soft folds and a raised turned-back corner; the texture
                # supplies actual yarn relief at its unchanged physical scale.
                z = .002 + .0045 * (1 + math.sin(2 * math.pi * (x / .067 + y / .24)))
                z += .003 * math.cos(y / depth * math.pi) ** 2
                z += .019 * math.exp(-((x - width * .40) / .031) ** 2) * ((y / depth + .5) ** 3)
                z += .004 * math.sin(y / .032 + x / .15) * (abs(x) / (width / 2)) ** 3
            verts.append((x, y, z))
    faces = [(j*n+i, j*n+i+1, (j+1)*n+i+1, (j+1)*n+i)
             for j in range(n - 1) for i in range(n - 1)]
    ob = assign_mesh(name, verts, faces, mat,
                     lambda k, p, l: (center[0] + verts[k][0] / spec['tile_m'],
                                      center[1] + verts[k][1] / spec.get('tile_y_m',spec['tile_m'])))
    ob.location = location
    relief(ob, spec)
    return ob


def lathe(name, profile, mat, spec, location, cylindrical=True, segments=256):
    displaced = next(n for n in mat.node_tree.nodes if n.type == 'GROUP').inputs['Displacement Mode'].default_value
    if displaced:
        dense = [profile[0]]
        for a, b in zip(profile, profile[1:]):
            count = max(1, math.ceil(math.dist(a, b) / .002))
            for k in range(1, count + 1):
                f = k / count
                dense.append((a[0] * (1-f) + b[0] * f, a[1] * (1-f) + b[1] * f))
        profile = dense
    verts = []
    for radius, z in profile:
        for i in range(segments + 1):
            a = 2 * math.pi * i / segments
            verts.append((radius * math.cos(a), radius * math.sin(a), z))
    row = segments + 1
    faces = [(j*row+i, j*row+i+1, (j+1)*row+i+1, (j+1)*row+i)
             for j in range(len(profile) - 1) for i in range(segments)]
    def uv_fn(k, p, l):
        if cylindrical:
            r, z = profile[k // row]
            return (k % row / segments * 2 * math.pi * max(r, .001) / spec['tile_m'] + .46,
                    z / spec.get('tile_y_m',spec['tile_m']) + .48)
        return (verts[k][0] / spec['tile_m'] + .35, verts[k][1] / spec.get('tile_y_m',spec['tile_m']) + .37)
    ob = assign_mesh(name, verts, faces, mat, uv_fn)
    ob.location = location
    sub = ob.modifiers.new('Continuous hand-turned profile', 'SUBSURF')
    sub.subdivision_type = 'CATMULL_CLARK'
    sub.levels = 2
    if displaced:
        relief(ob, spec)
    return ob


def ring(name, mat, spec, location, radius=.09, inner=.049, depth=.064):
    # A carved annulus with actual through-hole and rounded edges.
    profile = [(inner, -depth/2+.004), (inner+.002, -depth/2+.001),
               (radius-.003, -depth/2+.001), (radius, -depth/2+.004),
               (radius, depth/2-.004), (radius-.003, depth/2-.001),
               (inner+.002, depth/2-.001), (inner, depth/2-.004),
               (inner, -depth/2+.004)]
    ob = lathe(name, profile, mat, spec, location, cylindrical=False)
    ob.rotation_euler.x = math.pi/2
    return ob


def studio(name, width, height, target, camera_pos, focal=65):
    scene = configure(atelier.studio_scene(name, width, height, world_strength=.10), name)
    scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (.84, .88, .94, 1)
    scene.view_settings.look = 'AgX - Medium High Contrast'
    scene.view_settings.exposure = .1
    cam = atelier.camera(scene, camera_pos, target, .8)
    cam.data.type = 'PERSP'
    cam.data.lens = focal
    cam.data.clip_start = .001
    cam.data.clip_end = 100
    cam.data.dof.use_dof = False
    return scene


def scene_stone_timber():
    scene = studio('CYCLES / Still life / Stone and Timber', args.size, int(args.size*.75),
                   (0, .025, .17), (.54, -.93, .61), 56)
    plaster, ps = surface('09_lime_plaster', True)
    floor = patch('Lime plaster / continuous surface', plaster, ps, (3.2, 3.2), (0, 0, -.006), n=321)
    wall = patch('Lime plaster / grazing light wall', plaster, ps, (3.2, 1.8), (0, .51, .84), n=321)
    wall.rotation_euler.x = math.pi/2
    walnut, ws = surface('03_american_walnut')
    for i in range(3):
        ob = atelier.bevel_cube('Walnut / aged board %02d' % (i+1),
                                ((i-1)*.212, -.006, .043), (.21, .53, .086), walnut, .004)
        atelier.cube_uv_meters(ob, ws['tile_m'],ws.get('tile_y_m'))
        # Three physically sized cuts sample different radial parts of the
        # authored board volume, with slight finish-color differences.
        ob.data.materials[0]=walnut.copy()
        node=next(n for n in ob.data.materials[0].node_tree.nodes if n.type=='GROUP')
        node.inputs['Tint'].default_value=((.97 if i==0 else 1.025 if i==2 else 1),)*3+(1,)
        for loop in ob.data.uv_layers.active.data:
            loop.uv.x += (i-1)*.28
    marble, ms = surface('01_calacatta_oro')
    dish_profile = [(.001,.006), (.05,.006), (.10,.007), (.133,.010),
                    (.145,.019), (.1455,.022), (.143,.024), (.140,.023),
                    (.128,.016), (.099,.012), (.05,.011), (.001,.011)]
    lathe('Calacatta / shallow stone tray', dish_profile, marble, ms, (-.112,-.11,.088), False)
    trav, ts = surface('02_roman_travertine', True)
    ring('Travertine / carved circular study', trav, ts, (.105,.118,.217), .122, .067, .068)
    base_mat, _ = surface('02_roman_travertine')
    base = atelier.bevel_cube('Travertine / cut base', (.105,.115,.099), (.17,.12,.024), base_mat, .002)
    atelier.cube_uv_meters(base, ts['tile_m'])
    oak, os = surface('04_fumed_oak')
    profile=[(.001,0),(.050,0),(.052,.003),(.052,.072),(.049,.076),
             (.043,.078),(.041,.075),(.041,.070),(.001,.070)]
    lathe('Fumed oak / hand turned cup', profile, oak, os, (.19,-.135,.09))
    atelier.area(scene, 'Large window / side daylight', (-.60,-.28,.65), (0,.02,.16), 35, .52, .62, color=(1,.96,.9))
    atelier.area(scene, 'Softbox / open shadows', (.55,-.20,.75), (0,0,.15), 7, .45, .65, color=(.83,.9,1))
    atelier.area(scene, 'Long top reflection', (-.05,.37,.70), (0,0,.12), 10, .18, .50)
    scene['material_ids'] = ','.join([ps['id'], ws['id'], ms['id'], ts['id'], os['id']])
    return scene, '11_stone_and_timber'


def curve_line(name, points, mat, radius=.00025):
    data = bpy.data.curves.new(name, 'CURVE')
    data.dimensions = '3D'
    data.resolution_u = 12
    data.bevel_depth = radius
    data.bevel_resolution = 2
    spline = data.splines.new('POLY')
    spline.points.add(len(points)-1)
    for point, co in zip(spline.points, points):
        point.co = (*co, 1)
    ob = bpy.data.objects.new(name, data)
    bpy.context.scene.collection.objects.link(ob)
    ob.data.materials.append(mat)
    return ob


def scene_leather_linen():
    scene = studio('CYCLES / Still life / Leather and Linen', args.size, int(args.size*.75),
                   (0,.006,.026), (.27,-.34,.39), 63)
    oak, os = surface('04_fumed_oak')
    floor = atelier.bevel_cube('Fumed oak / workbench', (0,0,-.023), (.8,.8,.045), oak, .002)
    # This workbench top uses the finite oak cut once, with fitted UVs.
    atelier.cube_uv_meters(floor, .8, .8)
    linen, ls = surface('10_natural_linen', True)
    cloth = patch('Natural linen / folded woven swatch', linen, ls, (.185,.20),
                  (-.061,.019,.003), True, center=(.5,.5), n=769)
    cloth.rotation_euler.z = -.12
    solid = cloth.modifiers.new('Actual cloth edge thickness', 'SOLIDIFY')
    solid.thickness = .0006
    leather, les = surface('08_saddle_leather', True)
    swatch = patch('Saddle leather / cut hide sample', leather, les, (.148,.139),
                   (.043,-.018,.018), center=(.46,.46), n=641)
    swatch.rotation_euler.z = -.18
    thick = swatch.modifiers.new('Hide thickness', 'SOLIDIFY')
    thick.thickness = .0018
    # A raised leather fold gives curvature, contact shadow, and edge profile.
    verts=[]; faces=[]; nx=321; ny=97
    for j in range(ny):
        y=(j/(ny-1)-.5)*.116
        for i in range(nx):
            a=(i/(nx-1))*math.pi*1.10
            r=.021 + .0015*math.cos(y/.032)
            verts.append((r*math.cos(a), y, r*math.sin(a)))
    for j in range(ny-1):
        for i in range(nx-1):faces.append((j*nx+i,(j+1)*nx+i,(j+1)*nx+i+1,j*nx+i+1))
    roll=assign_mesh('Saddle leather / curled cut edge',verts,faces,leather,
                     lambda k,p,l:(.46+(k%nx)/(nx-1)*.021*math.pi*1.1/les['tile_m'],
                                   .46+verts[k][1]/les['tile_m']))
    roll.location=(.106,.056,.020);roll.rotation_euler.z=.12
    relief(roll,les)
    mod=roll.modifiers.new('Hide edge thickness','SOLIDIFY');mod.thickness=.0018
    thread = atelier.plain('Natural waxed flax / seam', (.43,.34,.23), .80)
    # Actual stitch geometry at 3.5mm spacing, each with a curved thread arc.
    for i in range(32):
        y=-.057+i*.0035
        x=.062
        local=[(x,y,.0009),(x+.0004,y+.0007,.0015),(x,y+.0025,.0009)]
        points=[]
        for px,py,pz in local:
            ca=math.cos(-.18);sa=math.sin(-.18)
            points.append((.043+ca*px-sa*py,-.018+sa*px+ca*py,.018+pz))
        curve_line('Hand stitched leather / %02d'%i,points,thread,.00023)
    brass, bs = surface('05_champagne_brass')
    for x,y in [(.087,-.062),(.104,.043)]:
        profile=[(.001,0),(.0045,0),(.005,.001),(.0048,.002),(.0035,.0028),(.001,.0029)]
        lathe('Aged brass / fastening stud',profile,brass,bs,(x,y,.020),False,64)
    scene.view_settings.exposure=-.30
    atelier.area(scene,'Raking daylight / fiber definition',(-.26,-.13,.12),(0,0,.017),2.5,.14,.25,color=(1,.965,.91))
    atelier.area(scene,'High soft fill',(.20,.15,.38),(0,0,.02),2.0,.35,.35,color=(.85,.92,1))
    atelier.area(scene,'Leather edge reflection',(.02,.28,.20),(.04,0,.025),1.1,.035,.25)
    scene['material_ids']=','.join([os['id'],ls['id'],les['id'],bs['id']])
    return scene,'12_leather_and_linen'


def scene_metal_ceramic():
    scene = studio('CYCLES / Still life / Metal and Ceramic', args.size, int(args.size*.75),
                   (0,-.005,.075), (.35,-.58,.37), 60)
    scene.view_settings.exposure = -.25
    plaster,ps=surface('09_lime_plaster',True)
    patch('Weathered lime / work surface',plaster,ps,(3.0,3.0),(0,0,-.009),n=321)
    ceramic,cs=surface('07_bone_porcelain')
    profile=[(.001,0),(.048,0),(.050,.003),(.050,.008),(.057,.025),(.065,.061),
             (.064,.092),(.057,.124),(.043,.150),(.040,.163),(.040,.170),
             (.038,.172),(.035,.171),(.035,.164),(.037,.151),(.050,.124),
             (.058,.092),(.059,.061),(.051,.029),(.044,.016),(.001,.012)]
    lathe('Crazed bone porcelain / open vessel',profile,ceramic,cs,(-.065,.067,0),True)
    steel, ss=surface('06_blackened_steel')
    plate=atelier.bevel_cube('Worn steel / solid sample plate',(.024,-.04,.006),(.255,.182,.012),steel,.003)
    plate.rotation_euler.z=-.22
    atelier.cube_uv_meters(plate,ss['tile_m'])
    exposed_edge=atelier.plain('Steel / rubbed machined edge',(.30,.33,.35),.34,1)
    plate.data.materials.append(exposed_edge)
    next(m for m in plate.modifiers if m.type=='BEVEL').material=1
    brass,bs=surface('05_champagne_brass')
    profile=[(.032,-.007),(.034,-.009),(.052,-.009),(.055,-.006),
             (.055,.007),(.052,.010),(.034,.010),(.032,.007),(.032,-.007)]
    lathe('Aged brass / machined annular fitting',profile,brass,bs,(.071,-.059,.023),False)
    small=[(.001,0),(.008,0),(.009,.002),(.009,.026),(.007,.028),(.001,.028)]
    lathe('Blackened steel / spacer',small,steel,ss,(.14,.036,0),True,96)
    atelier.area(scene,'Warm window / grazing glaze',(-.32,-.20,.45),(0,.01,.08),13,.30,.45,color=(1,.96,.90))
    atelier.area(scene,'Cool soft fill',(.35,-.02,.35),(0,0,.08),4.0,.30,.42,color=(.82,.91,1))
    atelier.area(scene,'Long strip / metal reflection',(-.01,.35,.42),(.04,0,.04),8.0,.045,.42)
    scene['material_ids']=','.join([ps['id'],cs['id'],ss['id'],bs['id']])
    return scene,'13_metal_and_ceramic'


def records():
    path = OUT / ('draft_manifest.json' if args.draft else 'render_manifest.json')
    return json.loads(path.read_text()) if path.exists() else {}


def release_pixels():
    for im in bpy.data.images:
        if im.source == 'FILE':
            im.buffers_free()


def render(scene, key, spec=None):
    folder = OUT / ('drafts' if args.draft else 'renders')
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / (key + '.png')
    scene.render.filepath = str(dest)
    print('PATH_TRACE_BEGIN ' + key, flush=True)
    start = time.monotonic()
    if not args.build_only:
        bpy.ops.render.render(write_still=True, scene=scene.name)
        data = records()
        data[key] = dict(
            scene=scene.name, file=str(dest.relative_to(OUT)),
            engine=scene.render.engine, blender=bpy.app.version_string,
            device='CPU', resolution=[scene.render.resolution_x, scene.render.resolution_y],
            maximum_samples=args.samples, minimum_samples=min(args.min_samples, args.samples),
            adaptive_threshold=args.threshold, denoising=False,
            postprocessing='Native AgX color transform only',
            color_depth='16-bit RGB PNG',
            view_transform=scene.view_settings.view_transform,
            look=scene.view_settings.look, exposure=scene.view_settings.exposure,
            material_ids=scene.get('material_ids', ''),
            elapsed_seconds=round(time.monotonic() - start, 2),
            provenance='Original deterministic procedural PBR; actual 3D path-traced rendering',
        )
        if spec:
            data[key].update(material_id=spec['id'], inspection_width_m=spec['preview_diameter_m'],
                             geometry_height='Height_Macro.png', residual_normal='Normal_Micro_OpenGL.png')
        meta = OUT / 'metadata'
        meta.mkdir(exist_ok=True)
        stage = '_draft' if args.draft else ''
        (meta / (key + stage + '.json')).write_text(json.dumps(data[key], indent=2) + '\n')
        (OUT / ('draft_manifest.json' if args.draft else 'render_manifest.json')).write_text(json.dumps(data, indent=2) + '\n')
    print('PATH_TRACE_READY ' + str(dest), flush=True)
    release_pixels()


def save_scenes(scenes, filename):
    if not scenes:
        return
    for image in bpy.data.images:
        if image.source == 'FILE':
            image.filepath = '//../materials/' + Path(image.filepath).parent.name + '/' + Path(image.filepath).name
    for sc in list(bpy.data.scenes):
        if not sc.objects and sc not in scenes:
            bpy.data.scenes.remove(sc)
    for sc in scenes:
        sc.render.filepath='//../path_traced/renders/' + Path(sc.render.filepath).name
    bpy.context.window.scene = scenes[0]
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                space = area.spaces.active
                space.region_3d.view_perspective = 'CAMERA'
                space.clip_start = .00001
                space.shading.type = 'MATERIAL'
                space.shading.use_scene_lights = True
                space.shading.use_scene_world = True
    notes = bpy.data.texts.new('READ ME / Path traced scenes')
    notes.write('CYBR MATERIAL 3 / Cycles path-traced studies\n'
                'Switch scenes to view the material studies. Every scene uses Cycles.\n'
                'Ten native 4K procedural surfaces; dimensions in meters; no source image inputs.\n'
                'Displaced objects pair Height_Macro geometry with the residual micro normal.\n'
                'Other objects use the full OpenGL normal.\n'
                'Rendering uses adaptive sampling without denoising or image smoothing.\n'
                'Texture paths are relative to the supplied suite materials directory.\n')
    bpy.ops.wm.save_as_mainfile(filepath=str(ROOT / 'blender' / filename), compress=True)


def main():
    OUT.mkdir(exist_ok=True)
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scenes = []
    if args.kind in ('macros', 'all'):
        for spec in atelier.SPECS:
            if args.only and spec['id'] not in args.only:
                continue
            scene = make_macro(spec)
            scenes.append(scene)
            render(scene, spec['id'] + '_detail', spec)
        if not args.draft:
            save_scenes(scenes, 'CYBR_Cycles_Details.blend')
    if args.kind in ('scenes', 'all'):
        scenes=[]
        builders=[('11_stone_and_timber',scene_stone_timber),
                  ('12_leather_and_linen',scene_leather_linen),
                  ('13_metal_and_ceramic',scene_metal_ceramic)]
        for key,builder in builders:
            if args.only and key not in args.only:continue
            scene,key=builder();scenes.append(scene)
            render(scene,key)
        if not args.draft:
            save_scenes(scenes,'CYBR_Cycles_Still_Lifes.blend')
    print('PATH_TRACING_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
