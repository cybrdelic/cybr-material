"""Build the CYBR asset library and render honest material previews.

blender -b --python source/build_blender.py -- --render --gallery-only --samples 64
Draft rendering uses --draft --no-save. UVs in library assets default to one tile.
"""
import argparse
import json
import math
import sys
import uuid
from pathlib import Path
import bpy
from mathutils import Vector

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'source'))
argv = sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else []
parser = argparse.ArgumentParser()
parser.add_argument('--draft',action='store_true')
parser.add_argument('--render',action='store_true')
parser.add_argument('--no-save',action='store_true')
parser.add_argument('--only',nargs='*',default=[])
parser.add_argument('--samples',type=int,default=64)
parser.add_argument('--size',type=int,default=1200)
parser.add_argument('--gallery-only',action='store_true')
args = parser.parse_args(argv)
MAP_ROOT = ROOT/('draft_materials' if args.draft else 'materials')
SPECS = json.loads((MAP_ROOT/'manifest.json').read_text())['materials']
if args.draft:
    SPECS=[s for s in SPECS if (MAP_ROOT/s['id']/'BaseColor.png').exists()]
CATALOGS={family:str(uuid.uuid5(uuid.NAMESPACE_DNS,'aurel.material.atelier.'+family))
          for family in sorted({s['family'] for s in SPECS})}


def socket(group, name, kind, direction='INPUT', default=None, lo=None, hi=None):
    s=group.interface.new_socket(name=name,in_out=direction,socket_type=kind)
    if default is not None:
        s.default_value=default
    if lo is not None:
        s.min_value=lo
    if hi is not None:
        s.max_value=hi
    return s


def material(spec):
    mat=bpy.data.materials.new('CYBR · '+spec['name'])
    mat.use_nodes=True
    mat.use_fake_user=True
    mat.diffuse_color=(.3,.3,.3,1)
    mat['aurel_id']=spec['id']
    mat['tile_m']=spec['tile_m']
    mat['tile_y_m']=spec.get('tile_y_m',spec['tile_m'])
    mat['tiling']=spec.get('tiling','periodic')
    mat['height_scale_m']=spec['height_scale_m']
    mat['normal_convention']='OpenGL +Y'
    mat['surface']=spec['description']
    mat.asset_mark()
    mat.asset_data.catalog_id=CATALOGS[spec['family']]
    mat.asset_data.description=spec['description']+' Native '+str(spec.get('resolution',4096))+' px; '+spec.get('tiling','periodic')+'.'
    for tag in ('CYBR','PBR',str(spec.get('resolution',4096))+'px',spec.get('tiling','periodic'),spec['family']):
        mat.asset_data.tags.new(tag)
    group=bpy.data.node_groups.new('CYBR / '+spec['name'],'ShaderNodeTree')
    socket(group,'Surface','NodeSocketShader','OUTPUT')
    socket(group,'Height','NodeSocketFloat','OUTPUT')
    socket(group,'Base Color','NodeSocketColor','OUTPUT')
    socket(group,'Macro Height','NodeSocketFloat','OUTPUT')
    socket(group,'Wear Mask','NodeSocketFloat','OUTPUT')
    socket(group,'Tint','NodeSocketColor',default=(1,1,1,1))
    socket(group,'Roughness Offset','NodeSocketFloat',default=0,lo=-1,hi=1)
    socket(group,'Normal Strength','NodeSocketFloat',default=1,lo=0,hi=2)
    socket(group,'UV Scale','NodeSocketVector',default=(1,1,1),lo=.001,hi=1000)
    socket(group,'Displacement Mode','NodeSocketBool',default=False)
    socket(group,'Wear Tint','NodeSocketColor',default=(1,1,1,1))
    if spec['family']=='Wood':socket(group,'Timber Axis Rotation','NodeSocketVector',default=(0,0,0))
    nodes=group.nodes; links=group.links
    gi=nodes.new('NodeGroupInput');gi.location=(-1100,550)
    go=nodes.new('NodeGroupOutput');go.location=(720,300)
    uv=nodes.new('ShaderNodeTexCoord');uv.location=(-1100,100)
    uv.label='Use UV unwrap; one UV tile = '+str(spec['tile_m'])+' m'
    scale=nodes.new('ShaderNodeVectorMath');scale.operation='MULTIPLY';scale.location=(-900,100)
    links.new(uv.outputs['UV'],scale.inputs[0]);links.new(gi.outputs['UV Scale'],scale.inputs[1])
    textures={}
    for i,filename in enumerate(['BaseColor.png','Roughness.png','Metallic.png','Normal_OpenGL.png','Height.png','AO.png','WearMask.png','Height_Macro.png','Normal_Micro_OpenGL.png','Opacity.png']):
        node=nodes.new('ShaderNodeTexImage')
        node.image=bpy.data.images.load(str(MAP_ROOT/spec['id']/filename),check_existing=True)
        node.image.colorspace_settings.name='sRGB' if filename=='BaseColor.png' else 'Non-Color'
        node.label=filename.removesuffix('.png')
        node.interpolation='Linear';node.extension='EXTEND' if spec.get('tiling')=='finite-cut' else 'REPEAT'
        node.location=(-640,610-i*260)
        links.new(scale.outputs[0],node.inputs['Vector'])
        textures[filename]=node
    tex=textures
    mul=nodes.new('ShaderNodeMixRGB');mul.blend_type='MULTIPLY';mul.inputs[0].default_value=1
    mul.location=(-280,550);mul.label='Base color tint'
    links.new(tex['BaseColor.png'].outputs['Color'],mul.inputs[1])
    links.new(gi.outputs['Tint'],mul.inputs[2])
    wear_tint=nodes.new('ShaderNodeMixRGB');wear_tint.blend_type='MIX';wear_tint.location=(-280,780)
    wear_tint.label='Tint only the worn surface features';wear_tint.inputs[1].default_value=(1,1,1,1)
    links.new(tex['WearMask.png'].outputs['Color'],wear_tint.inputs[0])
    links.new(gi.outputs['Wear Tint'],wear_tint.inputs[2])
    worn_color=nodes.new('ShaderNodeMixRGB');worn_color.blend_type='MULTIPLY';worn_color.inputs[0].default_value=1
    worn_color.location=(-40,740)
    links.new(mul.outputs[0],worn_color.inputs[1]);links.new(wear_tint.outputs[0],worn_color.inputs[2])
    rough=nodes.new('ShaderNodeMath');rough.operation='ADD';rough.use_clamp=True
    rough.location=(-250,320);rough.label='Calibrated roughness + offset'
    links.new(tex['Roughness.png'].outputs['Color'],rough.inputs[0])
    links.new(gi.outputs['Roughness Offset'],rough.inputs[1])
    normal=nodes.new('ShaderNodeNormalMap');normal.location=(-260,-140)
    normal.label='OpenGL / tangent space / calibrated relief'
    normals=nodes.new('ShaderNodeMixRGB');normals.blend_type='MIX';normals.location=(-460,-180)
    normals.label='Full normal / residual normal when geometry is displaced'
    links.new(gi.outputs['Displacement Mode'],normals.inputs[0])
    links.new(tex['Normal_OpenGL.png'].outputs['Color'],normals.inputs[1])
    links.new(tex['Normal_Micro_OpenGL.png'].outputs['Color'],normals.inputs[2])
    links.new(normals.outputs[0],normal.inputs['Color'])
    links.new(gi.outputs['Normal Strength'],normal.inputs['Strength'])
    principled=nodes.new('ShaderNodeBsdfPrincipled');principled.location=(100,500)
    principled.inputs['IOR'].default_value=spec['ior']
    for input_name,key in [('Coat Weight','coat'),('Coat Roughness','coat_roughness'),
                           ('Anisotropic','anisotropy'),('Sheen Weight','sheen'),
                           ('Sheen Roughness','sheen_roughness')]:
        if key in spec:
            principled.inputs[input_name].default_value=spec[key]
    principled.inputs['Coat IOR'].default_value=spec['ior']
    tangent=nodes.new('ShaderNodeTangent');tangent.direction_type='UV_MAP';tangent.location=(-250,-380)
    links.new(tangent.outputs['Tangent'],principled.inputs['Tangent'])
    links.new(worn_color.outputs[0],principled.inputs['Base Color'])
    links.new(rough.outputs[0],principled.inputs['Roughness'])
    links.new(tex['Metallic.png'].outputs['Color'],principled.inputs['Metallic'])
    links.new(tex['Opacity.png'].outputs['Color'],principled.inputs['Alpha'])
    links.new(normal.outputs[0],principled.inputs['Normal'])
    links.new(normal.outputs[0],principled.inputs['Coat Normal'])
    if spec.get('tiling')=='finite-cut':
        # Finite procedural wood cuts may be fitted to a prop.
        # Height-derived bump computes slopes for that actual mapping;
        # displaced meshes use only the remaining height band.
        used_macro=nodes.new('ShaderNodeMath');used_macro.operation='MULTIPLY'
        links.new(tex['Height_Macro.png'].outputs['Color'],used_macro.inputs[0])
        links.new(gi.outputs['Displacement Mode'],used_macro.inputs[1])
        remainder=nodes.new('ShaderNodeMath');remainder.operation='SUBTRACT'
        links.new(tex['Height.png'].outputs['Color'],remainder.inputs[0])
        links.new(used_macro.outputs[0],remainder.inputs[1])
        bump=nodes.new('ShaderNodeBump');bump.label='Finite wood cut / physical height derivative'
        bump.inputs['Distance'].default_value=spec['height_scale_m']
        links.new(gi.outputs['Normal Strength'],bump.inputs['Strength'])
        links.new(remainder.outputs[0],bump.inputs['Height'])
        links.new(bump.outputs['Normal'],principled.inputs['Normal'])
        links.new(bump.outputs['Normal'],principled.inputs['Coat Normal'])
    if spec['family']=='Wood':
        from timber_volume_shader import install
        install(group,spec)
    links.new(principled.outputs[0],go.inputs['Surface'])
    links.new(worn_color.outputs[0],go.inputs['Base Color'])
    links.new(tex['Height.png'].outputs['Color'],go.inputs['Height'])
    links.new(tex['Height_Macro.png'].outputs['Color'],go.inputs['Macro Height'])
    links.new(tex['WearMask.png'].outputs['Color'],go.inputs['Wear Mask'])
    tex['Height.png'].label='Height · linear 16-bit · optional true displacement'
    tex['AO.png'].label='AO · cavity export · intentionally not baked into color'
    mat.node_tree.nodes.clear()
    shader=mat.node_tree.nodes.new('ShaderNodeGroup');shader.node_tree=group
    shader.location=(-340,180);shader.width=320
    shader.label=spec['name']+' / '+str(spec['tile_m'])+' m tile'
    output=mat.node_tree.nodes.new('ShaderNodeOutputMaterial');output.location=(160,180)
    mat.node_tree.links.new(shader.outputs['Surface'],output.inputs['Surface'])
    aov=mat.node_tree.nodes.new('ShaderNodeOutputAOV');aov.name='CYBR_Albedo'
    aov.aov_name='CYBR_Albedo'
    aov.location=(160,-380);aov.label='Optional first-hit base color AOV'
    mat.node_tree.links.new(shader.outputs['Base Color'],aov.inputs['Color'])
    disp=mat.node_tree.nodes.new('ShaderNodeDisplacement');disp.location=(-80,-150)
    disp.inputs['Scale'].default_value=spec['height_scale_m']
    disp.inputs['Midlevel'].default_value=.5
    disp.label='Optional macro displacement: enable Displacement Mode and subdivide'
    mat.node_tree.links.new(shader.outputs['Macro Height'],disp.inputs['Height'])
    return mat


def plain(name,rgb,rough=.5,metal=0):
    mat=bpy.data.materials.new(name);mat.use_nodes=True
    p=mat.node_tree.nodes.get('Principled BSDF')
    p.inputs['Base Color'].default_value=(*rgb,1)
    p.inputs['Roughness'].default_value=rough;p.inputs['Metallic'].default_value=metal
    return mat


def look_at(obj,target):
    obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()


def studio_scene(name, width, height, world_strength=.24):
    scene=bpy.data.scenes.new(name)
    bpy.context.window.scene=scene
    scene.render.engine='BLENDER_EEVEE_NEXT'
    scene.eevee.taa_render_samples=args.samples
    scene.eevee.use_raytracing=False
    if hasattr(scene.eevee,'use_shadow_jitter'):scene.eevee.use_shadow_jitter=False
    scene.cycles.device='CPU';scene.cycles.samples=args.samples
    scene.cycles.use_adaptive_sampling=True
    scene.cycles.adaptive_threshold=.006
    # Some Linux distribution builds omit OIDN. Remain usable in those builds.
    try:
        scene.cycles.denoiser='OPENIMAGEDENOISE'
        scene.cycles.use_denoising=True
    except TypeError:
        scene.cycles.use_denoising=False
        scene.cycles.adaptive_threshold=.004
    scene.cycles.max_bounces=8
    scene.cycles.diffuse_bounces=3
    scene.cycles.glossy_bounces=5
    scene.cycles.sample_clamp_indirect=4
    scene.render.threads_mode='FIXED';scene.render.threads=5
    scene.render.resolution_x=width;scene.render.resolution_y=height
    scene.render.resolution_percentage=100
    scene.render.image_settings.file_format='PNG'
    scene.render.image_settings.color_mode='RGB'
    scene.render.image_settings.color_depth='8'
    scene.render.film_transparent=False
    scene.view_settings.view_transform='AgX'
    scene.view_settings.look='AgX - Medium High Contrast'
    scene.view_settings.exposure=-.1
    world=bpy.data.worlds.new(name+' / environment');world.use_nodes=True
    world.node_tree.nodes['Background'].inputs['Color'].default_value=(.78,.84,1,1)
    world.node_tree.nodes['Background'].inputs['Strength'].default_value=world_strength
    scene.world=world
    scene.unit_settings.system='METRIC'
    scene.unit_settings.length_unit='METERS'
    return scene


def area(scene,name,location,target,energy,size,size_y=None,color=(1,1,1)):
    data=bpy.data.lights.new(name,'AREA');data.energy=energy;data.color=color
    if size_y:
        data.shape='RECTANGLE';data.size=size;data.size_y=size_y
    else:
        data.shape='DISK';data.size=size
    ob=bpy.data.objects.new(name,data);scene.collection.objects.link(ob)
    ob.location=location;look_at(ob,target)
    return ob


def plane(name,location,size,mat):
    bpy.ops.mesh.primitive_plane_add(size=size,location=location)
    ob=bpy.context.object;ob.name=name;ob.data.materials.append(mat)
    return ob


def bevel_cube(name,location,dimensions,mat,bevel=.012):
    # Direct mesh construction avoids repeatedly evaluating every existing
    # object through transform_apply while assembling thousands of books.
    x,y,z=[v*.5 for v in dimensions]
    verts=[(-x,-y,-z),(x,-y,-z),(x,y,-z),(-x,y,-z),(-x,-y,z),(x,-y,z),(x,y,z),(-x,y,z)]
    faces=[(3,2,1,0),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(verts,[],faces);mesh.update()
    ob=bpy.data.objects.new(name,mesh);bpy.context.scene.collection.objects.link(ob);ob.location=location
    mesh.uv_layers.new(name='UVMap')
    ob.data.materials.append(mat)
    modifier=ob.modifiers.new('Soft machined edges','BEVEL');modifier.width=bevel;modifier.segments=5
    normal=ob.modifiers.new('Weighted face normals','WEIGHTED_NORMAL');normal.keep_sharp=True
    return ob


def cube_uv_meters(ob,tile,tile_y=None):
    uv=ob.data.uv_layers.active or ob.data.uv_layers.new()
    for poly in ob.data.polygons:
        axis=max(range(3),key=lambda i:abs(poly.normal[i]))
        axes=(1,2) if axis==0 else ((0,2) if axis==1 else (0,1))
        for idx in poly.loop_indices:
            co=ob.data.vertices[ob.data.loops[idx].vertex_index].co
            uv.data[idx].uv=(co[axes[0]]/tile+.5,co[axes[1]]/(tile_y or tile)+.5)


def ball(name,location,radius,mat,tile):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=192,ring_count=128,radius=radius,location=location)
    ob=bpy.context.object;ob.name=name
    for p in ob.data.polygons:p.use_smooth=True
    ob.data.materials.append(mat)
    uv=ob.data.uv_layers.active
    base=[v.uv.copy() for v in uv.data]
    for v,b in zip(uv.data,base):
        v.uv=(b.x*2*math.pi*radius/tile,b.y*math.pi*radius/tile)
    return ob,base


def pedestal(name,location,radius,mat):
    bpy.ops.mesh.primitive_cylinder_add(vertices=128,radius=radius,depth=.04,location=location)
    ob=bpy.context.object;ob.name=name;ob.data.materials.append(mat)
    for p in ob.data.polygons:p.use_smooth=abs(p.normal.z)<.5
    mod=ob.modifiers.new('Rounded rim','BEVEL');mod.width=.007;mod.segments=4
    ob.modifiers.new('Weighted normals','WEIGHTED_NORMAL')
    return ob


def camera(scene,location,target,ortho):
    data=bpy.data.cameras.new(scene.name+' / camera');data.type='ORTHO';data.ortho_scale=ortho
    ob=bpy.data.objects.new(scene.name+' / camera',data);scene.collection.objects.link(ob)
    ob.location=location;look_at(ob,target);scene.camera=ob
    return ob


def render_clean(scene,dest):
    scene.render.filepath=str(dest)
    bpy.ops.render.render(write_still=True)


def main():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    mats={s['id']:material(s) for s in SPECS}
    floor=plain('Studio / warm limestone',(.24,.225,.196),.57)
    stand=plain('Studio / graphite plinth',(.026,.031,.033),.3)
    scene=studio_scene('CYBR / Inspection Studio',args.size,args.size)
    plane('Infinite matte cyclorama',(0,0,-.005),200,floor)
    area(scene,'Key / broad silk',(-.55,-.8,1.5),(0,0,.25),115,1.25,color=(1,.94,.84))
    area(scene,'Fill / daylight',(1.1,-.1,.9),(0,0,.25),55,.8,color=(.8,.89,1))
    area(scene,'Rim / vertical strip',(-.75,.6,.95),(0,0,.3),85,.25,1.3,color=(1,.96,.88))
    camera(scene,(1.15,-2.25,1.04),(.035,0,.31),.985)
    pedestal('Graphite sample plinth',(-.135,-.045,.02),.267,stand)
    sphere,base_uv=ball('Curved material sample',(-.135,-.045,.28),.24,mats[SPECS[0]['id']],SPECS[0]['tile_m'])
    slab=bevel_cube('Planar material sample',(.23,.13,.315),(.35,.075,.59),mats[SPECS[0]['id']],.012)
    slab.rotation_euler.z=-.17
    cube_uv_meters(slab,SPECS[0]['tile_m'])
    studio_bases={ob.name:dict(location=ob.location.copy(),scale=ob.scale.copy(),
                  energy=ob.data.energy if ob.type=='LIGHT' else None,
                  size=ob.data.size if ob.type=='LIGHT' else None,
                  size_y=ob.data.size_y if ob.type=='LIGHT' else None)
                  for ob in scene.objects}
    gallery=studio_scene('CYBR / Material Gallery',1800,1260,world_strength=.2)
    plane('Gallery / floor',(0,0,-.005),200,floor)
    area(gallery,'Gallery / silk key',(-2.0,-2.3,5.2),(0,0,.2),1100,4.0,color=(1,.93,.83))
    area(gallery,'Gallery / cool fill',(3.2,.4,3.5),(0,0,.3),650,3.0,color=(.8,.9,1))
    area(gallery,'Gallery / back strip',(-2.4,3,2.7),(0,0,.3),950,1.0,4.0)
    camera(gallery,(0,-4.8,6.3),(0,.1,.15),5.0)
    for i,s in enumerate(SPECS):
        pos=((i%5-2)*.89,(.82 if i<5 else -.55),.3)
        if s.get('tiling')=='finite-cut':
            ob=bevel_cube(s['name'],pos,(.21,.06,.53),mats[s['id']],.005)
            cube_uv_meters(ob,s['tile_m'],s['tile_y_m'])
        else:
            ob,_=ball(s['name'],pos,.255,mats[s['id']],s['tile_m'])
        pedestal(s['name']+' / plinth',(pos[0],pos[1],.023),.294,stand)
        ob.location.z=.295
        sample=bevel_cube(s['name']+' / flat swatch',(pos[0],pos[1]-.41,.023),(.54,.20,.035),mats[s['id']],.006)
        cube_uv_meters(sample,s['tile_m'],s.get('tile_y_m'))
    if not args.no_save:
        bpy.context.window.scene=gallery
        for image in bpy.data.images:
            if image.source=='FILE':
                image.filepath='//../materials/'+Path(image.filepath).parent.name+'/'+Path(image.filepath).name
        bpy.ops.wm.save_as_mainfile(filepath=str(ROOT/'blender'/'CYBR_Material_Atelier.blend'))
        print('CYBR_LIBRARY_READY',flush=True)
    if args.render:
        if not args.only:
            bpy.context.window.scene=gallery
            render_clean(gallery,ROOT/'previews'/('gallery-draft.png' if args.draft else 'gallery.png'))
        if not args.gallery_only:
            bpy.context.window.scene=scene
            for s in SPECS:
                if args.only and s['id'] not in args.only:continue
                print('CYBR_RENDER '+s['id'],flush=True)
                factor=s.get('preview_diameter_m',.48)/.48
                for ob in scene.objects:
                    base=studio_bases[ob.name]
                    ob.location=base['location']*factor
                    ob.scale=base['scale']*(factor if ob.type not in ('LIGHT','CAMERA') else 1)
                    if ob.type=='LIGHT':
                        ob.data.energy=base['energy']*factor*factor
                        ob.data.size=base['size']*factor
                        ob.data.size_y=base['size_y']*factor
                scene.camera.data.ortho_scale=.985*factor
                sphere.data.materials[0]=mats[s['id']]
                slab.data.materials[0]=mats[s['id']]
                for v,b in zip(sphere.data.uv_layers.active.data,base_uv):
                    v.uv=(b.x*2*math.pi*.24*factor/s['tile_m'],b.y*math.pi*.24*factor/s['tile_m'])
                cube_uv_meters(slab,s['tile_m']/factor)
                suffix='-draft' if args.draft else ''
                render_clean(scene,ROOT/'previews'/('studio_'+s['id']+suffix+'.png'))
    if not args.no_save:
        # Portable external paths avoid duplicating all 4K maps in the .blend.
        bpy.context.window.scene=gallery
        for image in bpy.data.images:
            if image.source=='FILE':
                image.filepath='//../materials/'+Path(image.filepath).parent.name+'/'+Path(image.filepath).name
        notes=bpy.data.texts.new('READ ME / CYBR')
        notes.write('CYBR / Material Atelier\nTen native 4K procedural materials; the wood maps are finite board cuts.\n'
                    'Every material is marked as an Asset. Add this blender folder as an Asset Library.\n'
                    'Keep ../materials beside this folder. Image textures use portable relative paths.\n'
                    'Shader controls: Tint, Wear Tint, Roughness Offset, Normal Strength, UV Scale.\n'
                    'Default UV Scale = 1, one UV tile = physical tile_m in material.json.\n'
                    'Gallery samples are mapped at physical scale.\n'
                    'Full normal is the default. Macro height + residual micro normal are also supplied.\n'
                    'For true displacement: enable Displacement Mode, connect the prepared Displacement,\n'
                    'add subdivision and use Cycles Displacement Only. Consult docs/IMPORT_GUIDE.md.\n'
                    'Both OpenGL and DirectX normals and ORM channel packs are included.\n')
        filepath=ROOT/'blender'/'CYBR_Material_Atelier.blend'
        # Remove local capture paths from the portable deliverable. Standard
        # Blender builds support the default native OIDN setting on reopen.
        for s in (scene,gallery):
            s.use_nodes=False
            if s.node_tree:
                s.node_tree.nodes.clear()
            s.view_layers[0].cycles.denoising_store_passes=False
            s.view_layers[0].use_pass_normal=False
            for aov in list(s.view_layers[0].aovs):
                s.view_layers[0].aovs.remove(aov)
            s.cycles.use_denoising=False
            s.render.engine='CYCLES'
            s.cycles.samples=768
            s.cycles.adaptive_min_samples=128
            s.render.image_settings.color_depth='16'
        catalogs=['# CYBR / Material Atelier asset catalogs','VERSION 1']
        catalogs += [f'{ident}:CYBR/{family}:{family}' for family,ident in CATALOGS.items()]
        (ROOT/'blender'/'blender_assets.cats.txt').write_text('\n'.join(catalogs)+'\n')
        for s in SPECS:
            preview=ROOT/'path_traced/renders'/(s['id']+'_detail.png')
            if preview.exists():
                with bpy.context.temp_override(id=mats[s['id']]):
                    bpy.ops.ed.lib_id_load_custom_preview(filepath=str(preview))
        for screen in bpy.data.screens:
            for view_area in screen.areas:
                if view_area.type=='VIEW_3D':
                    space=view_area.spaces.active
                    space.shading.type='MATERIAL'
                    space.region_3d.view_perspective='CAMERA'
                    space.clip_start=.001
                    space.clip_end=1000
        bpy.ops.wm.save_as_mainfile(filepath=str(filepath),compress=True)
        print('CYBR_SAVED '+str(filepath),flush=True)


if __name__=='__main__':main()
