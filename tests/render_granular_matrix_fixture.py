"""Blender-only SI and equivalent-object-transform fixture, no scene dependency.
Run: blender -b -t 2 --python-exit-code 1 --python tests/render_granular_matrix_fixture.py -- /fresh/output
"""
from pathlib import Path
import sys,json,hashlib,math
import bpy,numpy as np
from mathutils import Vector
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'source'))
from granular_matrix import Nodes,load_profile,surface_footprint,aggregate_fields,attach_closure
O=Path(sys.argv[sys.argv.index('--')+1]).resolve();O.mkdir(parents=True,exist_ok=False)
p=load_profile(R/'source/granular_matrix_profiles.json')
bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene;s.unit_settings.system='METRIC';s.unit_settings.scale_length=1.;s.render.engine='CYCLES';s.cycles.device='CPU';s.cycles.samples=16;s.cycles.use_denoising=False;s.cycles.use_adaptive_sampling=False;s.render.resolution_x=640;s.render.resolution_y=480;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.render.image_settings.color_depth='16';s.render.image_settings.color_mode='RGB'
bpy.ops.mesh.primitive_plane_add(size=1);ob=bpy.context.object;ob.name='One metre optical substrate';m=bpy.data.materials.new('Neutral fixture');m.use_nodes=True;ob.data.materials.append(m);b=m.node_tree.nodes.get('Principled BSDF');b.inputs['Base Color'].default_value=(.27,.18,.095,1);b.inputs['Roughness'].default_value=.83
bpy.ops.object.camera_add(location=(0,-.63,.34));s.camera=bpy.context.object;s.camera.rotation_euler=(Vector((0,.015,0))-s.camera.location).to_track_quat('-Z','Y').to_euler();s.camera.data.lens=50;s.camera.data.clip_start=.001
bpy.ops.object.light_add(type='SUN',rotation=(math.radians(50),math.radians(-30),math.radians(-45)));bpy.context.object.data.energy=2.5;bpy.context.object.data.angle=math.radians(.53)
s.world=bpy.data.worlds.new('Neutral fixture ambient');s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.65,.76,1,1);s.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.18
builder=Nodes(m.node_tree);geom=builder.node('ShaderNodeNewGeometry','world metre coordinates');cam=builder.node('ShaderNodeCameraData','ray pixel footprint');footprint=surface_footprint(builder,geom.outputs['Normal'],geom.outputs['Incoming'],cam.outputs['View Distance'],s.camera.data.sensor_width/s.camera.data.lens/s.render.resolution_x);fields=aggregate_fields(builder,geom.outputs['Position'],footprint,1.,p);attach_closure(builder,b,fields['height_m'],builder.calc('ADD',1.,fields['pigment_delta']),fields['roughness_proxy'],p)
original=np.array([v.co[:]for v in ob.data.vertices]);cases=[]
for name,scale in [('identity',(1,1,1)),('nonuniform',(.25,4,2))]:
 ob.scale=scale;ob.data.vertices.foreach_set('co',(original/scale).astype('f4').ravel());ob.data.update();bpy.context.view_layer.update();assert np.array_equal(np.array([ob.matrix_world@v.co for v in ob.data.vertices]),original)
 s.render.filepath=str(O/(name+'.png'));bpy.ops.render.render(write_still=True);image=bpy.data.images.load(s.render.filepath,check_existing=False);pixels=np.array(image.pixels[:],dtype='f4');cases.append({'name':name,'matrix_world':list(map(list,ob.matrix_world)),'decoded_linear_RGBA_sha256':hashlib.sha256(pixels.tobytes()).hexdigest()})
 if name=='identity':reference=pixels.copy()
 else:assert np.array_equal(reference,pixels)
report={'renderer':bpy.app.version_string+' Cycles CPU','unit_scale':1.,'world_plane_size_m':[1.,1.],'cases':cases,'same_world_geometry':True,'decoded_RGB_exact':True,'source_sha256':hashlib.sha256((R/'source/granular_matrix.py').read_bytes()).hexdigest(),'profile_sha256':hashlib.sha256((R/'source/granular_matrix_profiles.json').read_bytes()).hexdigest(),'qualification':'Metric/transform diagnostic, not visual acceptance of a scene.'};(O/'VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
