"""Actual Blender/Cycles CPU review of metric native ground maps."""
import bpy,sys,json,math,time,hashlib
from pathlib import Path
import numpy as np
from mathutils import Vector
args=sys.argv[sys.argv.index('--')+1:];root=Path(args[0]).resolve();out=Path(args[1]).resolve();out.mkdir(parents=True,exist_ok=True)
view=args[2] if len(args)>2 else 'macro';samples=int(args[3]) if len(args)>3 else 24;res=int(args[4]) if len(args)>4 else 640
meta=json.loads((root/'material.json').read_text());tile=meta['tile_m'];height=meta['height_scale_m'];g=np.load(root/'SoilGeometryHeight.npy')[::-1];ng=g.shape[0]
bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene;s.render.engine='CYCLES';s.cycles.device='CPU';s.render.threads_mode='FIXED';s.render.threads=2;s.cycles.samples=samples;s.cycles.seed=1810;s.cycles.use_adaptive_sampling=False;s.cycles.use_denoising=False;s.cycles.use_light_tree=False
s.cycles.max_bounces=12;s.cycles.diffuse_bounces=4;s.cycles.glossy_bounces=4;s.cycles.transmission_bounces=12;s.cycles.transparent_max_bounces=12;s.cycles.volume_bounces=0;s.cycles.sample_clamp_direct=0;s.cycles.sample_clamp_indirect=10;s.cycles.blur_glossy=1
s.render.resolution_x=s.render.resolution_y=res;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.render.image_settings.color_depth='16';s.render.image_settings.color_mode='RGB';s.render.dither_intensity=0
s.view_settings.view_transform='AgX';s.view_settings.look='AgX - Medium High Contrast';s.view_settings.exposure=0;s.view_settings.gamma=1
s.world=bpy.data.worlds.new('Neutral outdoor fill');s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs[0].default_value=(.18,.21,.25,1);s.world.node_tree.nodes['Background'].inputs[1].default_value=.04
m=bpy.data.materials.new(meta['name']);m.use_nodes=True;nodes=m.node_tree.nodes;links=m.node_tree.links;p=nodes.get('Principled BSDF');p.inputs['IOR'].default_value=meta['ior'];p.inputs['Specular IOR Level'].default_value=.35
uv=nodes.new('ShaderNodeTexCoord')
for name,space,sock in [('BaseColor','sRGB','Base Color'),('Roughness','Non-Color','Roughness'),('Metallic','Non-Color','Metallic')]:
 tex=nodes.new('ShaderNodeTexImage');tex.name=name;tex.image=bpy.data.images.load(str(root/(name+'.png')));tex.image.colorspace_settings.name=space;tex.image.pack();tex.extension='REPEAT';links.new(uv.outputs['UV'],tex.inputs[0]);links.new(tex.outputs['Color'],p.inputs[sock])
tex=nodes.new('ShaderNodeTexImage');tex.name='Full metric normal';tex.image=bpy.data.images.load(str(root/'Soil_Normal_Object_RGB16.png'));tex.image.colorspace_settings.name='Non-Color';tex.image.pack();links.new(uv.outputs['UV'],tex.inputs[0]);nm=nodes.new('ShaderNodeNormalMap');nm.space='OBJECT';nm.inputs['Strength'].default_value=1.;links.new(tex.outputs['Color'],nm.inputs['Color']);links.new(nm.outputs[0],p.inputs['Normal'])
a=np.linspace(-tile/2,tile/2,ng,dtype='f4');y,x=np.meshgrid(a,a,indexing='ij');verts=np.stack((x,y,(g-.5)*height),axis=-1).reshape(-1,3)
ii=np.arange((ng-1)*(ng-1));j=ii//(ng-1);i=ii%(ng-1);c=j*ng+i;faces=np.stack((c,c+1,c+ng+1,c+ng),axis=-1)
gravel=np.load(root/'GravelGeometry.npz');gv=gravel['vertices'];gf=gravel['triangles'];offset=len(verts);mixed=faces.tolist()+((gf+offset).tolist());verts=np.concatenate((verts,gv));me=bpy.data.meshes.new('Metric soil and explicit convex gravel');me.from_pydata(verts,[],mixed);me.update();o=bpy.data.objects.new('Victorville ground / displaced surface',me);s.collection.objects.link(o);o.data.materials.append(m)
rock=m.copy();rock.name='Victorville / actual convex gravel';rn=rock.node_tree.nodes;rl=rock.node_tree.links;rp=rn.get('Principled BSDF')
for link in list(rl):
 if link.to_node==rp and link.to_socket==rp.inputs['Normal']:rl.remove(link)
micro=rn.new('ShaderNodeTexImage');micro.name='Coupled mica quartz microrelief';micro.image=bpy.data.images.load(str(root/'GravelMicroHeight.png'));micro.image.colorspace_settings.name='Non-Color';micro.image.pack();rl.new(rn.get('Texture Coordinate').outputs['UV'],micro.inputs['Vector']);bump=rn.new('ShaderNodeBump');bump.inputs['Distance'].default_value=.00012;bump.inputs['Strength'].default_value=1.;rl.new(micro.outputs['Color'],bump.inputs['Height']);rl.new(bump.outputs['Normal'],rp.inputs['Normal']);o.data.materials.append(rock)
material_ids=np.zeros(len(me.polygons),'i4');material_ids[len(faces):]=1;me.polygons.foreach_set('material_index',material_ids)

uvlayer=me.uv_layers.new(name='Physical ground UV');vi=np.empty(len(me.loops),dtype='i4');me.loops.foreach_get('vertex_index',vi);coords=np.stack((verts[vi,0]/tile+.5,verts[vi,1]/tile+.5),axis=-1).astype('f4');uvlayer.data.foreach_set('uv',coords.ravel());me.polygons.foreach_set('use_smooth',np.concatenate((np.ones(len(faces),'i4'),np.zeros(len(gf),'i4'))));me.update()
# Irregular physical edge weathering: bounded bevel width, gravel edges only.
edge_verts=np.empty((len(me.edges),2),'i4');me.edges.foreach_get('vertices',edge_verts.ravel());mask=(edge_verts[:,0]>=offset)&(edge_verts[:,1]>=offset)
weights=np.zeros(len(me.edges),'f4');rng=np.random.default_rng(meta['seed']+720);weights[mask]=rng.uniform(.20,1.0,int(mask.sum()))
attr=me.attributes.new(name='bevel_weight_edge',type='FLOAT',domain='EDGE');attr.data.foreach_set('value',weights)
bev=o.modifiers.new('Physical gravel edge weathering / 0.044–0.22 mm','BEVEL');bev.width=.00022;bev.segments=2;bev.limit_method='WEIGHT';bev.use_clamp_overlap=True
# Sun and a sky fill give genuine inter-grain shadowing; no image environment map.
# User-selected low-grazing asphalt disk rig, scaled in SI units for each crop.
factor=(.46 if view=='macro' else .13)/.135
focus=Vector((0,0,.0025) if view=='macro' else (.014,-.031,.0025))
ld=bpy.data.lights.new('Grazing disk / asphalt reference','AREA');ld.shape='DISK';ld.size=.055*factor;ld.energy=3.24*factor*factor;ld.color=(1,.92,.82)
so=bpy.data.objects.new(ld.name,ld);s.collection.objects.link(so);so.location=focus+Vector((-.22*factor,-.16*factor,.027*factor));so.rotation_euler=(focus-so.location).to_track_quat('-Z','Y').to_euler()
cam=bpy.data.cameras.new('Ground review camera');co=bpy.data.objects.new('Ground review camera',cam);s.collection.objects.link(co);s.camera=co;cam.type='ORTHO'
if view=='macro':target=Vector((0,0,.001));co.location=(0,-.45,.75);cam.ortho_scale=.46
elif view=='detail':target=Vector((.014,-.031,.001));co.location=target+Vector((.080,-.120,.15));cam.ortho_scale=.13
elif view=='tile':target=Vector((0,0,0));co.location=(0,0,1);cam.ortho_scale=.6
else:raise ValueError(view)
co.rotation_euler=(target-co.location).to_track_quat('-Z','Y').to_euler();cam.lens=50;s.camera.data.clip_start=.001
# Pack both review cameras into every portable scene.
co.name='Review / '+view
alt=bpy.data.objects.new('Review / '+('detail' if view=='macro' else 'macro'),cam.copy());s.collection.objects.link(alt)
if view=='macro':at=Vector((.014,-.031,.001));alt.location=at+Vector((.080,-.120,.15));alt.data.ortho_scale=.13
else:at=Vector((0,0,.001));alt.location=(0,-.45,.75);alt.data.ortho_scale=.46
alt.rotation_euler=(at-alt.location).to_track_quat('-Z','Y').to_euler()
s.unit_settings.system='METRIC';s.unit_settings.scale_length=1
s['CYBR_lighting']='User-selected grazing disk; relative to .135 m reference crop; position/diameter scaled by f, power by f squared; world .04';s['CYBR_material_id']=meta['id'];s['CYBR_material_status']='unselected candidate';s['CYBR_normal_contract']='Soil-only metric normal; physically weathered gravel with shared mica/quartz micro-relief, color and roughness';s['CYBR_physical_measurement']=False
bpy.context.preferences.filepaths.save_version=0;blend=out/f'CYBR_Victorville_{view}.blend';bpy.ops.wm.save_as_mainfile(filepath=str(blend),compress=True)
if len(args)>5 and args[5]=='build_only':print('SCENE_BUILT',blend);sys.exit(0)
s.render.filepath=str(out/f'CYBR_Victorville_{view}.png');start=time.monotonic();bpy.ops.render.render(write_still=True)
receipt=dict(material_id=meta['id'],source_maps=str(root),resolution=[res,res],samples=samples,render_seconds=time.monotonic()-start,engine='Cycles',device='CPU',threads=2,denoising='none',view=view,lighting='scaled low-grazing disk / user-selected asphalt reference',light_scale=factor,light_energy_w=ld.energy,light_disk_diameter_m=ld.size,world_strength=.04,physical_tile_m=tile,view_width_m=cam.ortho_scale,geometry_vertices=len(me.vertices),geometry_faces=len(me.polygons),geometry_model='Soil mesh plus 9389 independently clipped convex gravel bodies',all_images_packed=all(im.packed_file for im in bpy.data.images if im.source=='FILE'),map_resolution=meta['resolution'],selected=False)
receipt['image_sha256']=hashlib.sha256(Path(s.render.filepath).read_bytes()).hexdigest();(out/f'{view}_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
