"""Blender recovery of preserved selected32 construction; never renders.

Replays the frozen geometry.py slab and fix_specimen_edges.py algorithms.
Native image and normalized geometry hashes are checked before use. This is a
new binary scene, not a claim that the vanished .blend bytes were recovered.
"""
import bpy,sys,json,hashlib,math,time
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from materials.core import ROOT,read,sha,write_new,workspace,material

cfg=read(sys.argv[sys.argv.index('--')+1]);out=Path(cfg['out']);native=Path(cfg['native'])
assert bpy.app.version==(4,3,2)
for path,digest in cfg['input_hashes'].items():assert sha(path)==digest,path
meta=read(native/'native_recovery.json');maps=native/'maps'
assert meta['all_pixels_exact'] and meta['all_png_bytes_exact']
for channel,row in meta['channels'].items():assert sha(maps/(channel+'.png'))==row['sha256']
a=np.load(maps/'GeometryHeight.npy')[::-1].copy()
assert sha(maps/'GeometryHeight.npy')=='7cbbbea96e7268d94298dae2386738063c8ca8610449a73b5dfba93450781572'
recipe=read(ROOT/'vendor/aggregate-reconstruction/selected_recipe_1536.json')
bpy.ops.wm.read_factory_settings(use_empty=True)
scene=bpy.context.scene;scene.unit_settings.system='METRIC';scene.unit_settings.scale_length=1
mat=bpy.data.materials.new(recipe['name']);mat.use_nodes=True
p=mat.node_tree.nodes.get('Principled BSDF');nodes=mat.node_tree.nodes;links=mat.node_tree.links
p.inputs['Base Color'].default_value=(*recipe['color'],1)
p.inputs['Roughness'].default_value=recipe['roughness'];p.inputs['Metallic'].default_value=recipe['metallic'];p.inputs['IOR'].default_value=recipe['ior']
for channel,socket in [('BaseColor','Base Color'),('Roughness','Roughness'),('Metallic','Metallic')]:
 t=nodes.new('ShaderNodeTexImage');t.image=bpy.data.images.load(str(maps/(channel+'.png')));t.image.colorspace_settings.name='sRGB' if channel=='BaseColor' else 'Non-Color';t.label=channel;links.new(t.outputs['Color'],p.inputs[socket])
t=nodes.new('ShaderNodeTexImage');t.image=bpy.data.images.load(str(maps/'Normal_Micro_OpenGL.png'));t.image.colorspace_settings.name='Non-Color'
normal=nodes.new('ShaderNodeNormalMap');links.new(t.outputs['Color'],normal.inputs['Color']);links.new(normal.outputs[0],p.inputs['Normal'])
n=a.shape[0]-1;width=recipe['tile_m'];verts=[];uv=[];faces=[]
for j in range(n+1):
 for i in range(n+1):
  x=(i/n-.5)*width;y=(j/n-.5)*width
  verts.append((x,y,.008+(a[j,i]-.5)*recipe['height_scale_m']));uv.append((i/n,j/n))
for j in range(n):
 for i in range(n):
  q=j*(n+1)+i;faces.append((q,q+1,q+n+2,q+n+1))
top_count=len(faces);top_vertices=len(verts)
boundary=list(range(n+1))+[j*(n+1)+n for j in range(1,n+1)]+[n*(n+1)+i for i in range(n-1,-1,-1)]+[j*(n+1) for j in range(n-1,0,-1)]
bottom=[]
for idx in boundary:
 bottom.append(len(verts));x,y,z=verts[idx];verts.append((x,y,-.00095));uv.append(uv[idx])
for i,aidx in enumerate(boundary):
 k=(i+1)%len(boundary);faces.append((aidx,bottom[i],bottom[k],boundary[k]))
faces.append(tuple(reversed(bottom)))
mesh=bpy.data.meshes.new('Concrete / selected reconstructed mesh');mesh.from_pydata(verts,[],faces);mesh.update()
ob=bpy.data.objects.new('Concrete / real displaced face',mesh);scene.collection.objects.link(ob);mesh.materials.append(mat)
core=bpy.data.materials.new('Specimen / clean cut core');core.use_nodes=True;cp=core.node_tree.nodes.get('Principled BSDF');cp.inputs['Base Color'].default_value=(.19,.185,.165,1);cp.inputs['Roughness'].default_value=.91;mesh.materials.append(core)
layer=mesh.uv_layers.new(name='UVMap')
for poly in mesh.polygons:
 poly.use_smooth=poly.index<top_count;poly.material_index=0 if poly.index<top_count else 1
 for li in poly.loop_indices:layer.data[li].uv=uv[mesh.loops[li].vertex_index]
sys.path.insert(0,str(ROOT.parent/'material-slabs/source'));from studio import configure
studio=configure('32_concrete_polished',[ob]);bpy.context.view_layer.update()
fixture=read(ROOT/'tests/fixtures/concrete_inspection.json');target=fixture['objects'][0]
assert len(mesh.vertices)==target['vertices'] and len(mesh.polygons)==target['polygons']
assert np.max(np.abs(np.asarray(ob.dimensions)-target['dimensions_m']))<2e-9
co=np.empty(len(mesh.vertices)*3,np.float32);mesh.vertices.foreach_get('co',co);co=co.reshape(-1,3)
assert np.array_equal(co,np.asarray(verts,dtype=np.float32))
assert len(mesh.materials)==2 and not ob.modifiers
assert all(np.array_equal(np.asarray(layer.data[li].uv),np.asarray(uv[mesh.loops[li].vertex_index],dtype=np.float32)) for poly in mesh.polygons for li in poly.loop_indices)
edge_counts=np.zeros(len(mesh.edges),dtype=np.int32)
for poly in mesh.polygons:
 for li in poly.loop_indices:edge_counts[mesh.loops[li].edge_index]+=1
assert np.all(edge_counts==2),'Nonmanifold recovery'
scene['recovery_status']='Binary scene reconstructed from preserved source; exact native maps and GeometryHeight; no visual promotion'
scene['recovery_old_scene_sha256']=material('32')['selected']['scene_sha256']
path=out/'32_selected_reconstructed.blend';bpy.ops.wm.save_as_mainfile(filepath=str(path))
report={'source':str(path),'source_sha256':sha(path),'binary_scene_reconstructed':True,'original_binary_sha_match_claimed':False,'exact_geometry_state':True,'native_maps_all_bytes_exact':True,'vertices':len(mesh.vertices),'polygons':len(mesh.polygons),'top_vertices':top_vertices,'top_faces':top_count,'float32_coordinates_match_preserved_algorithm':True,'all_uvs_match_preserved_algorithm':True,'fixture_dimensions_match':True,'every_edge_has_two_faces':True,'studio':studio,'input_hashes':cfg['input_hashes'],'rendered':False,'selected_registry_changed':False,'known_limits':['Retained authored material remains visually unqualified','Prior residual-normal adapter is preserved, not independently repaired during recovery','Scene byte identity and unused historical datablocks cannot be recovered from source alone']}
write_new(out/'scene_reconstruction.json',report)
capture=material('32')['selected']['capture'];deps=[{'path':str((maps/(c+'.png')).relative_to(workspace())),'sha256':r['sha256']} for c,r in meta['channels'].items() if c in {'BaseColor','Roughness','Metallic','Normal_Micro_OpenGL'}]
manifest={'schema_version':1,'selected':False,'quality_acceptance':False,'experiment_id':'selected32_recovery','material_id':'32_concrete_polished','source_scene':str(path.relative_to(workspace())),'source_sha256':sha(path),'dependencies':deps,'capture':capture,'purpose':'Recovery of retained selected32 authored baseline. Exact maps/geometry state, semantically reconstructed scene. No new optical model.','pass_gates':report,'contract':{'map_resolution':[4096,4096]}}
write_new(out/'experiment.json',manifest)
print('RECONSTRUCTED32_READY',json.dumps(report),flush=True)
