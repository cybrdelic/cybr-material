"""Build-only reconstructed complete part. Never edits selected source or registry.
Requires explicit resource admission and the renderer-owned analytic-shader check.
"""
from pathlib import Path
import bpy,numpy as np,json,hashlib,sys,gc
from mathutils import Matrix,Vector
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'source'));from analytic_shader import build_group
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan=json.loads((R/'receipts/native_coverage.json').read_text());assert plan['passes_numerical_lod'];detail=json.loads((R/'receipts/molding_details.json').read_text());assert detail['passes_sampled_details'];lod_detail=json.loads((R/'receipts/detail_lod_verification.json').read_text());assert lod_detail['passes'];assert lod_detail['new_mesh_sha256']==detail['arrays_sha256'];assert detail['triangles']<2000000
plan={**plan,'original_lod_arrays_sha256':plan['arrays_sha256'],'detail_preservation':detail,'arrays_path':detail['arrays_path'],'arrays_sha256':detail['arrays_sha256'],'vertices':detail['vertices'],'triangles':detail['triangles']}
admission=json.loads((R/'receipts/full_build_admission.json').read_text());assert admission['admitted']
shader_gate=json.loads((R/'receipts/shader_gate_acceptance.json').read_text());assert shader_gate['accepted_by_parent'];assert shader_gate['shader_source_sha256']==sha(R/'source/analytic_shader.py')
array=Path(plan['arrays_path']);assert sha(array)==plan['arrays_sha256'];a=np.load(array);q=a['q_m'];pos=a['position_m'];tri=a['triangles'];chart=a['chart_ids'];assert pos.dtype==np.float32
provpath=R.parent/'plastic_native_relief_transfer/inputs/production_scene.json';prov=json.loads(provpath.read_text());sig=prov['studio_signature'];rec=json.loads((R.parent/'plastic_native_relief_transfer/receipts/native_recovery.json').read_text())
src=Path('/workspace/scratch/b4387906eb93/technical-recovery-20261006/restored/CYBR_technical_cloud_review/scenes/21_plastic_texture_r5_hero.blend');source_sha=sha(src);bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene
with bpy.data.libraries.load(str(src),link=False) as (old,new):new.materials=['Plastic / Molded stipple']
mat=new.materials[0];assert mat;nodes=mat.node_tree.nodes;links=mat.node_tree.links;principled=nodes['Principled BSDF'];defaults={i.name:list(i.default_value) if hasattr(i.default_value,'__len__') else i.default_value for i in principled.inputs if hasattr(i,'default_value')};assert defaults['Base Color']==prov['base_color_linear_rgba']
# Preserve the frozen r2 shutoff material construction before stripping the field.
split=mat.copy();split.name='Mold shutoff contact / retained r2 recipe';sp=split.node_tree.nodes['Principled BSDF']
for link in list(sp.inputs['Roughness'].links):split.node_tree.links.remove(link)
sp.inputs['Roughness'].default_value=.44
split_normal_links=[{'from_node':l.from_node.name,'from_socket':l.from_socket.name,'to_socket':l.to_socket.name} for l in sp.inputs['Normal'].links]
for node in list(nodes):
    if node.type not in {'OUTPUT_MATERIAL','BSDF_PRINCIPLED'}:nodes.remove(node)
principled.inputs['Roughness'].default_value=.475
for im in list(bpy.data.images):
    if im.users==0:bpy.data.images.remove(im)
height=R.parent/'plastic_native_relief_transfer/maps/production_retained_height_native4096.png';assert sha(height)==rec['native_height_sha256'];image=bpy.data.images.load(str(height));image.colorspace_settings.name='Non-Color';group=build_group(image)
mesh=bpy.data.meshes.new('Original molded form / camera-qualified retained surface');mesh.vertices.add(len(pos));mesh.vertices.foreach_set('co',pos.ravel());mesh.loops.add(tri.size);mesh.loops.foreach_set('vertex_index',tri.ravel());mesh.polygons.add(len(tri));mesh.polygons.foreach_set('loop_start',np.arange(len(tri),dtype='i4')*3);mesh.polygons.foreach_set('loop_total',np.full(len(tri),3,'i4'));mesh.polygons.foreach_set('use_smooth',np.ones(len(tri),bool));mesh.update();assert not mesh.validate(verbose=False)
ob=bpy.data.objects.new('21 / complete molded plastic / retained analytic native finish',mesh);s.collection.objects.link(ob);at=mesh.attributes.new('OriginalChartQ','FLOAT_VECTOR','POINT');at.data.foreach_set('vector',q.astype('f4').ravel());mesh.polygons.foreach_set('material_index',chart.astype('i4'))
tangents=[((1,0,0),(0,1,0)),((0,1,0),(1,0,0)),((0,1,0),(0,0,1)),((0,0,1),(0,1,0)),((1,0,0),(0,0,1)),((0,0,1),(1,0,0))]
for c,(u,v) in enumerate(tangents):
    m=mat.copy();m.name='Plastic / retained native field / chart '+str(c);n=m.node_tree.nodes;l=m.node_tree.links;attr=n.new('ShaderNodeAttribute');attr.attribute_name='OriginalChartQ';g=n.new('ShaderNodeGroup');g.node_tree=group;l.new(attr.outputs['Vector'],g.inputs['Original q']);g.inputs['Tangent U'].default_value=u;g.inputs['Tangent V'].default_value=v;tr=n.new('ShaderNodeVectorTransform');tr.vector_type='NORMAL';tr.convert_from='OBJECT';tr.convert_to='WORLD';l.new(g.outputs['Normal object'],tr.inputs[0]);norm=n.new('ShaderNodeVectorMath');norm.operation='NORMALIZE';l.new(tr.outputs[0],norm.inputs[0]);l.new(norm.outputs[0],n['Principled BSDF'].inputs['Normal']);
    if 'Coat Normal' in n['Principled BSDF'].inputs:l.new(norm.outputs[0],n['Principled BSDF'].inputs['Coat Normal'])
    m['intrinsic_roughness_retained']=.475;m['limits']='Assumed generic intrinsic optical roughness; no material-specific optical calibration. Explicit mesh is camera LOD of the same height surface; no global20nm geometry claim.';mesh.materials.append(m)
mesh.materials.append(split);material_ids=chart.astype('i4');shutoff=json.loads((R/'receipts/shutoff_contract.json').read_text());band_bounds=shutoff['retained_qz_bounds_m'];center_z=q[tri,2].mean(axis=1);band=(chart>=2)&(center_z>band_bounds[0])&(center_z<band_bounds[1]);assert np.min(q[tri[band],2])>=band_bounds[0]-1e-15 and np.max(q[tri[band],2])<=band_bounds[1]+1e-15;material_ids[band]=6;mesh.polygons.foreach_set('material_index',material_ids);shutoff_triangles=int(band.sum());bpy.data.materials.remove(mat);del q,pos,tri,chart,a;gc.collect()
# Reconstruct exact historical camera/light transforms and properties, not a new studio fit.
for name,row in sig['objects'].items():
    if row['type']=='CAMERA':
        data=bpy.data.cameras.new(name);data.type=row['data']['type'];data.ortho_scale=row['data']['ortho_scale'];data.lens=row['data']['lens'];data.clip_start=1e-5;data.clip_end=100;obj=bpy.data.objects.new(name,data);s.collection.objects.link(obj)
    elif row['type']=='LIGHT':
        data=bpy.data.lights.new(name,'AREA')
        for key,val in row['data'].items():setattr(data,key,val)
        obj=bpy.data.objects.new(name,data);s.collection.objects.link(obj)
    elif name=='Slab / neutral floor':
        bpy.ops.mesh.primitive_plane_add(size=16);obj=bpy.context.object;obj.name=name;m=bpy.data.materials.new(name);m.use_nodes=True;p=m.node_tree.nodes['Principled BSDF'];p.inputs['Base Color'].default_value=(.18,.18,.18,1);p.inputs['Roughness'].default_value=.85;obj.data.materials.append(m)
    else:raise ValueError(name)
    obj.matrix_world=Matrix(np.array(row['matrix']).reshape(4,4).tolist())
    if row['type'] in {'CAMERA','LIGHT'}:
        obj.location=Vector(np.array(row['matrix']).reshape(4,4)[:3,3]);obj.scale=(1,1,1);target=(.023,-.033,.0048) if name=='Plastic / authored molding detail camera' else (0,0,.004);obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()
s.world=bpy.data.worlds.new('Slab / retained neutral world');s.world.use_nodes=True;bg=s.world.node_tree.nodes['Background'];bg.inputs[0].default_value=sig['world'][0];bg.inputs[1].default_value=sig['world'][1]
s.view_settings.view_transform=sig['view']['transform'];s.view_settings.look=sig['view']['look'];s.view_settings.exposure=sig['view']['exposure'];s.view_settings.gamma=sig['view']['gamma'];s.camera=s.objects[sig['active_camera']]
s.render.engine='CYCLES';s.cycles.device='CPU';s.cycles.samples=128;s.cycles.use_denoising=False;s.cycles.use_light_tree=False;s.cycles.use_adaptive_sampling=False;s.cycles.max_bounces=12;s.cycles.transmission_bounces=12;s.render.resolution_x=960;s.render.resolution_y=960;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.render.image_settings.color_mode='RGB';s.render.image_settings.color_depth='16';s.render.dither_intensity=0;s.unit_settings.system='METRIC';s.unit_settings.scale_length=1
s['selected']=False;s['historical_binary_recovered']=False;s['scope']='Complete reconstruction from frozen original form, exact retained native state and historical rig contracts; camera-dependent explicit geometry approximation. Actual appearance pending renderer.';ob['native_field_width_m']=.08;ob['physical_state_sha256']=rec['production_state_sha256'];ob['height_sha256']=rec['native_height_sha256'];ob['macro_formula_sha256']='31280c3ec450d251fa5da8a77f4bf4485fdccfbe0fe805571f8052e67e29698f'
def current_signature():
    return {'objects':{o.name:{'type':o.type,'matrix':[v for row in o.matrix_world for v in row],'data':({'energy':o.data.energy,'size':o.data.size,'color':list(o.data.color),'shape':o.data.shape} if o.type=='LIGHT' else {'type':o.data.type,'ortho_scale':o.data.ortho_scale,'lens':o.data.lens} if o.type=='CAMERA' else {})} for o in s.objects if o.type in {'LIGHT','CAMERA'} or o.name=='Slab / neutral floor'},'world':[list(bg.inputs[0].default_value),bg.inputs[1].default_value],'view':{'transform':s.view_settings.view_transform,'look':s.view_settings.look,'exposure':s.view_settings.exposure,'gamma':s.view_settings.gamma},'active_camera':s.camera.name}
bpy.context.view_layer.update();actual=current_signature();(R/'receipts/rig_reconstruction_comparison.json').write_text(json.dumps({'actual':actual,'expected':sig},indent=2)+'\n');assert actual==sig,'Reconstructed historical studio differs.'
bpy.context.preferences.filepaths.save_version=0;files=[]
for tag,cam in [('full','Slab / shared camera'),('macro','Plastic / authored molding detail camera')]:
    s.camera=s.objects[cam];out=R/'scenes'/('plastic_complete_'+tag+'.blend');assert not out.exists();bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert out.stat().st_size<40*2**20;files.append({'case':tag,'path':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'camera':cam})
# A reproducible second incident-light direction: rotate all light transforms90deg
# about original plaque center, preserving energy, source shape/size and world.
center=Vector((0,0,.004));rot=Matrix.Translation(center)@Matrix.Rotation(np.pi/2,4,'Z')@Matrix.Translation(-center)
for obj in s.objects:
    if obj.type=='LIGHT':obj.matrix_world=rot@obj.matrix_world
s.camera=s.objects['Slab / shared camera'];out=R/'scenes/plastic_complete_second_light.blend';assert not out.exists();bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True);assert out.stat().st_size<40*2**20;files.append({'case':'second_light','path':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'camera':s.camera.name,'rig':current_signature()})
assert sum(x['bytes'] for x in files)<120*2**20;assert sha(src)==source_sha;assert sha(height)==rec['native_height_sha256'];receipt={'reconstruction':True,'historical_binary_identity_claim':False,'files':files,'source_material_scene':str(src),'source_material_scene_sha256':source_sha,'source_material_defaults':defaults,'intrinsic_roughness':.475,'legacy_variance_broadening_removed':True,'retained_shutoff_material':{'roughness':.44,'original_r2_qz_band_m':band_bounds,'source_contract':shutoff,'triangles':shutoff_triangles,'source_normal_links':split_normal_links,'source_graph_preserved_except_original_r2_roughness_override':True},'native_height_sha256':rec['native_height_sha256'],'native_roughness_archived_sha256':rec['native_roughness_sha256'],'physical_state_sha256':rec['production_state_sha256'],'process_parameters':prov['same_other_process_parameters'],'pressure_Pa':prov['pressure_Pa'],'exact_rig_signature_match':True,'source_rig_signature_sha256':prov['studio_signature_sha256'],'analytic_normal_lobe_sockets':['Normal','Coat Normal when present'],'shader_group_nodes':len(group.nodes),'shader_source_sha256':sha(R/'source/analytic_shader.py'),'arrays_sha256':plan['arrays_sha256'],'vertices':len(mesh.vertices),'triangles':len(mesh.polygons),'geometry_lod':plan,'second_light_contract':'All area light object transforms rotated+90deg around world Z about(0,0,.004)m. Energy,size,shape,color,world,geometry and material stay identical.','calibration_limits':rec['model_limitations']+['Intrinsic roughness0.475 and source dielectric defaults are retained assumptions, not calibrated polymer optics.','Finite-precision native knot ownership has explicit conditioning uncertainty; see knot_conditioning.json.','Shared full-part mesh is a sampled camera LOD, not the20nm bounded seam representation.'],'appearance_qualified':False,'selected':False};(R/'receipts/complete_build.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k not in {'geometry_lod','process_parameters','source_material_defaults'}}),flush=True)
