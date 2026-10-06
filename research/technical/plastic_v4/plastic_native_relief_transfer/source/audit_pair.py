import bpy,json,hashlib,numpy as np
from pathlib import Path
R=Path(__file__).resolve().parents[1];WORK=R.parents[3];EXP=R.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();fields=np.load(R/'state/interior_geometry_fields.npz');recovery=json.loads((R/'receipts/native_recovery.json').read_text());sigs=[];rows=[];positions={}
def val(x):
 try:return list(x)
 except TypeError:return x
for case,key in [('coarse','coarse'),('native','fine')]:
 report=json.loads((R/f'receipts/scene_{case}.json').read_text());src=Path(report['scene']);assert sha(src)==report['scene_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene;bpy.context.view_layer.update();ob=s.objects['21 / flat production interior / '+case];a=fields[key+'_axis_m'];n=len(a);h=fields[key+'_height_m'];verts=np.empty((len(ob.data.vertices),3),'f4');ob.data.vertices.foreach_get('co',verts.ravel());top=verts[:n*n].reshape(n,n,3);positions[case]=top.copy();expected=np.stack([np.tile(a,(n,1)),np.tile(a[:,None],(1,n)),.008+h],axis=2).astype('f4');assert np.array_equal(top,expected);assert tuple(ob.scale)==(1.,1.,1.);assert not ob.modifiers;assert not ob.data.uv_layers # World metric position is deliberately the only mapping.
 sig={'rig':{},'materials':{},'world':{},'view':{}}
 for o in s.objects:
  if o==ob:continue
  d={'type':o.type,'matrix':list(map(list,o.matrix_world))}
  if o.type=='CAMERA':d.update(type=o.data.type,ortho_scale=o.data.ortho_scale,lens=o.data.lens)
  if o.type=='LIGHT':d.update(type=o.data.type,energy=o.data.energy,color=list(o.data.color),size=o.data.size,shape=o.data.shape);assert tuple(o.rotation_euler)==(0,0,0)
  if o.type=='MESH':
   xyz=np.empty(len(o.data.vertices)*3,'f4');o.data.vertices.foreach_get('co',xyz);d['geometry_sha256']=hashlib.sha256(xyz.tobytes()).hexdigest()
  sig['rig'][o.name]=d
 for m in bpy.data.materials:
  if m.use_nodes and m.users:
   sig['materials'][m.name]={'nodes':{q.name:{'type':q.bl_idname,'inputs':[(i.identifier,val(i.default_value)) for i in q.inputs if hasattr(i,'default_value')],'image':sha(Path(bpy.path.abspath(q.image.filepath))) if q.type=='TEX_IMAGE' and q.image else None} for q in m.node_tree.nodes},'links':sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links)}
 bg=s.world.node_tree.nodes['Background'];sig['world']={'color':list(bg.inputs[0].default_value),'strength':bg.inputs[1].default_value};sig['view']={k:getattr(s.view_settings,k) for k in ['view_transform','look','exposure','gamma']};sigs.append(sig)
 p=ob.data.materials[0].node_tree.nodes['Principled BSDF'];assert not p.inputs['Normal'].links and not p.inputs['Coat Normal'].links;assert len(p.inputs['Roughness'].links)==1
 images=[]
 for im in bpy.data.images:
  if im.source=='FILE' and im.users:
   path=Path(bpy.path.abspath(im.filepath)).resolve();assert path.is_file() and list(im.size)==[4096,4096] and im.colorspace_settings.name=='Non-Color';assert path.read_bytes()[24]==16;images.append({'path':str(path.relative_to(WORK)),'sha256':sha(path)})
 assert len(images)==1 and images[0]['sha256']==recovery['native_roughness_sha256']
 eid='plastic_relief_'+case;eroot=EXP/eid;eroot.mkdir(exist_ok=True);manifest={'schema_version':1,'experiment_id':eid,'material_id':'21_plastic_texture','selected':False,'quality_acceptance':False,'source_scene':str(src.relative_to(WORK)),'source_sha256':sha(src),'dependencies':images,'contract':{'map_resolution':[4096,4096]},'capture':{'camera':'Plastic / matched physical interior','engine':'CYCLES','width':960,'height':960,'samples':128,'save_guides':True,'light_tree':False,'color_depth':16,'dither_intensity':0,'rss_cap_mib':3072,'time_limit_seconds':600},'purpose':'Fixed-state plastic flat-interior geometry sampling test; all optical coefficients and native roughness held equal','pass_gates':{'native_tool_state_bytes_exact':True,'production_state_bytes_exact':True,'fresh_saved_top_positions_exact':True,'closed_manifold_build_test':True,'material_realism':False},'physical_dimensions_m':[report['coupon_xy_extent_m'],report['coupon_xy_extent_m'],.003],'print_on_image':'Molded plastic · 19.972 mm interior study · '+('156 µm mesh control' if case=='coarse' else '19.5 µm native relief')}
 (eroot/'experiment.json').write_text(json.dumps(manifest,indent=2));rows.append({'case':case,'scene':str(src),'scene_sha256':sha(src),'experiment_id':eid,'experiment_manifest':str(eroot/'experiment.json'),'dependencies':images,'mesh_pitch_m':report['geometry_pitch_m'],'native_state_sha256':recovery['production_state_sha256']})
assert sigs[0]==sigs[1],'Rig, physical position mapping, shader or view changed';assert np.array_equal(positions['coarse'],positions['native'][::8,::8]);audit={'fresh_reopen':True,'saved_top_vertices_match_reconstructed_fields_exactly':True,'all_shared_coarse_native_samples_bit_exact':True,'rig_optical_coefficients_image_bytes_and_metric_mapping_equal':True,'native_roughness_4096_16bit':True,'roughness_partition_unqualified':True,'pair':rows};(R/'receipts/fresh_pair_audit.json').write_text(json.dumps(audit,indent=2));handoff={'ready_for_capture':True,'capture_released':True,'owner':'/root/improve_cloud_technical_materials','renderer':'/root/finish_clean_material_panels','capture_limit':2,'requested_views':['Plastic / matched physical interior'],'cases':rows,'scope':'One matched 19.972mm flat-interior coupon pair; no whole-part or curved-fillet claim','hypothesis':'Legacy coarse geometry under-resolves the unchanged native imprint state. Compare explicit native-pitch geometry without changing optics, height amplitude, palette, force, cooling or native roughness.','source_recovery':'New flat studies from byte-exact recovered original tool/state/maps; neither source is an old complete blend recovery','known_limitation':recovery['roughness_limitation'],'all_cases_unselected':True,'visual_gate_pending':True};(R/'capture_handoff.json').write_text(json.dumps(handoff,indent=2));print('PLASTIC_PAIR_READY',json.dumps(rows),flush=True)
