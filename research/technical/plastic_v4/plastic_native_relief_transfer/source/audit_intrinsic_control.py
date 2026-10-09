from pathlib import Path
import bpy,json,hashlib,numpy as np
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=R.parents[3];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();b=json.loads((R/'receipts/intrinsic_build.json').read_text());sig=[]
def val(x):
 try:return list(x)
 except TypeError:return x
for path,key in [(b['source'],'control'),(b['scene'],'candidate')]:
 assert sha(path)==b['source_sha256' if key=='control' else 'scene_sha256'];bpy.ops.wm.open_mainfile(filepath=path);s=bpy.context.scene;bpy.context.view_layer.update();d={'objects':{},'principled':{},'world':{},'view':{}}
 for o in s.objects:
  q={'type':o.type,'matrix':list(map(list,o.matrix_world))}
  if o.type=='MESH':
   a=np.empty(len(o.data.vertices)*3,'f4');o.data.vertices.foreach_get('co',a);l=np.empty(len(o.data.loops),'i4');o.data.loops.foreach_get('vertex_index',l);smooth=np.empty(len(o.data.polygons),bool);o.data.polygons.foreach_get('use_smooth',smooth);q.update(vertices=hashlib.sha256(a.tobytes()).hexdigest(),topology=hashlib.sha256(l.tobytes()).hexdigest(),smooth=hashlib.sha256(smooth.tobytes()).hexdigest());assert not o.modifiers
  if o.type=='CAMERA':q.update(ortho_scale=o.data.ortho_scale,type=o.data.type,lens=o.data.lens)
  if o.type=='LIGHT':q.update(type=o.data.type,energy=o.data.energy,size=o.data.size,color=list(o.data.color),shape=o.data.shape)
  d['objects'][o.name]=q
 for m in bpy.data.materials:
  if m.use_nodes and m.users:
   p=m.node_tree.nodes.get('Principled BSDF');assert p
   v={i.name:val(i.default_value) for i in p.inputs if hasattr(i,'default_value')}
   if m.name.startswith('Plastic /'):
    assert p.inputs['Normal'].is_linked is False and p.inputs['Coat Normal'].is_linked is False
    if key=='candidate':assert abs(p.inputs['Roughness'].default_value-.475)<1e-7 and not p.inputs['Roughness'].links
    else:assert len(p.inputs['Roughness'].links)==1
    v.pop('Roughness')
   d['principled'][m.name]=v
 d['world']={'color':list(s.world.node_tree.nodes['Background'].inputs[0].default_value),'strength':s.world.node_tree.nodes['Background'].inputs[1].default_value};d['view']={k:getattr(s.view_settings,k) for k in ['view_transform','look','exposure','gamma']};sig.append(d)
assert sig[0]==sig[1],'Unexpected non-roughness change';report={'fresh_reopen':True,'all_geometry_topology_smooth_flags_objects_lights_cameras_world_view_and_other_optical_coefficients_exact':True,'only_material_change':'Old source-height broadening removed; original0.475 intrinsic roughness retained','numeric_height_sampling_receipt':str(R/'receipts/roughness_overlap_audit.json'),'geometric_representation_unchanged':True,'material_qualified':False};(R/'receipts/intrinsic_fresh_audit.json').write_text(json.dumps(report,indent=2))
base=json.loads((EXP/'plastic_relief_native/experiment.json').read_text());eid='plastic_relief_intrinsic';eroot=EXP/eid;eroot.mkdir(exist_ok=True);base.update(experiment_id=eid,source_scene=str(Path(b['scene']).relative_to(WORK)),source_sha256=b['scene_sha256'],dependencies=[],contract={},purpose='Single fixed-native-geometry roughness overlap control, using original assumed intrinsic0.475 only',pass_gates=report,print_on_image='Molded plastic · 28 mm field · intrinsic-only roughness control');(eroot/'experiment.json').write_text(json.dumps(base,indent=2));h={'ready_for_capture':True,'capture_released':True,'capture_limit':1,'experiment_id':eid,'experiment_manifest':str(eroot/'experiment.json'),'scene':b['scene'],'scene_sha256':b['scene_sha256'],'matched_baseline_image':'plastic_relief_native_capture_20261006_0845/CLEAN_OIDN_PANEL.png','scope':'Explicit roughness-budget diagnostic; no model/realism promotion','native_source_state_resolution':[4096,4096],'active_bitmap_images':0,'native_source_geometry_unchanged':True};(R/'intrinsic_capture_handoff.json').write_text(json.dumps(h,indent=2));print('INTRINSIC_READY',json.dumps(h))
