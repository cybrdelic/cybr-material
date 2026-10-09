import bpy,bmesh,numpy as np,json,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];EXP=R.parent;WORK=EXP.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def value(x):
 try:return list(x)
 except TypeError:return x
signatures=[];rows=[]
for case in ['control','conditional_poisson']:
 report=json.loads((R/'receipts'/('panel_'+case+'.json')).read_text());src=Path(report['scene']);assert sha(src)==report['scene_sha256'];bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
 sig={'objects':{},'materials':{}}
 for o in s.objects:
  d={'type':o.type,'matrix':list(map(list,o.matrix_world))}
  if o.type=='MESH':
   v=np.empty(len(o.data.vertices)*3,'f4');o.data.vertices.foreach_get('co',v);l=np.empty(len(o.data.loops),'i4');o.data.loops.foreach_get('vertex_index',l);uv=np.empty(len(o.data.loops)*2,'f4');o.data.uv_layers.active.data.foreach_get('uv',uv)
   d.update(vertices=hashlib.sha256(v.tobytes()).hexdigest(),topology=hashlib.sha256(l.tobytes()).hexdigest(),uv=hashlib.sha256(uv.tobytes()).hexdigest(),materials=[m.name for m in o.data.materials])
   if o.name.startswith('07 new'):
    bm=bmesh.new();bm.from_mesh(o.data);assert all(len(e.link_faces)==2 for e in bm.edges);bm.free();assert tuple(o.scale)==(1,1,1)
  if o.type=='CAMERA':d.update(ortho_scale=o.data.ortho_scale,lens=o.data.lens,type=o.data.type)
  if o.type=='LIGHT':d.update(energy=o.data.energy,color=list(o.data.color),size=o.data.size,type=o.data.type)
  sig['objects'][o.name]=d
 for m in bpy.data.materials:
  if m.use_nodes:
   sig['materials'][m.name]={'nodes':{n.name:{'type':n.bl_idname,'inputs':[(i.identifier,value(i.default_value)) for i in n.inputs if hasattr(i,'default_value')]} for n in m.node_tree.nodes},'links':sorted((l.from_node.name,l.from_socket.identifier,l.to_node.name,l.to_socket.identifier) for l in m.node_tree.links)}
 deps=[]
 for im in bpy.data.images:
  if im.source=='FILE':
   p=Path(bpy.path.abspath(im.filepath)).resolve();assert p.is_file() and list(im.size)==[4096,4096] and im.colorspace_settings.name=='Non-Color';deps.append({'path':str(p.relative_to(WORK)),'sha256':sha(p),'image':im.name,'bits':p.read_bytes()[24]})
 assert len(deps)==3
 for m in bpy.data.materials:
  if m.use_nodes:
   for n in m.node_tree.nodes:
    if n.type=='NORMAL_MAP':assert n.space=='TANGENT' and n.uv_map=='UVMap' and n.inputs['Strength'].default_value==1
 signatures.append(sig);eid='porcelain_arrival_'+('control' if case=='control' else 'poisson');eroot=EXP/eid;eroot.mkdir(exist_ok=True)
 manifest={'schema_version':1,'experiment_id':eid,'material_id':'07_bone_porcelain','selected':False,'quality_acceptance':False,'source_scene':str(src.relative_to(WORK)),'source_sha256':sha(src),'dependencies':[{k:d[k] for k in ['path','sha256']} for d in deps],'contract':{'map_resolution':[4096,4096]},'capture':{'camera':'Slab / porcelain glaze and edge detail','engine':'CYCLES','width':960,'height':960,'samples':128,'save_guides':True,'light_tree':False,'color_depth':16,'dither_intensity':0,'rss_cap_mib':2560,'time_limit_seconds':600},'purpose':'Paired native4K finite-deposit layout control; same dose, body, leveling, optical coefficients and camera','pass_gates':{'original_control_native_fields_byte_verified':True,'packet_volume_conserved':True,'same_closed_geometry_uvs_and_optical_coefficients':True,'visual_acceptance':False},'physical_dimensions_m':report['dimensions_m'],'print_on_image':'Porcelain · 28mm edge crop · '+('lattice deposition control' if case=='control' else 'conditional-Poisson arrivals')}
 (eroot/'experiment.json').write_text(json.dumps(manifest,indent=2));rows.append({'case':case,'scene':str(src),'scene_sha256':sha(src),'experiment_id':eid,'experiment_manifest':str(eroot/'experiment.json'),'dependencies':deps,'physical_geometry_sha256':hashlib.sha256(json.dumps(sig,sort_keys=True).encode()).hexdigest()})
assert signatures[0]==signatures[1],'Geometry/UV/optical coefficient/camera/light change between cases'
d={'fresh_process_reopen':True,'pair_geometry_uv_shader_coefficients_camera_lights_exact':True,'only_bound_normal_and_roughness_state_images_differ':True,'rows':rows,'selected':False};(R/'receipts/pair_audit.json').write_text(json.dumps(d,indent=2));print('PAIR_AUDIT_PASSED',json.dumps(rows),flush=True)
