"""Metric molded panels, independent explicit process-state r2.
A geometric manufacturing approximation, not injection-flow or thermal simulation.
Selected panel camera, light transforms, materials' base colors and principal gloss are frozen.
"""
import bpy, math, json, hashlib, sys, numpy as np
from pathlib import Path
from mathutils import Vector
R=Path('/workspace/shared/material-native4k/plastic-detail'); OUT=Path('/workspace/shared/material-slabs/scenes')

def aim(o,t):o.rotation_euler=(Vector(t)-o.location).to_track_quat('-Z','Y').to_euler()
def sig(s):
 return {'objects':{o.name:{'matrix':[v for row in o.matrix_world for v in row], 'type':o.type,'data':({'energy':o.data.energy,'size':o.data.size,'color':list(o.data.color),'shape':o.data.shape} if o.type=='LIGHT' else {'ortho_scale':o.data.ortho_scale,'type':o.data.type} if o.type=='CAMERA' else {})} for o in s.objects if o.type in {'LIGHT','CAMERA'} or o.name=='Slab / neutral floor'},'view':{'transform':s.view_settings.view_transform,'look':s.view_settings.look,'exposure':s.view_settings.exposure,'gamma':s.view_settings.gamma},'world':[list(s.world.node_tree.nodes['Background'].inputs[0].default_value),s.world.node_tree.nodes['Background'].inputs[1].default_value]}
def sample(h,x,y):
 # Float source map is native4096, sampled at the same global projected XY as the shader.
 qx=max(0,min(4095,(x/.08+.5)*4096-.5));qy=max(0,min(4095,(y/.08+.5)*4096-.5));ix=int(qx);iy=int(qy);fx=qx-ix;fy=qy-iy
 return float(h[iy,ix]*(1-fx)*(1-fy)+h[iy,min(ix+1,4095)]*fx*(1-fy)+h[min(iy+1,4095),ix]*(1-fx)*fy+h[min(iy+1,4095),min(ix+1,4095)]*fx*fy)
def set_input(p,name,value):
 if name in p.inputs:p.inputs[name].default_value=value

def run(cid):
 textured=cid.startswith('21_');src=OUT/(cid+'_selected_analytic_panel.blend');bpy.ops.wm.open_mainfile(filepath=str(src));s=bpy.context.scene
 original_sig=sig(s);old=next(o for o in s.objects if o.name.startswith(cid));mat=old.data.materials[0].copy();mat.name+=' / process-coherent r2'
 palette=list(mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value)
 bpy.data.objects.remove(old,do_unlink=True)
 height=np.load(R/'maps'/'r2'/'replicated_relief_native4096_m.npy') if textured else None
 # Shared integer lattice yields one closed manifold. 80mm panel / 3mm Minkowski edge fillet.
 N=512 if textured else 160; axis=np.linspace(-.04,.04,N+1)
 zs=np.unique(np.concatenate((np.linspace(0,.003,17),np.array([.0034,.0038,.00388,.00393,.003965,.004,.004035,.00407,.00412,.0042,.0046]),np.linspace(.005,.008,17))))
 K=len(zs)-1;vertices=[];keys={};faces=[];tags=[]
 draft=math.tan(math.radians(.8)); seamh=.000022; seamsigma=.000027
 def v(i,j,k):
  key=(i,j,k)
  if key in keys:return keys[key]
  q=np.array((axis[i],axis[j],zs[k]-.004),dtype=np.float64); core=np.minimum(np.maximum(q,[-.037,-.037,-.001]),[.037,.037,.001]);dv=q-core;normal=dv/np.linalg.norm(dv);pos=core+normal*.003;pos[2]+=.004
  # Both mold halves release away from their continuous z=4mm split. .8 degree outward draft.
  pos[:2]*=1-draft*abs(pos[2]-.004)/.04
  sidefactor=max(0,1-normal[2]**2); seam=seamh*math.exp(-.5*((pos[2]-.004)/seamsigma)**2)
  pos[:2]+=normal[:2]*seam*sidefactor
  # Textured upper mold cavity: replicated relief on the top and its fillet, bounded to32-36um.
  if textured and normal[2]>0:
   amp=normal[2]**6; pos+=normal*(sample(height,pos[0],pos[1])*amp)
  # Ejector witness solely on underside: four shallow6um recess rings in core-supported locations.
  if normal[2]<-.999:
   for cx in (-.023,.023):
    for cy in (-.023,.023):
     d=math.hypot(pos[0]-cx,pos[1]-cy);pos[2]+=.000006*math.exp(-.5*((d-.0022)/.000095)**2)+.000002/(1+math.exp(min(50,(d-.0022)/.00006)))
  idx=len(vertices);keys[key]=idx;vertices.append(tuple(pos));return idx
 def face(a,tag=0):faces.append(a);tags.append(tag)
 # Winding is outward on all six faces.
 for j in range(N):
  for i in range(N):
   face((v(i,j,K),v(i+1,j,K),v(i+1,j+1,K),v(i,j+1,K)))
   face((v(i,j,0),v(i,j+1,0),v(i+1,j+1,0),v(i+1,j,0)))
 for k in range(K):
  tag=int(abs((zs[k]+zs[k+1])*.5-.004)<.000065)
  for i in range(N):
   face((v(i,0,k),v(i+1,0,k),v(i+1,0,k+1),v(i,0,k+1)),tag)
   face((v(i,N,k),v(i,N,k+1),v(i+1,N,k+1),v(i+1,N,k)),tag)
  for j in range(N):
   face((v(0,j,k),v(0,j,k+1),v(0,j+1,k+1),v(0,j+1,k)),tag)
   face((v(N,j,k),v(N,j+1,k),v(N,j+1,k+1),v(N,j,k+1)),tag)
 mesh=bpy.data.meshes.new(cid+' / continuous injection molded body');mesh.from_pydata(vertices,[],faces);mesh.update();o=bpy.data.objects.new(mesh.name,mesh);s.collection.objects.link(o)
 mesh.materials.append(mat);split=mat.copy();split.name='Mold shutoff contact / same polymer palette';p=split.node_tree.nodes.get('Principled BSDF')
 for l in list(p.inputs['Roughness'].links):split.node_tree.links.remove(l)
 p.inputs['Roughness'].default_value=.44 if textured else .26;mesh.materials.append(split)
 for poly,tag in zip(mesh.polygons,tags):poly.use_smooth=True;poly.material_index=tag
 # Real surface normal is supplied by explicit geometry. Subgrid residual stays diagnostic only.
 if textured:
  nodes=mat.node_tree.nodes;links=mat.node_tree.links;p=nodes.get('Principled BSDF')
  for n in list(nodes):
   if n.type not in {'OUTPUT_MATERIAL','BSDF_PRINCIPLED'}:nodes.remove(n)
  tex=nodes.new('ShaderNodeTexImage');tex.name='Rounded tool impressions → polymer surface roughness';tex.label='Native4096 / same process field as geometry';img=bpy.data.images.load(str(R/'maps'/'r2'/'tool_roughness_native4096.png'),check_existing=True);img.colorspace_settings.name='Non-Color';img.pack();tex.image=img;tex.interpolation='Linear';tex.extension='EXTEND'
  coord=nodes.new('ShaderNodeNewGeometry');coord.label='World metric position; fixed 80mm physical domain'
  vect=nodes.new('ShaderNodeVectorMath');vect.operation='SCALE';vect.inputs[3].default_value=12.5;links.new(coord.outputs['Position'],vect.inputs[0]);add=nodes.new('ShaderNodeVectorMath');add.operation='ADD';add.inputs[1].default_value=(.5,.5,0);links.new(vect.outputs['Vector'],add.inputs[0]);links.new(add.outputs['Vector'],tex.inputs['Vector']);links.new(tex.outputs['Color'],p.inputs['Roughness'])
  # Coat must see the same geometric normals; no second artificial bump path.
  mat['process_field']='tool_state_native4096.npy';mat['geometry_relief_m']=json.loads((R/'receipts'/'native_tool_state_r2.json').read_text())['replicated_relief_peak_to_valley_m'];mat['approximation']='Rounded blast/etch tool impressions and90% replication; no flow or heat solve.'
 o['catalog_id']=cid;o['draft_angle_degrees']=.8;o['parting_plane_z_m']=.004;o['parting_witness_peak_radial_m']=seamh;o['parting_witness_fwhm_m']=2.35482*seamsigma;o['ejector_witness_count']=4;o['ejector_witness_face']='Underside only';o['ejector_witness_radius_m']=.0022;o['tool_fillet_radius_m']=.003;o['manufacturing_model']='Closed drafted plaque, continuous real parting witness, underside ejector contact. Geometric approximation; no melt-flow/heat solver.'
 # Authored detail exposes upper surface, adjacent side wall, and actual split witness in one frame.
 camdata=s.camera.data.copy();camdata.ortho_scale=.024;cam=bpy.data.objects.new('Plastic / authored molding detail camera',camdata);s.collection.objects.link(cam);cam.location=(.070,-.085,.043);aim(cam,(.023,-.033,.0048));cam['detail_frame_width_m']=.024;cam['detail_intent']='Molded finish and side shutoff witness at physical scale';cam['preserves_matched_studio_lighting']=True
 # Audits, computed from actual data (not a predeclared pass).
 current_sig=sig(s);current_sig['objects'].pop(cam.name)
 assert original_sig==current_sig, 'Selected camera/light/floor/view changed'
 coords=np.asarray(vertices);bounds=[coords.min(axis=0).tolist(),coords.max(axis=0).tolist()]
 edge_counts={}
 for f in faces:
  for a,b in zip(f,f[1:]+f[:1]):e=(a,b) if a<b else (b,a);edge_counts[e]=edge_counts.get(e,0)+1
 boundary=sum(v==1 for v in edge_counts.values());nonmanifold=sum(v!=2 for v in edge_counts.values());assert nonmanifold==0
 normals=np.array([tuple(p.normal) for p in mesh.polygons]);normlength=np.linalg.norm(normals,axis=1);assert np.all(np.isfinite(normals));assert float(normlength.min())>.999
 volume=sum(float(Vector(vertices[f[0]]).dot(Vector(vertices[f[k]]).cross(Vector(vertices[f[k+1]]))))/6 for f in faces for k in range(1,len(f)-1));assert volume>0
 floor=s.objects['Slab / neutral floor'].location.z;clearance=float(coords[:,2].min()-floor);assert clearance>=0 and clearance<.000002
 p=mat.node_tree.nodes['Principled BSDF'];assert list(p.inputs['Base Color'].default_value)==palette
 s['plastic_process_revision']='r2';s['plastic_process_approximation']='Geometry and native tool-state replication approximation, not an injection-flow or thermal solver';s['plastic_authored_detail_camera']=cam.name;s['review_status']='Geometry/numeric checked; rendering and pixel review pending'
 out=OUT/(cid+'_molding_detail_r2_panel.blend');bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True)
 receipt={'id':cid,'revision':'r2','source_scene':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scene':str(out),'scene_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'base_color_linear_rgba':palette,'frozen_studio_verified':True,'active_matched_camera':s.camera.name,'authored_detail_camera':cam.name,'actual_bounds_m':bounds,'actual_dimensions_m':(coords.max(axis=0)-coords.min(axis=0)).tolist(),'nominal_dimensions_m':[.08,.08,.008],'fillet_radius_m':.003,'draft_degrees':.8,'parting_plane_z_m':.004,'parting_witness_radial_peak_m':seamh,'parting_witness_fwhm_m':2.35482*seamsigma,'ejector_witness_face':'Underside only','ejector_count':4,'ejector_relief_max_m':.000008,'floor_clearance_m':clearance,'vertex_count':len(vertices),'face_count':len(faces),'boundary_edges':boundary,'nonmanifold_edges':nonmanifold,'normal_length_range':[float(normlength.min()),float(normlength.max())],'positive_oriented_volume_m3':volume,'textured_surface_geometry_grid':[N+1,N+1] if textured else None,'native_process_state_size':[4096,4096] if textured else None,'tool_state_couples':['geometric_relief','roughness'] if textured else [],'mapping_note':'Top and upper fillet relief replicated from a single native4096 tool field; lower half remains smooth. No top-face ejectors. No pigment grunge.','known_approximation':'Molding geometry and empirical90% relief replication. Does not implement heat transfer, injection flow, shrink solver or subsurface polymer structure.','pixel_review':'Pending render-owner review'}
 (R/'receipts'/(cid+'_molding_detail_r2.json')).write_text(json.dumps(receipt,indent=2));print('READY',str(out),receipt['scene_sha256'],flush=True)
for cid in sys.argv[sys.argv.index('--')+1:]:run(cid)
