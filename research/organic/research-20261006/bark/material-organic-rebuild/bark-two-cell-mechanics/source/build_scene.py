"""Save a geometry-only Blender comparison. Never invokes a render operation."""
from pathlib import Path
import json,hashlib,math
import numpy as np
import bpy,bmesh
from mathutils import Vector

P=Path(__file__).resolve().parents[1]
data_path=P/'data/zero_shear_largest_admissible_96.npz'
receipt_path=P/'receipts/zero_shear_contact_limit_96.json'
record=json.loads(receipt_path.read_text());assert record['terminal']
assert record['qualified_accepted_states']
topology=record['topology'];par=topology['parameters']
data=np.load(data_path);strain=float(data['engineering_strain'])
ref=data['reference_positions_m'];deformed=data['positions_m']
assert strain==record['largest_accepted_strain']
state_sha=hashlib.sha256(data_path.read_bytes()).hexdigest()
bpy.ops.wm.read_factory_settings(use_empty=True)
sc=bpy.context.scene;sc.name='Planar radial cell-wall mechanism'
root=bpy.data.collections.new('BARK_TWO_CELL_MECHANISM');sc.collection.children.link(root)
annotations=bpy.data.collections.new('FLAT_LABELS');sc.collection.children.link(annotations)
stats=[]

def strip(points,name,collection,offset,wall):
    delta=np.diff(points,axis=0);tangent=delta/np.linalg.norm(delta,axis=1)[:,None]
    normals=np.c_[-tangent[:,1],tangent[:,0]]
    miters=np.empty_like(points);miters[0]=normals[0];miters[-1]=normals[-1]
    for i in range(1,len(points)-1):
        bisector=normals[i-1]+normals[i];bisector/=np.linalg.norm(bisector)
        miters[i]=bisector/(bisector@normals[i])
    edge=miters*par['wall_thickness_m']/2;depth=par['strip_width_m']/2
    verts=[]
    for p,e in zip(points,edge):
        verts.extend([(p[0]-e[0],-depth,p[1]-e[1]),
                      (p[0]+e[0],-depth,p[1]+e[1]),
                      (p[0]+e[0],depth,p[1]+e[1]),
                      (p[0]-e[0],depth,p[1]-e[1])])
    faces=[(3,2,1,0)]
    for i in range(len(points)-1):
        for k in range(4):faces.append((4*i+k,4*(i+1)+k,4*(i+1)+(k+1)%4,4*i+(k+1)%4))
    last=4*(len(points)-1);faces.append(tuple(last+k for k in range(4)))
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(verts,[],faces);mesh.update()
    bm=bmesh.new();bm.from_mesh(mesh);bmesh.ops.recalc_face_normals(bm,faces=list(bm.faces))
    manifold=all(e.is_manifold for e in bm.edges);volume=bm.calc_volume(signed=False)
    bm.to_mesh(mesh);bm.free();assert manifold and volume>0
    obj=bpy.data.objects.new(name,mesh);collection.objects.link(obj);obj.location=(offset,0,0)
    obj['wall_name']=wall['name'];obj['shared_wall']=wall['name']=='shared'
    obj['material_node_ids']=wall['nodes'];obj['mechanism_units']='metres; radians'
    obj['geometry_scope']='Piecewise-linear strip offset; unresolved fused junction overlap is retained'
    obj['source_state_sha256']=state_sha
    actual=np.array(verts).reshape(-1,4,3).mean(axis=1)[:,[0,2]]
    error=float(np.max(np.linalg.norm(actual-points,axis=1)))
    assert error<1e-15
    stats.append(dict(name=name,vertices=len(verts),faces=len(faces),closed_manifold=manifold,
                      strip_volume_m3=volume,material_centerline_error_m=error,shared_wall=obj['shared_wall']))

for name,positions,offset in [('REFERENCE',ref,-70e-6),('ZERO_MEAN_SHEAR',deformed,35e-6)]:
    coll=bpy.data.collections.new(name);root.children.link(coll)
    for wall in topology['walls']:
        strip(positions[np.array(wall['nodes'])],name+' / '+wall['name'],coll,offset,wall)

def label(name,text,xy,size):
    curve=bpy.data.curves.new(name,'FONT');curve.body=text;curve.size=size
    curve.extrude=0.;curve.bevel_depth=0.;curve.align_x='LEFT'
    obj=bpy.data.objects.new(name,curve);annotations.objects.link(obj)
    obj.location=(xy[0],-2e-6,xy[1]);obj.rotation_euler=(math.pi/2,0,0)
    obj['annotation_only']=True;return obj

label('Scope','Planar radial cell-wall mechanism',(-72e-6,69e-6),4.0e-6)
label('Limits','Uncalibrated elastic network; no PBR/fracture claim',(-72e-6,62e-6),2.7e-6)
label('Reference title','REFERENCE',(-70e-6,50e-6),3.0e-6)
label('Compressed title',f'{abs(strain)*100:g}% RADIAL SHORTENING',(35e-6,50e-6),3.0e-6)
label('Reference dimensions','43 um height; 1 um wall',(-70e-6,-8e-6),2.4e-6)
label('Shear boundary','Mean transverse shear = 0',(35e-6,-8e-6),2.4e-6)
label('Support boundary','Six prescribed radial junction coordinates; interiors and rotations free',(-72e-6,-18e-6),2.1e-6)
label('Axes','x: transverse section     z: radial     axial depth: 1 um',(-72e-6,-25e-6),2.1e-6)
cam_data=bpy.data.cameras.new('TWO_CELL_ORTHOGRAPHIC');cam=bpy.data.objects.new('TWO_CELL_ORTHOGRAPHIC',cam_data)
sc.collection.objects.link(cam);target=Vector((15e-6,0,22e-6))
cam.location=target+Vector((0,-300e-6,0));cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler()
cam.data.type='ORTHO';cam.data.ortho_scale=132e-6;cam.data.clip_start=1e-7;cam.data.clip_end=.1;sc.camera=cam
sc.unit_settings.system='METRIC';sc.unit_settings.scale_length=1.;sc.unit_settings.length_unit='MICROMETERS'
sc.render.engine='BLENDER_WORKBENCH';sc.display.shading.light='STUDIO'
sc.display.shading.color_type='SINGLE';sc.display.shading.single_color=(.55,.55,.55)
sc.display.shading.show_shadows=False;sc.display.shading.show_cavity=False
sc.display.shading.show_specular_highlight=False
sc.display.shading.background_type='WORLD'
sc.world=bpy.data.worlds.new('Neutral background');sc.world.color=(.92,.92,.92)
sc.render.resolution_x=1200;sc.render.resolution_y=720;sc.render.resolution_percentage=100
sc.render.image_settings.file_format='PNG';sc.render.image_settings.color_mode='RGB';sc.render.image_settings.color_depth='16'
sc.render.dither_intensity=0.;sc.view_settings.view_transform='Standard';sc.view_settings.look='None'
frame=np.array([tuple(v) for v in cam.data.view_frame(scene=sc)])
cam.data.ortho_scale*=132e-6/np.ptp(frame[:,1])
frame=np.array([tuple(v) for v in cam.data.view_frame(scene=sc)])
field_width=float(np.ptp(frame[:,0]));field_height=float(np.ptp(frame[:,1]))
sc['source_state_sha256']=state_sha;sc['source_scope']='Planar radial mechanism; zero-average-shear point-support control; not PBR or fracture'
sc['rendered_by_builder']=False;sc['wall_stock_count_per_specimen']=7;sc['unresolved_junction_radius_m']=2.5e-6
(P/'scenes').mkdir(exist_ok=True);out=P/'scenes/Bark_two_cell_zero_shear.blend'
bpy.ops.wm.save_as_mainfile(filepath=str(out))
handoff=dict(scene=str(out),scene_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_state_sha256=state_sha,
    source_state=str(data_path),append_collection='BARK_TWO_CELL_MECHANISM',labels_collection='FLAT_LABELS',
    camera='TWO_CELL_ORTHOGRAPHIC',field_width_m=field_width,field_height_m=field_height,
    resolution=[1200,720],render_engine='BLENDER_WORKBENCH',color_type='SINGLE',color=[.55,.55,.55],
    color_mode='RGB',color_depth='16',dither=0.,engineering_strain=strain,
    shared_wall_meshes_per_specimen=1,mesh_checks=stats,
    no_continuous_platen_solids=True,no_render_performed=True,visual_acceptance='pending central renderer capture')
(P/'receipts/render_handoff.json').write_text(json.dumps(handoff,indent=2));print(json.dumps(handoff,indent=2))
