"""New local-source porcelain surface candidate; never a recovered original material."""
import bpy,bmesh,sys,json,hashlib,math
from pathlib import Path
from mathutils import Vector
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parents[1];S=R/'inputs/material-slabs';sys.path.insert(0,str(S/'source'));from studio import configure
state=json.load(open(R/'receipts/generation.json'));CASE=sys.argv[sys.argv.index('--')+1];assert CASE in ['control','conditional_poisson'];case=next(x for x in state['cases'] if x['case']==CASE);bpy.ops.wm.read_factory_settings(use_empty=True);s=bpy.context.scene
base_linear=(.66556436,.632415,.53975433)
def shader(name,c,rough):
 m=bpy.data.materials.new(name);m.use_nodes=True;p=m.node_tree.nodes['Principled BSDF'];p.inputs['Base Color'].default_value=(*c,1);p.inputs['Roughness'].default_value=rough;p.inputs['IOR'].default_value=1.49;return m
def tex(mat,name):
 n=mat.node_tree.nodes.new('ShaderNodeTexImage');n.image=bpy.data.images.load(str(R/'maps'/'Substrate_Height.png') if name=='Substrate_Height' else str(R/'maps'/CASE/('Glaze_Normal_OpenGL_RGB16.png' if name=='Glaze_Normal_OpenGL' else name+'.png')),check_existing=True);assert list(n.image.size)==[4096,4096];n.image.colorspace_settings.name='Non-Color';n.label=name+' / solved glaze state';return n
bodymat=shader('07 new / unglazed porcelain body',[x*.90 for x in base_linear],.67);bn=tex(bodymat,'Substrate_Height');bump=bodymat.node_tree.nodes.new('ShaderNodeBump');bump.inputs['Distance'].default_value=20e-6;bump.inputs['Strength'].default_value=1;bodymat.node_tree.links.new(bn.outputs['Color'],bump.inputs['Height']);bodymat.node_tree.links.new(bump.outputs['Normal'],bodymat.node_tree.nodes['Principled BSDF'].inputs['Normal'])
glaze=shader('07 new / mean-calibrated pigmented clearcoat surface',base_linear,.17984575);p=glaze.node_tree.nodes['Principled BSDF'];p.inputs['Coat Weight'].default_value=.32;p.inputs['Coat Roughness'].default_value=.12
normal=glaze.node_tree.nodes.new('ShaderNodeNormalMap');normal.uv_map='UVMap';ni=tex(glaze,'Glaze_Normal_OpenGL');rr=tex(glaze,'Glaze_Roughness');glaze.node_tree.links.new(ni.outputs['Color'],normal.inputs['Color']);glaze.node_tree.links.new(normal.outputs['Normal'],p.inputs['Normal']);glaze.node_tree.links.new(normal.outputs['Normal'],p.inputs['Coat Normal']);glaze.node_tree.links.new(rr.outputs['Color'],p.inputs['Roughness'])
W=.096;D=.096;H=.0036;film=state['retained_mean_film_m'];edge=.0006;period=state['parameters']['tile_m']
bpy.ops.mesh.primitive_cube_add(size=1,location=(0,0,H/2));body=bpy.context.object;body.name='07 new / fired porcelain coupon body';body.dimensions=(W,D,H);bpy.ops.object.transform_apply(location=False,rotation=False,scale=True);body.data.materials.append(bodymat);bev=body.modifiers.new('Fired eased600micron edges','BEVEL');bev.width=edge;bev.segments=8;bpy.context.view_layer.objects.active=body;bpy.ops.object.modifier_apply(modifier=bev.name)
# Retain the actually curved upper body surface as the exact inner boundary of the coating.
coat=body.copy();coat.data=body.data.copy();s.collection.objects.link(coat);coat.name='07 new / finite60micron mean glaze shell';coat.data.materials.clear();coat.data.materials.append(glaze)
bm=bmesh.new();bm.from_mesh(coat.data);remove=[f for f in bm.faces if f.normal.z<=.02 or f.calc_center_median().z<H/2-edge-.00001];bmesh.ops.delete(bm,geom=remove,context='FACES');bm.to_mesh(coat.data);bm.free();assert len(coat.data.polygons)>1
solid=coat.modifiers.new('Mean-thickness film / outward from exact body surface','SOLIDIFY');solid.thickness=film;solid.offset=1;solid.use_even_offset=True;bpy.context.view_layer.objects.active=coat;bpy.ops.object.modifier_apply(modifier=solid.name)
for ob in [body,coat]:
 while ob.data.uv_layers:ob.data.uv_layers.remove(ob.data.uv_layers[0])
 uv=ob.data.uv_layers.new(name='UVMap');uv.active_render=True
 for f in ob.data.polygons:
  f.use_smooth=abs(f.normal.z)<.999
  axis=max(range(3),key=lambda i:abs(f.normal[i]));coords=(0,1) if axis==2 else (0,2) if axis==1 else (1,2)
  for k in f.loop_indices:
   v=ob.data.vertices[ob.data.loops[k].vertex_index].co;uv.data[k].uv=(v[coords[0]]/period+.5,v[coords[1]]/period+.5)
 ob['catalog_id']='07_bone_porcelain';bm=bmesh.new();bm.from_mesh(ob.data);assert all(len(e.link_faces)==2 for e in bm.edges),ob.name;bm.free()
# Normal-proxy field is zero-mode normalized to the stated mean film thickness.
a=np.array(Image.open(R/'maps'/CASE/'Glaze_Height.png'),dtype='f4')/65535;g=(a-.5)*case['height_scale_m'];del a;a=np.array(Image.open(R/'maps/Substrate_Height.png'),dtype='f4')/65535;thickness=film+g-(a-.5)*20e-6;assert float(thickness.min())>0;tr=[float(thickness.min()),float(thickness.max()),float(thickness.mean())];del a,g,thickness
report=configure('07_bone_porcelain',[body,coat]);wide=s.camera;target=Vector((.033,-.033,H));direction=(wide.location-Vector((0,0,H/2))).normalized();cam=wide.copy();cam.data=wide.data.copy();s.collection.objects.link(cam);cam.name='Slab / porcelain glaze and edge detail';cam.location=target+direction*.12;cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler();cam.data.ortho_scale=.028;s.camera=wide;s.render.dither_intensity=0;s.render.image_settings.color_depth='16';s.name='Panel / new07 porcelain / solved glaze surface'
report.update(provenance='NEW local-source process candidate; no original native-file recovery or512px texture enlargement',scope='Actual linear capillary surface-leveling solver on prescribed particles. Porcelain body chemistry, bulk forming, firing and glaze spectral optics are not simulated.',palette_calibration={'source':'Existing local512px BaseColor/roughness, reduced to mean scalar values only','mean_linear_rgb':base_linear,'mean_roughness':.17984575,'native_image_detail_reused':False},coating_state=str(R/'receipts/generation.json'),periodic_representative_surface_cell_m=period,nominal_geometry_film_m=film,shared_state_thickness_range_and_mean_m=tr,microgeometry_approximation='Coarse curved shell follows exact substrate; subpixel solved coating variation is represented by native4096 tangent normals',all_meshes_closed=True,detail_camera=cam.name)
s['deposition_layout_control']=CASE
for ob in s.objects:
 if ob.type=='LIGHT' and ob.data.type=='AREA' and ob.name.startswith('Slab /'):ob.rotation_euler=(0,0,0)
out=R/'scenes'/('07_porcelain_arrival_'+CASE+'.blend');bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(out),compress=True)
report.update(scene=str(out),scene_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),arrival_case=CASE,source_recovery='Reconstructed from preserved exact builder and byte-verified native original fields; not byte-exact old blend recovery',optical_coefficients_unchanged=True,coating_mean_preserved_m=film,leveled_RMS_m=case['leveled_rms_m'],quality_acceptance=False)
(R/'receipts'/('panel_'+CASE+'.json')).write_text(json.dumps(report,indent=2));print('PORCELAIN_ARRIVAL_SOURCE_READY',CASE,report['scene_sha256'],flush=True)
