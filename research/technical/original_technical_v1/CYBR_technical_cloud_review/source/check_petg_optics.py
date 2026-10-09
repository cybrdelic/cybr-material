import bpy,bmesh,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'scenes/23_petg_transparent_r6_sideproof.blend';bpy.ops.wm.open_mainfile(filepath=str(p));s=bpy.context.scene;o=s.objects['PETG / unified printed wall rim and solid floor']
mat=o.data.materials[0];bs=mat.node_tree.nodes['Principled BSDF'];assert abs(bs.inputs['IOR'].default_value-1.57)<1e-5;assert bs.inputs['Transmission Weight'].default_value>.95;assert bs.inputs['Alpha'].default_value==1
assert abs(o['wall_thickness_m']-.0012)<1e-9 and abs(o['base_thickness_m']-.0012)<1e-9
nt=160;nz=round(.026/.0002)*8;ni=round((.026-.0012)/.0002)*8
outer=[f.normal.x*f.center.x+f.normal.y*f.center.y for f in o.data.polygons[:nt*nz]]
inner=[f.normal.x*f.center.x+f.normal.y*f.center.y for f in o.data.polygons[nt*nz:nt*(nz+ni)]]
assert min(outer)>0 and max(inner)<0
bm=bmesh.new();bm.from_mesh(o.data);assert all(e.is_manifold for e in bm.edges);vol=bm.calc_volume(signed=True);assert vol>0;bm.free()
bottom=min((o.matrix_world@v.co).z for v in o.data.vertices);floor=s.objects['Studio / floor'];top=floor.location.z+floor.dimensions.z/2;assert abs(bottom-top)<1e-7
report={'scene':str(p),'closed_manifold_boundary':True,'positive_volume_m3':vol,'outer_normals_radial_out':True,'inner_normals_radial_in':True,'ground_contact_error_m':bottom-top,'ior':bs.inputs['IOR'].default_value,'transmission':bs.inputs['Transmission Weight'].default_value,'alpha':bs.inputs['Alpha'].default_value,'wall_thickness_m':o['wall_thickness_m'],'layer_pitch_m':o['layer_height_m'],'status':'Physical setup checks only, rendered optical outcome still separate'}
(ROOT/'tests/clear_petg_r6_checks.json').write_text(json.dumps(report,indent=2)+'\n');print('PETG_OPTICS_SETUP_PASS',report,flush=True)
