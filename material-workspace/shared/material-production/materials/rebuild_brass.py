"""Driver for preserved brass scene stages; no selected pointer mutation."""
import json,sys,math
from pathlib import Path
from .core import ROOT,read,sha,write_new,ContractError
from .__main__ import destination,run

def prepare(native,out):
 native=Path(native).resolve();out=destination(out);frozen=ROOT/'vendor/reconstruction-source'
 paths={'builder':frozen/'build_original_native4k_panels.py','tangent':frozen/'persist_brass_uv_tangent.py','normal_patch':frozen/'guarded_top_only_patch.py','golden_candidate_receipt':frozen/'candidate_receipt.json'}
 config={'stage':'upstream','out':str(out),'native':str(native),'commit':'de4146ac759164e29323b152d37f64417c00c3fb','metadata_sha256':'f3105ef3feb8df007177e3a861cdc8a8e5fd78ee88c5baa34f67893763ca8937','studio_source':str(ROOT.parent/'material-slabs/source'),**{k:str(p) for k,p in paths.items()},'source_sha256':{str(p):sha(p) for p in paths.values()}}
 write_new(out/'upstream_config.json',config)
 run(['blender','-b','--threads','1','--python-exit-code','1','--python',str(ROOT/'materials/brass_rebuild_stage.py'),'--',str(out/'upstream_config.json')],out/'upstream_build.log')
 return read(out/'upstream_ready.json')

def finish(out,admission):
 out=Path(out).resolve();config=read(out/'upstream_config.json');config.update(stage='top_normal',admission=str(Path(admission).resolve()));write_new(out/'top_normal_config.json',config)
 run(['blender','-b','--threads','1','--python-exit-code','1','--python',str(ROOT/'materials/brass_rebuild_stage.py'),'--',str(out/'top_normal_config.json')],out/'top_normal_build.log')
 return read(out/'reconstruction_result.json')


RECIPE_FIELDS=['dimensions_m','bounds_min_m','bounds_max_m','physical_height_scale_m','height_midlevel','physical_thickness_m','physical_edge_bevel_m','bevel_segments','bevel_limit_degrees','bevel_overlap_clamp','manifold_two_faces_per_edge','normals_consistent_outward','evaluated_geometry_counts','top_uv_domain','side_uvs','macro_geometry_grid','macro_geometry_sample_pitch_m','map_pixel_pitch_m','builder_script_sha256']
def compare_recipe_fields(old,new):
 checks={k:old.get(k)==new.get(k) for k in RECIPE_FIELDS}
 for k in ['signed_volume_m3','beveled_signed_volume_m3']:
  checks[k]=math.isclose(old[k],new[k],rel_tol=1e-10,abs_tol=1e-18)
 return checks

def qualify_reconstruction(out):
 out=Path(out).resolve();old=read(ROOT/'vendor/reconstruction-source/original_brass_builder_receipt.json');new=read(out/'receipts/05_champagne_brass_original_native4k_slab.json');candidate=read(out/'reconstruction_result.json');cfg=read(out/'upstream_config.json');native=read(Path(cfg['native'])/'native_recovery.json')
 checks=compare_recipe_fields(old,new)
 checks['all_native_integer_pixels_exact']=all(r['pixels_exact'] for r in native['channels'])
 checks['patch_preserves_reconstructed_geometry_and_controls']=not candidate['camera_light_world_geometry_changed']
 checks['top_face_partition']=candidate['corrected_evaluated_faces']==65536 and candidate['unchanged_control_evaluated_faces']==7237
 report={'passed':all(checks.values()),'checks':checks,'classification':'Reconstructed source-recipe conformance; NOT exact archived mesh-array or.blend identity','volume_relative_tolerance':1e-10,'volume_tolerance_purpose':'Allows floating-point reduction-order differences without changing physical parameters','archived_mesh_array_hashes_match':all(candidate['matches_archived_structure'].values()),'selected_registry_changed':False,'visual_acceptance':False,'generated_bevel_uv_limitation':'Automatic bevel edge/UV ordering can vary with Blender threads; top charts and source pixels are preserved','candidate':candidate['candidate'],'candidate_sha256':candidate['candidate_sha256']}
 write_new(out/'source_recipe_qualification.json',report)
 if not report['passed']:raise ContractError('Reconstructed source does not match the archived physical recipe')
 return report
