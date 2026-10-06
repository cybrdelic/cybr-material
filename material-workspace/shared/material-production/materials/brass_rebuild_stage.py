"""Blender recovery adapter. Executes preserved geometry/shader stages unchanged."""
import ast,bpy,bmesh,sys,json,hashlib,math,resource,importlib.util
from pathlib import Path
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from materials.core import read,sha,check_hash,write_new,ROOT

args=sys.argv[sys.argv.index('--')+1:];config=read(args[0]);out=Path(config['out']);native=Path(config['native']);cid='05_champagne_brass'
for p,h in config['source_sha256'].items():check_hash(p,h)
recovery=read(native/'native_recovery.json')
for row in recovery['channels']:check_hash(row['path'],row['sha256']);assert row['pixels_exact']
P=native/'materials'/cid;metadata=read(P/'material.json');assert sha(P/'material.json')==config['metadata_sha256']
assert tuple(bpy.app.version)==(4,3,2)

def script(path,overrides):
 tree=ast.parse(Path(path).read_text());seen=set()
 for n in tree.body:
  if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in overrides:
   name=n.targets[0].id;n.value=ast.Call(func=ast.Name(id='Path',ctx=ast.Load()),args=[ast.Constant(str(overrides[name]))],keywords=[]);seen.add(name)
 assert seen==set(overrides),seen
 ast.fix_missing_locations(tree);scope={'__name__':'recovered_verified_stage','__file__':str(path)};exec(compile(tree,str(path),'exec'),scope);return scope

if config['stage']=='upstream':
 base=Path(config['builder']);tree=ast.parse(base.read_text());indices=[i for i,n in enumerate(tree.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='bpy.ops.wm.read_factory_settings'];assert len(indices)==1
 # Only input verification/root plumbing and provenance labels are adapted. The
 # entire mesh/UV/bevel/material/studio algorithm below is the preserved source.
 body=ast.Module(body=tree.body[indices[0]:],type_ignores=[])
 class Provenance(ast.NodeTransformer):
  def visit_Constant(self,node):
   replacements={'Byte-verified original native-4096 maps; no upsampling or regeneration':'Source-regenerated native4096 maps; all integer pixels match saved canonical hashes; PNG compression differs','Exact recovered original native-4096 material':'Reconstructed original native4096 material; selected registry unchanged'}
   if isinstance(node.value,str) and node.value in replacements:return ast.copy_location(ast.Constant(replacements[node.value]),node)
   return node
 body=Provenance().visit(body);ast.fix_missing_locations(body)
 sys.path.insert(0,config['studio_source']);from studio import configure
 sources={row['channel']+'.png':{'sha256':row['sha256'],'pixels_exact':True,'canonical_pixel_sha256':row['canonical_integer_pixel_sha256'],'classification':'source replay; not original compressed bytes'} for row in recovery['channels']}
 sources['material.json']={'sha256':sha(P/'material.json'),'classification':'exact original metadata bytes'}
 scope=globals().copy();scope.update(OUT=out,CID=cid,COMMIT=config['commit'],sources=sources,__file__=str(base),configure=configure)
 exec(compile(body,str(base),'exec'),scope)
 # The preserved tangent stage reads the exact base panel we just constructed.
 script(config['tangent'],{'R':out})
 upstream=out/'scenes/05_champagne_brass_original_native4k_slab_r2_uv_tangent.blend'
 write_new(out/'upstream_ready.json',{'source':str(upstream),'source_sha256':sha(upstream),'reconstructed':True,'source_maps_integer_exact':True,'geometry_source_sha256':sha(base),'tangent_source_sha256':sha(config['tangent']),'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'selected_registry_changed':False})
 print('UPSTREAM_READY',upstream,flush=True)
elif config['stage']=='top_normal':
 spec=importlib.util.spec_from_file_location('preserved_top_normal_patch',config['normal_patch']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 module.ROOT=out/'normal_transfer';module.ROOT.mkdir(exist_ok=True);module.SOURCE_ROOT=out/'scenes';module.MAP_ROOT=native/'materials'
 module.main(['--material',cid,'--admission-file',config['admission']])
 report=read(module.ROOT/'candidates'/cid/'candidate_receipt.json');gold=read(config['golden_candidate_receipt'])
 comparisons={key:report[key]==gold[key] for key in ['base_geometry_uv_normal_sha256','evaluated_geometry_uv_normal_sha256','corrected_evaluated_faces','unchanged_control_evaluated_faces']}
 report.update(recovery_classification='Reconstructed derivative; original binary scene SHA is not claimed',matches_archived_structure=comparisons,selected_registry_changed=False)
 write_new(out/'reconstruction_result.json',report)
 if not all(comparisons.values()):print('RAW_ARRAY_IDENTITY_UNMATCHED_REQUIRE_SOURCE_RECIPE_QUALIFICATION',json.dumps(comparisons),flush=True)
 else:print('RECONSTRUCTION_STRUCTURE_EXACT',json.dumps(comparisons),flush=True)
else:raise RuntimeError('Unknown stage')
