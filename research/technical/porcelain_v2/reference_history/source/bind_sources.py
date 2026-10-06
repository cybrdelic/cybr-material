"""Bind admitted history-only scene and a separately HELD intrinsic-finish variant.
The effective layer control, geometry, optics, lighting and capture settings remain frozen.
"""
import ast
import copy
import hashlib
import json
from pathlib import Path
import sys
import bpy

HERE=Path(__file__).resolve().parents[1]
PRODUCTION=HERE.parents[1]
AUDIT=PRODUCTION/'experiments/porcelain_glaze_transfer_audit'
LAYER=PRODUCTION/'experiments/porcelain_glaze_layer_assignment'
OUT=HERE/'scenes'
GLAZE='07 new / mean-calibrated pigmented clearcoat surface'
R0=.17984575

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()

AUDIT_HELPERS_SHA='0ea5f1500b25e28ab7c5b3b3c505054711273f4c4df40ce822e709a68fd44087'
helper=AUDIT/'build_control_pair.py'
assert sha(helper)==AUDIT_HELPERS_SHA
for path,names in [(helper,{'scalar_rna','tree','snapshot','normal_state'}),
                   (LAYER/'source/build_layer_candidate.py',{'id_properties','frozen_state'})]:
    tree_ast=ast.parse(path.read_text())
    definitions=[n for n in tree_ast.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(definitions)==len(names)
    exec(compile(ast.Module(body=definitions,type_ignores=[]),str(path),'exec'))

def verify_history(before,after,normal_name,new_hash):
    normalized=copy.deepcopy(after)
    old_image=before['materials'][GLAZE]['tree']['nodes']['Image Texture']['image']
    new_image=normalized['materials'][GLAZE]['tree']['nodes']['Image Texture']['image']
    assert old_image['sha256']!=new_hash and new_image['sha256']==new_hash
    new_image['sha256']=old_image['sha256']
    assert normalized['images'][normal_name]['sha256']==new_hash
    normalized['images'][normal_name]=copy.deepcopy(before['images'][normal_name])
    assert normalized==before, 'Undeclared history-only source change'
    return {'normal_image_changed':True,'old_normal_sha256':old_image['sha256'],
            'new_normal_sha256':new_hash,'all_other_observed_data_equal':True,
            'roughness_image_and_shader_inputs_exact':True,'height_map_is_companion_not_second_bump':True}

def verify_intrinsic(before,after):
    normalized=copy.deepcopy(after)
    old=before['materials'][GLAZE]['tree']
    new=normalized['materials'][GLAZE]['tree']
    old_links=set(map(tuple,old['links']))
    new_links=set(map(tuple,new['links']))
    removed=('Image Texture.001','Color','Principled BSDF','Coat Roughness')
    assert old_links-new_links=={removed} and new_links-old_links==set()
    new['links']=copy.deepcopy(old['links'])
    changes=[]
    for i,(a,b) in enumerate(zip(old['nodes']['Principled BSDF']['inputs'],new['nodes']['Principled BSDF']['inputs'])):
        if a!=b:
            assert a['identifier']=='Coat Roughness'
            changes.append({'before':a,'after':b})
            new['nodes']['Principled BSDF']['inputs'][i]=copy.deepcopy(a)
    assert len(changes)==1 and normalized==before, 'Undeclared held optical-accounting source change'
    return {'removed_link':removed,'coat_roughness':changes,'intrinsic_R0':R0,
            'all_other_observed_data_equal':True,'status':'HELD; rendering not authorized'}

args=sys.argv[sys.argv.index('--')+1:]
assert args==['--admitted']
assert bpy.app.version[:3]==(4,3,2)
OUT.mkdir(exist_ok=True)
generation=json.loads((HERE/'receipts/generation.json').read_text())
resources=json.loads((HERE/'receipts/generation_resources.json').read_text())
assert resources['returncode']==0 and resources['stop_reason'] is None
assert generation['all_sensitivity_films_positive'] and generation['mean_and_volume_conservation_pass']
assert all(generation['checks'].values())
layer=json.loads((LAYER/'scenes/build_receipt.json').read_text())
control=next(s for s in layer['scenes'] if s['rig']=='main')
source=Path(control['scene'])
assert sha(source)==control['sha256']
frozen=json.loads((LAYER/'frozen_handoff.json').read_text())
assert frozen['first_capture']['source_sha256']==control['sha256']
normal=generation['maps']['Normal']
height=generation['maps']['Height']
assert sha(normal['path'])==normal['sha256'] and sha(height['path'])==height['sha256']
rough_hash=generation['first_binding_scope']['roughness_source_sha256']
bpy.ops.wm.open_mainfile(filepath=str(source),load_ui=False)
bpy.context.preferences.filepaths.save_version=0
before=frozen_state()
glaze=bpy.data.materials[GLAZE]
nt=glaze.node_tree
principled=nt.nodes['Principled BSDF']
normal_image=nt.nodes['Image Texture'].image
normal_name=normal_image.name
assert normal_image.colorspace_settings.name=='Non-Color'
assert sha(bpy.path.abspath(nt.nodes['Image Texture.001'].image.filepath))==rough_hash
assert [(l.from_node.name,l.from_socket.identifier) for l in principled.inputs['Coat Roughness'].links]==[('Image Texture.001','Color')]
assert [(l.from_node.name,l.from_socket.identifier) for l in principled.inputs['Coat Normal'].links]==[('Normal Map','Normal')]
assert not normal_image.packed_file
normal_image.filepath=normal['path']
normal_image.reload()
assert tuple(normal_image.size)==(4096,4096)
after=frozen_state()
history_delta=verify_history(before,after,normal_name,normal['sha256'])
history_path=OUT/'07_porcelain_reference_history_main.blend'
assert not history_path.exists()
bpy.ops.wm.save_as_mainfile(filepath=str(history_path),compress=True)
bpy.ops.wm.open_mainfile(filepath=str(history_path),load_ui=False)
reopened=frozen_state()
assert reopened==after
assert verify_history(before,reopened,normal_name,normal['sha256'])==history_delta
history_receipt={'variant':'history_only','status':'Source ready; no render performed',
                 'source':str(history_path),'sha256':sha(history_path),'bytes':history_path.stat().st_size,
                 'control':str(source),'control_sha256':control['sha256'],
                 'delta':history_delta,'before_state_sha256':digest(before),'after_state_sha256':digest(after),
                 'fresh_reopen_equal':True,'height_companion':height,
                 'coat_roughness_sha256':rough_hash,'shared_state_PBR_qualification':False}
# A separate source is saved for later optical accounting, held out of this first comparison.
nt=bpy.data.materials[GLAZE].node_tree
principled=nt.nodes['Principled BSDF']
for link in list(principled.inputs['Coat Roughness'].links): nt.links.remove(link)
principled.inputs['Coat Roughness'].default_value=R0
intrinsic=frozen_state()
intrinsic_delta=verify_intrinsic(after,intrinsic)
held_path=OUT/'07_porcelain_reference_history_intrinsic_HELD.blend'
assert not held_path.exists()
bpy.ops.wm.save_as_mainfile(filepath=str(held_path),compress=True)
bpy.ops.wm.open_mainfile(filepath=str(held_path),load_ui=False)
assert frozen_state()==intrinsic
assert verify_intrinsic(after,frozen_state())==intrinsic_delta
held_receipt={'variant':'intrinsic_only_finish','status':'HELD; do not render without a later release',
              'source':str(held_path),'sha256':sha(held_path),'bytes':held_path.stat().st_size,
              'control':str(history_path),'control_sha256':sha(history_path),
              'delta':intrinsic_delta,'fresh_reopen_equal':True,'after_state_sha256':digest(intrinsic)}
assert sha(source)==control['sha256']
assert sha(normal['path'])==normal['sha256'] and sha(height['path'])==height['sha256']
assert sha(generation['first_binding_scope']['roughness_source_path'])==rough_hash
report={'experiment':'porcelain_reference_history','source_only':True,'rendered':False,
        'selected':False,'visual_acceptance':False,'geometry_normals_uv_optics_lights_capture_unchanged':True,
        'original_layer_control_unchanged':True,'all_numeric_checks_pass':True,
        'generation_receipt':str(HERE/'receipts/generation.json'),
        'generation_receipt_sha256':sha(HERE/'receipts/generation.json'),
        'source_builder_sha256':sha(__file__),'observation_helpers_sha256':AUDIT_HELPERS_SHA,
        'scenes':[history_receipt,held_receipt],
        'first_comparison':'History-only, holding exact old coat roughness; not final shared-state PBR qualification',
        'held_variant':'Same new-history common state and retained intrinsic source roughness R0, without re-adding resolved normal variance',
        'strong_complete_baseline_library_file_id':'libfile_edc83cc981548191906982cc4fa7763b'}
(OUT/'build_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print('SOURCES_READY',json.dumps([{'variant':s['variant'],'source':s['source'],'sha256':s['sha256'],'status':s['status']} for s in report['scenes']]),flush=True)
