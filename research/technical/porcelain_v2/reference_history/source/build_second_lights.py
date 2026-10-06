"""Copy the already-frozen +90 degree world-Z light rig; never refit or render."""
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
OUT=HERE/'second_light_sources'
GLAZE='07 new / mean-calibrated pigmented clearcoat surface'
R0=.17984575

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()

helper=AUDIT/'build_control_pair.py'
assert sha(helper)=='0ea5f1500b25e28ab7c5b3b3c505054711273f4c4df40ce822e709a68fd44087'
for path,names in [(helper,{'scalar_rna','tree','snapshot','normal_state'}),
                  (LAYER/'source/build_layer_candidate.py',{'id_properties','frozen_state'}),
                  (HERE/'source/bind_sources.py',{'verify_history','verify_intrinsic'})]:
    parsed=ast.parse(path.read_text())
    defs=[n for n in parsed.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(defs)==len(names)
    exec(compile(ast.Module(body=defs,type_ignores=[]),str(path),'exec'))

assert sys.argv[sys.argv.index('--')+1:]==['--admitted']
assert bpy.app.version[:3]==(4,3,2)
OUT.mkdir(exist_ok=True)
TRANSFORMS=('rotation_mode','location','rotation_euler','rotation_quaternion','rotation_axis_angle',
            'scale','delta_location','delta_rotation_euler','delta_rotation_quaternion','delta_scale')
control=LAYER/'scenes/07_porcelain_effective_layer_second_light.blend'
control_hash='bd26cc4d693f2196e5fd5d39c12f02b5f3357848c540d90d5aed67727409a785'
assert sha(control)==control_hash
bpy.ops.wm.open_mainfile(filepath=str(control),load_ui=False)
bpy.context.preferences.filepaths.save_version=0
control_state=frozen_state()
transforms={}
for ob in bpy.context.scene.objects:
    if ob.type=='LIGHT' and ob.data.type=='AREA' and ob.name.startswith('Slab /'):
        assert ob.parent is None
        transforms[ob.name]={key:getattr(ob,key) if isinstance(getattr(ob,key),str) else tuple(getattr(ob,key)) for key in TRANSFORMS}
assert transforms
build=json.loads((HERE/'scenes/build_receipt.json').read_text())
normal=json.loads((HERE/'receipts/generation.json').read_text())['maps']['Normal']
receipts=[]
main_history_state=None
second_history_state=None
for index,source_rec in enumerate(build['scenes']):
    variant='history_only' if index==0 else 'intrinsic_R0'
    source=Path(source_rec['source'])
    assert sha(source)==source_rec['sha256']
    assert source_rec['sha256']==('6e54d1c375d421bd2ffcc1e70b5e503b16d8d677566357dcf807176f6468ddf8' if index==0 else '6385b8dcc70876cab8b4030fbc16d89d769df47e0b6aaec9778f409f437011a3')
    bpy.ops.wm.open_mainfile(filepath=str(source),load_ui=False)
    bpy.context.preferences.filepaths.save_version=0
    before=frozen_state()
    if index==0:main_history_state=before
    else:verify_intrinsic(main_history_state,before)
    expected=copy.deepcopy(before)
    for name,values in transforms.items():
        ob=bpy.context.scene.objects[name]
        assert ob.type=='LIGHT' and ob.parent is None
        old_light=before['objects'][name]
        saved_light=control_state['objects'][name]
        # Verify saved rig is the prior +90 world-Z rig, without computing replacement values.
        a=old_light['matrix']; b=saved_light['matrix']
        assert abs(b[0][3]+a[1][3])<1e-7 and abs(b[1][3]-a[0][3])<1e-7 and abs(b[2][3]-a[2][3])<1e-7
        normalized=copy.deepcopy(saved_light)
        normalized['matrix']=copy.deepcopy(old_light['matrix'])
        for key in TRANSFORMS:
            if key in old_light['properties']:
                normalized['properties'][key]=copy.deepcopy(old_light['properties'][key])
        assert normalized==old_light,'Saved second rig changes more than transforms: '+name
        for key in TRANSFORMS:setattr(ob,key,values[key])
        expected['objects'][name]=copy.deepcopy(saved_light)
    assert before['scene_properties']['light_rig_variant']=='main'
    assert control_state['scene_properties']['light_rig_variant']=='second_light'
    bpy.context.scene['light_rig_variant']='second_light'
    expected['scene_properties']['light_rig_variant']='second_light'
    after=frozen_state()
    assert after==expected,'Copy differs from exact saved control transforms'
    normal_name=bpy.data.materials[GLAZE].node_tree.nodes['Image Texture'].image.name
    if index==0:
        control_delta=verify_history(control_state,after,normal_name,normal['sha256'])
        second_history_state=after
    else:
        control_delta=verify_intrinsic(second_history_state,after)
    destination=OUT/('07_porcelain_reference_history_second_light.blend' if index==0 else '07_porcelain_reference_intrinsic_second_light.blend')
    assert not destination.exists()
    bpy.ops.wm.save_as_mainfile(filepath=str(destination),compress=True)
    bpy.ops.wm.open_mainfile(filepath=str(destination),load_ui=False)
    assert frozen_state()==after,'Fresh reopen differs'
    assert sha(source)==source_rec['sha256'] and sha(control)==control_hash
    receipts.append({'variant':variant,'source_main':str(source),'source_main_sha256':source_rec['sha256'],
                     'scene':str(destination),'sha256':sha(destination),'bytes':destination.stat().st_size,
                     'light_control':str(control),'light_control_sha256':control_hash,
                     'copied_exact_saved_transform_components':transforms,
                     'only_physical_change_is_light_transforms':True,'all_other_observed_data_equal':True,
                     'metadata_change':{'light_rig_variant':['main','second_light']},
                     'fresh_reopen_equal':True,'before_state_sha256':digest(before),'after_state_sha256':digest(after),
                     'comparison_delta':control_delta,'rendered':False})
report={'experiment':'porcelain_reference_history_second_lights','scenes':receipts,
        'source_builder_sha256':sha(__file__),'rendered':False,'selected':False,
        'retained_mean_state_and_images':True,'new_roughness_or_shader_changes':False,
        'scope':'Exact second-light derivatives for both previously frozen main sources, in one admitted source-only save'}
(OUT/'build_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print('SECOND_LIGHT_SOURCES_READY',json.dumps([{k:r[k] for k in ('variant','scene','sha256')} for r in receipts]),flush=True)
