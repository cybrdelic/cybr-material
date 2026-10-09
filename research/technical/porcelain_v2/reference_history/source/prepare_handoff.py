"""Freeze the single macro handoff; never start a render."""
import copy
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
P=R.parents[1]
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
build=json.loads((R/'scenes/build_receipt.json').read_text())
resources=json.loads((R/'scenes/resources.json').read_text())
generation=json.loads((R/'receipts/generation.json').read_text())
assert resources['returncode']==0 and resources['stop_reason'] is None
scene,held=build['scenes']
assert scene['variant']=='history_only' and held['variant']=='intrinsic_only_finish'
assert scene['fresh_reopen_equal'] and sha(scene['source'])==scene['sha256']
control=P/'runs/porcelain_layer_macro_capture_20261006_1836/CLEAN_OIDN_PANEL.png'
control_receipt=control.with_name(control.name+'.json')
capture=json.loads(control_receipt.read_text())
assert capture['source_sha256']==scene['control_sha256']
template=json.loads((P/'experiments/porcelain_effective_layer_candidate_macro/experiment.json').read_text())
manifest=copy.deepcopy(template)
manifest['experiment_id']='porcelain_reference_history_macro'
manifest['source_scene']='shared/material-production/'+str(Path(scene['source']).relative_to(P))
manifest['source_sha256']=scene['sha256']
for d in manifest['dependencies']:
    if d['path'].endswith('Glaze_Normal_OpenGL_RGB16.png'):
        d['path']='shared/material-production/'+str(Path(generation['maps']['Normal']['path']).relative_to(P))
        d['sha256']=generation['maps']['Normal']['sha256']
    assert sha(P/Path(d['path']).relative_to('shared/material-production'))==d['sha256']
manifest['purpose']='Imposed Q92 reference history, exact old coat roughness held for attribution; unselected'
manifest['print_on_image']='Porcelain · reference history · macro'
manifest['pass_gates']={'native4096_same_initial_packet_state':True,'positive_film':True,
                        'volume_conservation':True,'fresh_reopen_equal':True,
                        'geometry_uv_custom_normals_lights_cameras_world_exposure_equal_to_control':True,
                        'roughness_and_other_optical_inputs_equal_to_control':True,
                        'shared_state_PBR_qualification':False,'visual_acceptance':False}
assert manifest['capture']==template['capture']
assert all(capture['job'][key]==value for key,value in manifest['capture'].items())
(R/'capture').mkdir(exist_ok=True)
mp=R/'capture/experiment.json'
assert not mp.exists()
mp.write_text(json.dumps(manifest,indent=2)+'\n')
reference=json.loads((P/'experiments/porcelain_glaze_transfer_audit/baseline_reference.json').read_text())
handoff={'created_utc':datetime.now(timezone.utc).isoformat(),'experiment':'porcelain_reference_history',
         'source_only':True,'rendered':False,'selected':False,'quality_acceptance':False,
         'first_capture':{'view':'macro','manifest':str(mp),'sha256':sha(mp),
                          'source':scene['source'],'source_sha256':scene['sha256'],
                          'capture_gate':'One matched macro only. Full/second-light require later release.'},
         'matched_macro_control':{'source':scene['control'],'source_sha256':scene['control_sha256'],
                                  'image':str(control),'image_sha256':sha(control),
                                  'receipt':str(control_receipt),'receipt_sha256':sha(control_receipt),
                                  'capture_profile_bit_equal':True},
         'generation_receipt':str(R/'receipts/generation.json'),
         'generation_receipt_sha256':sha(R/'receipts/generation.json'),
         'build_receipt':str(R/'scenes/build_receipt.json'),'build_receipt_sha256':sha(R/'scenes/build_receipt.json'),
         'binding_resources':resources,'held_variant':held,
         'height_companion':generation['maps']['Height'],
         'protected_historical_appearance':reference['baseline'],
         'scope':'Reference history changes resolved height/normal only. Exact old coat roughness retained, so this is not final shared-state PBR qualification.',
         'renderer_instruction':'Reuse verified layer-control macro; capture only history-only source once with exact matched profile. Keep intrinsic-only source HELD. No promotion or upload.'}
hp=R/'frozen_handoff.json'
assert not hp.exists()
hp.write_text(json.dumps(handoff,indent=2)+'\n')
print(json.dumps({'handoff':str(hp),'sha256':sha(hp),'manifest':str(mp),'manifest_sha256':sha(mp)},indent=2))
