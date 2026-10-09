"""Freeze exact shared second-light manifests; no capture is started."""
import copy,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
R=Path(__file__).resolve().parents[1];P=R.parents[1]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
b=json.loads((R/'second_light_sources/build_receipt.json').read_text())
r=json.loads((R/'second_light_sources/resources.json').read_text())
assert r['returncode']==0 and r['stop_reason'] is None
template=json.loads((P/'experiments/porcelain_effective_layer_candidate_second_light/experiment.json').read_text())
normal=json.loads((R/'receipts/generation.json').read_text())['maps']['Normal']
records=[]
for name,scene in [('control',None),('history',b['scenes'][0]),('intrinsic',b['scenes'][1])]:
    m=copy.deepcopy(template);m['experiment_id']='porcelain_reference_'+name+'_second_light'
    if scene:
        assert sha(scene['scene'])==scene['sha256']
        assert scene['fresh_reopen_equal'] and scene['only_physical_change_is_light_transforms']
        m['source_scene']='shared/material-production/'+str(Path(scene['scene']).relative_to(P));m['source_sha256']=scene['sha256']
        for d in m['dependencies']:
            if d['path'].endswith('Glaze_Normal_OpenGL_RGB16.png'):
                d['path']='shared/material-production/'+str(Path(normal['path']).relative_to(P));d['sha256']=normal['sha256']
    m['purpose']='Matched +90 degree saved light rig; '+name+' comparison; unselected'
    m['print_on_image']='Porcelain · '+name+' · second light'
    m['pass_gates']={'saved_control_light_transforms_exact':True,'fresh_reopen_equal':True,'frozen_main_material_retained':True,'visual_acceptance':False}
    assert m['capture']==template['capture']
    assert sha(P/Path(m['source_scene']).relative_to('shared/material-production'))==m['source_sha256']
    for d in m['dependencies']:assert sha(P/Path(d['path']).relative_to('shared/material-production'))==d['sha256']
    out=R/'capture'/f'{name}_second_light.json';assert not out.exists();out.write_text(json.dumps(m,indent=2)+'\n')
    records.append({'variant':name,'manifest':str(out),'manifest_sha256':sha(out),'source':str(P/Path(m['source_scene']).relative_to('shared/material-production')),'source_sha256':m['source_sha256']})
h={'created_utc':datetime.now(timezone.utc).isoformat(),'experiment':'porcelain_reference_second_lights','manifest_records':records,
   'shared_capture_profile':template['capture'],'all_capture_profiles_bit_equal':True,
   'build_receipt':str(R/'second_light_sources/build_receipt.json'),'build_receipt_sha256':sha(R/'second_light_sources/build_receipt.json'),
   'resources':r,'exact_saved_light_control_sha256':'bd26cc4d693f2196e5fd5d39c12f02b5f3357848c540d90d5aed67727409a785',
   'copied_transforms_without_refitting':True,'metadata_change':'Descriptive light_rig_variant main -> second_light only; all physical changes limited to exact copied light transforms.',
   'rendered_here':False,'selected':False,'scope':'Both exact derivatives prepared; renderer/operator chooses released captures. Intrinsic macro uses unchanged6385… main under separate release.'}
out=R/'second_light_handoff.json';assert not out.exists();out.write_text(json.dumps(h,indent=2)+'\n')
print(json.dumps({'handoff':str(out),'sha256':sha(out),'records':records},indent=2))
