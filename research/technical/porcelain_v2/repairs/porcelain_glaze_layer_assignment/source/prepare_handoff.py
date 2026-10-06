"""Freeze source receipts and renderer inputs without starting a capture."""
import copy
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parents[1]
production = root.parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


receipt = json.loads((root / 'scenes/build_receipt.json').read_text())
resources = json.loads((root / 'scenes/resources.json').read_text())
assert resources['returncode'] == 0 and resources['stop_reason'] is None
assert receipt['coefficient_identical'] is False and receipt['all_map_bytes_unchanged']
scenes = {s['rig']: s for s in receipt['scenes']}
for s in scenes.values():
    assert sha(s['scene']) == s['sha256']
    assert sha(s['control']) == s['control_sha256']
    assert s['fresh_reopen_equal']
for name, expected in receipt['original_map_files_sha256'].items():
    assert sha(production / name) == expected

reference_path = production / 'experiments/porcelain_glaze_transfer_audit/baseline_reference.json'
reference = json.loads(reference_path.read_text())
assert sha(reference['baseline']['labeled_image']) == reference['baseline']['labeled_image_sha256']
control_capture = production / 'runs/porcelain_shell_candidate_macro_capture_20261006_1817/CLEAN_OIDN_PANEL.png'
control_receipt = control_capture.with_name(control_capture.name + '.json')
capture = json.loads(control_receipt.read_text())
assert capture['source_sha256'] == scenes['main']['control_sha256']
records = []
for view, rig in [('macro', 'main'), ('full', 'main'), ('second_light', 'second_light')]:
    template_path = production / f'experiments/porcelain_shell_normal_candidate_{view}/experiment.json'
    template = json.loads(template_path.read_text())
    manifest = copy.deepcopy(template)
    experiment_id = f'porcelain_effective_layer_candidate_{view}'
    manifest['experiment_id'] = experiment_id
    manifest['source_scene'] = 'shared/material-production/' + str(Path(scenes[rig]['scene']).relative_to(production))
    manifest['source_sha256'] = scenes[rig]['sha256']
    manifest['purpose'] = 'Effective substrate/glaze field-role candidate; optical response deliberately changed; unselected'
    manifest['print_on_image'] = f'Porcelain · effective layer · {view.replace("_", " ")}'
    manifest['pass_gates'] = {
        'geometry_uv_custom_normals_lights_cameras_world_exposure_equal_to_control': True,
        'all_map_bytes_unchanged': True, 'fresh_reopen_equal': True,
        'optical_response_changed': True, 'coefficient_identical': False,
        'visual_acceptance': False,
    }
    assert manifest['capture'] == template['capture']
    if view == 'macro':
        assert all(capture['job'][key] == value for key, value in manifest['capture'].items())
    output = production / 'experiments' / experiment_id / 'experiment.json'
    assert not output.exists(), f'Refusing overwrite: {output}'
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2))
    records.append({'view': view, 'manifest': str(output), 'sha256': sha(output),
                    'source_sha256': manifest['source_sha256'],
                    'capture_gate': 'macro-first; full and second-light await useful macro visual result'})

handoff = {
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'experiment': 'porcelain_glaze_layer_assignment',
    'build_receipt': str(root / 'scenes/build_receipt.json'),
    'build_receipt_sha256': sha(root / 'scenes/build_receipt.json'),
    'resources': resources,
    'source_only': True, 'rendered': False, 'selected': False, 'quality_acceptance': False,
    'first_capture': records[0], 'later_views': records[1:],
    'matched_macro_control': {'source': scenes['main']['control'],
                             'source_sha256': scenes['main']['control_sha256'],
                             'image': str(control_capture), 'image_sha256': sha(control_capture),
                             'receipt': str(control_receipt), 'receipt_sha256': sha(control_receipt),
                             'capture_profile_bit_equal': True,
                             'noise_acceptance': False, 'visual_acceptance': False},
    'protected_historical_appearance': reference['baseline'],
    'protected_reference_receipt': str(reference_path),
    'declared_optical_changes': scenes['main']['field_roles'],
    'limitation': 'Effective opaque layered closure; not explicit refracted paths through a finite glass volume or calibrated material optics.',
    'renderer_instruction': 'Reuse the exact matched macro control if its frozen source/profile/pipeline receipts validate. Capture the main candidate macro first. Compare actual pixels before admitting full-panel or second-light captures. Do not promote or publish.',
}
path = root / 'frozen_handoff.json'
assert not path.exists()
path.write_text(json.dumps(handoff, indent=2))
print(json.dumps({'handoff': str(path), 'sha256': sha(path), 'macro_manifest': records[0]}, indent=2))
