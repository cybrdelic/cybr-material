#!/usr/bin/env python3
"""Check preparation receipt against numeric plan and immutable file hashes. No render."""
import argparse
import json
import math
from pathlib import Path
from grazing_preview import sha256


def equivalent(actual, expected):
    if isinstance(expected, dict):
        return isinstance(actual, dict) and all(k in actual and equivalent(actual[k], v) for k, v in expected.items())
    if isinstance(expected, (tuple, list)):
        return isinstance(actual, (tuple, list)) and len(actual) == len(expected) and all(equivalent(a,b) for a,b in zip(actual, expected))
    if isinstance(expected, (float, int)) and not isinstance(expected, bool):
        return math.isclose(actual, expected, rel_tol=2e-6, abs_tol=2e-7)
    return actual == expected


def verify(path):
    receipt = json.loads(Path(path).read_text())
    actual, plan = receipt['inspection'], receipt['plan']
    expected = {'camera_location_m': plan['camera']['location_m'],
                'camera_rotation_euler': plan['camera']['rotation_euler'],
                'camera_matrix_world': plan['camera']['matrix_world'],
                'camera_ortho_scale_m': plan['camera']['ortho_scale'],
                'key_location_m': plan['key']['location_m'],
                'key_energy_w': plan['key']['energy'], 'key_diameter_m': plan['key']['size'],
                'key_color': plan['key']['color'], 'world_rgba': plan['world']['color_rgba'],
                'world_strength': plan['world']['strength'], 'cycles': plan['cycles'],
                'display': {k:plan['display'][k] for k in ('view_transform','look','exposure','gamma')},
                'png': {'depth':'16','mode':'RGB','dither':0.}, 'rendered':False}
    errors = [key for key, value in expected.items() if not equivalent(actual.get(key), value)]
    if errors:
        raise ValueError('Saved derivative settings mismatch: ' + ', '.join(errors))
    if any(actual['other_light_energies'].values()):
        raise ValueError('Other explicit lights are still emitting')
    if sha256(receipt['immutable_source']) != receipt['source_sha256_before']:
        raise ValueError('Immutable source hash mismatch')
    if sha256(receipt['derivative']) != receipt['derivative_sha256']:
        raise ValueError('Derivative hash mismatch')
    for path, digest in receipt['protected_files_sha256'].items():
        if sha256(path) != digest:
            raise ValueError('Protected file hash mismatch: ' + path)
    return {'valid': True, 'preset_id': plan['preset_id'], 'saved_settings_match_plan': True,
            'source_unchanged': True, 'protected_files_unchanged': True, 'rendered': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt')
    args = parser.parse_args()
    print(json.dumps(verify(args.receipt), indent=2))
