"""Reconstruct pinned R5D inputs/geometry in a new directory; never render."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify_source():
    lock = json.loads((ROOT / 'SOURCE_LOCK.json').read_text())
    for row in lock['files']:
        path = ROOT / row['path']
        if not path.is_file() or digest(path) != row['sha256']:
            raise ValueError('Pinned source differs: ' + row['path'])
    return lock

def verify_inputs(directory):
    import numpy as np
    directory = Path(directory)
    expected = json.loads((ROOT / 'expected_inputs.json').read_text())
    if digest(directory / 'prototypes/fracture_solids.json') != expected['fracture_solids.json']['sha256']:
        raise ValueError('Fracture prototype bytes differ from the approved input')
    with np.load(directory / 'prototypes/microbed.npz', allow_pickle=False) as data:
        rows = expected['microbed.npz']['arrays']
        if set(data.files) != set(rows):
            raise ValueError('Microbed array identities differ')
        for name, row in rows.items():
            array = data[name]
            if list(array.shape) != row['shape'] or str(array.dtype) != row['dtype']:
                raise ValueError('Microbed array contract differs: ' + name)
            if hashlib.sha256(array.tobytes()).hexdigest() != row['sha256']:
                raise ValueError('Microbed array content differs: ' + name)
    return {'fracture_prototype_bytes_match': True, 'microbed_array_contents_match': True,
            'npz_container_byte_parity_claimed': False}

def reconstruct(out, blender='blender', prepare_only=False):
    verify_source()
    out = Path(out).resolve()
    if out.exists() or out.is_relative_to(ROOT):
        raise ValueError('Use a new output directory outside the pinned source bundle')
    out.mkdir(parents=True)
    for name in ['src', 'prototypes', 'preview', 'receipts', 'deliverables']:
        (out / name).mkdir()
    for name in ['fracture_prototypes.py', 'compacted_microbed.py', 'build_diagnostic.py']:
        shutil.copyfile(ROOT / 'src' / name, out / 'src' / name)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    # Convex-hull face/vertex ordering can differ across hosts. The accepted
    # normalized fracture geometry is an explicit text source recipe.
    expected = json.loads((ROOT / 'expected_inputs.json').read_text())
    prototype = expected['fracture_solids.json']['canonical_json'].encode('utf-8')
    if hashlib.sha256(prototype).hexdigest() != expected['fracture_solids.json']['sha256']:
        raise ValueError('Pinned fracture geometry source is corrupt')
    (out / 'prototypes/fracture_solids.json').write_bytes(prototype)
    subprocess.run([sys.executable, str(out / 'src/compacted_microbed.py')], check=True, env=env)
    result = verify_inputs(out)
    result.update(version='R5D', fracture_prototypes_mode='pinned normalized geometry source',
                  geometry_rebuilt=False, render_performed=False,
                  selected_scene_binary_replacement=False)
    if not prepare_only:
        subprocess.run([blender, '-b', '-t', '2', '--python-exit-code', '1',
                        '--python', str(out / 'src/build_diagnostic.py')], check=True, env=env)
        scene = out / 'preview/diagnostic_d.blend'
        result.update(geometry_rebuilt=True, scene_sha256=digest(scene))
    (out / 'reconstruction_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    parser.add_argument('--blender', default='blender')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.out, args.blender, args.prepare_only), indent=2))
