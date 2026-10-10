"""Validate final input hashes and package only the verified reference files."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    reference = json.loads((ROOT/'receipts/reference_receipt.json').read_text())
    tests = json.loads((ROOT/'receipts/test_receipt.json').read_text())
    assert reference['status'] == 'PASS_NUMERICAL_REFERENCE_ONLY'
    assert reference['resource_gates_passed'] and tests['passed']
    assert tests['tests_run'] == 6 and not tests['failures'] and not tests['errors']
    gates = reference['protocol']['gates']
    assert tests['wall_seconds'] < gates['process_wall_seconds']
    assert tests['maximum_resident_memory_KiB'] < gates['process_maximum_resident_KiB']
    for key, expected in reference['source_sha256'].items():
        path = (ROOT.parent/key) if key.startswith('bark-constitutive-upgrade/') else ROOT/key
        assert digest(path) == expected, f'Stale receipt: {key}'
    delivered = [ROOT/'README.md', ROOT/'SOURCES.md', ROOT/'PROTOCOL.json', ROOT/'run_checks.sh',
                 *sorted((ROOT/'source').glob('*.py')), *sorted((ROOT/'tests').glob('*.py')),
                 ROOT/'receipts/reference_receipt.json', ROOT/'receipts/test_receipt.json',
                 ROOT/'receipts/reference.log', ROOT/'receipts/tests.log']
    manifest = {'status': 'VERIFIED_NUMERICAL_REFERENCE_ONLY',
                'all_recorded_source_hashes_match': True,
                'inherited_evidence_hashes_match': True,
                'sha256': {str(p.relative_to(ROOT)): digest(p) for p in delivered},
                'limits': reference['limitations']}
    (ROOT/'receipts/manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('Verified all source/evidence hashes and both passing receipts; wrote manifest.')


if __name__ == '__main__':
    main()
