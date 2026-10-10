"""Crash-between-files regression without launching another mechanical solve."""
from pathlib import Path
import json,tempfile
import numpy as np
from checkpoint_io import write_checkpoint,read_checkpoint
from model import ROOT

checks=[]
with tempfile.TemporaryDirectory() as name:
    path=Path(name)/'restart.npz';side=Path(name)/'restart.json'
    u=np.arange(30,dtype=float).reshape(10,3)*1e-6
    original=dict(parameters={'identity':'frozen'},trace=[{'station':3,'work_J':.123}],status='running')
    sha=write_checkpoint(path,dict(displacement=u),original)
    side.write_text(json.dumps(dict(status='stale',trace=[{'station':2}])))
    a,m=read_checkpoint(path,side)
    assert m==original and np.array_equal(a['displacement'],u)
    checks.append('Embedded state and ledger survive stale companion JSON')
    # An interrupted later write cannot replace the accepted checkpoint.
    path.with_name('restart.tmp.npz').write_bytes(b'incomplete write')
    a,m=read_checkpoint(path,side)
    assert m==original and np.array_equal(a['displacement'],u)
    checks.append('Incomplete temporary checkpoint leaves accepted generation intact')
    legacy=Path(name)/'legacy.npz';np.savez_compressed(legacy,displacement=u)
    rejected=False
    try:read_checkpoint(legacy,side)
    except (KeyError,ValueError):rejected=True
    assert rejected;checks.append('Mixed legacy state/sidecar pair cannot resume')
    states=[]
    for label in ('coarse16','coarse16_half','medium16_half','fine16_half','medium32_half'):
        a,m=read_checkpoint(ROOT/'state'/label/'restart.npz',ROOT/'state'/label/'restart.json')
        assert np.isfinite(a['displacement']).all() and m['trace'][-1]['station']>0
        states.append(dict(label=label,station=m['trace'][-1]['station'],load=m['trace'][-1]['load'],status=m['status']))
    checks.append('All five retained legacy states verify against their own sidecars')
result=dict(status='PASS',checks=checks,retained_states=states,
    scope='Checkpoint-only correctness patch after the mechanical solves. The archived original runner reproduces their exact source hashes; the new runner has identical mechanics and requires its own new label/source identity.')
(ROOT/'receipts/checkpoint_integrity.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
