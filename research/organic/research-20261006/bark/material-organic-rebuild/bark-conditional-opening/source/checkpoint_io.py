"""One authoritative atomic checkpoint; JSON sidecars are only readable copies."""
from pathlib import Path
import hashlib,json,os
import numpy as np


def write_checkpoint(path,arrays,metadata):
    path=Path(path);temp=path.with_name(path.stem+'.tmp.npz')
    np.savez_compressed(temp,metadata_json=np.array(json.dumps(metadata)),**arrays)
    with temp.open('rb') as handle:os.fsync(handle.fileno())
    os.replace(temp,path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_checkpoint(path,legacy_sidecar=None):
    path=Path(path)
    with np.load(path,allow_pickle=False) as data:
        arrays={key:data[key].copy() for key in data.files if key!='metadata_json'}
        if 'metadata_json' in data.files:
            metadata=json.loads(str(data['metadata_json']))
        else:
            if legacy_sidecar is None:raise ValueError('Legacy checkpoint requires its matching sidecar')
            metadata=json.loads(Path(legacy_sidecar).read_text())
            if hashlib.sha256(path.read_bytes()).hexdigest()!=metadata['state_sha256']:
                raise ValueError('Legacy checkpoint and sidecar disagree; do not resume a mixed pair')
    return arrays,metadata
