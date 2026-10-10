"""Small same-directory atomic numerical checkpoints."""
from pathlib import Path
import json,os
import numpy as np

def save_npz(path,**arrays):
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,**arrays);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)

def save_json(path,value):
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value,stream,indent=2);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)
