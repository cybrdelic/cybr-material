import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from incremental_envelope import IncrementalEnvelope
import audit_candidate
original_factory=audit_candidate.synthetic_candidate

def incremental_factory(axes=None):
    old=original_factory(axes)
    return IncrementalEnvelope(old.tension,old.compression,old.G,old.axes,old.minimum_stretch,old.maximum_stretch)

audit_candidate.synthetic_candidate=incremental_factory
audit_candidate.RECEIPTS=Path(__file__).resolve().parents[1]/'receipts/full_prism'
audit_candidate.RECEIPTS.mkdir(parents=True,exist_ok=True)
audit_candidate.main()
