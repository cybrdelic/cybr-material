"""Explicit native dependency closure for incremental ring experiments."""
from pathlib import Path
import hashlib
from run_ring_r2 import hashes as original_hashes
R=Path(__file__).resolve().parent
NATIVE=('batch_volume.cpp','libvolume_native.so','cohesion_batch.cpp','libcohesion_native.so','batch_volume_incremental.cpp','libvolume_incremental.so','cohesion_incremental.cpp','libcohesion_incremental.so')
def hashes():
 out=original_hashes()
 for name in NATIVE:out['volume-core/source/'+name]=hashlib.sha256((R/name).read_bytes()).hexdigest()
 return out
