from pathlib import Path
import ctypes,numpy as np
from batch_volume import BatchPrisms
class NativePrisms(BatchPrisms):
 def __init__(self,cells):
  super().__init__(cells);self.grad=np.ascontiguousarray(self.grad);self.weights=np.ascontiguousarray(self.weights);self.mu=np.ascontiguousarray(self.mu.ravel());self.lam=np.ascontiguousarray(self.lam.ravel());self.lib=ctypes.CDLL(str(Path(__file__).with_name('libvolume_native.so')));self.fun=self.lib.volume_batch;P=ctypes.POINTER(ctypes.c_double);self.fun.argtypes=[ctypes.c_int,ctypes.c_int]+[P]*9;self.fun.restype=ctypes.c_int
 def evaluate(self,nodes):
  x=np.ascontiguousarray(nodes,dtype=np.float64)
  if x.shape!=(self.count,18,3) or not np.isfinite(x).all():raise ValueError('Invalid complete material-cell array')
  grad=np.zeros_like(x);energy=np.zeros(self.count);vol=np.zeros(self.count);minimum=np.zeros(1);P=ctypes.POINTER(ctypes.c_double);ptr=lambda a:a.ctypes.data_as(P)
  rc=self.fun(self.count,self.weights.shape[1],ptr(x),ptr(self.grad),ptr(self.weights),ptr(self.mu),ptr(self.lam),ptr(grad),ptr(energy),ptr(vol),ptr(minimum))
  if rc:raise ValueError('Invalid material cell geometry or input, native status '+str(rc))
  if not all(np.isfinite(a).all() for a in (grad,energy,vol,minimum)):raise FloatingPointError('Nonfinite native constitutive result')
  return float(energy.sum()),grad,{'cell_energy_J':energy,'cell_volume_m3':vol,'minimum_det_F':float(minimum[0]),'cell_material_mass_kg':self.mass.copy()}
