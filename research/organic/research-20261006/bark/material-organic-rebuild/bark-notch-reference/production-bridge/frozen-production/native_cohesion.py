from pathlib import Path
import numpy as np,ctypes
class NativeCohesion:
 def __init__(self,cells,operators):
  self.cells=int(cells);self.ops={k:np.ascontiguousarray(operators[k],dtype=np.int32 if k=='pairs' else np.float64) for k in ('pairs','J','T1','T2','area','kt','kn','normal_sign')};self.N=len(self.ops['area']);self.damage=np.zeros(self.N);self.maximum_opening=np.zeros(self.N);self.strength=np.asarray(operators['strength']);self.Gc=np.asarray(operators['Gc']);self.delta0=self.strength/self.ops['kt'];self.deltaf=2*self.Gc/self.strength
  if self.cells<1 or self.N<1 or self.ops['pairs'].shape!=(self.N,2) or self.ops['pairs'].min()<0 or self.ops['pairs'].max()>=self.cells or any(self.ops[k].shape!=(self.N,36) for k in ('J','T1','T2')) or any(self.ops[k].shape!=(self.N,) for k in ('area','kt','kn','normal_sign')):raise ValueError('Invalid complete cohesive operator shapes')
  if any(not np.isfinite(a).all() for a in self.ops.values()) or any(np.any(self.ops[k]<=0) for k in ('area','kt','kn')) or not np.isin(self.ops['normal_sign'],[-1.,1.]).all() or np.any(self.deltaf<=self.delta0):raise ValueError('Invalid cohesive operator material values')
  self.lib=ctypes.CDLL(str(Path(__file__).with_name('libcohesion_native.so')));self.fun=self.lib.cohesive_batch;P=ctypes.POINTER(ctypes.c_double);self.fun.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]+[P]*12;self.fun.restype=ctypes.c_int
 def evaluate(self,q):
  q=np.ascontiguousarray(q,dtype=np.float64)
  if q.shape!=(self.cells,18,3) or not np.isfinite(q).all():raise ValueError('Invalid material position field')
  g=np.zeros_like(q);E=np.zeros(1);opening=np.zeros(self.N);P=ctypes.POINTER(ctypes.c_double);ptr=lambda a:a.ctypes.data_as(P);o=self.ops
  code=self.fun(self.cells,self.N,o['pairs'].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ptr(q),ptr(o['J']),ptr(o['T1']),ptr(o['T2']),ptr(o['area']),ptr(o['kt']),ptr(o['kn']),ptr(self.damage),ptr(o['normal_sign']),ptr(E),ptr(g),ptr(opening))
  if code or not all(np.isfinite(a).all() for a in (E,g,opening)):raise ValueError('Cohesive evaluation failed: '+str(code))
  return float(E[0]),g,opening
 def evolve(self,opening):
  self.maximum_opening=np.maximum(self.maximum_opening,opening);k=self.maximum_opening;d=np.where(k<=self.delta0,0.,np.where(k>=self.deltaf,1.,1-self.delta0*(self.deltaf-k)/(np.maximum(k,1e-30)*(self.deltaf-self.delta0))));self.damage=np.maximum(self.damage,d)
 def dissipation(self):
  k=np.minimum(self.maximum_opening,self.deltaf);work=np.where(k<=self.delta0,.5*self.ops['kt']*k*k,.5*self.strength*self.delta0+self.strength*(k-self.delta0)-.5*self.strength*(k-self.delta0)**2/(self.deltaf-self.delta0));return float(self.ops['area']@np.maximum(work-.5*(1-self.damage)*self.ops['kt']*k*k,0.))
