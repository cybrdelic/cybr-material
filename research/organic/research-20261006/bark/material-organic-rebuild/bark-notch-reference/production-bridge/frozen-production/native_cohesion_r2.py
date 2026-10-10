"""Validated native cohesive state; same constitutive law as the frozen wrapper.

Material arrays are private owned copies. History restoration is checked against
its irreversible bilinear law before any native call. This changes validation,
not the traction law or the fracture-energy calibration.
"""
import numpy as np
from native_cohesion import NativeCohesion as _FrozenNativeCohesion

class NativeCohesion(_FrozenNativeCohesion):
 def __init__(self,cells,operators):
  if isinstance(cells,(bool,np.bool_)) or not np.isscalar(cells) or not np.isfinite(cells) or int(cells)!=cells or cells<1:
   raise ValueError('Cell count must be a positive integer')
  p=np.asarray(operators.get('pairs'))
  if p.ndim!=2 or p.shape[1]!=2 or p.shape[0]<1 or not np.issubdtype(p.dtype,np.number) or not np.isfinite(p).all() or np.any(p!=np.floor(p)) or p.min()<0 or p.max()>=cells:
   raise ValueError('Cohesive cell indices must be exact integers in range')
  N=p.shape[0];owned={}
  for key in ('pairs','J','T1','T2','area','kt','kn','normal_sign','strength','Gc'):
   a=np.array(operators.get(key),dtype=np.int32 if key=='pairs' else np.float64,copy=True,order='C')
   expected=(N,2) if key=='pairs' else ((N,36) if key in ('J','T1','T2') else (N,))
   if a.shape!=expected or not np.isfinite(a).all():raise ValueError('Invalid cohesive '+key+' shape or finite-value contract')
   if key in ('area','kt','kn','strength','Gc') and np.any(a<=0):raise ValueError('Cohesive '+key+' must be positive')
   owned[key]=a
  if not np.isin(owned['normal_sign'],(-1.,1.)).all():raise ValueError('Cohesive orientation must be +1 or -1')
  # Gradient coefficients must annihilate translations. The C++ evaluator also
  # centers coordinates, but rejecting a malformed operator avoids relying on
  # that centering to silently change the supplied interpolation contract.
  for key in ('J','T1','T2'):
   a=owned[key];tol=2e-13*np.maximum(1.,np.sum(np.abs(a),axis=1))
   if np.any(np.abs(a.sum(axis=1))>tol):raise ValueError('Cohesive '+key+' does not annihilate translations')
  super().__init__(int(cells),owned)
  for a in self.ops.values():a.flags.writeable=False
  self.strength.flags.writeable=False;self.Gc.flags.writeable=False
  self.delta0.flags.writeable=False;self.deltaf.flags.writeable=False
 def _law(self,k):
  d=np.zeros_like(k);soft=(k>self.delta0)&(k<self.deltaf);d[k>=self.deltaf]=1.
  d[soft]=1-self.delta0[soft]*(self.deltaf[soft]-k[soft])/(k[soft]*(self.deltaf[soft]-self.delta0[soft]))
  return d
 def _check_state(self):
  for key in ('damage','maximum_opening'):
   a=getattr(self,key)
   if not isinstance(a,np.ndarray) or a.dtype!=np.float64 or a.shape!=(self.N,) or not a.flags.c_contiguous or not np.isfinite(a).all():raise ValueError('Invalid cohesive history array '+key)
  if np.any(self.maximum_opening<0) or np.any(self.damage<0) or np.any(self.damage>1):raise ValueError('Cohesive history is outside its physical domain')
  if not np.allclose(self.damage,self._law(self.maximum_opening),rtol=0.,atol=2e-12):raise ValueError('Cohesive damage is inconsistent with retained maximum opening')
 def restore(self,maximum_opening,damage):
  k=np.array(maximum_opening,dtype=np.float64,copy=True,order='C');d=np.array(damage,dtype=np.float64,copy=True,order='C')
  old=(self.maximum_opening,self.damage);self.maximum_opening=k;self.damage=d
  try:self._check_state()
  except Exception:self.maximum_opening,self.damage=old;raise
 def evaluate(self,q):
  self._check_state();return super().evaluate(q)
 def evolve(self,opening):
  self._check_state();a=np.asarray(opening,dtype=np.float64)
  if a.shape!=(self.N,) or not np.isfinite(a).all() or np.any(a<0):raise ValueError('Opening must be a finite nonnegative per-point vector')
  self.maximum_opening=np.maximum(self.maximum_opening,a);self.damage=self._law(self.maximum_opening)
 def dissipation(self):
  self._check_state();return super().dissipation()
