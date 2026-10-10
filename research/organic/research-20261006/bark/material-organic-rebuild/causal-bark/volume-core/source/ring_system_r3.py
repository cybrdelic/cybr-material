"""Post-checkpoint ring with complete moving-support mass/force accounting.

Constrained positions are supplied by a separate physical loading history.
Damping is the declared mass-proportional drag alpha*M*v, with dimensions 1/s.
The kinetic metric and material mass retain all free/constrained cross terms.
"""
import numpy as np
from ring_system_r2 import RingSystem as _RingSystem
from native_cohesion_r2 import NativeCohesion

class RingSystem(_RingSystem):
 def __init__(self,sectors=8,stem_increment_m=0.):
  super().__init__(sectors,stem_increment_m)
  self.bonds=NativeCohesion(self.N,self.ops)
 def mass_product(self,x):
  x=np.asarray(x,dtype=float)
  if x.shape!=self.q.shape or not np.isfinite(x).all():raise ValueError('Invalid full kinetic vector')
  return np.einsum('cij,cja->cia',self.M,x)
 def moving_force(self,gradient,constrained_velocity,constrained_acceleration,alpha):
  """Return F in vdot_free=F-alpha*v_free.

  Input boundary arrays contain only constrained entries. No mass term can be
  discarded merely because a positional degree of freedom is prescribed.
  """
  g=np.asarray(gradient,dtype=float);vc=np.asarray(constrained_velocity,dtype=float);ac=np.asarray(constrained_acceleration,dtype=float)
  if any(a.shape!=self.q.shape or not np.isfinite(a).all() for a in (g,vc,ac)) or not np.isfinite(alpha) or alpha<0:raise ValueError('Invalid moving-boundary state')
  if np.any(vc[self.free_mask]!=0) or np.any(ac[self.free_mask]!=0):raise ValueError('Boundary arrays must vanish on free coordinates')
  return -np.einsum('cij,cja->cia',self.Minv,g+self.mass_product(ac+alpha*vc))
 def moving_acceleration_reaction(self,gradient,velocity,constrained_acceleration,alpha):
  vc=np.zeros_like(velocity);vc[~self.free_mask]=velocity[~self.free_mask]
  F=self.moving_force(gradient,vc,constrained_acceleration,alpha)
  a=F-alpha*np.where(self.free_mask,velocity,0.)+constrained_acceleration
  residual=self.mass_product(a)+alpha*self.mass_product(velocity)+gradient
  return a,residual
