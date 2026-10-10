"""Origin-relative coordinates with stable constitutive displacement gradients.

No rest metric is reset. H0 stores the retained reference deformation relative
 to each cell's birth geometry. The chosen displacement origin is a geometrical
 state whose intended cohesive faces coincide; it need not be stress free.
"""
from pathlib import Path
import ctypes,numpy as np
from ring_system_r3 import RingSystem
P=ctypes.POINTER(ctypes.c_double)
def ptr(a):return a.ctypes.data_as(P)
class IncrementalRing(RingSystem):
 def __init__(self,sectors=8,stem_increment_m=0.):
  super().__init__(sectors,stem_increment_m);self.origin=self.q0.copy();natural=np.stack([c.X for c in self.cells]);u0=self.origin-natural
  self.H0=np.ascontiguousarray(np.einsum('cia,cqib->cqab',u0-u0.mean(axis=1,keepdims=True),self.bulk.grad))
  o=self.bonds.ops;x=np.stack([np.vstack((self.origin[i],self.origin[j])) for i,j in o['pairs']]);xc=x-x[:,:1,:];initial_gap=np.einsum('pi,pia->pa',o['J'],xc)
  if np.max(np.linalg.norm(initial_gap,axis=1))>1e-14:raise ValueError('Displacement origin requires intended coincident cohesive faces')
  self.reference_interface_gap_roundoff_m=float(np.max(np.linalg.norm(initial_gap,axis=1)));self.tref=np.ascontiguousarray(np.column_stack((np.einsum('pi,pia->pa',o['T1'],xc),np.einsum('pi,pia->pa',o['T2'],xc))))
  self.vlib=ctypes.CDLL(str(Path(__file__).with_name('libvolume_incremental.so')));self.vfun=self.vlib.volume_incremental;self.vfun.argtypes=[ctypes.c_int,ctypes.c_int]+[P]*10;self.vfun.restype=ctypes.c_int
  self.clib=ctypes.CDLL(str(Path(__file__).with_name('libcohesion_incremental.so')));self.cfun=self.clib.cohesive_incremental;self.cfun.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]+[P]*13;self.cfun.restype=ctypes.c_int
 def evaluate_displacement(self,displacement):
  u=np.ascontiguousarray(displacement,dtype=np.float64)
  if u.shape!=self.q.shape or not np.isfinite(u).all():raise ValueError('Invalid origin-relative displacement field')
  self.bonds._check_state();g=np.zeros_like(u);E=np.zeros(self.N);V=np.zeros(self.N);J=np.zeros(1)
  rc=self.vfun(self.N,self.bulk.weights.shape[1],ptr(u),ptr(self.H0),ptr(self.bulk.grad),ptr(self.bulk.weights),ptr(self.bulk.mu),ptr(self.bulk.lam),ptr(g),ptr(E),ptr(V),ptr(J))
  if rc:raise ValueError('Incremental volume status '+str(rc))
  gc=np.zeros_like(u);Ec=np.zeros(1);opening=np.zeros(self.bonds.N);o=self.bonds.ops
  rc=self.cfun(self.N,self.bonds.N,o['pairs'].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ptr(u),ptr(self.tref),ptr(o['J']),ptr(o['T1']),ptr(o['T2']),ptr(o['area']),ptr(o['kt']),ptr(o['kn']),ptr(self.bonds.damage),ptr(o['normal_sign']),ptr(Ec),ptr(gc),ptr(opening))
  if rc or not all(np.isfinite(a).all() for a in (g,E,V,J,gc,Ec,opening)):raise ValueError('Invalid incremental force/energy result; cohesive status '+str(rc))
  # Retain the original potential datum. Gravity is affine, so no full q array
  # is formed to evaluate its very small physical displacement contribution.
  energy=float(E.sum()+Ec[0]-np.sum(self.gravity*self.origin)-np.sum(self.gravity*u))
  return energy,g+gc-self.gravity,opening,{'cell_energy_J':E,'cell_volume_m3':V,'minimum_det_F':float(J[0]),'cell_material_mass_kg':self.bulk.mass.copy()}
