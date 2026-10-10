"""Closed radial-growth ring, finite-volume body and objective fracture bonds."""
import ctypes,numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from ring_geometry import build
from ring_interfaces import build_interfaces
from native_volume import NativePrisms
from native_cohesion import NativeCohesion
class RingSystem:
 def __init__(self,sectors=8,stem_increment_m=.006):
  self.history,self.cells,self.q,self.metadata,self.triangles=build(sectors,stem_increment_m=stem_increment_m);self.q0=self.q.copy();self.sectors=sectors;self.N=len(self.cells);self.bulk=NativePrisms(self.cells);self.ops=build_interfaces(self.cells,self.q,self.metadata,sectors);self.bonds=NativeCohesion(self.N,self.ops);self.free_mask=np.ones_like(self.q,dtype=bool);self.free_mask[:sectors*2,:6]=False;self.free=np.flatnonzero(self.free_mask.ravel());self.fixed=np.flatnonzero(~self.free_mask.ravel());self.M=np.stack([c.mass_scalar for c in self.cells]);self.Minv=np.zeros_like(self.M);self.Mhalf=np.zeros_like(self.M)
  for c in range(self.N):
   ids=np.arange(6,18) if c<sectors*2 else np.arange(18);a=self.M[c][np.ix_(ids,ids)];e,V=np.linalg.eigh(a);self.Minv[c][np.ix_(ids,ids)]=np.linalg.inv(a);self.Mhalf[c][np.ix_(ids,ids)]=(V*(1/np.sqrt(e)))@V.T
  self.gravity=np.zeros_like(self.q);self.gravity[:,:,1]=-9.81*self.M.sum(axis=2)
  G=self.bulk.grad;W=self.bulk.weights;A=np.einsum('cqia,cqjb,cq->cijab',G,G,W);mu=self.bulk.mu[:,0];la=self.bulk.lam[:,0];K=la[:,None,None,None,None]*A+mu[:,None,None,None,None]*(A.swapaxes(-1,-2)+np.einsum('cij,ab->cijab',np.trace(A,axis1=-2,axis2=-1),np.eye(3)));self.K=K.transpose(0,1,3,2,4).reshape(self.N,54,54)
  P=ctypes.POINTER(ctypes.c_double);self.hvp=self.bonds.lib.cohesive_linear_hvp;self.hvp.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]+[P]*5;self.hvp.restype=ctypes.c_int
 def evaluate(self,q=None):
  q=self.q if q is None else q;E,g,report=self.bulk.evaluate(q);Ec,gc,opening=self.bonds.evaluate(q);return E+Ec-float((self.gravity*q).sum()),g+gc-self.gravity,opening,report
 def acceleration(self,gradient):return -np.einsum('cij,cja->cia',self.Minv,gradient)
 def kinetic(self,v):return float(.5*np.einsum('cia,cij,cja->',v,self.M,v))
 def linear_operator(self,v):
  x=np.ascontiguousarray(v.reshape(self.N,18,3));out=np.zeros_like(x);P=ctypes.POINTER(ctypes.c_double);ptr=lambda a:a.ctypes.data_as(P);o=self.ops
  rc=self.hvp(self.N,self.bonds.N,self.bonds.ops['pairs'].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ptr(x),ptr(self.bonds.ops['J']),ptr(self.bonds.ops['area']),ptr(self.bonds.ops['kt']),ptr(out))
  if rc:raise ValueError('Invalid linear cohesive operation')
  return np.einsum('cij,cj->ci',self.K,x.reshape(self.N,54)).reshape(self.N,18,3)+out
 def maximum_reference_frequency(self):
  def mat(v):
   x=np.zeros(self.q.size);x[self.free]=v;x=np.einsum('cij,cja->cia',self.Mhalf,x.reshape(self.q.shape));a=self.linear_operator(x);a=np.einsum('cij,cja->cia',self.Mhalf,a);return a.ravel()[self.free]
  A=LinearOperator((len(self.free),len(self.free)),matvec=mat,dtype=float);eig=eigsh(A,k=1,which='LM',tol=1e-8,return_eigenvectors=False,v0=np.random.default_rng(825).normal(size=len(self.free)))[0];return float(np.sqrt(eig))
