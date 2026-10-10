"""Vectorized exact same finite material-cell law; no material coarse-graining."""
import numpy as np
from quadratic_prism import positive_small_energy

class BatchPrisms:
 def __init__(self,cells):
  self.cells=tuple(cells);self.count=len(cells)
  if not self.count or len({len(c.weights) for c in cells})!=1:raise ValueError('A common quadrature rule is required')
  self.grad=np.stack([c.grad for c in cells]);self.weights=np.stack([c.weights for c in cells]);self.mu=np.array([c.mu for c in cells])[:,None];self.lam=np.array([c.lam for c in cells])[:,None];self.mass=np.array([c.mass for c in cells])
 def evaluate(self,nodes):
  x=np.asarray(nodes,dtype=float).reshape(self.count,18,3);F=np.einsum('cia,cqib->cqab',x-x.mean(axis=1)[:,None,:],self.grad);J=np.linalg.det(F)
  if np.any(J<=0) or not np.isfinite(J).all():raise ValueError('Invalid current material-cell volume')
  delta=np.linalg.svd(F,compute_uv=False)-1;lnJ=np.log1p(delta).sum(-1);psi=.5*self.mu*np.sum(delta*delta+2*positive_small_energy(delta),axis=-1)+.5*self.lam*lnJ*lnJ;FiT=np.linalg.inv(F).swapaxes(-1,-2);P=self.mu[:,:,None,None]*(F-FiT)+(self.lam*lnJ)[:,:,None,None]*FiT
  energies=np.sum(self.weights*psi,axis=1);gradient=np.einsum('cqib,cqab,cq->cia',self.grad,P,self.weights);volumes=np.sum(self.weights*J,axis=1)
  return float(energies.sum()),gradient,{'cell_energy_J':energies,'cell_volume_m3':volumes,'minimum_det_F':float(J.min()),'cell_material_mass_kg':self.mass.copy()}
