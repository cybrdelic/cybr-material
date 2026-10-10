"""Finite-volume material cell for growth-state research: P2 triangle x P2 thickness.
18 geometric nodes,54 positional DOFs. Through-thickness deformation is free,
so normal contraction is represented rather than hidden in a unit director.
Compressible neo-Hookean law is an uncalibrated benchmark, not a cork law.
Reference constitutive form: https://docs.fenicsproject.org/dolfinx/main/cpp/demos/demo_hyperelasticity.html
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

def shapes(r,s,z):
 L=np.array([1-r-s,r,s]);dL=np.array([[-1.,1.,0.],[-1.,0.,1.]])
 H=np.r_[L*(2*L-1),4*L[0]*L[1],4*L[1]*L[2],4*L[2]*L[0]];dH=np.zeros((2,6));dH[:,:3]=dL*(4*L-1)
 for j,(a,b) in enumerate(((0,1),(1,2),(2,0)),3):dH[:,j]=4*(dL[:,a]*L[b]+L[a]*dL[:,b])
 Z=np.array([.5*z*(z-1),1-z*z,.5*z*(z+1)]);dZ=np.array([z-.5,-2*z,z+.5]);N=np.outer(Z,H).ravel();dN=np.column_stack((np.outer(Z,dH[0]).ravel(),np.outer(Z,dH[1]).ravel(),np.outer(dZ,H).ravel()));return N,dN

def flat_reference(xy,thickness):
 xy=np.asarray(xy);six=np.vstack((xy,(xy[0]+xy[1])/2,(xy[1]+xy[2])/2,(xy[2]+xy[0])/2));return np.concatenate([np.column_stack((six,np.full(6,z))) for z in (-thickness/2,0,thickness/2)])

def positive_small_energy(delta):
 # delta-log(1+delta), evaluated without small-strain cancellation.
 out=delta-np.log1p(delta);small=abs(delta)<.001;x=delta[small];v=np.zeros_like(x)
 for k in range(2,11):v+=((-1)**k)*x**k/k
 out[small]=v;return out

class QuadraticPrism:
 def __init__(self,reference_nodes,E,nu,density,nq=5):
  self.X=np.asarray(reference_nodes,dtype=float).copy()
  if self.X.shape!=(18,3) or not np.isfinite(self.X).all() or not np.isfinite([E,nu,density]).all() or E<=0 or density<=0 or not -1<nu<.5:raise ValueError('Invalid material cell input')
  self.E=float(E);self.nu=float(nu);self.rho=float(density);self.mu=E/(2*(1+nu));self.lam=E*nu/((1+nu)*(1-2*nu));u,w=leggauss(nq);u01=(u+1)/2;w01=w/2;Ns=[];grads=[];weights=[]
  for r,wr in zip(u01,w01):
   for v,wv in zip(u01,w01):
    s=(1-r)*v
    for z,wz in zip(u,w):
     N,dN=shapes(r,s,z);J=(self.X-self.X.mean(axis=0)).T@dN;det=np.linalg.det(J)
     if det<=0:raise ValueError('Inverted reference material cell')
     material_grad=dN@np.linalg.inv(J);material_grad-=material_grad.mean(axis=0)
     Ns.append(N);grads.append(material_grad);weights.append(wr*wv*(1-r)*wz*det)
  self.N=np.asarray(Ns);self.grad=np.asarray(grads);self.weights=np.asarray(weights);self.reference_volume=float(self.weights.sum());self.mass=self.rho*self.reference_volume;self.mass_scalar=self.rho*np.einsum('qi,qj,q->ij',self.N,self.N,self.weights)
  self.X.setflags(write=False)
 def evaluate(self,nodes):
  x=np.asarray(nodes,dtype=float).reshape(18,3);F=np.einsum('ia,qib->qab',x-x.mean(axis=0),self.grad);J=np.linalg.det(F)
  if np.any(J<=0) or not np.isfinite(J).all():raise ValueError('Inverted or nonfinite material deformation')
  stretch=np.linalg.svd(F,compute_uv=False);delta=stretch-1;logs=np.log1p(delta);lnJ=logs.sum(axis=1);psi=.5*self.mu*np.sum(delta*delta+2*positive_small_energy(delta),axis=1)+.5*self.lam*lnJ*lnJ
  FiT=np.linalg.inv(F).transpose(0,2,1);P=self.mu*(F-FiT)+(self.lam*lnJ)[:,None,None]*FiT;energy=float(self.weights@psi);gradient=np.einsum('qib,qab,q->ia',self.grad,P,self.weights);volume=float(self.weights@J)
  return energy,gradient,{'mapped_volume_m3':volume,'minimum_det_F':float(J.min()),'maximum_det_F':float(J.max()),'material_mass_kg':self.mass,'mapped_mean_density_kg_m3':self.mass/volume}
