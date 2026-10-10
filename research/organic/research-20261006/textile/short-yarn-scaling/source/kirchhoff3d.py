"""Objective, material-labelled, inextensible/unshearable discrete director rods.

R columns are material directors in world coordinates; edge = ds0*R[:,2].
Bend/twist curvature at an internal vertex is log(R_i.T R_{i+1})/dual_ds.
Energy = 1/2 sum dual_ds*(kappa-kappa0).B.(kappa-kappa0).
The discrete constitutive law is geometrically objective, with linear elastic
material moments; it is not an exact continuum solution or a manufacturing law.
Boundary half-segments have no curvature cell, as in the usual vertex/Voronoi
rod discretization. Continuum comparisons must include boundary refinement.
Whole fixed material fibres are represented. There is no mass-only feed reserve.
"""
from dataclasses import dataclass
import numpy as np
from scipy.spatial.transform import Rotation

def skew(v):
 v=np.asarray(v);a=np.zeros(v.shape[:-1]+(3,3));a[...,0,1]=-v[...,2];a[...,0,2]=v[...,1];a[...,1,0]=v[...,2];a[...,1,2]=-v[...,0];a[...,2,0]=-v[...,1];a[...,2,1]=v[...,0];return a

def exp(v):
 v=np.asarray(v,dtype=float);return Rotation.from_rotvec(v.reshape(-1,3)).as_matrix().reshape(v.shape[:-1]+(3,3))

def log(R):
 R=np.asarray(R,dtype=float);return Rotation.from_matrix(R.reshape(-1,3,3)).as_rotvec().reshape(R.shape[:-2]+(3,))

def right_jacobian(v,inverse=False):
 v=np.asarray(v);t2=np.sum(v*v,axis=-1);small=t2<1e-8;t=np.sqrt(np.maximum(t2,1e-30));K=skew(v);K2=K@K;I=np.broadcast_to(np.eye(3),K.shape)
 if inverse:
  a=np.where(small,1/12+t2/720+t2*t2/30240,(1-.5*t/np.tan(.5*t))/np.maximum(t2,1e-30));return I+.5*K+a[...,None,None]*K2
 a=np.where(small,.5-t2/24+t2*t2/720,(1-np.cos(t))/np.maximum(t2,1e-30));b=np.where(small,1/6-t2/120+t2*t2/5040,(t-np.sin(t))/np.maximum(t*t2,1e-30));return I-a[...,None,None]*K+b[...,None,None]*K2

def mv(A,v):return np.einsum('...ij,...j->...i',A,v)

@dataclass(frozen=True)
class FiberMaterial:
 radius_m:float=1.597930424803307e-5
 young_modulus_Pa:float=1e9
 poisson_ratio:float=.30 # new torsional assumption; no measured wool claim
 density_kg_m3:float=1310.
 @property
 def B(self):
  I=np.pi*self.radius_m**4/4;G=self.young_modulus_Pa/(2*(1+self.poisson_ratio));return np.array([self.young_modulus_Pa*I,self.young_modulus_Pa*I,G*2*I])

class KirchhoffRod:
 def __init__(self,material_arc_m,material=FiberMaterial(),intrinsic_curvature_per_m=None):
  self.s=np.asarray(material_arc_m,dtype=float).copy();self.mat=material
  if self.s.ndim!=1 or not 3<=self.s.size<=1026 or not np.isfinite(self.s).all() or np.any(np.diff(self.s)<=0):raise ValueError('Strictly increasing finite material arc labels required.')
  if not np.isfinite([material.radius_m,material.young_modulus_Pa,material.poisson_ratio,material.density_kg_m3]).all() or min(material.radius_m,material.young_modulus_Pa,material.density_kg_m3)<=0 or not -1<material.poisson_ratio<.5:raise ValueError('Invalid circular-fibre material.')
  self.ds=np.diff(self.s);self.dual=.5*(self.ds[:-1]+self.ds[1:]);self.M=len(self.ds);self.s.setflags(write=False);self.ds.setflags(write=False);self.dual.setflags(write=False)
  self.k0=np.zeros((self.M-1,3)) if intrinsic_curvature_per_m is None else np.array(intrinsic_curvature_per_m,dtype=float)
  if self.k0.shape!=(self.M-1,3) or not np.isfinite(self.k0).all():raise ValueError('Intrinsic curvature must be labelled at each material vertex.')
  self.k0.setflags(write=False)
 def validate(self,R):
  R=np.asarray(R,dtype=float)
  if R.shape!=(self.M,3,3) or not np.isfinite(R).all() or np.max(abs(np.swapaxes(R,-1,-2)@R-np.eye(3)))>1e-9 or np.min(np.linalg.det(R))<1-1e-9:raise ValueError('Each material frame must be a proper SO(3) rotation.')
  return R
 def points(self,R,origin=(0.,0.,0.)):
  R=self.validate(R);origin=np.asarray(origin,dtype=float)
  if origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('Invalid physical root position.')
  return np.vstack([origin,origin+np.cumsum(self.ds[:,None]*R[:,:,2],axis=0)])
 def boundary_terms(self,R,end_frames=(None,None),intrinsic=None):
  R=self.validate(R);intrinsic=np.zeros((2,3)) if intrinsic is None else np.asarray(intrinsic,dtype=float)
  if intrinsic.shape!=(2,3) or not np.isfinite(intrinsic).all():raise ValueError("Two material-labelled boundary intrinsic curvatures required.")
  E=0.;g=np.zeros((self.M,3));external=np.zeros((2,3))
  for end,Q in enumerate(end_frames):
   if Q is None:continue
   Q=np.asarray(Q,dtype=float)
   if Q.shape!=(3,3) or not np.isfinite(Q).all() or np.max(abs(Q.T@Q-np.eye(3)))>1e-9 or np.linalg.det(Q)<1-1e-9:raise ValueError("Invalid endpoint director.")
   i=0 if end==0 else self.M-1;h=self.ds[i]/2;A=Q if end==0 else R[i];B=R[i] if end==0 else Q;rel=A.T@B;w=log(rel)
   if np.linalg.norm(w)>=np.pi-1e-5:raise ValueError("Endpoint half-cell rotation needs refinement.")
   d=w/h-intrinsic[end];tau=self.mat.B*d;v=right_jacobian(w,inverse=True).T@tau;left=-rel@v;right=v;E+=.5*h*float(d@tau)
   g[i]+=right if end==0 else left;external[end]=left if end==0 else right
  return E,g,external
 def energy_gradient(self,R,end_frames=(None,None),boundary_intrinsic=None):
  """Return energy and body torque-gradient for right variations R Exp(dtheta)."""
  R=self.validate(R);rel=np.swapaxes(R[:-1],-1,-2)@R[1:];w=log(rel)
  if np.max(np.linalg.norm(w,axis=1))>=np.pi-1e-5:raise ValueError('Adjacent-frame log branch requires spatial refinement.')
  k=w/self.dual[:,None];d=k-self.k0;tau=d*self.mat.B[None,:];E=.5*np.sum(self.dual[:,None]*d*tau)
  # d log(Q) = Jr(log(Q))^-1 (dtheta_next - Q.T dtheta_prev).
  v=mv(np.swapaxes(right_jacobian(w,inverse=True),-1,-2),tau);g=np.zeros((self.M,3));g[:-1]-=mv(rel,v);g[1:]+=v
  Eb,gb,_=self.boundary_terms(R,end_frames,boundary_intrinsic);return float(E+Eb),g+gb
 def pullback_position_gradient(self,R,gp):
  """Exact chain rule for positions integrated from fixed-length directors."""
  R=self.validate(R);gp=np.asarray(gp,dtype=float)
  if gp.shape!=(self.M+1,3) or not np.isfinite(gp).all():raise ValueError('Nodal force-gradient shape mismatch.')
  cumulative=np.cumsum(gp[::-1],axis=0)[::-1][1:];body=mv(np.swapaxes(R,-1,-2),cumulative);return self.ds[:,None]*np.cross(np.array([0.,0.,1.]),body)
 def diagnostics(self,R):
  p=self.points(R);length=np.linalg.norm(np.diff(p,axis=0),axis=1);w=log(np.swapaxes(R[:-1],-1,-2)@R[1:]);vol=np.pi*self.mat.radius_m**2*float(self.ds.sum())
  return dict(material_arc_length_m=float(self.ds.sum()),current_polygon_arc_length_m=float(length.sum()),maximum_relative_length_error=float(np.max(abs(length-self.ds)/self.ds)),reference_and_current_volume_m3=vol,mass_kg=vol*self.mat.density_kg_m3,max_adjacent_rotation_rad=float(np.linalg.norm(w,axis=1).max()),max_radius_curvature=float((self.mat.radius_m*np.linalg.norm(w[:,:2],axis=1)/self.dual).max()),curvature_integration_length_m=float(self.dual.sum()),intrinsic_curvature_labels='immutable material vertex arc coordinates',axial_feed='none; whole fibre remains represented')
 def coordinates_energy_gradient(self,w,origin=(0.,0.,0.),nodal_load_N=None):
  w=np.asarray(w,dtype=float).reshape(self.M,3);R=exp(w);E,g=self.energy_gradient(R)
  if nodal_load_N is not None:
   F=np.asarray(nodal_load_N,dtype=float);p=self.points(R,origin)
   if F.shape!=p.shape or not np.isfinite(F).all():raise ValueError('Invalid dead-load nodal array.')
   E-=np.sum(F*p);g+=self.pullback_position_gradient(R,-F)
  return float(E),mv(np.swapaxes(right_jacobian(w),-1,-2),g).ravel()
