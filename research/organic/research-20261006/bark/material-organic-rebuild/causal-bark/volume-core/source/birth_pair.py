"""Two material-born curved cells coupled by an objective cohesive interface.
A bounded assembly test. Its coefficients are inherited audit inputs; no weak
annual-ring adhesion or calibrated cork identity is inferred from a birth date.
"""
from pathlib import Path
import sys,numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'accretion/source'))
from phellem_history import AccretionHistory
from birth_cell import make_cohort_cell
from quadratic_prism import shapes
from surface_cohesion import evaluate as cohesive_evaluate

def softening(k,k0,kf):
 return np.where(k<=k0,0.,np.where(k>=kf,1.,1-k0*(kf-k)/(np.maximum(k,1e-30)*(kf-k0))))

class BirthPair:
 def __init__(self,cohesive_order=4):
  self.history=AccretionHistory(.105,1.,.1);self.history.advance(1,new_phellem_thickness_m=.0003,new_density_kg_m3=240,identity='old');self.history.advance(2,.001,new_phellem_thickness_m=.00025,new_density_kg_m3=240,identity='young');triangle=np.array([[0.,0.],[.02,0.],[.002,.006]])
  young,qy,my=make_cohort_cell(self.history,'young',triangle,2.8e6,.28);old,qo,mo=make_cohort_cell(self.history,'old',triangle,2.8e6,.28);self.cells=[young,old];self.material_metadata=[my,mo];self.q0=np.r_[qy.ravel(),qo.ravel()];self.q=self.q0.copy();self.free=np.arange(18,108);self.fixed=np.arange(18)
  self.M=np.zeros((108,108))
  for l,cell in enumerate(self.cells):self.M[l*54:(l+1)*54,l*54:(l+1)*54]=np.kron(cell.mass_scalar,np.eye(3))
  self.Mf=self.M[np.ix_(self.free,self.free)];self.Minv=np.linalg.inv(self.Mf);self.kt=self.kn=2.8e6/(.01*.002);self.strength=100000.;self.Gc=10.;self.delta0=self.strength/self.kt;self.deltaf=2*self.Gc/self.strength;self.points=[];u,w=leggauss(cohesive_order);u=(u+1)/2;w=w/2
  for r,wr in zip(u,w):
   for v,wv in zip(u,w):
    s=(1-r)*v;Ny,dNy=shapes(r,s,1.);No,dNo=shapes(r,s,-1.);J=np.zeros((3,108));T1=J.copy();T2=J.copy()
    for offset,N,dN,sign in ((0,Ny,dNy,-1),(54,No,dNo,1)):
     for k in range(18):J[:,offset+3*k:offset+3*k+3]=sign*N[k]*np.eye(3);T1[:,offset+3*k:offset+3*k+3]=.5*dN[k,0]*np.eye(3);T2[:,offset+3*k:offset+3*k+3]=.5*dN[k,1]*np.eye(3)
    area=wr*wv*(1-r)*np.linalg.norm(np.cross(T1@self.q0,T2@self.q0));self.points.append((J,T1,T2,area))
  self.damage=np.zeros(len(self.points));self.maximum_opening=np.zeros(len(self.points))
 def evaluate(self,q=None):
  q=self.q if q is None else q;E=0.;gradient=np.zeros(108);volumes=[]
  for l,cell in enumerate(self.cells):
   e,g,v=cell.evaluate(q[l*54:(l+1)*54]);E+=e;gradient[l*54:(l+1)*54]+=g.ravel();volumes.append(v)
  openings=[]
  for i,(J,T1,T2,a) in enumerate(self.points):
   e,g,d=cohesive_evaluate(q,J,T1,T2,a,self.kt,self.kn,self.damage[i]);E+=e;gradient+=g;openings.append(d['equivalent_opening_m'])
  return float(E),gradient,np.asarray(openings),volumes
 def evolve_damage(self,openings):
  self.maximum_opening=np.maximum(self.maximum_opening,openings);self.damage=np.maximum(self.damage,softening(self.maximum_opening,self.delta0,self.deltaf))
 def dissipation(self):
  k=np.minimum(self.maximum_opening,self.deltaf);work=np.where(k<=self.delta0,.5*self.kt*k*k,.5*self.strength*self.delta0+self.strength*(k-self.delta0)-.5*self.strength*(k-self.delta0)**2/(self.deltaf-self.delta0));loss=np.maximum(work-.5*(1-self.damage)*self.kt*k*k,0.);return float(loss@np.array([p[3] for p in self.points]))
