"""Finite-pose material-pair kinematics with the existing Coulomb return map.
Kinematic constitutive probe only: constant normal load, fixed material labels,
no moving quadrature, rod equilibrium, changing-load or contact-birth claim.
"""
from dataclasses import dataclass
from pathlib import Path
import sys,numpy as np
from scipy.spatial.transform import Rotation
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P.parent/'compact-friction/source'))
from friction_compact import return_map

def exp(w):return Rotation.from_rotvec(w).as_matrix()
def log(R):return Rotation.from_matrix(R).as_rotvec()
def align(a,b):
 v=np.cross(a,b);s=np.linalg.norm(v);c=float(np.clip(a@b,-1,1))
 if s<1e-14:
  if c<0:raise ValueError('Antipodal contact normal requires substepping')
  return np.eye(3)
 return exp(v*(np.arctan2(s,c)/s))
def pose(p,R):
 p=np.asarray(p,dtype=float);R=np.asarray(R,dtype=float)
 if p.shape!=(3,) or R.shape!=(3,3) or not np.isfinite(p).all() or not np.isfinite(R).all() or np.linalg.norm(R.T@R-np.eye(3))>1e-10 or np.linalg.det(R)<0:raise ValueError('Invalid material pose')
 return p.copy(),R.copy()
@dataclass(frozen=True)
class State:
 labels:tuple
 pi:np.ndarray
 pj:np.ndarray
 Ri:np.ndarray
 Rj:np.ndarray
 normal:float
 mu:float
 delta0:float
 elastic:np.ndarray
 heat:float=0.
 algorithmic:float=0.

def make(labels,pi,pj,Ri,Rj,normal=1e-6,mu=.3,delta0=1e-8,elastic=None):
 pi,Ri=pose(pi,Ri);pj,Rj=pose(pj,Rj)
 if not isinstance(labels,tuple) or len(labels)!=2 or labels[0]==labels[1] or normal<=0 or mu<0 or delta0<=0 or not np.isfinite([normal,mu,delta0]).all() or np.linalg.norm(pi-pj)<1e-15:raise ValueError('Invalid labelled contact')
 n=(pi-pj)/np.linalg.norm(pi-pj);e=np.zeros(3) if elastic is None else np.asarray(elastic,dtype=float).copy()
 if e.shape!=(3,) or not np.isfinite(e).all() or abs(e@n)>1e-15 or np.linalg.norm(e)>mu*delta0*(1+1e-12):raise ValueError('Invalid accepted tangential spring')
 return State(labels,pi,pj,Ri,Rj,float(normal),float(mu),float(delta0),e)

def advance(old,pi,pj,Ri,Rj,*,labels,normal):
 if labels!=old.labels:raise ValueError('Material labels changed; no anchor reassignment')
 if normal!=old.normal:raise NotImplementedError('Changing normal load needs its separate stiffness/normal-work coupling gate')
 pi,Ri=pose(pi,Ri);pj,Rj=pose(pj,Rj);d0=old.pi-old.pj;d1=pi-pj
 if np.linalg.norm(d1)<1e-15:raise ValueError('Degenerate material-pair normal')
 n0=d0/np.linalg.norm(d0);n1=d1/np.linalg.norm(d1);Qi=Ri@old.Ri.T;Qj=Rj@old.Rj.T;relative=log(Qi.T@Qj)
 if np.linalg.norm(relative)>np.pi-1e-7:raise ValueError('Ambiguous mean spin requires substepping')
 Qmean=Qi@exp(.5*relative);Q=align(Qmean@n0,n1)@Qmean
 # Advect the old coincident midpoint surface samples with each actual body.
 # Centre translations and finite rotations both contribute. Common rigid
 # motion produces exactly coincident advected points, including about n.
 midpoint=.5*(old.pi+old.pj);ai=midpoint-old.pi;aj=midpoint-old.pj
 relative_motion=(pi-pj)+Qi@ai-Qj@aj;du=relative_motion-n1*(n1@relative_motion);e0=Q@old.elastic;trial=e0+du
 rm=return_map(trial[None],np.zeros((1,3)),np.array([normal]),old.mu,old.delta0);e=rm['elastic'][0];dp=rm['plastic'][0];traction=rm['traction'][0];k=normal/old.delta0;stored0=.5*k*(old.elastic@old.elastic);stored1=.5*k*(e@e);heat=rm['heat_J'];A=.5*k*np.sum((e-e0)**2);work=float(traction@du);balance=work-(stored1-stored0)-heat-A
 new=State(old.labels,pi,pj,Ri,Rj,old.normal,old.mu,old.delta0,e,old.heat+heat,old.algorithmic+A)
 return new,dict(surface_increment_m=du,transport=Q,traction_N=traction,plastic_increment_m=dp,stored_J=stored1,stored_change_J=stored1-stored0,friction_heat_J=heat,endpoint_discrete_work_J=work,algorithmic_energy_J=float(A),work_balance_residual_J=float(balance),sliding=bool(rm['slip'][0]))
