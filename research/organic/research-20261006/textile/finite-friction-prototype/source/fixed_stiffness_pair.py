"""Explicit alternative regularization: k_pair = k_density * dS0_i dS0_j.
Tangential stored energy is independent of current normal force. The Coulomb
limit alone is mu*N. This is an uncalibrated constitutive probe, not a silent
replacement of the older N/delta0 law or an integrated spinning solver.
"""
from dataclasses import dataclass
import numpy as np
from objective_pair import pose,exp,log,align
@dataclass(frozen=True)
class FixedState:
 labels:tuple
 pi:np.ndarray
 pj:np.ndarray
 Ri:np.ndarray
 Rj:np.ndarray
 normal:float
 mu:float
 stiffness_density_N_m3:float
 reference_pair_measure_m2:float
 elastic:np.ndarray
 heat:float=0.
 numerical_loss:float=0.
 @property
 def stiffness_N_m(self):return self.stiffness_density_N_m3*self.reference_pair_measure_m2

def make(labels,pi,pj,Ri,Rj,*,normal,mu=.3,stiffness_density_N_m3=2.5e11,reference_pair_measure_m2=4e-10,elastic=None):
 pi,Ri=pose(pi,Ri);pj,Rj=pose(pj,Rj);k=stiffness_density_N_m3*reference_pair_measure_m2
 if not isinstance(labels,tuple) or len(labels)!=2 or labels[0]==labels[1] or normal<0 or mu<0 or stiffness_density_N_m3<=0 or reference_pair_measure_m2<=0 or not np.isfinite([normal,mu,k]).all() or np.linalg.norm(pi-pj)<1e-15:raise ValueError('Invalid persistent material contact')
 n=(pi-pj)/np.linalg.norm(pi-pj);e=np.zeros(3) if elastic is None else np.asarray(elastic,dtype=float).copy()
 if e.shape!=(3,) or not np.isfinite(e).all() or abs(e@n)>1e-15 or np.linalg.norm(e)*k>mu*normal*(1+1e-12):raise ValueError('Invalid accepted tangential spring')
 return FixedState(labels,pi,pj,Ri,Rj,float(normal),float(mu),float(stiffness_density_N_m3),float(reference_pair_measure_m2),e)

def advance(old,pi,pj,Ri,Rj,*,labels,normal):
 if labels!=old.labels:raise ValueError('Material labels changed; anchor reassignment forbidden')
 if not np.isfinite(normal) or normal<0:raise ValueError('Normal force must be finite and compressive')
 pi,Ri=pose(pi,Ri);pj,Rj=pose(pj,Rj);d0=old.pi-old.pj;d1=pi-pj
 if np.linalg.norm(d1)<1e-15:raise ValueError('Degenerate normal')
 n0=d0/np.linalg.norm(d0);n1=d1/np.linalg.norm(d1);Qi=Ri@old.Ri.T;Qj=Rj@old.Rj.T;relative=log(Qi.T@Qj)
 if np.linalg.norm(relative)>np.pi-1e-7:raise ValueError('Ambiguous material spin needs substepping')
 Qmean=Qi@exp(.5*relative);Q=align(Qmean@n0,n1)@Qmean;midpoint=.5*(old.pi+old.pj);ai=midpoint-old.pi;aj=midpoint-old.pj;motion=(pi-pj)+Qi@ai-Qj@aj;du=motion-n1*(n1@motion);e0=Q@old.elastic;trial=e0+du;k=old.stiffness_N_m;cap=old.mu*normal;length=np.linalg.norm(trial);ratio=min(1.,cap/(k*length)) if length else 1.;e=ratio*trial;dp=trial-e;traction=k*e;stored0=.5*k*(old.elastic@old.elastic);stored1=.5*k*(e@e);heat=float(cap*np.linalg.norm(dp));loss=float(.5*k*np.sum((e-e0)**2));work=float(traction@du);residual=work-(stored1-stored0)-heat-loss
 new=FixedState(old.labels,pi,pj,Ri,Rj,float(normal),old.mu,old.stiffness_density_N_m3,old.reference_pair_measure_m2,e,old.heat+heat,old.numerical_loss+loss)
 return new,dict(surface_increment_m=du,traction_N=traction,plastic_increment_m=dp,stored_J=float(stored1),stored_change_J=float(stored1-stored0),friction_heat_J=heat,numerical_loss_J=loss,endpoint_work_J=work,work_balance_residual_J=float(residual),contact_active=normal>0,normal_force_N=float(normal),stiffness_N_m=k)
