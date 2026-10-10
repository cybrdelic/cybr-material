"""Objective finite-surface cohesive point with conservative closure contact.
All stiffnesses are traction/separation (N/m³), area is reference bond area m².
Reference-normal approximations are avoided: tangent-derived current normals
and their full position derivatives participate in the force gradient.
"""
import numpy as np

def skew(v):return np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])

def evaluate(q,Jjump,Jt1,Jt2,area,kt,kn,damage,normal_sign=1.):
 if not np.isfinite([area,kt,kn,damage,normal_sign]).all() or area<=0 or kt<=0 or kn<=0 or not 0<=damage<=1 or normal_sign not in (-1.,1.):raise ValueError('Invalid cohesive point parameters or orientation sign')
 q=np.asarray(q,dtype=float);Jjump=np.asarray(Jjump,dtype=float);Jt1=np.asarray(Jt1,dtype=float);Jt2=np.asarray(Jt2,dtype=float)
 if q.ndim!=1 or not np.isfinite(q).all() or any(J.shape!=(3,len(q)) or not np.isfinite(J).all() for J in (Jjump,Jt1,Jt2)):raise ValueError('Invalid cohesive state/operator contract')
 jump=Jjump@q;t1=Jt1@q;t2=Jt2@q;v=np.cross(t1,t2);length=np.linalg.norm(v)
 if length<1e-18:raise ValueError('Singular current cohesive surface')
 n=normal_sign*v/length;Jn=normal_sign*(np.eye(3)-np.outer(n,n))/length@(-skew(t2)@Jt1+skew(t1)@Jt2)
 dn=float(jump@n);Jdn=n@Jjump+jump@Jn;closed=min(dn,0.);w=1-damage
 energy=.5*area*(w*kt*(jump@jump)+w*(kn-kt)*dn*dn+damage*kn*closed*closed)
 grad=area*(w*kt*(Jjump.T@jump)+(w*(kn-kt)*dn+damage*kn*closed)*Jdn)
 eq2=float(jump@jump-dn*dn+(kn/kt)*max(dn,0.)**2)
 return float(energy),grad,{'normal':n,'normal_gap_m':dn,'equivalent_opening_m':np.sqrt(max(eq2,0))}
