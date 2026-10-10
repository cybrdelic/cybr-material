"""Unique-projection moving-contact gate for two prescribed rigid fibre spans.
Fixed primary reference quadrature; partner material coordinate is state, not key.
This one-sided line-contact energy is not the prior double material-pair law.
"""
import numpy as np
import hashlib
from scipy.spatial.transform import Rotation
from scipy.optimize import brentq
L=.0024;RADIUS=16e-6;KN=1e8;KT_LINE=2.5e9;MU=.3;EVENT_EPS=1e-15

def exp(v):return Rotation.from_rotvec(v).as_matrix()
def log(R):return Rotation.from_matrix(R).as_rotvec()
def mv(A,b):return np.einsum('...ij,...j->...i',A,b)
def skew(v):
 v=np.asarray(v);a=np.zeros(v.shape[:-1]+(3,3));a[...,0,1]=-v[...,2];a[...,0,2]=v[...,1];a[...,1,0]=v[...,2];a[...,1,2]=-v[...,0];a[...,2,0]=-v[...,1];a[...,2,1]=v[...,0];return a

def quadrature(order=9,cells=128):
 x,w=np.polynomial.legendre.leggauss(order);edges=np.linspace(-L/2,L/2,cells+1);s=(.5*(edges[:-1,None]+edges[1:,None])+.5*(L/cells)*x).ravel();weights=np.tile(w*.5*L/cells,cells);labels=np.c_[np.zeros(len(s),int),np.full(len(s),cells),np.full(len(s),order),np.repeat(np.arange(cells),order),np.tile(np.arange(order),cells),np.ones(len(s),int)];return s,weights,labels

def controls(t,case):
 if case=='rolling':
  R=2*RADIUS-.25e-6;Rd=0.;theta=.8*t;td=.8;C=L/2+5e-6-.5*R*theta;Cd=-.5*R*td
 elif case=='sliding':
  R=2*RADIUS-.25e-6;Rd=0.;theta=0.;td=0.;C=L/2-5e-6+20e-6*t;Cd=20e-6
 elif case=='birth_loss':
  R=2*RADIUS+.5e-6*np.cos(2*np.pi*t);Rd=-np.pi*1e-6*np.sin(2*np.pi*t);theta=.4*t;td=.4;C=L/2-1e-6-.5*td*(2*RADIUS*t+.5e-6*np.sin(2*np.pi*t)/(2*np.pi))+10e-6*t;Cd=-.5*td*R+10e-6
 else:raise ValueError(case)
 return R,Rd,theta,td,C,Cd

def geometry(s,weights,t,case,superpose=False):
 R,Rd,theta,td,C,Cd=controls(t,case);Ry=exp([0,theta,0]);n=Ry@np.array([0.,0.,1.]);axis=Ry@np.array([1.,0.,0.]);nodes_s=np.array([0.,L/2,L]);nodes=-R*n+(nodes_s-C)[:,None]*axis;pa=np.c_[np.zeros(len(s)),s,np.zeros(len(s))]
 G=exp(t*np.array([.7,-.4,.6])) if superpose else np.eye(3);omegaG=np.array([.7,-.4,.6]) if superpose else np.zeros(3);Tdot=np.array([.002,-.001,.0007]) if superpose else np.zeros(3);T=t*Tdot
 pa=pa@G.T+T;nodes=nodes@G.T+T;v=nodes[1:]-nodes[:-1];rel=pa[:,None,:]-nodes[:-1][None];u=np.clip(np.sum(rel*v[None],axis=2)/np.sum(v*v,axis=1)[None],0,1);feet=nodes[:-1][None]+u[:,:,None]*v[None];dist2=np.sum((pa[:,None]-feet)**2,axis=2);segment=np.argmin(dist2,axis=1);ids=np.arange(len(s));pb=feet[ids,segment];sb=nodes_s[segment]+u[ids,segment]*(L/2);difference=pa-pb;distance=np.linalg.norm(difference,axis=1);normal=difference/distance[:,None];penetration=np.maximum(2*RADIUS-distance,0.);penetration[np.abs(2*RADIUS-distance)<EVENT_EPS]=0.;N=KN*weights*penetration;Un=.5*KN*np.sum(weights*penetration**2)
 # Partial time derivative at the CURRENT material coordinate sb, not the
 # derivative of the geometric closest point (which includes Cdot advection).
 vb0=-(Rd+(sb-C)*td)[:,None]*n-(R*td+Cd)*axis;pb0=-R*n+(sb-C)[:,None]*axis;va=np.cross(omegaG,pa-T)+Tdot;vb=vb0@G.T+np.cross(omegaG,pb0@G.T)+Tdot;wa=omegaG;wb=omegaG+G@np.array([0.,td,0.]);mid=.5*(pa+pb);surface_relative=va+np.cross(wa,mid-pa)-vb-np.cross(wb,mid-pb);vt=surface_relative-normal*np.sum(surface_relative*normal,axis=1)[:,None]
 # Deliberately wrong control: derivative of closest GEOMETRIC position.
 vb_closest0=-Rd*n-R*td*axis;vb_wrong=np.broadcast_to(vb_closest0,(len(s),3))@G.T+np.cross(omegaG,pb0@G.T)+Tdot;wrong=va+np.cross(wa,mid-pa)-vb_wrong-np.cross(wb,mid-pb);wrong-=normal*np.sum(wrong*normal,axis=1)[:,None]
 return dict(t=t,partner_nodes=nodes,pa=pa,pb=pb,sb=sb,segment=segment,n=normal,N=N,Un=float(Un),Ra=G,Rb=G@Ry,va=va,vb=vb,vt=vt,wrong_vt=wrong)

def transport(a,b):
 Qi=b['Ra']@a['Ra'].T;Qj=b['Rb']@a['Rb'].T;Qmean=Qi@exp(.5*log(Qi.T@Qj));pred=a['n']@Qmean.T;cross=np.cross(pred,b['n']);c=np.sum(pred*b['n'],axis=1)
 if np.min(c)<-.99:raise ValueError('Normal turn requires substepping')
 K=skew(cross);Qcorr=np.eye(3)+K+(K@K)/(1+c)[:,None,None];return Qcorr@Qmean

def times(steps,s,case):
 values=list(np.linspace(0,1,steps+1));values.append(.5)
 if case=='birth_loss':
  for y in s:
   if abs(y)>=2*RADIUS:continue
   threshold=np.sqrt((2*RADIUS)**2-y*y);c=(threshold-2*RADIUS)/(.5e-6)
   if -1<c<1:
    t=np.arccos(c)/(2*np.pi);values.extend([t,1-t])
 if (controls(0,case)[4]-L/2)*(controls(1,case)[4]-L/2)<0:values.append(brentq(lambda t:controls(t,case)[4]-L/2,0,1,xtol=1e-14))
 return np.unique(values)

def run(case,steps=128,order=9,superpose=False):
 s,w,labels=quadrature(order);grid=times(steps,s,case);k=KT_LINE*w;g=geometry(s,w,0.,case,superpose);e=np.zeros((len(s),3));active=g['N']>0;births=active.astype(int);epochs=births.copy();changes=0;resets_at_segment_crossing=0
 if case=='rolling':e[:,0]=.2*MU*g['N']/k
 if superpose:e=e@g['Ra'].T
 U0=g['Un']+.5*np.sum(k*np.sum(e*e,axis=1));Ut0=.5*np.sum(k*np.sum(e*e,axis=1));Wn=Wt=D=A=0.;max_local=0.;slip=np.zeros(len(s));wrong_slip=np.zeros(len(s));trace=[];loading_error=None
 for t1 in grid[1:]:
  t0=g['t'];dt=t1-t0;mid=geometry(s,w,.5*(t0+t1),case,superpose);new=geometry(s,w,t1,case,superpose);Q=transport(g,new);Qmid=transport(mid,new);e0=mv(Q,e);du=dt*mv(Qmid,mid['vt']);trial=e0+du;length=np.linalg.norm(trial,axis=1);cap=MU*new['N'];ratio=np.minimum(1.,np.divide(cap,k*length,out=np.ones_like(length),where=length>0));enew=ratio[:,None]*trial;dp=trial-enew;tau=k[:,None]*enew;Ut=.5*np.sum(k*np.sum(enew*enew,axis=1));oldUt=.5*np.sum(k*np.sum(e*e,axis=1));heat=float(np.sum(cap*np.linalg.norm(dp,axis=1)));loss=float(.5*np.sum(k*np.sum((enew-e0)**2,axis=1)));work=float(np.sum(tau*du));Wnormal=float(-dt*np.sum(mid['N']*np.sum(mid['n']*(mid['va']-mid['vb']),axis=1)));max_local=max(max_local,abs(work-(Ut-oldUt)-heat-loss));Wn+=Wnormal;Wt+=work;D+=heat;A+=loss;slip+=np.linalg.norm(du,axis=1)*(mid['N']>0);wrong_slip+=dt*np.linalg.norm(mid['wrong_vt'],axis=1)*(mid['N']>0)
  newactive=new['N']>0;born=(~active)&newactive;epochs+=born;births+=born;cross=(new['segment']!=g['segment'])&active&newactive;changes+=int(cross.sum());resets_at_segment_crossing+=int(np.sum(cross&born));active=newactive;e=enew;g=new
  if abs(t1-.5)<1e-12 and case=='birth_loss':loading_error=abs(Wn-g['Un'])/max(g['Un'],1e-30)
  if np.any(cross):trace.append(dict(time=float(t1),counterpart_arc_m=float(g['sb'][len(s)//2]),active_crossings=int(cross.sum()),total_spring_J=Ut,epochs_changed_at_boundary=int(np.sum(cross&born))))
 U=g['Un']+.5*np.sum(k*np.sum(e*e,axis=1));balance=Wn+Wt-(U-U0)-D-A;active_ever=births>0;ww=w*active_ever;mean=lambda x:float(np.sum(ww*x)/max(np.sum(ww),1e-30))
 return dict(case=case,steps=steps,actual_intervals=len(grid)-1,order=order,primary_samples=len(s),material_identity_schema='(primary fibre, immutable reference rule, cell, Gauss-node index, partner fibre); closest coordinate and segment are state only',material_identity_sha256=hashlib.sha256(labels.tobytes()+s.tobytes()+w.tobytes()).hexdigest(),maximum_births_per_label=int(births.max()),total_births=int(births.sum()),active_segment_boundary_crossings=changes,history_resets_at_segment_boundary=resets_at_segment_crossing,mean_integrated_material_slip_m=mean(slip),mean_wrong_geometric_point_slip_m=mean(wrong_slip),normal_work_J=Wn,tangential_endpoint_work_J=Wt,initial_energy_J=U0,final_energy_J=U,initial_tangential_energy_J=Ut0,final_tangential_energy_J=float(.5*np.sum(k*np.sum(e*e,axis=1))),friction_heat_J=D,numerical_loss_J=A,numerical_loss_fraction=A/max(abs(Wt),D,Ut0,1e-30),work_balance_residual_J=balance,relative_work_balance_residual=abs(balance)/max(abs(Wn)+abs(Wt),D,U0,U,1e-30),maximum_local_tangential_balance_error_J=max_local,normal_loading_work_relative_error=loading_error,final_active_contacts=int(active.sum()),boundary_trace=trace,final_counterpart_arc_m=g['sb'].tolist(),final_elastic_m=e.tolist())
