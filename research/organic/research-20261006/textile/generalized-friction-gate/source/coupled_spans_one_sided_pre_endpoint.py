"""Twelve generalized rigid-span DOFs with the original double normal energy.
Friction history is immutable during trials. Its full residual tangent may be
nonsymmetric. Fixed reference quadrature and a unique connected component only.
"""
from pathlib import Path
import sys,hashlib,numpy as np
from scipy.spatial.transform import Rotation
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P.parent/'suspended-spinning/source'))
from kirchhoff3d import right_jacobian
from segment_contact_sine import evaluate as native_pair
L=.0024;r=16e-6;RS=2*r;S=RS;KC=1e8;KT=2.5e9;MU=.3;UNIT=1e-4*S

def exp(v):return Rotation.from_rotvec(v).as_matrix()
def log(R):return Rotation.from_matrix(R).as_rotvec()
def mv(A,b):return np.einsum('...ij,...j->...i',A,b)
def skew(v):
 a=np.zeros(v.shape[:-1]+(3,3));a[...,0,1]=-v[...,2];a[...,0,2]=v[...,1];a[...,1,0]=v[...,2];a[...,1,2]=-v[...,0];a[...,2,0]=-v[...,1];a[...,2,1]=v[...,0];return a
class Spans:
 def __init__(self,primary_order=65,partner_order=32):
  x,w=np.polynomial.legendre.leggauss(primary_order);edge=np.linspace(-L/2,L/2,129);self.s=(.5*(edge[:-1,None]+edge[1:,None])+.5*L/128*x).ravel();self.w=np.tile(w*.5*L/128,128);self.k=KT*self.w;self.localA=np.c_[np.zeros(len(self.s)),self.s,np.zeros(len(self.s))];self.gx,self.gw=np.polynomial.legendre.leggauss(partner_order);self.identity=hashlib.sha256(self.s.tobytes()+self.w.tobytes()).hexdigest();self.order=primary_order
 def geometry(self,q,normal=True):
  Ta=q[:3]*S;Ra=exp(q[3:6]);Tb=q[6:9]*S;Rb=exp(q[9:12]);axis=Rb[:,0];pa=self.localA@Ra.T+Ta;sb=np.clip((pa-Tb)@axis,-L/2,L/2);pb=Tb+sb[:,None]*axis;delta=pa-pb;d=np.linalg.norm(delta,axis=1);n=delta/d[:,None];ids=np.flatnonzero(d<RS);N=np.zeros(len(d));ga=np.zeros_like(pa);Fa=np.zeros(3);Fb=np.zeros(3);Ma=np.zeros(3);Mb=np.zeros(3);E=0.
  if normal and len(ids):
   half=np.sqrt(RS*RS-d[ids]**2);center=sb[ids]
   if np.any(center-half<=-L/2) or np.any(center+half>=L/2):raise ValueError('Endpoint-truncated pressure mapping is outside this unique-interior-component gate')
   theta=.5*np.pi*self.gx;sv=center[:,None]+half[:,None]*np.sin(theta);weights=self.w[ids,None]*half[:,None]*np.cos(theta)*(.5*np.pi*self.gw);pbt=Tb+sv[:,:,None]*axis;diff=pa[ids,None]-pbt;rr=np.linalg.norm(diff,axis=2);root=np.sqrt(np.maximum(RS*RS-rr*rr,0));z=root/RS;term=np.empty_like(z);small=z<.02;z2=z*z;term[small]=RS*z[small]*z2[small]*(1/3+z2[small]*(1/5+z2[small]*(1/7+z2[small]/9)));term[~small]=RS*np.arccosh(RS/rr[~small])-root[~small];E=float(np.sum(KC/np.pi*term*weights));co=-KC/np.pi*root/(rr*rr)*weights;gp=co[:,:,None]*diff;ga[ids]=gp.sum(axis=1);Fa=ga.sum(axis=0);Fb=-Fa;Ma=np.cross(pa-Ta,ga).sum(axis=0);Mb=np.cross(pbt-Tb,-gp).sum(axis=(0,1));N[ids]=-np.sum(ga[ids]*n[ids],axis=1)
  return dict(q=np.array(q).copy(),Ta=Ta,Tb=Tb,Ra=Ra,Rb=Rb,axis=axis,pa=pa,pb=pb,sb=sb,n=n,ids=ids,N=N,E=E,normal_wrenches=np.r_[Fa,Ma,Fb,Mb],segment=(sb>=0).astype(int))
 def state(self,q,elastic=None):
  g=self.geometry(q);e=np.zeros((len(self.s),3)) if elastic is None else np.asarray(elastic).copy();return dict(q=np.array(q).copy(),geometry=g,elastic=e,identity=self.identity)
 def Q(self,a,b,ids):
  Qi=b['Ra']@a['Ra'].T;Qj=b['Rb']@a['Rb'].T;mean=Qi@exp(.5*log(Qi.T@Qj));pred=a['n'][ids]@mean.T;cross=np.cross(pred,b['n'][ids]);c=np.sum(pred*b['n'][ids],axis=1)
  if np.min(c,initial=1)<-.99:raise ValueError('Normal rotation requires subdivision')
  K=skew(cross);return (np.eye(3)+K+K@K/(1+c)[:,None,None])@mean
 def trial(self,q,old):
  if old['identity']!=self.identity:raise ValueError('Material quadrature history mismatch')
  a=old['geometry'];b=self.geometry(q);ids=np.flatnonzero((b['N']>0)|np.any(old['elastic']!=0,axis=1));e=np.zeros_like(old['elastic']);sliding=np.zeros(len(self.s),dtype=bool);heat=loss=work=0.;fric=np.zeros(12);du=np.empty((0,3));tau=np.empty((0,3));Q=np.empty((0,3,3))
  if len(ids):
   Qi=b['Ra']@a['Ra'].T;Qj=b['Rb']@a['Rb'].T;Q=self.Q(a,b,ids);m0=.5*(a['pa'][ids]+a['pb'][ids]);m1=.5*(b['pa'][ids]+b['pb'][ids]);pb1_old=a['sb'][ids,None]*b['axis']+b['Tb'];pb0_new=b['sb'][ids,None]*a['axis']+a['Tb'];forward=b['pa'][ids]+(m0-a['pa'][ids])@Qi.T-pb1_old-(m0-a['pb'][ids])@Qj.T;back=a['pa'][ids]+(m1-b['pa'][ids])@Qi-pb0_new-(m1-b['pb'][ids])@Qj;forward-=b['n'][ids]*np.sum(forward*b['n'][ids],axis=1)[:,None];back-=a['n'][ids]*np.sum(back*a['n'][ids],axis=1)[:,None];du=.5*(forward-mv(Q,back));e0=mv(Q,old['elastic'][ids]);trial=e0+du;mag=np.linalg.norm(trial,axis=1);cap=MU*b['N'][ids];ratio=np.minimum(1.,np.divide(cap,self.k[ids]*mag,out=np.ones_like(mag),where=mag>0));e[ids]=ratio[:,None]*trial;sliding[ids]=ratio<1.;dp=trial-e[ids];tau=self.k[ids,None]*e[ids];heat=float(np.sum(cap*np.linalg.norm(dp,axis=1)));loss=float(.5*np.sum(self.k[ids]*np.sum((e[ids]-e0)**2,axis=1)));work=float(np.sum(tau*du));F=tau.sum(axis=0);Ma=np.cross(m1-b['Ta'],tau).sum(axis=0);Mb=-np.cross(m1-b['Tb'],tau).sum(axis=0);fric=np.r_[F,Ma,-F,Mb]
  wrench=b['normal_wrenches']+fric;A=b['Ra']@right_jacobian(q[3:6]);B=b['Rb']@right_jacobian(q[9:12]);res=np.r_[S*wrench[:3],A.T@wrench[3:6],S*wrench[6:9],B.T@wrench[9:12]]/UNIT;Ut=float(.5*np.sum(self.k*np.sum(e*e,axis=1)));Uold=float(.5*np.sum(self.k*np.sum(old['elastic']**2,axis=1)));new=dict(q=np.array(q).copy(),geometry=b,elastic=e,identity=self.identity)
  return dict(residual=res,state=new,wrenches=wrench,normal_wrenches=b['normal_wrenches'],friction_wrenches=fric,normal_energy_J=b['E'],tangential_energy_J=Ut,friction_heat_J=heat,numerical_loss_J=loss,tangential_endpoint_work_J=work,tangential_balance_J=work-(Ut-Uold)-heat-loss,ids=ids,tau=tau,du=du,contact_pa=b['pa'][ids],contact_pb=b['pb'][ids],active_set=b['N']>0,sliding_set=sliding)
 def jacobian(self,q,old,h=2e-6,columns=None):
  columns=np.arange(12) if columns is None else np.array(columns);out=np.empty((12,len(columns)));baseline=self.trial(q,old);audits=[]
  for k,j in enumerate(columns):
   v=np.zeros(12);v[j]=h;trials=[self.trial(q+c*v,old) for c in [-2,-1,1,2]];out[:,k]=(trials[0]['residual']-8*trials[1]['residual']+8*trials[2]['residual']-trials[3]['residual'])/(12*h)
   audits.append(dict(column=int(j),step=h,changes=[dict(offset=c,active=int(np.count_nonzero(z['active_set']!=baseline['active_set'])),sliding=int(np.count_nonzero(z['sliding_set']!=baseline['sliding_set']))) for c,z in zip([-2,-1,1,2],trials)]))
  self.last_jacobian_audit=audits;return out
 def native(self,q):
  g=self.geometry(q);arc=np.array([-L/2,0.,L/2]);pa=g['Ta']+arc[:,None]*g['Ra'][:,1];pb=g['Tb']+arc[:,None]*g['Rb'][:,0];p=np.array([pa,pb]);E,gp,stats=native_pair(p,np.array([r,r]),np.full((2,2),L/2),KC,order=48,self_contact=False);wrench=np.r_[gp[0].sum(axis=0),np.cross(pa-g['Ta'],gp[0]).sum(axis=0),gp[1].sum(axis=0),np.cross(pb-g['Tb'],gp[1]).sum(axis=0)];return E,wrench
