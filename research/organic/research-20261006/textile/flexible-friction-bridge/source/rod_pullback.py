"""Exact material-point/edge-director virtual-work map for existing rods."""
from pathlib import Path
import sys,numpy as np
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P.parent/'suspended-spinning/source'))
from kirchhoff3d import KirchhoffRod,exp,log,right_jacobian,mv,skew
L=.0024;S=32e-6
class RodCoordinates:
 def __init__(self,n=3,segments=16):
  self.n=n;self.M=segments;self.rod=KirchhoffRod(np.linspace(0,L,segments+1));self.radius=self.rod.mat.radius_m;self.size=n*(3+3*segments)
 def unpack(self,q):
  a=np.asarray(q).reshape(self.n,3+3*self.M);o=a[:,:3]*S;w=a[:,3:].reshape(self.n,self.M,3);R=exp(w);p=np.array([self.rod.points(R[i],o[i]) for i in range(self.n)]);return dict(q=np.array(q).copy(),o=o,w=w,R=R,p=p)
 def pack(self,o,R):return np.c_[np.asarray(o)/S,log(R).reshape(self.n,-1)].ravel()
 def surface(self,g,i,s):
  """Objective C0 contact-spin field with quadratic endpoint half-cells.
  Interior uses two-edge geodesics. Endpoint reconstruction uses three internal
  directors in a local SO(3) chart; all exact derivative maps are returned.
  This is a convergent contact discretization, not pointwise Kirchhoff framing.
  """
  s=np.asarray(s);sm=.5*(self.rod.s[:-1]+self.rod.s[1:]);k=np.searchsorted(sm,s,side='right')-1;left=np.clip(k,0,self.M-2);right=left+1;den=sm[right]-sm[left];a=(s-sm[left])/den;Rl=g['R'][i,left];Rr=g['R'][i,right];Q=np.swapaxes(Rl,-1,-2)@Rr;u=log(Q);C=exp(a[:,None]*u);D=Rl@C;T=(a[:,None,None]*right_jacobian(a[:,None]*u))@right_jacobian(u,inverse=True);A=np.swapaxes(C,-1,-2)-T@np.swapaxes(Q,-1,-2);ids=np.stack([left,right,right],axis=1);maps=np.stack([A,T,np.zeros_like(T)],axis=1)
  for end,mask in [(0,s<sm[0]),(1,s>sm[-1])]:
   at=np.flatnonzero(mask)
   if not len(at):continue
   js=np.array([0,1,2]) if end==0 else np.array([self.M-1,self.M-2,self.M-3]);xi=(s[at]-sm[js[0]])/(sm[js[1]]-sm[js[0]]);b1=xi*(2-xi);b2=.5*xi*(xi-1);R0=g['R'][i,js[0]];Q1=R0.T@g['R'][i,js[1]];Q2=R0.T@g['R'][i,js[2]];u1=log(Q1);u2=log(Q2);v=b1[:,None]*u1+b2[:,None]*u2;Cv=exp(v);Jv=right_jacobian(v);T1=(b1[:,None,None]*Jv)@right_jacobian(u1,inverse=True);T2=(b2[:,None,None]*Jv)@right_jacobian(u2,inverse=True);D[at]=R0@Cv;ids[at]=js;maps[at]=np.stack([np.swapaxes(Cv,-1,-2)-T1@Q1.T-T2@Q2.T,T1,T2],axis=1)
  return D,ids,maps
 def sample(self,g,i,s):
  s=np.asarray(s);j=np.minimum(np.searchsorted(self.rod.s,s,side='right')-1,self.M-1);j=np.maximum(j,0);u=(s-self.rod.s[j])/self.rod.ds[j];x=(1-u[:,None])*g['p'][i,j]+u[:,None]*g['p'][i,j+1];D,_,_=self.surface(g,i,s);return x,D,j,u
 def point_load(self,g,i,s,F,M,gp,gm):
  _,_,j,u=self.sample(g,i,s);np.add.at(gp[i],j,(1-u[:,None])*F);np.add.at(gp[i],j+1,u[:,None]*F);D,ids,maps=self.surface(g,i,s);body=mv(np.swapaxes(D,-1,-2),M)
  for k in range(3):np.add.at(gm[i],ids[:,k],mv(g['R'][i,ids[:,k]],mv(np.swapaxes(maps[:,k],-1,-2),body)))
 def pullback(self,g,gp,gm):
  out=np.zeros((self.n,3+3*self.M))
  for i in range(self.n):
   body=self.rod.pullback_position_gradient(g['R'][i],gp[i])+mv(np.swapaxes(g['R'][i],-1,-2),gm[i]);out[i,:3]=S*gp[i].sum(axis=0);out[i,3:]=mv(np.swapaxes(right_jacobian(g['w'][i]),-1,-2),body).ravel()
  return out.ravel()
 def velocity(self,g,v,i,s):
  a=np.asarray(v).reshape(self.n,3+3*self.M);body=mv(right_jacobian(g['w'][i]),a[i,3:].reshape(self.M,3));omega=mv(g['R'][i],body);dp=np.vstack([a[i,:3]*S,a[i,:3]*S+np.cumsum(self.rod.ds[:,None]*np.cross(omega,g['R'][i,:,:,2]),axis=0)]);_,_,j,u=self.sample(g,i,s);D,ids,maps=self.surface(g,i,s);surface_spin=mv(D,sum((mv(maps[:,k],body[ids[:,k]]) for k in range(3)),np.zeros((len(s),3))));return (1-u[:,None])*dp[j]+u[:,None]*dp[j+1],surface_spin
