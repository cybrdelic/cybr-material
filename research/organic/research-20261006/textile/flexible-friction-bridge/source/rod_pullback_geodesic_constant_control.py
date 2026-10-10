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
  """C0 geodesic material-spin interpolation between edge-centre directors.
  Returns exact body-spin maps to the neighbouring rod directors. Centreline
  positions and the original rod constitutive energy remain unchanged.
  """
  s=np.asarray(s);sm=.5*(self.rod.s[:-1]+self.rod.s[1:]);k=np.searchsorted(sm,s,side='right')-1;left=np.clip(k,0,self.M-1);right=np.clip(k+1,0,self.M-1);den=sm[right]-sm[left];a=np.divide(s-sm[left],den,out=np.zeros_like(s),where=den>0);a=np.clip(a,0,1);Rl=g['R'][i,left];Rr=g['R'][i,right];Q=np.swapaxes(Rl,-1,-2)@Rr;u=log(Q);C=exp(a[:,None]*u);D=Rl@C;T=(a[:,None,None]*right_jacobian(a[:,None]*u))@right_jacobian(u,inverse=True);A=np.swapaxes(C,-1,-2)-T@np.swapaxes(Q,-1,-2);return D,left,right,A,T
 def sample(self,g,i,s):
  s=np.asarray(s);j=np.minimum(np.searchsorted(self.rod.s,s,side='right')-1,self.M-1);j=np.maximum(j,0);u=(s-self.rod.s[j])/self.rod.ds[j];x=(1-u[:,None])*g['p'][i,j]+u[:,None]*g['p'][i,j+1];D,_,_,_,_=self.surface(g,i,s);return x,D,j,u
 def point_load(self,g,i,s,F,M,gp,gm):
  _,_,j,u=self.sample(g,i,s);np.add.at(gp[i],j,(1-u[:,None])*F);np.add.at(gp[i],j+1,u[:,None]*F);D,left,right,A,B=self.surface(g,i,s);body=mv(np.swapaxes(D,-1,-2),M);np.add.at(gm[i],left,mv(g['R'][i,left],mv(np.swapaxes(A,-1,-2),body)));np.add.at(gm[i],right,mv(g['R'][i,right],mv(np.swapaxes(B,-1,-2),body)))
 def pullback(self,g,gp,gm):
  out=np.zeros((self.n,3+3*self.M))
  for i in range(self.n):
   body=self.rod.pullback_position_gradient(g['R'][i],gp[i])+mv(np.swapaxes(g['R'][i],-1,-2),gm[i]);out[i,:3]=S*gp[i].sum(axis=0);out[i,3:]=mv(np.swapaxes(right_jacobian(g['w'][i]),-1,-2),body).ravel()
  return out.ravel()
 def velocity(self,g,v,i,s):
  a=np.asarray(v).reshape(self.n,3+3*self.M);body=mv(right_jacobian(g['w'][i]),a[i,3:].reshape(self.M,3));omega=mv(g['R'][i],body);dp=np.vstack([a[i,:3]*S,a[i,:3]*S+np.cumsum(self.rod.ds[:,None]*np.cross(omega,g['R'][i,:,:,2]),axis=0)]);_,_,j,u=self.sample(g,i,s);D,left,right,A,B=self.surface(g,i,s);surface_spin=mv(D,mv(A,body[left])+mv(B,body[right]));return (1-u[:,None])*dp[j]+u[:,None]*dp[j+1],surface_spin
