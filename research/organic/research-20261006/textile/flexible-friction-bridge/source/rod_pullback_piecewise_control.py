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
 def sample(self,g,i,s):
  s=np.asarray(s);j=np.minimum(np.searchsorted(self.rod.s,s,side='right')-1,self.M-1);j=np.maximum(j,0);u=(s-self.rod.s[j])/self.rod.ds[j];x=(1-u[:,None])*g['p'][i,j]+u[:,None]*g['p'][i,j+1];return x,g['R'][i,j],j,u
 def point_load(self,g,i,s,F,M,gp,gm):
  _,_,j,u=self.sample(g,i,s);np.add.at(gp[i],j,(1-u[:,None])*F);np.add.at(gp[i],j+1,u[:,None]*F);np.add.at(gm[i],j,M)
 def pullback(self,g,gp,gm):
  out=np.zeros((self.n,3+3*self.M))
  for i in range(self.n):
   body=self.rod.pullback_position_gradient(g['R'][i],gp[i])+mv(np.swapaxes(g['R'][i],-1,-2),gm[i]);out[i,:3]=S*gp[i].sum(axis=0);out[i,3:]=mv(np.swapaxes(right_jacobian(g['w'][i]),-1,-2),body).ravel()
  return out.ravel()
 def velocity(self,g,v,i,s):
  a=np.asarray(v).reshape(self.n,3+3*self.M);omega=mv(g['R'][i]@right_jacobian(g['w'][i]),a[i,3:].reshape(self.M,3));dp=np.vstack([a[i,:3]*S,a[i,:3]*S+np.cumsum(self.rod.ds[:,None]*np.cross(omega,g['R'][i,:,:,2]),axis=0)]);_,_,j,u=self.sample(g,i,s);return (1-u[:,None])*dp[j]+u[:,None]*dp[j+1],omega[j]
