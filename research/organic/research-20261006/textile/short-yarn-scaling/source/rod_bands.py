"""Local rod Hessian/positive search bands; no per-fibre dense square arrays."""
import numpy as np
from scipy.sparse import coo_matrix
from loop_bundle import exp,log,right_jacobian,mv

def physical_bands(rod,w,ends,intrinsic,eps=1e-5):
 M=rod.M;diag=np.zeros((M,3,3));upper=np.zeros((M-1,3,3));lower=np.zeros_like(upper)
 def gradient(a):
  _,g=rod.energy_gradient(exp(a),ends,intrinsic);return mv(np.swapaxes(right_jacobian(a),-1,-2),g)
 for color in range(3):
  ids=np.arange(color,M,3)
  for component in range(3):
   d=np.zeros_like(w);d[ids,component]=eps;derivative=(gradient(w+d)-gradient(w-d))/(2*eps)
   diag[ids,:,component]=derivative[ids]
   js=ids[ids>0];upper[js-1,:,component]=derivative[js-1]
   js=ids[ids<M-1];lower[js,:,component]=derivative[js+1]
 return .5*(diag+np.swapaxes(diag,-1,-2)),.5*(upper+np.swapaxes(lower,-1,-2))

def positive_bands(rod,w,ends):
 R=exp(w);J=right_jacobian(w);M=rod.M;diag=np.zeros((M,3,3));rel=np.swapaxes(R[:-1],-1,-2)@R[1:];Q=right_jacobian(log(rel),inverse=True);a=-Q@np.swapaxes(rel,-1,-2)@J[:-1];b=Q@J[1:];at=np.swapaxes(a,-1,-2);bt=np.swapaxes(b,-1,-2);weights=rod.mat.B[None,:,None]/rod.dual[:,None,None];diag[:-1]+=at@(weights*a);diag[1:]+=bt@(weights*b);off=at@(weights*b)
 for end,Q in enumerate(ends):
  if Q is None:continue
  i=0 if end==0 else M-1;rel=Q.T@R[i] if end==0 else R[i].T@Q;iv=right_jacobian(log(rel),inverse=True);a=iv@J[i] if end==0 else -iv@rel.T@J[i];diag[i]+=a.T@(rod.mat.B[:,None]*a)/(.5*rod.ds[i])
 return diag,off

def sparse_bands(diag,off):
 M=len(diag);i=np.arange(M);rows=(3*i[:,None,None]+np.arange(3)[None,:,None])+np.zeros((M,1,3),dtype=int);cols=(3*i[:,None,None]+np.arange(3)[None,None,:])+np.zeros((M,3,1),dtype=int);u=rows[:-1];v=cols[:-1]+3
 return coo_matrix((np.r_[diag.ravel(),off.ravel(),off.ravel()],(np.r_[rows.ravel(),u.ravel(),v.ravel()],np.r_[cols.ravel(),v.ravel(),u.ravel()])),shape=(3*M,3*M)).tocsr()
