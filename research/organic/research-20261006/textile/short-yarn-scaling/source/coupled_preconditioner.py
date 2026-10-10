"""Exact sparse inverse of the positive contact/elastic search metric (N<=7).
All inter-fibre contact blocks are retained. This changes numerical conditioning,
not the physical energy or material coefficients. No dense global Hessian exists.
"""
import numpy as np
from scipy.sparse import coo_matrix,bmat,block_diag,diags,eye,kron
from scipy.sparse.linalg import splu
from loop_bundle import mv
from contact_blocks import UPPER
class CoupledPreconditioner:
 def __init__(self,op):
  if op.N>37 or op.M>128:raise ValueError('Scaling preconditioner bounded to 37 fibres ×128 segments.')
  self.op=op;n=op.size;self.n=n;C=op.contact;ai=C.nodes[:,UPPER[:,0]];bi=C.nodes[:,UPPER[:,1]];B=C.gn;rows=np.broadcast_to(ai[:,:,None,None]*3+np.arange(3)[None,None,:,None],B.shape).ravel();cols=np.broadcast_to(bi[:,:,None,None]*3+np.arange(3)[None,None,None,:],B.shape).ravel();off=UPPER[:,0]!=UPPER[:,1];Bo=B[:,off];row2=np.broadcast_to(bi[:,off,None,None]*3+np.arange(3)[None,None,None,:],Bo.shape).ravel();col2=np.broadcast_to(ai[:,off,None,None]*3+np.arange(3)[None,None,:,None],Bo.shape).ravel();H=coo_matrix((np.r_[B.ravel(),Bo.ravel()],(np.r_[rows,row2],np.r_[cols,col2])),shape=(op.N*op.K*3,op.N*op.K*3)).tocsr()+op.surface;keep=np.ones(op.N*op.K*3,dtype=bool);keep[np.arange(op.N)[:,None]*op.K*3+np.arange(3)[None,:]]=False;Dp=H[keep][:,keep]*(op.L**2/op.unit);D=block_diag(list(op.D.reshape(-1,3,3)/op.L),format='csc');diff=kron(eye(op.N),kron(diags([np.ones(op.M),-np.ones(op.M-1)],[0,-1],shape=(op.M,op.M)),eye(3)),format='csc');B=block_diag(op.positive_rods,format='csc');A=bmat([[B,None,-D.T],[None,Dp,diff.T],[-D,diff,None]],format='csc');self.perm=np.arange(3*n).reshape(3,op.N,op.M,3).transpose(2,1,0,3).ravel();self.inv=np.argsort(self.perm);self.A=A[self.perm][:,self.perm];self.factor=splu(self.A,permc_spec='NATURAL');self.factor_nonzeros=self.factor.L.nnz+self.factor.U.nnz;J=np.column_stack([op.root_tv(np.eye(op.N*3)[i]) for i in range(op.N*3)]).T;response=self.rotate_inverse(J.T);S=J@response;free=np.repeat(~op.c.endpoint_fixed,3);S[free,free]=1.;self.Sinv=np.linalg.inv(.5*(S+S.T))
 def rotate_inverse(self,v):
  v=np.asarray(v);rhs=np.zeros((3*self.n,)+v.shape[1:]);rhs[:self.n]=v;rhs=rhs[self.perm];x=self.factor.solve(rhs)
  for _ in range(2):x+=self.factor.solve(rhs-self.A@x)
  return x[self.inv][:self.n]
 def kkt_inverse(self,v):return np.r_[self.rotate_inverse(v[:self.n]),self.Sinv@v[self.n:]]
