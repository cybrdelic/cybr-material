"""Exact material-rod Hessian/Jacobian products without dense global matrices.
Per-fibre sparse mixed solves apply the positive elastic/contact block inverse.
No same-area effective solid rod or mass/stiffness rescaling is used.
"""
import numpy as np
from scipy.sparse import bmat,block_diag,coo_matrix,diags,eye,kron,csr_matrix
from scipy.sparse.linalg import splu,LinearOperator
from contact_blocks import ContactBlocks
from loop_bundle import exp,mv,skew,right_jacobian,planes
from geometry_hessian import geometry_hessian

from rod_bands import physical_bands,positive_bands,sparse_bands

def bands(A,M):
 return np.array([A[3*i:3*i+3,3*i:3*i+3] for i in range(M)]),np.array([A[3*i:3*i+3,3*i+3:3*i+6] for i in range(M-1)])
def band_apply(diag,off,v):
 out=mv(diag,v);out[:,:-1]+=mv(off,v[:,1:]);out[:,1:]+=mv(np.swapaxes(off,-1,-2),v[:,:-1]);return out

def plane_csr(c):
 N,M=c.N,c.M;K=M+1;p=c.points();rows=[];cols=[];values=[]
 for i in range(N):
  for j in range(M):
   for va,vb in [(c.radius-p[i,j,2],c.radius-p[i,j+1,2]),(p[i,j,2]+c.radius-.02,p[i,j+1,2]+c.radius-.02)]:
    if max(va,vb)<=0:continue
    lo,hi=0.,1.
    if min(va,vb)<0:
     t=-va/(vb-va)
     if va<0:lo=t
     else:hi=t
    def F(t):return np.array([[t-t*t+t**3/3,t*t/2-t**3/3],[t*t/2-t**3/3,t**3/3]])
    H=c.kc*c.ref[i,j]*(F(hi)-F(lo));ids=[(i*K+j)*3+2,(i*K+j+1)*3+2]
    for a in range(2):
     for b in range(2):rows.append(ids[a]);cols.append(ids[b]);values.append(H[a,b])
 return coo_matrix((values,(rows,cols)),shape=(N*K*3,N*K*3)).tocsr()

class RodHessianOperator:
 def __init__(self,c,lam,max_cache_MiB=200):
  self.c=c;self.N=c.N;self.M=c.M;self.K=c.M+1;self.size=c.N*c.M*3;self.unit=c.unit;self.L=.0024
  R=exp(c.w);self.D=-c.rod.ds[None,:,None,None]*(R@skew(np.array([0.,0.,1.]))@right_jacobian(c.w));self.JJ=np.einsum('nmij,nmkj->nik',self.D,self.D)/self.L**2;self.JJ[~c.endpoint_fixed]=np.eye(3);self.JJinv=np.linalg.inv(self.JJ)
  p=c.points();self.contact=ContactBlocks(p,c.radius,c.ref,c.kc,c.contact_order,max_cache_MiB);self.surface=plane_csr(c);_,gp,_,_=planes(p,c.radius,c.kc,.02,0.,reference_lengths=c.ref,current_measure=False);force=self.contact.gradient+gp;self.diag=[];self.off=[];self.gndiag=[];self.gnoff=[];self.geom=[];self.positive_rods=[]
  for i in range(c.N):
   dd,oo=physical_bands(c.rods[i],c.w[i],c.frames_for(i),c.boundary_k0[i]);self.diag.append(dd);self.off.append(oo)
   dd,oo=positive_bands(c.rod,c.w[i],c.frames_for(i));self.gndiag.append(dd);self.gnoff.append(oo);self.positive_rods.append(sparse_bands(dd,oo)/c.unit)
   C=np.cumsum(force[i][::-1],axis=0)[::-1][1:]+lam.reshape(c.N,3)[i]*c.unit/self.L;self.geom.append(geometry_hessian(c.rod,c.w[i],C))
  self.diag=np.array(self.diag);self.off=np.array(self.off);self.gndiag=np.array(self.gndiag);self.gnoff=np.array(self.gnoff);self.geom=np.array(self.geom)
 def jvp(self,v):
  v=np.asarray(v).reshape(self.N,self.M,3);p=np.zeros((self.N,self.K,3));p[:,1:]=np.cumsum(mv(self.D,v),axis=1);return p
 def jtv(self,g):
  g=np.asarray(g).reshape(self.N,self.K,3);suffix=np.cumsum(g[:,::-1],axis=1)[:,::-1][:,1:];return mv(np.swapaxes(self.D,-1,-2),suffix)
 def root_mv(self,v):return (np.einsum('nmij,nmj->ni',self.D,np.asarray(v).reshape(self.N,self.M,3))*self.c.endpoint_fixed[:,None]).ravel()/self.L
 def root_tv(self,v):return (mv(np.swapaxes(self.D,-1,-2),(np.asarray(v).reshape(self.N,1,3)*self.c.endpoint_fixed[:,None,None]))/self.L).ravel()
 def project(self,v):return np.asarray(v)-self.root_tv(mv(self.JJinv,self.root_mv(v).reshape(self.N,3)))
 def matvec(self,v,positive=False):
  vv=np.asarray(v).reshape(self.N,self.M,3);p=self.jvp(vv);hp=self.contact.matvec(p,positive=positive).ravel()+self.surface@p.ravel();out=self.jtv(hp)
  if positive:out+=band_apply(self.gndiag,self.gnoff,vv)
  else:out+=band_apply(self.diag,self.off,vv)+mv(self.geom,vv)
  return out.ravel()/self.unit
 def linear_operator(self,positive=False):return LinearOperator((self.size,self.size),matvec=lambda v:self.matvec(v,positive),dtype=float)

class MixedFiberPreconditioner:
 def __init__(self,op):
  self.op=op;N,M,K=op.N,op.M,op.K;n=M*3;self.n=n;self.factors=[];self.mixed_matrices=[];self.perm=np.arange(9*M).reshape(3,M,3).transpose(1,0,2).ravel();self.invperm=np.argsort(self.perm);self.roots=[];self.schur_inv=[];self.factor_nonzeros=0
  difference=kron(diags([np.ones(M),-np.ones(M-1)],[0,-1],shape=(M,M)),eye(3),format='csc')
  for i in range(N):
   ids=slice(i*K*3,(i+1)*K*3);Dp=(op.contact.diagonal_fibre_csr(i)+op.surface[ids,ids])[3:,3:]*(op.L**2/op.unit);D=block_diag(list(op.D[i]/op.L),format='csc');B=op.positive_rods[i]
   mixed=bmat([[B,None,-D.T],[None,Dp,difference.T],[-D,difference,None]],format='csc');permuted=mixed[self.perm][:,self.perm];factor=splu(permuted,permc_spec='NATURAL');self.factors.append(factor);self.mixed_matrices.append(permuted);self.factor_nonzeros+=factor.L.nnz+factor.U.nnz
   J=np.transpose(op.D[i],(1,0,2)).reshape(3,n)/op.L*op.c.endpoint_fixed[i];self.roots.append(J);response=self.solve_fibre(i,J.T);S=J@response;self.schur_inv.append(np.linalg.inv((S+S.T)*.5) if op.c.endpoint_fixed[i] else np.eye(3))
  self.schur_inv=np.array(self.schur_inv)
 def solve_fibre(self,i,v):
  v=np.asarray(v);rhs=np.zeros((3*self.n,)+v.shape[1:]);rhs[:self.n]=v;permuted=rhs[self.perm];sol=self.factors[i].solve(permuted)
  for _ in range(2):sol+=self.factors[i].solve(permuted-self.mixed_matrices[i]@sol)
  return sol[self.invperm][:self.n]
 def rotate_inverse(self,v):
  v=np.asarray(v);vv=v.reshape((self.op.N,self.n)+v.shape[1:]);out=np.empty_like(vv)
  for i in range(self.op.N):out[i]=self.solve_fibre(i,vv[i])
  return out.reshape(v.shape)
 def kkt_inverse(self,v):
  n=self.op.size;return np.r_[self.rotate_inverse(v[:n]),mv(self.schur_inv,v[n:].reshape(self.op.N,3)).ravel()]
