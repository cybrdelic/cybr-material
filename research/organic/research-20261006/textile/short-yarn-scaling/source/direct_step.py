"""Bounded sparse mixed Newton solve, never a dense global rod Hessian.
Physical energy, all rods and root constraints are unchanged. This diagnostic
factor route is admitted only up to 37 x128; report actual fill and memory.
"""
import numpy as np
from scipy.sparse import coo_matrix,bmat,block_diag,diags,eye,kron
from scipy.sparse.linalg import splu
from contact_blocks import UPPER
from rod_bands import sparse_bands

def step(op,g,con,positive=False):
 c=op.c;N,M,K=c.N,c.M,c.M+1;n=op.size;C=op.contact;ai=C.nodes[:,UPPER[:,0]];bi=C.nodes[:,UPPER[:,1]];B=C.gn if positive else C.full
 rows=np.broadcast_to(ai[:,:,None,None]*3+np.arange(3)[None,None,:,None],B.shape).ravel();cols=np.broadcast_to(bi[:,:,None,None]*3+np.arange(3)[None,None,None,:],B.shape).ravel();off=UPPER[:,0]!=UPPER[:,1];Bo=B[:,off];row2=np.broadcast_to(bi[:,off,None,None]*3+np.arange(3)[None,None,None,:],Bo.shape).ravel();col2=np.broadcast_to(ai[:,off,None,None]*3+np.arange(3)[None,None,:,None],Bo.shape).ravel();H=coo_matrix((np.r_[B.ravel(),Bo.ravel()],(np.r_[rows,row2],np.r_[cols,col2])),shape=(N*K*3,N*K*3)).tocsr()+op.surface;keep=np.ones(N*K*3,dtype=bool);keep[np.arange(N)[:,None]*K*3+np.arange(3)[None,:]]=False;Hp=H[keep][:,keep]*(op.L**2/op.unit);Hq=block_diag(op.positive_rods,format='csc') if positive else block_diag([sparse_bands(op.diag[i]+op.geom[i],op.off[i])/op.unit for i in range(N)],format='csc');D=block_diag(list(op.D.reshape(-1,3,3)/op.L),format='csc');diff=kron(eye(N),kron(diags([np.ones(M),-np.ones(M-1)],[0,-1],shape=(M,M)),eye(3)),format='csc');rootids=np.array([(i*M+M-1)*3+a for i in range(N) if c.endpoint_fixed[i] for a in range(3)]);Jp=coo_matrix((np.ones(len(rootids)),(np.arange(len(rootids)),rootids)),shape=(len(rootids),n)).tocsc();A=bmat([[Hq,None,-D.T,None],[None,Hp,diff.T,Jp.T],[-D,diff,None,None],[None,Jp,None,None]],format='csc');perm=np.r_[np.arange(3*n).reshape(3,N,M,3).transpose(2,1,0,3).ravel(),np.arange(3*n,A.shape[0])];inv=np.argsort(perm);Ap=A[perm][:,perm];factor=splu(Ap,permc_spec='COLAMD');rhs=np.r_[-g,np.zeros(2*n),-con[np.repeat(c.endpoint_fixed,3)]];xp=factor.solve(rhs[perm])
 for _ in range(2):xp+=factor.solve(rhs[perm]-Ap@xp)
 x=xp[inv];dq=x[:n];lam=np.zeros(3*N);lam[np.repeat(c.endpoint_fixed,3)]=x[3*n:];res=np.r_[op.matvec(dq,positive)+op.root_tv(lam)+g,op.root_mv(dq)+con];rel=float(np.linalg.norm(res)/max(np.linalg.norm(np.r_[g,con]),1e-30));return dq,lam,dict(iterations=1,actual_relative_linear_residual=rel,actual_max_linear_residual=float(abs(res).max()),positive_search_metric=positive,sparse_factor_nonzeros=factor.L.nnz+factor.U.nnz)
