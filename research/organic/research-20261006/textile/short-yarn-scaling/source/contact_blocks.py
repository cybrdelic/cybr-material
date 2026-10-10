"""Bounded two-pass sparse contact cache: exact same potential and quadrature."""
from pathlib import Path
import ctypes,numpy as np
P=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS');lib=ctypes.CDLL(str(Path(__file__).with_name('_contact_blocks.so')));f=lib.segment_pair_contact_blocks;f.argtypes=[ctypes.c_int,ctypes.c_int,P,P,P,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double,P,P,ctypes.c_int,I,P,P];f.restype=ctypes.c_double;apply=lib.apply_contact_blocks;apply.argtypes=[ctypes.c_int,ctypes.c_int,I,P,P,P];apply.restype=None
UPPER=np.array([(i,j) for i in range(4) for j in range(i,4)])
class ContactBlocks:
 def __init__(self,p,radius,ref,kc,order=16,max_cache_MiB=200):
  p=np.ascontiguousarray(p,dtype=np.float64);ref=np.ascontiguousarray(ref,dtype=np.float64)
  if p.ndim!=3 or p.shape[2]!=3:raise ValueError('Exact (fibres,nodes,3) shape required.')
  N,K,_=p.shape
  if not 1<=N<=700 or not 2<=K<=1025 or N*(K-1)>50000 or ref.shape!=(N,K-1) or not np.isfinite(p).all() or not np.isfinite(ref).all() or np.min(ref)<=0 or np.any(np.linalg.norm(np.diff(p,axis=1),axis=2)<=1e-14):raise ValueError('Invalid bounded material geometry.')
  if not np.isfinite([radius,kc,max_cache_MiB]).all() or radius<1e-10 or kc<0 or max_cache_MiB<=0 or not isinstance(order,int) or not 4<=order<=48:raise ValueError('Invalid contact/cache inputs.')
  self.N=N;self.K=K;self.nodes_total=N*K;r=np.full(N,radius);g=np.zeros_like(p);s=np.zeros(7);dummy=np.zeros(1);ids=np.zeros(1,dtype=np.int32)
  E=f(N,K,p,r,ref,kc,order,1,4.,g,s,0,ids,dummy,dummy)
  if not np.isfinite(E):raise FloatingPointError('Invalid contact in cache counting pass.')
  self.count=int(s[5]);required=self.count*(4*4+2*90*8)
  if required>max_cache_MiB*1024**2:raise MemoryError('Contact cache admission exceeded: '+str(required)+' bytes.')
  self.nodes=np.empty((self.count,4),dtype=np.int32);self.full=np.empty((self.count,10,3,3));self.gn=np.empty_like(self.full)
  if self.count:
   E=f(N,K,p,r,ref,kc,order,1,4.,g,s,self.count,self.nodes,self.full,self.gn)
   if not np.isfinite(E) or int(s[5])!=self.count:raise FloatingPointError('Contact cache count/fill mismatch.')
  self.energy=float(E);self.gradient=g;self.bytes=required;self.stats=s.copy()
 def matvec(self,v,positive=False):
  v=np.ascontiguousarray(v,dtype=float)
  if v.size!=self.nodes_total*3 or not np.isfinite(v).all():raise ValueError('Contact matvec shape/finiteness mismatch.')
  out=np.zeros_like(v);apply(self.count,self.nodes_total,self.nodes,self.gn if positive else self.full,v,out);return out
 def diagonal_fibre_csr(self,fibre):
  from scipy.sparse import coo_matrix
  ai=self.nodes[:,UPPER[:,0]];bi=self.nodes[:,UPPER[:,1]];mask=(ai//self.K==fibre)&(bi//self.K==fibre);a=(ai[mask]%self.K).ravel();b=(bi[mask]%self.K).ravel();B=self.gn[mask];pair_i=np.broadcast_to(UPPER[:,0],ai.shape)[mask];pair_j=np.broadcast_to(UPPER[:,1],ai.shape)[mask]
  row=np.broadcast_to(a[:,None,None]*3+np.arange(3)[None,:,None],B.shape).ravel();col=np.broadcast_to(b[:,None,None]*3+np.arange(3)[None,None,:],B.shape).ravel();values=B.ravel();off=pair_i!=pair_j;Bo=B[off];row2=np.broadcast_to(b[off,None,None]*3+np.arange(3)[None,None,:],Bo.shape).ravel();col2=np.broadcast_to(a[off,None,None]*3+np.arange(3)[None,:,None],Bo.shape).ravel();return coo_matrix((np.r_[values,Bo.ravel()],(np.r_[row,row2],np.r_[col,col2])),shape=(self.K*3,self.K*3)).tocsr()
