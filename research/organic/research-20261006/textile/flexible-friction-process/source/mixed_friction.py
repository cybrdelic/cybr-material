"""Compact positive mixed-coordinate friction metric, never the physical law."""
from pathlib import Path
import sys,numpy as np
from scipy.sparse import coo_matrix
P=Path(__file__).resolve().parents[1];BRIDGE=P.parent/'flexible-friction-bridge';sys.path.insert(0,str(BRIDGE/'source'))
from rod_pullback import L,right_jacobian,skew

def mixed_metric(f,t):
 g=t['state']['g'];n=f.n*f.M*3;rows=[];cols=[];values=[];weights=[];offset=0
 def block(col,mat,mask=None):
  count=len(col);ids=np.arange(count)
  if mask is not None:ids=ids[mask]
  if not len(ids):return
  mat=mat[ids];row=np.broadcast_to(offset+ids[:,None,None]*3+np.arange(3)[None,:,None],mat.shape);column=np.broadcast_to(col[ids,None,None]+np.arange(3)[None,None,:],mat.shape);rows.append(row.ravel());cols.append(column.ravel());values.append(mat.ravel())
 for a in t['surface']:
  count=len(a['s']);side=f.pairs.index((a['i'],a['j']));ids=np.searchsorted(f.s,a['s']);normal=g['sides'][side]['n'][ids];ki=.5*f.k[ids]*(g['sides'][side]['N'][ids]>0);weights.extend(np.sqrt(ki)[:,None,None]*(np.eye(3)-normal[:,:,None]*normal[:,None,:]))
  for i,s,arm,sign in [(a['i'],a['s'],a['armA'],1),(a['j'],a['sb'],a['armB'],-1)]:
   _,_,edge,u=f.sample(g,i,s)
   for node,weight in [(edge,1-u),(edge+1,u)]:block(n+(i*f.M+node-1)*3,sign*L*weight[:,None,None]*np.broadcast_to(np.eye(3),(count,3,3)),node>0)
   D,neighbours,maps=f.surface(g,i,s)
   for k in range(3):
    j=neighbours[:,k];W=D@maps[:,k]@right_jacobian(g['w'][i,j]);block((i*f.M+j)*3,-sign*skew(arm)@W)
  offset+=3*count
 B=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(offset,2*n)).tocsr();W=np.array(weights);q=np.arange(len(W))*3;wr=np.broadcast_to(q[:,None,None]+np.arange(3)[None,:,None],W.shape);wc=np.broadcast_to(q[:,None,None]+np.arange(3)[None,None,:],W.shape);C=coo_matrix((W.ravel(),(wr.ravel(),wc.ravel())),shape=(offset,offset)).tocsr();A=C@B;H=(A.T@A).tocsc();return H,B,dict(contact_samples=offset//3,mixed_size=2*n,velocity_nnz=B.nnz,metric_nnz=H.nnz,storage_bytes=B.data.nbytes+B.indices.nbytes+B.indptr.nbytes+H.data.nbytes+H.indices.nbytes+H.indptr.nbytes)
