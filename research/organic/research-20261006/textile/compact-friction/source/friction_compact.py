"""Associative tangential spring/Coulomb slider on persistent material pairs.
N is a frozen compressive normal force, kt=N/delta0. mu and delta0 are declared
uncalibrated constitutive regularizations. Stored spring energy and plastic slip
heat are separate; rejected nonlinear iterates never mutate accepted history.
"""
import numpy as np
from scipy.sparse import coo_matrix,bmat
from scipy.sparse.linalg import splu

def return_map(relative,plastic,normal,mu,delta0):
 if mu<0 or delta0<=0 or np.min(normal)<0:raise ValueError('Invalid friction coefficients or tension load.')
 t=np.asarray(relative)-np.asarray(plastic);norm=np.linalg.norm(t,axis=1);kt=normal/delta0;cap=mu*normal;limit=mu*delta0;stick=norm<=limit
 ratio=np.minimum(1.,np.divide(limit,norm,out=np.ones_like(norm),where=norm>0));elastic=t*ratio[:,None];newplastic=np.asarray(relative)-elastic;traction=kt[:,None]*elastic
 H=kt[:,None,None]*np.broadcast_to(np.eye(3),(len(t),3,3)).copy();slip=~stick
 if np.any(slip):
  direction=t[slip]/norm[slip,None];H[slip]=((cap[slip]/norm[slip])[:,None,None]*(np.eye(3)-direction[:,:,None]*direction[:,None,:]))
 stored=.5*np.sum(kt*np.sum(elastic**2,axis=1));heat=np.sum(cap*np.linalg.norm(newplastic-plastic,axis=1))
 return dict(traction=traction,H=H,plastic=newplastic,elastic=elastic,stored_J=float(stored),heat_J=float(heat),incremental_potential_J=float(stored+heat),slip=slip,cap=cap)

def potential_change(t,dt,normal,mu,delta0):
 """Stable exact radial Huber-potential difference, including yield crossings."""
 a=np.linalg.norm(t,axis=1);b=np.linalg.norm(t+dt,axis=1);dn=np.divide(2*np.sum(t*dt,axis=1)+np.sum(dt*dt,axis=1),a+b,out=np.zeros_like(a),where=(a+b)>0);limit=mu*delta0;kt=normal/delta0;cap=mu*normal;bothstick=(a<=limit)&(b<=limit);bothslip=(a>=limit)&(b>=limit);change=np.empty_like(a);change[bothstick]=.5*kt[bothstick]*(a[bothstick]+b[bothstick])*dn[bothstick];change[bothslip]=cap[bothslip]*dn[bothslip];up=(a<limit)&(b>limit);down=(a>limit)&(b<limit);change[up]=.5*kt[up]*(limit+a[up])*(limit-a[up])+cap[up]*(b[up]-limit);change[down]=cap[down]*(limit-a[down])+.5*kt[down]*(b[down]+limit)*(b[down]-limit);return float(change.sum())

def blocks_csr(blocks):
 n=len(blocks);ids=np.arange(n)*3;rows=np.broadcast_to(ids[:,None,None]+np.arange(3)[None,:,None],blocks.shape).ravel();cols=np.broadcast_to(ids[:,None,None]+np.arange(3)[None,None,:],blocks.shape).ravel();return coo_matrix((blocks.ravel(),(rows,cols)),shape=(3*n,3*n)).tocsr()

def solve(f,force_N,mu,delta0,previous,maxiter=40):
 y=np.asarray(previous['y'],dtype=float).copy();plastic=np.asarray(previous['plastic'],dtype=float)
 if y.shape!=(2*f.n,) or plastic.shape!=(len(f.table),3) or not np.isfinite(y).all() or not np.isfinite(plastic).all():raise ValueError('Friction history shape/finiteness mismatch')
 identity=previous.get('contact_identity')
 if identity!=f.contact_identity and (identity is not None or np.any(y) or np.any(plastic)):raise ValueError('Contact material identity mismatch; refusing to reset/reassign anchors')
 if previous.get('mu',mu)!=mu or previous.get('delta0_m',delta0)!=delta0:raise ValueError('History coefficient mismatch; a separate material intervention needs its own initial history')
 history=[];nf=f.C.shape[0];ext=f.force*force_N
 def eval(y):
  r=f.friction.evaluate(y,plastic,mu,delta0);g=f.H@y+r['gradient']/f.unit-ext;value=.5*y@(f.H@y)+r['incremental_potential_J']/f.unit-ext@y;return r,g,value
 for it in range(maxiter):
  r,g,value=eval(y);B=f.H+r['H']/f.unit;K=bmat([[B,f.C.T],[f.C,None]],format='csc');factor=splu(K);rhs=-np.r_[g,f.C@y];step=factor.solve(rhs)
  for _ in range(2):step+=factor.solve(rhs-K@step)
  lag=step[len(y):];res=g+f.C.T@lag;resid=float(abs(res).max());linear=float(np.linalg.norm(K@step-rhs)/max(np.linalg.norm(rhs),1e-30));history.append(dict(iteration=it,residual=resid,linear_residual=linear))
  if resid<1e-10 and abs(f.C@y).max()<1e-12:break
  direction=step[:len(y)];projected_gradient=g+f.C.T@lag;slope=projected_gradient@direction;alpha=1.
  for ls in range(30):
   dy=alpha*direction;trial=y+dy
   du=f.friction.relative(dy);delta_phi=potential_change(r['relative']-plastic,du,f.normal,mu,delta0);change=float(dy@projected_gradient+.5*dy@(f.H@dy))+(delta_phi-float(dy@r['gradient']))/f.unit
   if change<=1e-4*alpha*slope:break
   alpha*=.5
  else:raise RuntimeError('Local friction line search failed: '+str(dict(residual=resid,slope=slope,change=change,alpha=alpha,history=history[-3:])))
  y=trial
 else:raise RuntimeError('Local friction equilibrium not converged: '+str(history[-6:]))
 r,g,value=eval(y);stored=float(.5*y@(f.H@y)*f.unit+r['stored_J']);tip=float(f.L*y[f.n+(f.c.N*f.c.M-1)*3+1]);cone=float(np.max(np.linalg.norm(r['elastic'],axis=1)*f.normal/delta0-mu*f.normal))
 return dict(contact_identity=f.contact_identity,mu=mu,delta0_m=delta0,y=y,plastic=r['plastic'],spring=r['elastic'],heat_increment_J=r['heat_J'],stored_J=stored,rod_stored_J=float(.5*y@(f.H@y)*f.unit),spring_stored_J=r['stored_J'],tip_y_m=tip,force_N=float(force_N),force_residual=resid,constraint_residual=float(abs(f.C@y).max()),cone_violation_N=cone,sliding_contacts=r['sliding_contacts'] if mu>0 else 0,iterations=it+1,iteration_trace=history)
