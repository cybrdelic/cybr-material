"""Matrix-free equality-constrained Newton with per-fibre physical preconditioning.
Only search directions use a positive metric fallback. Every acceptance and
terminal residual is evaluated with the unchanged physical energy/gradient.
"""
import time,warnings
import numpy as np
from scipy.sparse.linalg import LinearOperator,minres,lobpcg,eigsh
from loop_bundle import LoopBundle,exp,log,skew,right_jacobian,mv
from operators import RodHessianOperator,MixedFiberPreconditioner
from contact_path import certify_path
from types import SimpleNamespace
from direct_step import step as direct_step

def root_blocks(c):
 D=-c.rod.ds[None,:,None,None]*(exp(c.w)@skew(np.array([0.,0.,1.]))@right_jacobian(c.w))/.0024
 J=np.transpose(D,(0,2,1,3)).reshape(c.N,3,c.M*3)*c.endpoint_fixed[:,None,None];JJ=J@np.swapaxes(J,-1,-2);JJ[~c.endpoint_fixed]=np.eye(3);return J,np.linalg.inv(JJ)
def roots_mv(J,v):return np.einsum('nij,nj->ni',J,np.asarray(v).reshape(J.shape[0],-1)).ravel()
def roots_tv(J,v):return np.einsum('nji,nj->ni',J,np.asarray(v).reshape(J.shape[0],3)).ravel()
def multipliers(J,inv,g):return -mv(inv,roots_mv(J,g).reshape(J.shape[0],3)).ravel()

def kkt_step(op,pre,g,con,positive=False,rtol=1e-13,maxiter=3000):
 n=op.size;nr=3*op.N;count=[0]
 def product(x):return np.r_[op.matvec(x[:n],positive)+op.root_tv(x[n:]),op.root_mv(x[:n])+np.repeat(~op.c.endpoint_fixed,3)*x[n:]]
 A=LinearOperator((n+nr,n+nr),matvec=product,dtype=float);P=LinearOperator(A.shape,matvec=pre.kkt_inverse,dtype=float);rhs=-np.r_[g,con]
 def cb(x):count[0]+=1
 sol,info=minres(A,rhs,M=P,rtol=rtol,maxiter=maxiter,callback=cb,check=False)
 rawroot=float(abs(op.root_mv(sol[:n])+con).max())
 # Remove the iterative solve's equality error before using a merit derivative.
 # Near equilibrium an inaccurate normal component can swamp tangent descent.
 for _ in range(2):sol[:n]-=op.root_tv(mv(op.JJinv,(op.root_mv(sol[:n])+con).reshape(op.N,3)))
 residual=product(sol)-rhs;relative=np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-30)
 return sol[:n],sol[n:],dict(raw_max_linear_constraint_residual=rawroot,corrected_max_linear_constraint_residual=float(abs(op.root_mv(sol[:n])+con).max()),iterations=count[0],minres_info=int(info),actual_relative_linear_residual=float(relative),actual_max_linear_residual=float(abs(residual).max()),positive_search_metric=bool(positive))

def local_stability(op,pre,tol=2e-6,maxiter=350):
 """Lanczos on projected physical Hessian; normal-space search eigenvalue = 1.
 The normal shift is a numerical eigenproblem device, not a material energy.
 Actual physical eigenvectors must independently satisfy root tangency.
 """
 n=op.size;calls=[0]
 def product(v):
  calls[0]+=1;p=op.project(v);return op.project(op.matvec(p))+v-p
 A=LinearOperator((n,n),matvec=product,dtype=float)
 with warnings.catch_warnings(record=True) as notes:
  values,vectors=eigsh(A,k=3,which='SA',tol=1e-9,maxiter=10000,ncv=40,v0=np.random.default_rng(923).normal(size=n))
 order=np.argsort(values);values=values[order];vectors=vectors[:,order];res=[float(np.linalg.norm(op.project(op.matvec(v))-ev*v)) for ev,v in zip(values,vectors.T)];root=[float(np.linalg.norm(op.root_mv(v))) for v in vectors.T]
 return dict(minimum_projected_lagrangian_eigenvalue=float(values[0]),eigenvalues=values.tolist(),projected_eigen_residuals=res,root_eigen_residuals=root,matrix_vector_calls=calls[0],warning_count=len(notes),passed=bool(values[0]>-1e-5 and res[0]<1e-5 and root[0]<1e-10))

class MatrixFreeBundle(LoopBundle):
 def __init__(self,*args,**kw):
  super().__init__(*args,**kw)
  if self.N>37 or self.M>128:raise ValueError('Scaling ladder is bounded to at most 37 actual fibres ×128 segments.')
 def solve(self,d,alpha,wall_s=240,maxiter=100):
  roots=self.controls(d,alpha);start=time.monotonic();a=self.w.ravel().copy();history=[];message='Iteration budget';success=False;accepted=0;stability=None;maxbytes=0;maxnz=0;path_rejections=0;path_audits=0
  for it in range(maxiter):
   if time.monotonic()-start>wall_s:message='Wall-time budget';break
   E,g=self.evaluate(a);previous_points=self.last['points'].copy();con=self.endpoint_constraint();J,inv=root_blocks(self);lam=multipliers(J,inv,g);res=float(abs(g+roots_tv(J,lam)).max());cerr=float(abs(con).max()*.0024);op=RodHessianOperator(self,lam);pre=SimpleNamespace(factor_nonzeros=0);maxbytes=max(maxbytes,op.contact.bytes);maxnz=max(maxnz,pre.factor_nonzeros);row=dict(iteration=it,energy_J=E*self.unit,projected_scaled_force_residual=res,root_error_m=cerr,contact_cache_bytes=op.contact.bytes,mixed_factor_nonzeros=pre.factor_nonzeros);history.append(row)
   print('MATRIX_FREE_ITER',it,'res',res,'root',cerr,'elapsed',time.monotonic()-start,flush=True)
   if res<2e-5 and cerr<1e-9:
    stability=local_stability(op,pre);row['stability']=stability
    if stability['passed']:success=True;message='Actual projected force, root constraints and matrix-free local stability passed';break
    message='Local stability not established';break
   candidates=[]
   for positive in [False,True]:
    step,lagstep,diagnostic=direct_step(op,g,con,positive);maxnz=max(maxnz,diagnostic['sparse_factor_nonzeros']);rho=max(1.,1.5*max(float(abs(lam).max()),float(abs(lagstep).max()))+1);pred=float(op.project(g)@op.project(step)+lam@con-rho*np.linalg.norm(con,1));diagnostic['unscaled_predicted_merit_descent']=pred
    row.setdefault('linear_solves',[]).append(diagnostic)
    if pred>=0 or diagnostic['actual_relative_linear_residual']>=1e-5:continue
    angle=float(np.linalg.norm(step.reshape(self.N,self.M,3),axis=2).max());scale=min(1.,.25/max(angle,1e-30));merit=E+rho*np.linalg.norm(con,1)
    for ls in range(30):
     trial=a+scale*step
     try:
      for soc in range(4):
       self.w[:]=trial.reshape(self.w.shape);ct=self.endpoint_constraint();remainder=ct-(1-scale)*con
       if abs(remainder).max()<1e-12:break
       Jt,iv=root_blocks(self);trial-=roots_tv(Jt,mv(iv,remainder.reshape(self.N,3)))
      Et,_=self.evaluate(trial);ct=self.endpoint_constraint()
      if Et+rho*np.linalg.norm(ct,1)<=merit+1e-4*scale*pred+1e-13:
       certificate=certify_path(previous_points,self.last['points'],self.radius,self.ref);path_audits+=certificate['audits']
       if certificate['passed']:
        candidates.append((Et,trial.copy(),dict(step_scale=scale,merit_penalty=rho,line_search_steps=ls+1,accepted_path_certificate=certificate,chosen_positive_metric=positive,actual_energy_drop_J=(E-Et)*self.unit)));break
       path_rejections+=1;row.setdefault('rejected_path_certificates',[]).append(certificate)
     except (FloatingPointError,ValueError,np.linalg.LinAlgError):pass
     scale*=.5
    diagnostic.update(line_search_steps=ls+1,final_scale=scale)
   if not candidates:message='Both physical and positive exact-energy search candidates failed';break
   Et,trial,chosen=min(candidates,key=lambda item:item[0]);row.update(chosen);self.w[:]=trial.reshape(self.w.shape)
   a=log(exp(trial.reshape(self.N,self.M,3))).ravel();accepted+=1
   if hasattr(self,'checkpoint_path'):
    np.savez_compressed(self.checkpoint_path,w=a.reshape(self.N,self.M,3),points=self.points(),material_arc_m=self.rod.s,endpoint_fixed=self.endpoint_fixed,accepted_steps=accepted)
  row=self.record(a,d,alpha,roots,start,success,message,it+1);row.update(algorithm='Competing sparse physical and positive-metric search steps; choose lower actual energy after unchanged SOC/contact gates',stability=stability,accepted_steps=accepted,continuous_contact_path_gate=True,path_rejections=path_rejections,path_capsule_audits=path_audits,maximum_contact_cache_bytes=maxbytes,maximum_mixed_factor_nonzeros=maxnz,solver_iteration_trace=history);return row
