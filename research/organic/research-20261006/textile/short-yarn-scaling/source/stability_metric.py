"""Sign-equivalent generalized constrained physical-Hessian stability test.
B is the positive rod/contact search metric, not added physical energy.
The left material-frame clamps make its rod bending/twist form positive.
"""
import numpy as np,time
from scipy.sparse.linalg import LinearOperator,eigsh
from coupled_preconditioner import CoupledPreconditioner

def assess(op,pre_unused=None):
 start=time.monotonic();c=op.c;pre=CoupledPreconditioner(op);Q=op.project;n=op.size;J=np.column_stack([op.root_tv(np.eye(3*c.N)[i]) for i in range(3*c.N)]).T;response=pre.rotate_inverse(J.T);S=J@response;free=np.repeat(~c.endpoint_fixed,3);S[free,free]=1.;Sinv=np.linalg.inv(.5*(S+S.T))
 def bi(v):
  p=Q(v);y=pre.rotate_inverse(p);return y-response@(Sinv@op.root_mv(y))+v-p
 A=LinearOperator((n,n),matvec=lambda v:Q(op.matvec(Q(v)))+v-Q(v),dtype=float);B=LinearOperator((n,n),matvec=lambda v:Q(op.matvec(Q(v),True))+v-Q(v),dtype=float);BI=LinearOperator((n,n),matvec=bi,dtype=float);rng=np.random.default_rng(923);test=rng.normal(size=n);inverse_error=float(np.linalg.norm(B@bi(test)-test)/np.linalg.norm(test));ev,v=eigsh(A,M=B,Minv=BI,k=3,which='SA',tol=1e-9,maxiter=1000,ncv=25,v0=rng.normal(size=n));rr=[float(np.linalg.norm(A@v[:,i]-ev[i]*(B@v[:,i]))/np.linalg.norm(B@v[:,i])) for i in range(len(ev))];roots=[float(np.linalg.norm(op.root_mv(v[:,i]))) for i in range(len(ev))];rayleigh=[float(v[:,i]@(A@v[:,i])/(v[:,i]@v[:,i])) for i in range(len(ev))]
 return dict(method='Projected physical Hessian generalized against positive rod/contact metric; no residual stiffness added to energy',minimum_generalized_eigenvalue=float(ev[0]),generalized_eigenvalues=ev.tolist(),generalized_eigen_residuals=rr,root_eigen_residuals=roots,physical_Rayleigh_quotients=rayleigh,positive_metric_inverse_relative_error=inverse_error,wall_s=time.monotonic()-start,passed=bool(ev[0]>1e-6 and max(rr)<1e-6 and max(roots)<1e-9 and inverse_error<1e-7))
