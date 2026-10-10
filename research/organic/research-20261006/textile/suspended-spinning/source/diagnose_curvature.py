import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import time,json,signal,resource,numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from newton_globalized import MatrixFreeBundle,root_blocks,multipliers,roots_tv
from operators import RodHessianOperator
from coupled_preconditioner import CoupledPreconditioner
from direct_step import step
P=Path(__file__).resolve().parents[1];c=MatrixFreeBundle(19,128);c.controls(.00238,0);c.w[:]=np.load(P/'data/yarn_19_force.npz')['w'];E,g=c.evaluate(c.w.ravel());J,iv=root_blocks(c);lam=multipliers(J,iv,g);res=(g+roots_tv(J,lam)).reshape(c.w.shape);op=RodHessianOperator(c,lam);pre=CoupledPreconditioner(op);Q=op.project;n=op.size;JJ=np.column_stack([op.root_tv(np.eye(3*c.N)[i]) for i in range(3*c.N)]).T;response=pre.rotate_inverse(JJ.T);S=JJ@response;free=np.repeat(~c.endpoint_fixed,3);S[free,free]=1.;Sinv=np.linalg.inv(.5*(S+S.T))
def bi(v):
 p=Q(v);y=pre.rotate_inverse(p);return y-response@(Sinv@op.root_mv(y))+v-p
A=LinearOperator((n,n),matvec=lambda v:Q(op.matvec(Q(v)))+v-Q(v),dtype=float);B=LinearOperator((n,n),matvec=lambda v:Q(op.matvec(Q(v),True))+v-Q(v),dtype=float);BI=LinearOperator((n,n),matvec=bi,dtype=float);rng=np.random.default_rng(21);test=rng.normal(size=n);inv_error=float(np.linalg.norm(B@bi(test)-test)/np.linalg.norm(test));dx,ll,diag=step(op,g,c.endpoint_constraint(),False);out=dict(maximum_force_dof=np.unravel_index(np.argmax(abs(res)),res.shape),per_fibre_max_force=np.max(abs(res),axis=(1,2)).tolist(),raw_physical_newton_max_rotation_rad=float(np.linalg.norm(dx.reshape(c.w.shape),axis=2).max()),positive_metric_inverse_error=inv_error,direct_step=diag)
def timeout(*a):raise TimeoutError('Bounded curvature diagnostic timed out')
signal.signal(signal.SIGALRM,timeout);signal.alarm(35);start=time.time()
try:
 ev,v=eigsh(A,M=B,Minv=BI,k=2,which='SA',tol=1e-7,maxiter=600,ncv=20,v0=rng.normal(size=n));rr=[float(np.linalg.norm(A@v[:,i]-ev[i]*(B@v[:,i]))/np.linalg.norm(B@v[:,i])) for i in range(len(ev))];out.update(generalized_minimum_eigenvalues=ev.tolist(),eigen_residuals=rr,root_tangent_error=[float(np.linalg.norm(op.root_mv(v[:,i]))) for i in range(len(ev))]);np.savez_compressed(P/'data/curvature_modes.npz',eigenvalues=ev,vectors=v)
except Exception as e:out['diagnostic_error']=str(e)
finally:signal.alarm(0)
out.update(wall_s=time.time()-start,scope='Curvature of an unconverged 19-fibre equilibrium search. This is a solver diagnosis, not a stable-state qualification.');out['maximum_force_dof']=[int(x) for x in out['maximum_force_dof']];(P/'receipts/curvature_diagnosis.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
