import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import time,json,resource,numpy as np
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
from newton_globalized import MatrixFreeBundle,root_blocks,multipliers
from operators import RodHessianOperator
from direct_step import step
P=Path(__file__).resolve().parents[1];c=MatrixFreeBundle(19,128);c.controls(.00238,0);c.w[:]=np.load(P/'data/yarn_19_128.npz')['w'];e,g=c.evaluate(c.w.ravel());J,iv=root_blocks(c);lam=multipliers(J,iv,g);op=RodHessianOperator(c,lam);rows=[]
for positive in [False,True]:
 t=time.time();dx,l,r=step(op,g,c.endpoint_constraint(),positive);r.update(wall_s=time.time()-t,projected_descent=float(op.project(g)@op.project(dx)),max_rotation_increment=float(np.linalg.norm(dx.reshape(c.w.shape),axis=2).max()));rows.append(r);print(r,flush=True)
out=dict(tests=rows,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,passed=bool(max(r['actual_relative_linear_residual'] for r in rows)<1e-8));(P/'receipts/direct_step_operator.json').write_text(json.dumps(out,indent=2));assert out['passed']
