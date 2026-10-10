import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import json,numpy as np,time
from carriage import Carriage,RodHessianOperator,stability,P,direction
s=Carriage();z=np.load(P/'data/zero_equilibrium.npz');s.x=z['x'].copy();s.theta=float(z['theta_rad']);E,g=s.evaluate(s.x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);t=time.monotonic();st=stability(op,s,J,S)
rows=[]
for shift in [0.,.05,.5,1.,2.,4.,8.]:
 dx,ck=direction(op,s,g,s.constraints(),shift);rows.append(dict(shift=shift,projected_descent=float(pg@dx),maximum_rotation_rad=float(np.linalg.norm(dx[:-1].reshape(s.c.w.shape),axis=2).max()),linear=ck))
out=dict(force_residual=float(abs(pg).max()),stability=st,steps=rows,wall_s=time.monotonic()-t,scope='Numerical metric diagnostic at the unconverged suspended preparation, not an equilibrium qualification.')
(P/'receipts/pending_curvature.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
