import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import json,numpy as np,time
from carriage import Carriage,RodHessianOperator,stability,exp,P
rows={}
for label in ['spectral','force']:
 s=Carriage();z=np.load(P/'data'/f'{label}_equilibrium.npz');s.x=z['x'].copy();s.theta=float(z['theta_rad']);E,g=s.evaluate(s.x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);st=stability(op,s,J,S);curv=max(s.c.rod.diagnostics(exp(s.c.w[i]))['max_radius_curvature'] for i in range(s.c.N))
 rows[label]=dict(force_residual=float(abs(pg).max()),per_fibre_scaled_residual=np.max(abs(pg[:-1].reshape(s.c.w.shape)),axis=(1,2)).tolist(),carriage_scaled_residual=float(abs(pg[-1])),terminal_stability=st,maximum_radius_curvature=curv,boundary_planes_m=s.c.plane_bounds,qualified=False)
(P/'receipts/terminal_preparation_audit.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
