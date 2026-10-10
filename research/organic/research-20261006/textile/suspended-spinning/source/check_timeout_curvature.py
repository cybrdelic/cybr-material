import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,numpy as np
from carriage import Carriage,RodHessianOperator,stability,certify_path,P
accepted=np.load(P/'ramp08/data/accepted_007.npz');out={}
for name in ['rejected_004','rejected_008']:
 z=np.load(P/'ramp08/data'/f'{name}.npz');h=json.loads((P/'ramp08/receipts'/f'{name}.json').read_text());s=Carriage();s.x=z['x'].copy();s.theta=float(z['theta_rad']);_,g=s.evaluate(s.x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);st=stability(op,s,J,S);out[name]=dict(theta_rad=s.theta,force_residual=float(abs(pg).max()),last_recorded_residuals=[r['force_residual'] for r in h['trace'][-5:]],stability=st,path_from_latest_accepted=certify_path(accepted['points'],s.c.points(),s.c.radius,s.c.ref,plane_bounds=None))
(P/'receipts/timeout_curvature.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
