import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import json,numpy as np
from carriage import Carriage,solve,exp,P
s=Carriage();z=np.load(P/'data/force_equilibrium.npz');s.x=z['x'].copy();s.theta=float(z['theta_rad']);s.evaluate(s.x)
def save(x,theta,trace):
 target=P/'data/positive_pending.npz';tmp=target.with_suffix('.tmp.npz');np.savez_compressed(tmp,x=x,theta_rad=theta,tension_N=s.tension,qualified=False);os.replace(tmp,target)
 (P/'receipts/positive_progress.json').write_text(json.dumps(dict(last=trace[-1],qualified=False),indent=2))
r=solve(s,s.theta,wall_s=45,checkpoint=save);r['maximum_radius_curvature']=max(s.c.rod.diagnostics(exp(s.c.w[i]))['max_radius_curvature'] for i in range(s.c.N));r['numerical_path_contact']=all(t['accepted_path_certificate']['passed'] for t in r['trace'] if 'step_scale' in t and 'failure' not in t);r['passed']=bool(r['passed'] and r['relative_stock_error']<1e-12 and r['capsule_penetration_m']<1e-6 and r['maximum_radius_curvature']<.1 and r['numerical_path_contact']);r['scope']='Static suspended preparation only. Endpoint homotopy is diagnostic; the actual numerical iterate sequence has separately retained contact certificates. No physical unloading, finite frictional spinning or work-history qualification.';r['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
np.savez_compressed(P/'data/positive_equilibrium.npz',x=s.x,w=s.c.w,points=s.c.points(),stock_lengths_m=s.c.ref,roots=np.stack([s.c.origins,s.c.targets],axis=1),radius_m=s.c.radius,endpoint_fixed=s.c.endpoint_fixed,tension_N=s.tension,theta_rad=s.theta,carriage_m=s.d,qualified=r['passed'])
(P/'receipts/positive_equilibrium.json').write_text(json.dumps(r,indent=2));print('RESULT',json.dumps({k:v for k,v in r.items() if k not in ('trace','root_forces_N','root_torques_Nm')}),flush=True)
