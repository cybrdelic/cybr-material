import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import json,numpy as np,hashlib,time
from carriage import Carriage,solve,P
s=Carriage();x0=s.x.copy();s.evaluate(x0);initial=s.record(x0,False,[])
def save_trial(x,theta,trace):
 path=P/'data/zero_equilibrium_pending.npz';tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,x=x,theta_rad=theta,tension_N=s.tension,qualified=False);os.replace(tmp,path)
 (P/'receipts/zero_equilibrium_progress.json').write_text(json.dumps(dict(iterations=len(trace),last=trace[-1],qualified=False)))
r=solve(s,0.,wall_s=90,checkpoint=save_trial);r['initial_backed_state_released_to_air']=initial;r['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024;r['boundary_mode']='Suspended between guides, no backing planes';r['loading_assumption']='Constant tension inherited from prior backed benchmark, not calibrated spinning tension.';r['scope']='Equilibrated preparation for a separate spinning stage, not a finite frictional removal or manufacturing history.'
np.savez_compressed(P/'data/zero_equilibrium.npz',x=s.x,w=s.c.w,points=s.c.points(),stock_lengths_m=s.c.ref,roots=np.stack([s.c.origins,s.c.targets],axis=1),radius_m=s.c.radius,endpoint_fixed=s.c.endpoint_fixed,tension_N=s.tension,theta_rad=s.theta,carriage_m=s.d,qualified=r['passed'])
(P/'receipts/zero_equilibrium.json').write_text(json.dumps(r,indent=2));print('RESULT',json.dumps({k:v for k,v in r.items() if k not in ('trace','initial_backed_state_released_to_air','root_forces_N','root_torques_Nm')}),flush=True)
