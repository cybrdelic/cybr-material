"""Physical Newton residual correction of a suspended static preparation.
No claim that the numerical line search is a physical dissipative trajectory.
"""
import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import time,json,numpy as np
from carriage import Carriage,RodHessianOperator,stability,direction,certify_path,log,exp,P
s=Carriage();z=np.load(P/'data/spectral_equilibrium.npz');x=z['x'].copy();s.theta=float(z['theta_rad']);start=time.monotonic();trace=[];success=False;stable=None
for it in range(100):
 if time.monotonic()-start>45:break
 E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);con=s.constraints();op=RodHessianOperator(s.c,lam);res=float(abs(pg).max());merit=float(pg@pg);row=dict(iteration=it,force_residual=res,force_l2_merit=merit,potential_J=E*s.unit,carriage_m=s.d);trace.append(row);print('FORCE',it,res,'elapsed',time.monotonic()-start,flush=True)
 if res<2e-5 and abs(con).max()*s.L<1e-9:
  stable=stability(op,s,J,S);success=stable['passed'];break
 dx,linear=direction(op,s,g,con,0.);angle=np.linalg.norm(dx[:-1].reshape(s.c.w.shape),axis=2).max();scale=min(1.,.05/max(angle,1e-30));p0=s.c.points().copy();accepted=False
 for ls in range(24):
  trial=x+scale*dx
  try:
   for _ in range(5):
    s.evaluate(trial);ct=s.constraints();remainder=ct-(1-scale)*con
    if abs(remainder).max()<1e-12:break
    _,_,Jt,St=s.projection(np.zeros_like(g));trial-=Jt.T@np.linalg.solve(St,remainder)
   et,gt=s.evaluate(trial);pt,*_=s.projection(gt)
   if pt@pt<=(1-1e-4*scale)*merit:
    cert=certify_path(p0,s.c.points(),s.c.radius,s.c.ref,plane_bounds=None)
    if cert['passed']:accepted=True;break
  except (ValueError,FloatingPointError,np.linalg.LinAlgError):pass
  scale*=.5
 row.update(linear=linear,step_scale=scale,line_search_steps=ls+1,accepted=accepted)
 if not accepted:row['failure']='No force-decreasing contact-safe step';break
 row['contact_certificate']=cert;x=trial;x[:-1]=log(exp(x[:-1].reshape(s.c.w.shape))).ravel();s.evaluate(x)
 target=P/'data/force_pending.npz';tmp=target.with_suffix('.tmp.npz');np.savez_compressed(tmp,x=x,theta_rad=s.theta,tension_N=s.tension,qualified=False);os.replace(tmp,target)
 (P/'receipts/force_progress.json').write_text(json.dumps(dict(last=row,qualified=False),indent=2))
r=s.record(x,success,trace,stable);r.update(wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Suspended static preparation corrected using physical-Newton force merit; no unloading trajectory or finite frictional manufacturing claim.',numerical_path_contact=all(t.get('contact_certificate',{}).get('passed',False) for t in trace if t.get('accepted')))
r['passed']=bool(r['passed'] and r['relative_stock_error']<1e-12 and r['capsule_penetration_m']<1e-6 and r['numerical_path_contact'])
np.savez_compressed(P/'data/force_equilibrium.npz',x=s.x,w=s.c.w,points=s.c.points(),stock_lengths_m=s.c.ref,roots=np.stack([s.c.origins,s.c.targets],axis=1),radius_m=s.c.radius,endpoint_fixed=s.c.endpoint_fixed,tension_N=s.tension,theta_rad=s.theta,carriage_m=s.d,qualified=r['passed'])
(P/'receipts/force_equilibrium.json').write_text(json.dumps(r,indent=2));print('RESULT',json.dumps({k:v for k,v in r.items() if k not in ('trace','root_forces_N','root_torques_Nm')}),flush=True)
