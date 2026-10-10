"""Curvature-controlled search metric; unchanged physical potential and constraints.
This finds a suspended preparation. Numerical descent is not physical unloading.
"""
import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
import time,json,numpy as np
from prepare_three import ThreeCarriage, P
from carriage import RodHessianOperator,stability,direction,certify_path,log,exp

def prepare(s,x,wall_s=90,prefix='spectral'):
 start=time.monotonic();trace=[];success=False;stable=None
 for it in range(100):
  if time.monotonic()-start>wall_s:break
  E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);con=s.constraints();op=RodHessianOperator(s.c,lam);res=float(abs(pg).max());stable=stability(op,s,J,S);mineig=stable['minimum_generalized_eigenvalue'];row=dict(iteration=it,force_residual=res,minimum_generalized_eigenvalue=mineig,potential_J=E*s.unit,carriage_m=s.d);trace.append(row)
  print('SPECTRAL',it,'res',res,'eig',mineig,'elapsed',time.monotonic()-start,flush=True)
  if res<2e-5 and abs(con).max()*s.L<1e-9 and stable['passed']:success=True;break
  shift=max(0.,-mineig+.1);p0=s.c.points().copy();accepted=False
  for attempt in range(5):
   dx,linear=direction(op,s,g,con,shift);rho=max(1.,1.5*max(float(abs(lam).max()),linear['max_abs_multiplier'])+1);pred=float(pg@(dx-J.T@np.linalg.solve(S,J@dx))+lam@con-rho*np.linalg.norm(con,1))
   if pred>=0 or linear['relative_linear_residual']>1e-7:shift=max(.1,2*shift);continue
   angle=np.linalg.norm(dx[:-1].reshape(s.c.w.shape),axis=2).max();scale=min(1.,.05/max(angle,1e-30))
   for ls in range(20):
    trial=x+scale*dx
    try:
     for _ in range(5):
      s.evaluate(trial);ct=s.constraints();remainder=ct-(1-scale)*con
      if abs(remainder).max()<1e-12:break
      _,_,Jt,St=s.projection(np.zeros_like(g));trial-=Jt.T@np.linalg.solve(St,remainder)
     et,gt=s.evaluate(trial)
     if et+rho*np.linalg.norm(s.constraints(),1)<=E+rho*np.linalg.norm(con,1)+1e-4*scale*pred+1e-13:
      cert=certify_path(p0,s.c.points(),s.c.radius,s.c.ref,plane_bounds=s.c.plane_bounds)
      if cert['passed']:accepted=True;break
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):pass
    scale*=.5
   if accepted:break
   shift=max(.1,2*shift)
  row.update(metric_shift=shift,linear=linear,step_scale=scale,line_search_steps=ls+1,accepted=accepted)
  if not accepted:row['failure']='No physical-energy decreasing contact-safe step';break
  row['contact_certificate']=cert;x=trial;x[:-1]=log(exp(x[:-1].reshape(s.c.w.shape))).ravel();s.evaluate(x)
  target=P/'data'/f'{prefix}_pending.npz';tmp=target.with_suffix('.tmp.npz');np.savez_compressed(tmp,x=x,theta_rad=s.theta,tension_N=s.tension,qualified=False);os.replace(tmp,target)
  (P/'receipts'/f'{prefix}_progress.json').write_text(json.dumps(dict(iteration=it,last=row,qualified=False),indent=2))
 r=s.record(x,success,trace,stable);r.update(wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Suspended static preparation under unchanged benchmark force; no dynamic release or finite frictional manufacturing claim.',numerical_path_contact=all(t.get('contact_certificate',{}).get('passed',False) for t in trace if t.get('accepted')),boundary_planes_m=s.c.plane_bounds)
 # Equilibrium qualification concerns the actual terminal state; the numerical
 # descent is only a sequence of individually certified contact-safe steps.
 r['passed']=bool(r['passed'] and r['relative_stock_error']<1e-12 and r['capsule_penetration_m']<1e-6 and r['numerical_path_contact'])
 np.savez_compressed(P/'data'/f'{prefix}_equilibrium.npz',x=s.x,w=s.c.w,points=s.c.points(),stock_lengths_m=s.c.ref,roots=np.stack([s.c.origins,s.c.targets],axis=1),radius_m=s.c.radius,endpoint_fixed=s.c.endpoint_fixed,tension_N=s.tension,theta_rad=s.theta,carriage_m=s.d,qualified=r['passed'])
 (P/'receipts'/f'{prefix}_equilibrium.json').write_text(json.dumps(r,indent=2));print('RESULT',json.dumps({k:v for k,v in r.items() if k not in ('trace','root_forces_N','root_torques_Nm')}),flush=True)
 return r
if __name__=='__main__':
 s=ThreeCarriage();z=np.load(P/'data/three_preparation_terminal.npz');s.x=z['x'].copy();s.theta=0.;prepare(s,s.x,90,prefix='three_spectral')
