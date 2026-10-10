import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import time,json,numpy as np,argparse,hashlib
from tangent_suspended import TangentSuspended,solve,P
from carriage import exp
ap=argparse.ArgumentParser();ap.add_argument('--state',default=str(P/'data/positive_equilibrium.npz'));ap.add_argument('--receipt',default=str(P/'receipts/qualified_suspended_restart.json'));ap.add_argument('--batch',default='ramp01');ap.add_argument('--seconds',type=float,default=120);ap.add_argument('--step',type=float,default=.02);ap.add_argument('--max-step',type=float,default=.1);ap.add_argument('--stop-theta',type=float,default=None);args=ap.parse_args()
s=TangentSuspended();z=np.load(args.state);s.x=z['x'].copy();s.theta=float(z['theta_rad']);s.evaluate(s.x);previous=json.loads(__import__('pathlib').Path(args.receipt).read_text());BATCH=P/args.batch;(BATCH/'data').mkdir(parents=True,exist_ok=True);(BATCH/'receipts').mkdir(parents=True,exist_ok=True);start=time.monotonic();step=args.step;R=float(np.max(np.linalg.norm(s.c.offset[:,:2],axis=1)));reference_goal=.1*.0023797191653033645/R;goal=reference_goal if args.stop_theta is None else min(reference_goal,args.stop_theta);milestone=.8395831849955003;attempts=[];accepted=[];latest=args.state;status='running';initial_energy=previous['stored_energy_J'];total_work=0.
while s.theta<goal-1e-12:
 remaining=args.seconds-(time.monotonic()-start)
 if remaining<12:status='paused at runtime bound';break
 oldx=s.x.copy();oldtheta=s.theta;target=min(goal,oldtheta+step)
 if oldtheta<milestone-1e-12:target=min(target,milestone)
 def pending(x,theta,trace):
  f=BATCH/'data/pending_trial.npz';tmp=f.with_suffix('.tmp.npz');np.savez_compressed(tmp,x=x,theta_rad=theta,tension_N=s.tension,qualified=False);os.replace(tmp,f)
 r=solve(s,target,wall_s=min(20,remaining-3),checkpoint=pending)
 work=.5*(previous['torque_conjugate_Nm']+r['torque_conjugate_Nm'])*(r['twist_rad']-oldtheta)+s.tension*(r['carriage_m']-previous['carriage_m']);de=r['stored_energy_J']-previous['stored_energy_J'];gap=abs(work-de)/max(abs(work),abs(de),1e-18);curv=max(s.c.rod.diagnostics(exp(s.c.w[i]))['max_radius_curvature'] for i in range(s.c.N));good=bool(r['passed'] and r.get('endpoint_path_contact',{}).get('passed',False) and r['relative_stock_error']<1e-12 and r['capsule_penetration_m']<1e-6 and curv<.1 and gap<.02)
 pts=s.c.points();t=np.diff(pts,axis=1)/s.c.ref[:,:,None];angles=np.degrees(np.arccos(np.clip(t[:,:,0],-1,1)));r.update(step_work_J=work,stored_change_J=de,relative_work_gap=gap,history_accepted=good,resume_state=str(latest),requested_theta=target,bearing_turns=r['twist_rad']/(2*np.pi),maximum_radius_curvature=curv,actual_filament_angle_degrees=dict(median=float(np.median(angles)),p95=float(np.percentile(angles,95)),maximum=float(angles.max())),angle_definition='Actual local filament tangent angle to fixed carriage axis; nominal helix angle is separately derived from root omega R.',scope='Suspended normal-contact quasistatic twist after an independently equilibrated static preparation. Fixed stock and benchmark axial force. No finite friction, wool calibration or physical support-removal history.')
 name=('accepted' if good else 'rejected')+f'_{len(attempts)+1:03d}';file=BATCH/'data'/f'{name}.npz';np.savez_compressed(file,x=s.x,w=s.c.w,points=pts,stock_lengths_m=s.c.ref,roots=np.stack([s.c.origins,s.c.targets],axis=1),radius_m=s.c.radius,endpoint_fixed=s.c.endpoint_fixed,tension_N=s.tension,theta_rad=s.theta,carriage_m=s.d,qualified=good);r['state_sha256']=hashlib.sha256(file.read_bytes()).hexdigest();(BATCH/'receipts'/f'{name}.json').write_text(json.dumps(r,indent=2))
 attempts.append(dict(name=name,target=target,accepted=good,force=r['force_residual'],gap=gap,wall_s=r.get('wall_s'),stability=r['stability']))
 if good:
  total_work+=work;latest=str(file);previous=r;accepted.append(name)
  if gap<.004 and r.get('wall_s',99)<10:step=min(args.max_step,step*1.35)
  print('ACCEPTED',name,s.theta,'tpm',r['twist_turns_per_m'],'omegaR',r['dimensionless_root_twist'],'gap',gap,'elapsed',time.monotonic()-start,flush=True)
 else:
  s.x=oldx;s.theta=oldtheta;s.evaluate(oldx);step*=.5;print('REJECTED',name,target,'force',r['force_residual'],'gap',gap,'next_step',step,flush=True)
  if step<.001:status='paused at unresolved small control interval';break
 summary=dict(status=status,accepted_steps=len(accepted),latest_accepted_state=latest,achieved_theta_rad=s.theta,attempts=attempts,elapsed_s=time.monotonic()-start)
 tmp=BATCH/'receipts/progress.tmp.json';tmp.write_text(json.dumps(summary,indent=2));os.replace(tmp,BATCH/'receipts/progress.json')
else:status='requested bearing-angle endpoint reached'
change=previous['stored_energy_J']-initial_energy;summary=dict(status=status,accepted_steps=len(accepted),latest_accepted_state=latest,achieved_theta_rad=s.theta,achieved_turns_per_m=previous['twist_turns_per_m'],achieved_omega_R=previous['dimensionless_root_twist'],attempts=attempts,elapsed_s=time.monotonic()-start,cumulative_boundary_work_J=total_work,cumulative_stored_change_J=change,cumulative_work_relative_gap=abs(total_work-change)/max(abs(total_work),abs(change),1e-18),peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,initial_state=args.state,normal_only=True,goal_reference_theta=reference_goal,run_endpoint_theta=goal,maximum_numerical_increment_rad=args.max_step,next_numerical_step_rad=step)
(BATCH/'receipts/progress.json').write_text(json.dumps(summary,indent=2));print('RAMP_STATUS',json.dumps(summary),flush=True)
