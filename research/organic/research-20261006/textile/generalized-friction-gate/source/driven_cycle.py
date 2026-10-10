"""Actual scalar equilibrium with unchanged-history trials and explicit ledger."""
from pathlib import Path
import argparse,hashlib,json,time,resource,numpy as np
from scipy.optimize import brentq
from symmetric_spans import SymmetricSpans
from coupled_spans import RS,S,UNIT,KT,MU
P=Path(__file__).resolve().parents[1]
def digest(s):return hashlib.sha256(s['q'].tobytes()+s['a']['elastic'].tobytes()+s['b']['elastic'].tobytes()).hexdigest()
def atomic(path,obj):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,default=lambda v:v.item()));tmp.replace(path)
def run(n,order=65):
 start=time.monotonic();f=SymmetricSpans(primary_order=order);K=100.;amp=8e-6;z0=RS-.45e-6;q=np.zeros(12);q[2]=z0/S;q[6]=.4e-6/S
 def preparation(x):
  q[0]=x/S;return f.trial(q,f.state(q))['wrenches'][0]+K*x
 xinit=brentq(preparation,-1e-8,1e-8,xtol=1e-20,rtol=1e-14);q[0]=xinit/S;old=f.state(q);initial=f.trial(q,old);old=initial['state'];Uold=.5*K*xinit*xinit;E0=initial['normal_energy_J']+initial['tangential_energy_J']+Uold;uold=0.;Enold=initial['normal_energy_J'];Etold=initial['tangential_energy_J'];rows=[];mutations=0;status='running';evaluations=0;contacts_crossing=0
 controls=[]
 for ua,ub in [(0.,amp),(amp,-amp),(-amp,0.)]:
  controls += [(float(u),z0) for u in np.linspace(ua,ub,n+1)[1:]]
 controls += [(0.,float(z)) for z in np.linspace(z0,RS+.1e-6,n+1)[1:]]
 tag=f'driven_cycle_{n:03d}'+('' if order==65 else f'_q{order}');out=P/'receipts'/f'{tag}.json';statepath=P/'receipts'/f'{tag}_accepted.npz'
 for i,(u,z) in enumerate(controls,1):
  if time.monotonic()-start>112:status='bounded_pause';break
  olddigest=digest(old);qt=old['q'].copy();qt[2]=z/S;calls=0
  def fun(x):
   nonlocal calls,evaluations
   calls+=1;evaluations+=1;qt[0]=x/S;t=f.trial(qt,old);return t['wrenches'][0]+K*(x-u)
  width=20e-6;lo=u-width;hi=u+width;fl=fun(lo);fh=fun(hi)
  if fl*fh>0:status='failed_root_bracket';break
  x=brentq(fun,lo,hi,xtol=1e-17,rtol=1e-14,maxiter=60);qt[0]=x/S;t=f.trial(qt,old);fr=float(t['wrenches'][0]+K*(x-u));mutation=digest(old)!=olddigest;mutations+=int(mutation);h=1e-12;slopes=[float((fr-fun(x-h))/h),float((fun(x+h)-fr)/h)];qt[0]=x/S
  dx=x-old['q'][0]*S;dz=z-old['q'][2]*S;du=u-uold;Ud=.5*K*(x-u)**2;E=t['normal_energy_J']+t['tangential_energy_J']+Ud;dE=E-(Enold+Etold+Uold);W=K*(u-x)*du+t['wrenches'][2]*dz;An=t['normal_wrenches'][0]*dx+t['normal_wrenches'][2]*dz-(t['normal_energy_J']-Enold);Ad=.5*K*((x-u)-(old['q'][0]*S-uold))**2;D=t['friction_heat_J'];At=t['numerical_loss_J'];gap=W-dE-D-At-An-Ad;scale=max(abs(W),abs(dE),D+At+abs(An)+Ad,1e-18);rel=abs(gap)/scale
  crossed=0
  for side in ['a','b']:
   go=old[side]['geometry'];gn=t['state'][side]['geometry'];crossed+=int(np.count_nonzero((go['N']>0)&(gn['N']>0)&(go['segment']!=gn['segment'])))
  row=dict(step=i,accepted_history_digest=olddigest,trial_history_digest=digest(t['state']),one_sided_restoring_N_m=slopes,normal_wrenches=t['normal_wrenches'].tolist(),friction_wrenches=t['friction_wrenches'].tolist(),control_u_m=u,normal_z_m=z,x_m=x,force_residual_N=fr,force_residual_scaled=fr/1e-4,evaluations=calls,normal_energy_J=t['normal_energy_J'],tangential_energy_J=t['tangential_energy_J'],drive_energy_J=Ud,external_endpoint_work_J=W,stored_increment_J=dE,friction_heat_J=D,return_map_numerical_loss_J=At,normal_endpoint_defect_J=An,drive_endpoint_defect_J=Ad,finite_kinematic_ledger_gap_J=gap,ledger_relative=rel,active=int(t['active_set'].sum()),sliding=int(t['sliding_set'].sum()),material_segment_crossings=crossed,history_unchanged=not mutation)
  accepted=abs(fr)<1e-10 and abs(fr)/1e-4<1e-6 and rel<.02 and not mutation and min(slopes)>0
  row['accepted']=accepted;rows.append(row)
  if not accepted:status='rejected_force_or_work_gate';break
  contacts_crossing+=crossed;old=t['state'];uold=u;Uold=Ud;Enold=t['normal_energy_J'];Etold=t['tangential_energy_J'];np.savez_compressed(statepath,q=old['q'],elastic_a=old['a']['elastic'],elastic_b=old['b']['elastic'],identity=old['identity'],control_u=u,step=i)
  if i%n==0:np.savez_compressed(P/'receipts'/f'{tag}_leg_{i//n}.npz',q=old['q'],elastic_a=old['a']['elastic'],elastic_b=old['b']['elastic'],identity=old['identity'],control_u=u,step=i)
  if i%8==0:atomic(out,dict(status='running',n=n,rows=rows,elapsed_s=time.monotonic()-start))
 if status=='running':status='completed'
 acceptedrows=[z for z in rows if z['accepted']];W=sum(z['external_endpoint_work_J'] for z in acceptedrows);D=sum(z['friction_heat_J'] for z in acceptedrows);A=sum(z['return_map_numerical_loss_J']+abs(z['normal_endpoint_defect_J'])+z['drive_endpoint_defect_J'] for z in acceptedrows);G=sum(z['finite_kinematic_ledger_gap_J'] for z in acceptedrows);absW=sum(abs(z['external_endpoint_work_J']) for z in acceptedrows);den=max(absW,1e-18)
 result=dict(status=status,n=n,primary_order=order,initial_x_m=xinit,initial_force_residual_N=float(initial['wrenches'][0]+K*xinit),accepted_steps=len(acceptedrows),planned_steps=len(controls),elapsed_s=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,evaluations=evaluations,initial_stored_J=E0,external_work_J=W,absolute_external_work_J=absW,physical_friction_heat_J=D,absolute_numerical_terms_J=A,numerical_terms_relative=A/den,separated_ledger_gap_J=G,ledger_relative=abs(G)/den,material_segment_crossings=contacts_crossing,history_mutations=mutations,rows=rows)
 result['passed']=status=='completed' and abs(result['initial_force_residual_N'])<1e-10 and result['ledger_relative']<.02 and result['numerical_terms_relative']<.02 and rows[-1]['active']==0 and mutations==0
 atomic(out,result);print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2,default=lambda v:v.item()),flush=True);return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=16);ap.add_argument('--order',type=int,default=65);args=ap.parse_args();run(args.n,args.order)
