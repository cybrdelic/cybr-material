"""Bounded actual three-rod coordinate test, before any prepared equilibrium."""
from pathlib import Path
import argparse,hashlib,json,time,resource,numpy as np
from flexible_contact import FlexibleContact
from rod_pullback import exp,log,S,mv,right_jacobian,skew
P=Path(__file__).resolve().parents[1]
def fixture(f):
 angle=np.arange(3)*2*np.pi/3;rad=1.99*f.radius/np.sqrt(3);o=np.c_[rad*np.cos(angle),rad*np.sin(angle),np.zeros(3)];s=(np.arange(f.M)+.5)/f.M;w=np.zeros((3,f.M,3));w[:,:,0]=.004*np.sin(2*np.pi*s);w[:,:,1]=.003*np.cos(2*np.pi*s);w[:,:,0]+=2e-5*np.sin(2*np.pi*s+angle[:,None]);return f.pack(o,exp(w))
def run(order):
 start=time.monotonic();f=FlexibleContact(order=order);q0=fixture(f);old=f.state(q0);q=q0.copy();q.reshape(3,-1)[0,2]+=.1e-9/S;z=f.trial(q,old);digest=lambda s:hashlib.sha256(s['q'].tobytes()+s['elastic'].tobytes()).hexdigest();before=digest(old);g=z['state']['g'];rng=np.random.default_rng(527);directions=rng.normal(size=(3,len(q)));directions/=np.linalg.norm(directions,axis=1)[:,None];power=[]
 for v in directions:
  direct=0.
  for a in z['surface']:
   vi,oi=f.velocity(g,v,a['i'],a['s']);vj,oj=f.velocity(g,v,a['j'],a['sb']);vel=vi+np.cross(oi,a['armA'])-vj-np.cross(oj,a['armB']);direct+=np.sum(a['tau']*vel)
  pulled=float(z['friction_residual']@v);power.append(abs(direct-pulled)/max(abs(direct),abs(pulled),1e-18))
 gp=z['friction_gp'];gm=z['friction_gm'];force=gp.sum(axis=(0,1));moment=(np.cross(g['p'],gp).sum(axis=1)+gm.sum(axis=1)).sum(axis=0);en,gn,stats=f.native(q);normalE=abs(en/z['normal_energy_J']-1);normalF=np.linalg.norm(gn-g['normal_gp'])/np.linalg.norm(gn);grad=[];jac=[]
 for v in directions:
  row=[]
  for h in [1e-6,1e-7,1e-8]:
   ep=f.geometry(q+h*v)['normal_energy_J'];em=f.geometry(q-h*v)['normal_energy_J'];numeric=(ep-em)/(2*h);analytic=z['normal_residual']@v;row.append(dict(h=h,relative_error=float(abs(numeric-analytic)/max(abs(analytic),1e-18))))
  grad.append(row);jr=[]
  for h in [1e-7,5e-8,2.5e-8]:
   ts=[f.trial(q+c*h*v,old) for c in [-2,-1,1,2]];J=(ts[0]['residual']-8*ts[1]['residual']+8*ts[2]['residual']-ts[3]['residual'])/(12*h);diff=(ts[2]['residual']-ts[1]['residual'])/(2*h);changes=[dict(active=int(np.count_nonzero(t['active_set']!=z['active_set'])),slip=int(np.count_nonzero(t['sliding_set']!=z['sliding_set'])),projection=int(np.count_nonzero(t['partner_segments']!=z['partner_segments']))) for t in ts];jr.append(dict(h=h,J=J,central_relative=float(np.linalg.norm(J-diff)/max(np.linalg.norm(J),1e-18)),changes=changes))
  agreements=[float(np.linalg.norm(jr[k]['J']-jr[k+1]['J'])/np.linalg.norm(jr[k+1]['J'])) for k in range(len(jr)-1)]
  for r in jr:del r['J']
  jac.append(dict(steps=jr,agreements=agreements))
 U=exp(np.array([.7,-.4,.2]));shift=np.array([.0001,-.0002,.0003]);go=old['g'];qr=f.pack(go['o']@U.T+shift,U@go['R']);rigid=f.trial(qr,old);objold=f.state(qr,old['elastic']@U.T);qnew=f.pack(g['o']@U.T+shift,U@g['R']);obj=f.trial(qnew,objold);gp_error=np.linalg.norm(obj['friction_gp']-gp@U.T)/max(np.linalg.norm(gp),1e-18);gm_error=np.linalg.norm(obj['friction_gm']-gm@U.T)/max(np.linalg.norm(gm),1e-18);length=np.linalg.norm(np.diff(g['p'],axis=1),axis=2);stock=float(np.max(abs(length-f.rod.ds)/f.rod.ds));checks=dict(virtual_power=max(power)<1e-8,net_force=np.linalg.norm(force)<1e-12,net_moment=np.linalg.norm(moment)<1e-15,normal_native=max(normalE,normalF)<.02,normal_gradient=all(min(t['relative_error'] for t in r)<1e-5 for r in grad),jacobian_plateau=all(max(j['agreements'])<1e-4 for j in jac),objectivity=max(gp_error,gm_error)<1e-8,rigid_no_slip=abs(rigid['friction_heat_J'])<1e-20,stock=stock<1e-12,history_immutable=digest(old)==before)
 result=dict(scope='Three flexible rod-coordinate interface; arbitrary test configuration, not prepared equilibrium or spun yarn.',order=order,segments=f.M,checks={k:bool(v) for k,v in checks.items()},passed=bool(all(checks.values())),power_relative=power,net_force_N=float(np.linalg.norm(force)),net_moment_Nm=float(np.linalg.norm(moment)),normal_energy_native_relative=normalE,normal_force_native_relative=normalF,normal_gradient=grad,jacobian=jac,objectivity_force_relative=gp_error,objectivity_moment_relative=gm_error,common_rigid_heat_J=rigid['friction_heat_J'],common_rigid_max_slip_m=max(float(np.linalg.norm(s['du'],axis=1).max(initial=0)) for s in rigid['surface']),stock_relative=stock,elapsed_s=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
 (P/'receipts'/f'rod_coupling_q{order}.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,default=9);a=ap.parse_args();run(a.order)
