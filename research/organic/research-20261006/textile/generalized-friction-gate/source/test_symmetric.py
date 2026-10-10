import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
import json,time,numpy as np
from coupled_spans import exp,log,right_jacobian,S,RS,KT,L,UNIT,P
from symmetric_spans import SymmetricSpans,swap,reorder
start=time.monotonic();f=SymmetricSpans();q0=np.zeros(12);q0[2]=(RS-.45e-6)/S;q0[6]=.4e-6/S;old=f.state(q0);exchanged=f.exchange_state(old);exchange=[]
for dx,dz in [(10e-9,0),(-10e-9,0),(-10e-9,-.05e-6)]:
 q=q0.copy();q[0]+=dx/S;q[2]+=dz/S;a=f.trial(q,old);b=f.trial(swap(q),exchanged);exchange.append(dict(force=float(np.linalg.norm(a['wrenches']-reorder(b['wrenches']))/np.linalg.norm(a['wrenches'])),heat=abs(a['friction_heat_J']-b['friction_heat_J'])/max(abs(a['friction_heat_J']),1e-30),work=abs(a['tangential_endpoint_work_J']-b['tangential_endpoint_work_J'])/max(abs(a['tangential_endpoint_work_J']),1e-30)));old=a['state'];exchanged=b['state']
# Generalized power/Jacobian tests at a nonzero populated history.
q0[3:6]=[.03,.08,-.05];q0[9:12]=[.02,-.03,.01];old=f.state(q0);seed=q0.copy();seed[0]+=.5e-9/S;old=f.trial(seed,old)['state'];oldbytes=old['a']['elastic'].tobytes()+old['b']['elastic'].tobytes();q=seed.copy();q[0]+=20e-9/S;q[2]-=.04e-6/S;q[4]+=.001;r=f.trial(q,old);g=r['geometry'];w=r['wrenches'];A=g['Ra']@right_jacobian(q[3:6]);B=g['Rb']@right_jacobian(q[9:12]);rng=np.random.default_rng(790);virtual=[]
for spin_only in [False,True]:
 v=rng.normal(size=12)
 if spin_only:v[:3]=0;v[6:9]=0
 Va=S*v[:3];wa=A@v[3:6];Vb=S*v[6:9];wb=B@v[9:12];pa=r['contact_pa'];pb=r['contact_pb'];m=.5*(pa+pb);va=Va+np.cross(wa,pa-g['Ta']);vb=Vb+np.cross(wb,pb-g['Tb']);surface=va+np.cross(wa,m-pa)-vb-np.cross(wb,m-pb);power=np.sum(r['tau']*surface);wf=r['friction_wrenches'];general=wf[:3]@Va+wf[3:6]@wa+wf[6:9]@Vb+wf[9:12]@wb;virtual.append(float(abs(power-general)/max(abs(power),abs(general),1e-30)))
force_balance=float(np.linalg.norm(w[:3]+w[6:9]));moment_balance=float(np.linalg.norm(w[3:6]+np.cross(g['Ta'],w[:3])+w[9:12]+np.cross(g['Tb'],w[6:9])));ne,nw=f.native(q);normal_native=dict(energy=abs(r['normal_energy_J']-ne)/ne,force=np.linalg.norm(r['normal_wrenches']-nw)/np.linalg.norm(nw))
def transform(q,H,b):
 z=q.copy();z[:3]=(H@(q[:3]*S)+b)/S;z[3:6]=log(H@exp(q[3:6]));z[6:9]=(H@(q[6:9]*S)+b)/S;z[9:12]=log(H@exp(q[9:12]));return z
H0=exp([.7,-.5,.9]);H1=exp([-.4,.8,1.1]);ot=f.state(transform(old['q'],H0,np.array([.001,-.002,.001])),old['a']['elastic']@H0.T,old['b']['elastic']@H0.T);rt=f.trial(transform(q,H1,np.array([-.001,.001,.002])),ot);expected=np.concatenate([H1@w[i:i+3] for i in range(0,12,3)]);obj=float(np.linalg.norm(rt['wrenches']-expected)/np.linalg.norm(expected));heatobj=abs(rt['friction_heat_J']-r['friction_heat_J'])/max(r['friction_heat_J'],1e-30)
Js=[];audits=[]
for h in [2.5e-7,1.25e-7,6.25e-8,3.125e-8,1.5625e-8]:
 Js.append(f.jacobian(q,old,h=h));audits.append(dict(step=h,changes=f.last_jacobian_audit,switched_perturbations=sum(any(x[c]>0 for c in ['active','sliding']) for col in f.last_jacobian_audit for x in col['changes'])))
agreement=[float(np.linalg.norm(b-a)/np.linalg.norm(b)) for a,b in zip(Js[:-1],Js[1:])];directional=[]
for i in range(3):
 v=rng.normal(size=12);v/=np.linalg.norm(v);rows=[]
 for h in [1e-6,1e-7,1e-8]:
  plus=f.trial(q+h*v,old);minus=f.trial(q-h*v,old);fd=(plus['residual']-minus['residual'])/(2*h);rows.append(dict(step=h,error=float(np.linalg.norm(fd-Js[-1]@v)/np.linalg.norm(Js[-1]@v)),active_changes=int(np.count_nonzero(plus['active_set']!=minus['active_set'])),sliding_changes=int(np.count_nonzero(plus['sliding_set']!=minus['sliding_set']))))
 directional.append(rows)
# Coextensive parallel spans: the half-weighted closure must not double k or muN.
p=np.zeros(12);p[2]=(RS-.45e-6)/S;p[3:6]=[0,0,-np.pi/2];po=f.state(p);tiny=p.copy();tiny[0]+=1e-10/S;small=f.trial(tiny,po);stiffness_error=abs(small['friction_wrenches'][0]-KT*L*1e-10)/(KT*L*1e-10);large=p.copy();large[0]+=1e-7/S;slip=f.trial(large,po);capacity=.3*abs(slip['normal_wrenches'][2]);capacity_error=abs(slip['friction_wrenches'][0]-capacity)/capacity
roll=p.copy();roll[3:6]=log(exp([.8,0,0])@exp(p[3:6]));roll[9:12]=[-.8,0,0];rolled=f.trial(roll,po);rolling_slip=float(np.max(np.linalg.norm(rolled['du'],axis=1),initial=0));rolling_heat=rolled['friction_heat_J'];asym=float(np.linalg.norm(Js[-1]-Js[-1].T)/np.linalg.norm(Js[-1]))
checks=dict(exchange=max(max(x.values()) for x in exchange)<1e-8,normal_native=max(normal_native.values())<.005,virtual_power=max(virtual)<1e-8,action_reaction=force_balance<1e-12 and moment_balance<1e-15,objectivity=obj<1e-8 and heatobj<1e-8,stable_fine_jacobian=max(agreement[-2:])<1e-4 and all(a['switched_perturbations']==0 for a in audits[-3:]),directional_jacobian=max(min(z['error'] for z in row if z['active_changes']==0 and z['sliding_changes']==0) for row in directional)<1e-4,parallel_stiffness=stiffness_error<1e-5,parallel_capacity=capacity_error<1e-5,parallel_rolling=rolling_slip<1e-14 and rolling_heat<1e-22,history_immutable=oldbytes==old['a']['elastic'].tobytes()+old['b']['elastic'].tobytes(),nonsymmetric=asym>1e-5);checks={k:bool(v) for k,v in checks.items()};out=dict(passed=all(checks.values()),checks=checks,exchange_history=exchange,normal_native=normal_native,virtual_power_errors=virtual,force_balance_N=force_balance,moment_balance_Nm=moment_balance,objectivity_error=obj,heat_objectivity_error=heatobj,jacobian_agreement=agreement,jacobian_active_set_audit=audits,directional_jacobian=directional,jacobian_asymmetry=asym,parallel_stiffness_relative_error=float(stiffness_error),parallel_capacity_relative_error=float(capacity_error),parallel_rolling_surface_increment_m=rolling_slip,parallel_rolling_heat_J=rolling_heat,wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Explicit symmetric two-sided half-weighted constitutive benchmark, original normal double energy counted once; rigid spans only.');(P/'receipts/symmetric_tests.json').write_text(json.dumps(out,indent=2));print('SUMMARY',json.dumps({k:v for k,v in out.items() if k!='jacobian_active_set_audit'}),flush=True)
