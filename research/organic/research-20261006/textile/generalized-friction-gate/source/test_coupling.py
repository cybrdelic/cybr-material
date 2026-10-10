import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
import json,time,numpy as np
from coupled_spans import Spans,exp,log,right_jacobian,S,UNIT,L,RS,KC,KT,MU,P
start=time.monotonic();f=Spans();q0=np.zeros(12);q0[2]=(RS-.45e-6)/S;q0[6]=.4e-6/S;q0[3:6]=[.03,.08,-.05];q0[9:12]=[.02,-.03,.01];old=f.state(q0);n=old['geometry']['n'];v=np.broadcast_to(np.array([1.,0,0]),n.shape).copy();v-=n*np.sum(v*n,axis=1)[:,None];v/=np.linalg.norm(v,axis=1)[:,None];old['elastic']=.2*MU*old['geometry']['N'][:,None]/f.k[:,None]*v;oldhash=old['elastic'].tobytes();q=q0.copy();q[0]+=20e-9/S;q[2]-=.04e-6/S;q[4]+=.001;r=f.trial(q,old);g=r['state']['geometry'];nativeE,nativeW=f.native(q);normal_native=dict(energy=abs(g['E']-nativeE)/abs(nativeE),force_moment=np.linalg.norm(g['normal_wrenches']-nativeW)/np.linalg.norm(nativeW));rng=np.random.default_rng(790);direction=rng.normal(size=12);direction/=np.linalg.norm(direction);A=g['Ra']@right_jacobian(q[3:6]);B=g['Rb']@right_jacobian(q[9:12]);normal_g=np.r_[S*g['normal_wrenches'][:3],A.T@g['normal_wrenches'][3:6],S*g['normal_wrenches'][6:9],B.T@g['normal_wrenches'][9:12]]/UNIT;normal_fd=[]
for h in [1e-4,1e-5,1e-6]:
 fd=(f.geometry(q+h*direction)['E']-f.geometry(q-h*direction)['E'])/(2*h*UNIT);normal_fd.append(dict(step=h,relative_error=abs(fd-normal_g@direction)/max(abs(normal_g@direction),1e-30)))
# Surface-spin virtual power is evaluated directly at material centreline points.
virtual=[]
for translation in [True,False]:
 vq=rng.normal(size=12)
 if not translation:vq[:3]=0;vq[6:9]=0
 Va=S*vq[:3];wa=A@vq[3:6];Vb=S*vq[6:9];wb=B@vq[9:12];ids=r['ids'];m=.5*(g['pa'][ids]+g['pb'][ids]);va=Va+np.cross(wa,g['pa'][ids]-g['Ta']);vb=Vb+np.cross(wb,g['pb'][ids]-g['Tb']);surface=va+np.cross(wa,m-g['pa'][ids])-vb-np.cross(wb,m-g['pb'][ids]);power=float(np.sum(r['tau']*surface));FW=r['friction_wrenches'];general=float(FW[:3]@Va+FW[3:6]@wa+FW[6:9]@Vb+FW[9:12]@wb);virtual.append(dict(surface_spin_only=not translation,relative_error=abs(power-general)/max(abs(power),abs(general),1e-30)))
w=r['wrenches'];force_balance=float(np.linalg.norm(w[:3]+w[6:9]));moment_balance=float(np.linalg.norm(w[3:6]+np.cross(g['Ta'],w[:3])+w[9:12]+np.cross(g['Tb'],w[6:9])))
def transform(q,H,b):
 z=q.copy();z[:3]=(H@(q[:3]*S)+b)/S;z[3:6]=log(H@exp(q[3:6]));z[6:9]=(H@(q[6:9]*S)+b)/S;z[9:12]=log(H@exp(q[9:12]));return z
H0=exp([.7,-.5,.9]);H1=exp([-.4,.8,1.1]);ot=f.state(transform(q0,H0,np.array([.001,-.002,.001])),old['elastic']@H0.T);rt=f.trial(transform(q,H1,np.array([-.001,.001,.002])),ot);expected=np.concatenate([H1@w[i:i+3] for i in range(0,12,3)]);objectivity=np.linalg.norm(rt['wrenches']-expected)/np.linalg.norm(expected);heat_objectivity=abs(rt['friction_heat_J']-r['friction_heat_J'])/max(r['friction_heat_J'],1e-30)
rigid=f.trial(transform(q0,H1,np.array([.001,.002,-.001])),old);rigid_slip=float(np.max(np.linalg.norm(rigid['du'],axis=1)));rigid_heat=rigid['friction_heat_J']
J1=f.jacobian(q,old,h=1e-6);J2=f.jacobian(q,old,h=5e-7);J3=f.jacobian(q,old,h=2.5e-7);jac_convergence=[float(np.linalg.norm(J2-J1)/np.linalg.norm(J2)),float(np.linalg.norm(J3-J2)/np.linalg.norm(J3))];directions=[]
for i in range(3):
 v=rng.normal(size=12);v/=np.linalg.norm(v);rows=[]
 for h in [1e-4,1e-5,1e-6]:
  fd=(f.trial(q+h*v,old)['residual']-f.trial(q-h*v,old)['residual'])/(2*h);rows.append(dict(step=h,relative_error=float(np.linalg.norm(fd-J3@v)/np.linalg.norm(J3@v))))
 directions.append(rows)
asym=float(np.linalg.norm(J3-J3.T)/np.linalg.norm(J3))
# Independent straight-span translating limits use the analytic Abel reduction:
# each primary point has N_i = k_n*w_i*(Rsum-distance_to_infinite_partner).
straight=np.zeros(12);straight[2]=(RS-.45e-6)/S;straight[6]=.4e-6/S;control=f.state(straight);dist=np.linalg.norm(control['geometry']['pa']-control['geometry']['pb'],axis=1);analyticN=KC*f.w*np.maximum(RS-dist,0);limits=[]
for dx in [1e-10,1e-9,1e-8,1e-7]:
 trial=straight.copy();trial[0]+=dx/S;rr=f.trial(trial,control);expected=float(np.minimum(f.k*dx,MU*analyticN).sum());actual=float(rr['friction_wrenches'][0]);limits.append(dict(displacement_m=dx,expected_force_N=expected,actual_force_N=actual,relative_error=abs(actual-expected)/max(abs(expected),1e-30)))
checks=dict(normal_native=max(normal_native.values())<.005,normal_energy_gradient=min(z['relative_error'] for z in normal_fd)<1e-5,virtual_power=max(x['relative_error'] for x in virtual)<1e-8,action_reaction=force_balance<1e-12 and moment_balance<1e-15,objectivity=objectivity<1e-8 and heat_objectivity<1e-8,finite_rigid_motion=rigid_slip<1e-14 and rigid_heat<1e-22,jacobian_directional=max(min(a['relative_error'] for a in rows) for rows in directions)<1e-4,jacobian_step_agreement=max(jac_convergence)<1e-4,analytic_translating_limits=max(x['relative_error'] for x in limits)<1e-5,history_immutable=oldhash==old['elastic'].tobytes(),nonsymmetry_measured=asym>1e-5)
checks={k:bool(v) for k,v in checks.items()}
out=dict(passed=all(checks.values()),checks=checks,normal_native=normal_native,normal_energy_gradient=normal_fd,virtual_power=virtual,force_balance_N=force_balance,moment_balance_Nm=moment_balance,objectivity_relative_error=float(objectivity),heat_objectivity_relative_error=float(heat_objectivity),finite_rigid_surface_increment_m=rigid_slip,finite_rigid_heat_J=rigid_heat,jacobian_relative_asymmetry=asym,jacobian_step_agreement=jac_convergence,directional_jacobian= directions,analytic_translating_limits=limits,wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='12 rigid-span DOFs; numerical full residual Jacobian, not an independent analytic derivative. Directional step convergence plus analytic translating limits. No flexible-yarn equilibrium or appearance qualification.');(P/'receipts/coupling_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
