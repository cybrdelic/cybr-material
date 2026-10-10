import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np
from objective_pair import make as old_make,advance as old_advance,exp,P
from fixed_stiffness_pair import make,advance
I=np.eye(3);d=32e-6;pi=np.array([-d/2,0.,0.]);pj=-pi;labels=((0,'arc',.001),(1,'arc',.001));N=1e-6;mu=.3;k=100.;measure=4e-10;kd=k/measure;e=np.array([0.,mu*N/k,0.]);base=make(labels,pi,pj,I,I,normal=N,mu=mu,stiffness_density_N_m3=kd,reference_pair_measure_m2=measure,elastic=e);U0=.5*k*(e@e)
# Increasing N does not alter a sticking spring's energy or traction.
raised,ri=advance(base,pi,pj,I,I,labels=labels,normal=2*N)
ramps=[]
for steps in [1,8,16,32,64,128,256]:
 state=base;heat=loss=work=0.;balance=0.
 for j in range(1,steps+1):
  state,r=advance(state,pi,pj,I,I,labels=labels,normal=N*(1-j/steps));heat+=r['friction_heat_J'];loss+=r['numerical_loss_J'];work+=r['endpoint_work_J'];balance=max(balance,abs(r['work_balance_residual_J']))
 ramps.append(dict(steps=steps,initial_elastic_J=U0,final_elastic_J=r['stored_J'],friction_heat_J=heat,numerical_loss_J=loss,numerical_loss_fraction=loss/U0,endpoint_work_J=work,balance_error_J=abs(work-(r['stored_J']-U0)-heat-loss),maximum_local_residual_J=balance,contact_active=r['contact_active']))
# Pure open-contact motion leaves no fictitious tangential memory; recontact at
# zero incremental motion has no force. Crossing-time integration remains open.
openstate,_=advance(base,pi,pj,I,I,labels=labels,normal=0.)
openstate,openr=advance(openstate,pi+[0,2e-6,0],pj,I,I,labels=labels,normal=0.)
recontact,re=advance(openstate,openstate.pi,openstate.pj,I,I,labels=labels,normal=N)
# Quadrature subdivision: both normal load and stiffness use reference measure.
partitions=[]
for count in [1,2,7,32]:
 sums=np.zeros(5)
 for j in range(count):
  state=make(((0,j),(1,j)),pi,pj,I,I,normal=N/count,mu=mu,stiffness_density_N_m3=kd,reference_pair_measure_m2=measure/count)
  state,r=advance(state,pi+[0,12e-9,0],pj,I,I,labels=state.labels,normal=N/count);sums+=np.array([r['stored_J'],r['friction_heat_J'],r['numerical_loss_J'],r['endpoint_work_J'],np.linalg.norm(r['traction_N'])])
 partitions.append(sums)
partition_error=float(max(np.linalg.norm(v-partitions[0])/max(np.linalg.norm(partitions[0]),1e-30) for v in partitions))
# Old constant-N control agrees when the fixed k equals N/delta0 initially.
old=old_make(labels,pi,pj,I,I,N,mu,N/k,elastic=.5*e);new=make(labels,pi,pj,I,I,normal=N,mu=mu,elastic=.5*e);Ri=exp([.05,-.03,.1]);Rj=exp([-.04,.06,-.02]);o,ro=old_advance(old,pi+[0,2e-9,0],pj,Ri,Rj,labels=labels,normal=N);n,rn=advance(new,pi+[0,2e-9,0],pj,Ri,Rj,labels=labels,normal=N);control_error=float(np.linalg.norm(o.elastic-n.elastic)/max(np.linalg.norm(o.elastic),1e-30))
# Arbitrary independent superposed frames, also during normal unloading.
A=exp([.7,-.2,.4]);B=exp([-.8,.9,.3]);b0=np.array([.001,.002,-.001]);b1=np.array([-.002,.001,.002]);transformed=make(labels,A@pi+b0,A@pj+b0,A,A,normal=N,elastic=A@(.5*e));physical=make(labels,pi,pj,I,I,normal=N,elastic=.5*e);one,r1=advance(physical,pi+[0,2e-9,0],pj,Ri,Rj,labels=labels,normal=.4*N);two,r2=advance(transformed,B@(pi+[0,2e-9,0])+b1,B@pj+b1,B@Ri,B@Rj,labels=labels,normal=.4*N);objectivity=float(np.linalg.norm(two.elastic-B@one.elastic)/max(np.linalg.norm(one.elastic),1e-30));heat_objectivity=abs(r2['friction_heat_J']-r1['friction_heat_J'])/max(r1['friction_heat_J'],1e-30)
out=dict(scope='Explicit fixed-SI-stiffness alternative, k_pair=k_density*reference material-pair measure. Constitutive/kinematic tests only; moving quadrature, contact birth timing, normal-load mechanics and coupled rods remain unqualified.',parameters=dict(mu=mu,normal_reference_N=N,stiffness_density_N_m3=kd,reference_pair_measure_m2=measure,stiffness_N_m=k),normal_increase_stored_change_J=ri['stored_change_J'],normal_increase_heat_J=ri['friction_heat_J'],normal_unloading=ramps,abrupt_loss_interpretation='One-step complete unloading puts all initial spring energy into numerical loss. It is not friction heat or a qualified physical jump.',resolved_unloading_256_relative_numerical_loss=ramps[-1]['numerical_loss_fraction'],open_contact_memory_norm_m=float(np.linalg.norm(openstate.elastic)),recontact_force_N=float(np.linalg.norm(re['traction_N'])),reference_measure_partition_relative_error=partition_error,old_constant_normal_control_relative_error=control_error,superposed_motion_relative_error=objectivity,superposed_heat_relative_error=heat_objectivity)
out['passed']=bool(abs(ri['stored_change_J'])<1e-28 and ri['friction_heat_J']==0 and ramps[-1]['numerical_loss_fraction']<.02 and abs(ramps[0]['numerical_loss_fraction']-1)<1e-12 and all(a['balance_error_J']<1e-28 for a in ramps) and partition_error<1e-12 and control_error<1e-12 and objectivity<1e-8 and heat_objectivity<1e-8 and np.linalg.norm(openstate.elastic)==0 and np.linalg.norm(re['traction_N'])<1e-15)
(P/'receipts/fixed_stiffness_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
