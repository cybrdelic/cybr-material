import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from objective_pair import make,advance,exp,State,P
labels=((0,'arc',.001),(1,'arc',.001));d=32e-6;I=np.eye(3);pi=np.array([-d/2,0.,0.]);pj=-pi;N=1e-6;mu=.3;delta=1e-8;seed=make(labels,pi,pj,I,I,N,mu,delta,elastic=np.array([0.,.5*mu*delta,0.]));Q=exp(np.array([.8,1.2,-.7]));b=np.array([.001,-.002,.003]);rigid,rg=advance(seed,Q@pi+b,Q@pj+b,Q,Q,labels=labels,normal=N)
rigid_error=float(np.linalg.norm(rigid.elastic-Q@seed.elastic));rigid_du=float(np.linalg.norm(rg['surface_increment_m']));rolling,rr=advance(seed,pi,pj,exp([0,0,.4]),exp([0,0,-.4]),labels=labels,normal=N)
# Objectivity under independent superposed old/new rigid frames.
a=np.array([0.,3e-9,1e-9]);Ri=exp([.03,-.02,.07]);Rj=exp([-.01,.04,-.02]);base,rb=advance(seed,pi+a,pj,Ri,Rj,labels=labels,normal=N);A=exp([.7,-.8,.2]);B=exp([-.9,.3,1.1]);t0=np.array([.002,.001,-.002]);t1=np.array([-.001,.001,.001]);transformed=make(labels,A@pi+t0,A@pj+t0,A,A,N,mu,delta,elastic=A@seed.elastic);obj,ro=advance(transformed,B@(pi+a)+t1,B@pj+t1,B@Ri,B@Rj,labels=labels,normal=N);objectivity=float(np.linalg.norm(obj.elastic-B@base.elastic)/max(np.linalg.norm(base.elastic),1e-30));heat_error=abs(ro['friction_heat_J']-rb['friction_heat_J'])/max(abs(rb['friction_heat_J']),1e-30)
# A frozen first-order operator is a negative control for a finite rigid turn.
angle=.8;Z=exp([0,0,angle]);oldn=(pi-pj)/d;rel=(Z@pi-pi)-(Z@pj-pj)+np.cross([0,0,angle],-pi)-np.cross([0,0,angle],-pj);false_slip=rel-oldn*(oldn@rel)
cycles=[]
for count in [64,128,256,512]:
 state=make(labels,pi,pj,I,I,N,mu,delta);work=0.;maxbalance=0.;heat=0.;alg=0.
 for j in range(1,count+1):
  t=j/count;u=10*mu*delta*np.sin(2*np.pi*t);state,r=advance(state,pi+[0,u,0],pj,I,I,labels=labels,normal=N);work+=r['endpoint_discrete_work_J'];heat+=r['friction_heat_J'];alg+=r['algorithmic_energy_J'];maxbalance=max(maxbalance,abs(r['work_balance_residual_J']))
 cycles.append(dict(steps=count,endpoint_work_J=work,stored_J=r['stored_J'],heat_J=heat,algorithmic_energy_J=alg,algorithmic_fraction=alg/max(abs(work),1e-30),maximum_local_balance_error_J=maxbalance,global_balance_error_J=abs(work-r['stored_J']-heat-alg)))
rejections={}
for tag,kw in [('changed_label',dict(labels=((0,'arc',.0011),labels[1]),normal=N)),('changed_normal',dict(labels=labels,normal=1.01*N))]:
 try:advance(seed,pi,pj,I,I,**kw);rejections[tag]=False
 except (ValueError,NotImplementedError):rejections[tag]=True
out=dict(scope='Finite-pose constant-normal material-pair constitutive probe only. No rod solve, moving quadrature, normal-load coupling or spinning-friction qualification.',rigid_surface_increment_m=rigid_du,rigid_elastic_transport_error_m=rigid_error,rigid_heat_J=rg['friction_heat_J'],rolling_surface_increment_m=float(np.linalg.norm(rr['surface_increment_m'])),superposed_motion_relative_error=objectivity,superposed_heat_relative_error=heat_error,frozen_linear_operator_false_rigid_slip_m=float(np.linalg.norm(false_slip)),cycles=cycles,guards=rejections,passed=bool(rigid_du<1e-15 and rigid_error<1e-15 and rg['friction_heat_J']<1e-25 and np.linalg.norm(rr['surface_increment_m'])<1e-15 and objectivity<1e-8 and heat_error<1e-8 and all(rejections.values()) and all(c['global_balance_error_J']<1e-25 for c in cycles) and cycles[-1]['algorithmic_fraction']<.02 and all(cycles[i+1]['algorithmic_energy_J']<cycles[i]['algorithmic_energy_J'] for i in range(3))))
(P/'receipts/objective_pair_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
