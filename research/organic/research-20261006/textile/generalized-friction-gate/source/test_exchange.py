import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
import json,time,numpy as np
from coupled_spans import Spans,exp,log,S,RS,P
start=time.monotonic();C=exp([0,0,-np.pi/2]);D=C.T

def swap(q):
 z=np.r_[q[6:9],log(exp(q[9:12])@C),q[:3],log(exp(q[3:6])@D)];return z
q0=np.zeros(12);q0[2]=(RS-.45e-6)/S;q0[6]=.4e-6/S;qs=[]
for dx,dz in [(10e-9,0),(-10e-9,0),(-10e-9,-.05e-6)]:
 q=q0.copy();q[0]+=dx/S;q[2]+=dz/S;qs.append(q)
rows=[];saved={}
for order in [65,129,257]:
 f=Spans(primary_order=order);a=f.state(q0);b=f.state(swap(q0));Ha=Hb=Wa=Wb=0.;steps=[]
 for index,q in enumerate(qs):
  ra=f.trial(q,a);rb=f.trial(swap(q),b);back=np.r_[rb['friction_wrenches'][6:],rb['friction_wrenches'][:6]];scale=max(np.linalg.norm(ra['friction_wrenches']),np.linalg.norm(back),1e-30);difference=float(np.linalg.norm(ra['friction_wrenches']-back)/scale);steps.append(dict(step=index,friction_force_moment_exchange_difference=difference,normal_energy_exchange_difference=abs(ra['normal_energy_J']-rb['normal_energy_J'])/ra['normal_energy_J']));Ha+=ra['friction_heat_J'];Hb+=rb['friction_heat_J'];Wa+=ra['tangential_endpoint_work_J'];Wb+=rb['tangential_endpoint_work_J'];a=ra['state'];b=rb['state']
 rows.append(dict(order=order,steps=steps,heat_original_J=Ha,heat_exchanged_J=Hb,heat_relative_difference=abs(Ha-Hb)/max(abs(Ha),abs(Hb)),work_original_J=Wa,work_exchanged_J=Wb,work_relative_difference=abs(Wa-Wb)/max(abs(Wa),abs(Wb))));saved[order]=[ra['friction_wrenches'],back]
base=saved[257];prior=saved[129];error=sum(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-30) for x,y in zip(base,prior));exchange=rows[-1]['steps'][-1]['friction_force_moment_exchange_difference'];asym=bool(exchange>max(3*error,1e-5));out=dict(one_sided_exchange_gate_passed=not asym,resolved_constitutive_asymmetry=asym,rows=rows,summed_final_quadrature_change=float(error),resolved_exchange_difference=exchange,wall_s=time.monotonic()-start,scope='Zero initial history followed by nonzero displacement/reversal/compression histories. One-sided primary reduction; not a bundle qualification.');(P/'receipts/exchange_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
