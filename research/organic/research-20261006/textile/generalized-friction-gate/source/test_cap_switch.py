"""One-sided derivatives at a committed Coulomb cap; no ordinary Jacobian there.
Use actual full symmetric residual. Analytic infinite/coextensive-span limit:
forward dF/dx=0; reverse dF/dx=KT*L, allowing 0.1% finite-end correction.
"""
from pathlib import Path
import json,numpy as np
from symmetric_spans import SymmetricSpans
from coupled_spans import S,RS,KT,MU,L
P=Path(__file__).resolve().parents[1];f=SymmetricSpans();q=np.zeros(12);q[2]=(RS-.45e-6)/S;q[5]=-np.pi/2;old=f.state(q)
# Both histories refer to the same physical positive tangential traction on A.
for name,sign in [('a',1),('b',-1)]:
 old[name]['elastic'][:,0]=sign*MU*old[name]['geometry']['N']/f.k
base=f.trial(q,old);rows=[]
for h in [1e-11,1e-12,1e-13,1e-14]:
 qp=q.copy();qm=q.copy();qp[0]+=h/S;qm[0]-=h/S;a=f.trial(qp,old);b=f.trial(qm,old);forward=(a['friction_wrenches'][0]-base['friction_wrenches'][0])/h;backward=(base['friction_wrenches'][0]-b['friction_wrenches'][0])/h;rows.append(dict(h_m=h,forward_N_m=forward,backward_N_m=backward,forward_relative_stiffness=abs(forward)/(KT*L),backward_relative_error=abs(backward/(KT*L)-1),forward_sliding=int(a['sliding_set'].sum()),backward_sliding=int(b['sliding_set'].sum())))
r=dict(scope='Actual symmetric coupled residual at a cap. One-sided limits differ; bracketed driven solver handles switches without a symmetric Hessian.',analytic_reverse_N_m=KT*L,finite_end_relative_tolerance=.001,rows=rows,passed=all(z['forward_relative_stiffness']<.001 and z['backward_relative_error']<.001 for z in rows))
(P/'receipts/cap_switch_tests.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
