import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from zero_shear_network import ZeroAverageShearNetwork
P=Path(__file__).resolve().parents[1]
n=ZeroAverageShearNetwork(24,8);rng=np.random.default_rng(91)
dq=rng.normal(size=len(n.free))*1e-4
q=n.increment(n.q0,dq)
e,g,H,_=n.evaluate(q,True);gf=n.reduced_gradient(g);hf=n.reduced_hessian(H)
F=np.zeros((n.ndof,len(n.free)));F[n.free,np.arange(len(n.free))]=1;F[n.dependent]=n.coeff
v=rng.normal(size=len(n.free));v/=np.linalg.norm(v)
checks=[]
for h in [1e-5,1e-6,1e-7]:
    ep,gp,_=n.evaluate(n.increment(q,v,h));em,gm,_=n.evaluate(n.increment(q,v,-h))
    checks.append(dict(step=h,gradient_error=abs((ep-em)/(2*h)-gf@v)/max(abs(gf@v),1.),
                       Hessian_error=float(np.linalg.norm(n.reduced_gradient(gp-gm)/(2*h)-hf@v)/np.linalg.norm(hf@v))))
state,receipt=n.solve(-.005)
_,gs,_=n.evaluate(state);support=n.support_gradient(gs)
res=gs-support
out=dict(constraint_error=float(abs(n.shear_gradient@q)),
    elimination_gradient_error=float(np.max(abs(gf-F.T@g))),
    elimination_Hessian_error=float(np.max(abs(hf-F.T@H@F))),
    reduced_derivative_checks=checks,small_increment=receipt,
    KKT_position_residual_N=float(np.max(abs(res[:2*n.nn]))*n.unit/n.L),
    linear_constraint_nullspace_error=float(np.max(abs(n.shear_gradient@F))),
    independent_junction_x_motion_dofs=4,
    passed=bool(receipt['passed'] and abs(n.shear_gradient@q)<1e-12
        and np.max(abs(gf-F.T@g))<1e-9 and np.max(abs(hf-F.T@H@F))<1e-8
        and min(c['gradient_error'] for c in checks)<1e-7
        and min(c['Hessian_error'] for c in checks)<1e-7
        and np.linalg.norm(receipt['support_force_balance_N'])<1e-13
        and abs(receipt['support_moment_balance_Nm'])<1e-18))
(P/'receipts/zero_shear_primitive_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='small_increment'},indent=2));assert out['passed']
