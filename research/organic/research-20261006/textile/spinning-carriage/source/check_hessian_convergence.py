import json,time,numpy as np
from carriage import Carriage,RodHessianOperator,P
s=Carriage();x=s.x.copy();rng=np.random.default_rng(478);v=rng.normal(size=len(x));v/=np.linalg.norm(v);rows=[]
for order in [16,24,32]:
 s.c.contact_order=order;E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);diff_force=float(np.linalg.norm(op.contact.gradient+s.c.last['backing_gradient_N']-s.c.last['position_gradient_N']));expected=np.r_[op.matvec(v[:-1]),0.]
 for eps in [1e-4,1e-5,1e-6,1e-7]:
  _,gp=s.evaluate(x+eps*v);jp=s.jacobian();_,gm=s.evaluate(x-eps*v);jm=s.jacobian();fd=(gp+jp.T@lam-gm-jm.T@lam)/(2*eps);err=float(np.linalg.norm(fd-expected)/np.linalg.norm(expected));rows.append(dict(order=order,epsilon=eps,Hessian_relative_error=err,position_force_disagreement_N=diff_force,maximum_force_residual=float(abs(pg).max())));print(rows[-1],flush=True)
 s.evaluate(x)
(P/'receipts/Hessian_convergence.json').write_text(json.dumps(rows,indent=2))
