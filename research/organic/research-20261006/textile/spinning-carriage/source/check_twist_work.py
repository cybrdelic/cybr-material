import json,numpy as np
from carriage import Carriage,P
s=Carriage();z=np.load(P/'data/first_twist_step.npz');s.x=z['x'];s.theta=float(z['theta_rad']);E,g=s.evaluate(s.x);pg,lam,J,S=s.projection(g);r=s.record(s.x,True,[]);theta=s.theta;eps=1e-6;values=[]
for d in [eps,-eps]:
 s.theta=theta+d;e,_=s.evaluate(s.x);values.append((e+lam@s.constraints())*s.unit)
s.theta=theta;s.evaluate(s.x);fd=(values[0]-values[1])/(2*eps);analytic=r['torque_conjugate_Nm'];rel=abs(fd-analytic)/max(abs(fd),abs(analytic));out=dict(finite_Lagrangian_control_derivative_Nm=fd,analytic_root_force_moment_work_Nm=analytic,relative_error=rel,passed=bool(rel<1e-5),scope='Full root translation and frame rotation at fixed carriage coordinate and bearing height; the axial force-control term is independent of theta.');(P/'receipts/twist_work_derivative.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
