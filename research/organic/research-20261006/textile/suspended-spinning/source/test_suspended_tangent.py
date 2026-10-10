import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,json
from tangent_suspended import TangentSuspended,P
from contact_path import certify_path
s=TangentSuspended();z=np.load(P/'data/positive_equilibrium.npz');s.x=z['x'].copy();s.theta=float(z['theta_rad']);x=s.x.copy();E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);oldpoints=s.c.points().copy();r=s.record(x,True,[]);rows=[]
for dt in [.01,.005,.0025]:
 s.x=x.copy();s.theta=0.;s.evaluate(x);pred=s.predictor_boundary(dt);s.theta=dt;Et,gt=s.evaluate(pred);pt,*_=s.projection(gt);rows.append(dict(step_rad=dt,force_l2=float(np.linalg.norm(pt)),force_infinity=float(abs(pt).max()),root_error_m=float(abs(s.constraints()).max()*s.L),predictor=s.predictor_diagnostic,contact=certify_path(oldpoints,s.c.points(),s.c.radius,s.c.ref,plane_bounds=None)))
s.theta=0.;s.evaluate(x);eps=1e-6;values=[]
for sign in [1,-1]:
 s.theta=sign*eps;e,_=s.evaluate(x);values.append((e+lam@s.constraints())*s.unit)
s.theta=0.;s.evaluate(x);torque_fd=(values[0]-values[1])/(2*eps);torque_err=abs(torque_fd-r['torque_conjugate_Nm'])/max(abs(torque_fd),1e-30)
# Root constrained predictor errors should decrease quadratically until the
# starting equilibrium's finite force tolerance becomes the floor.
ratios=[rows[i]['force_l2']/rows[i+1]['force_l2'] for i in range(2)]
out=dict(steps=rows,force_error_reduction_ratios=ratios,torque_conjugate_Nm=r['torque_conjugate_Nm'],torque_fd_Nm=torque_fd,torque_derivative_relative_error=torque_err,passed=all(v>3 for v in ratios) and torque_err<1e-5 and all(a['contact']['passed'] for a in rows),boundary_planes_m=None)
(P/'receipts/suspended_tangent_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
