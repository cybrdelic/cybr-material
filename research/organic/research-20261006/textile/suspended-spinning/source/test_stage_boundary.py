"""Check a manufacturing boundary change across every mechanical operator."""
import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
import json,numpy as np
from carriage import Carriage,RodHessianOperator,P
from contact_path import certify_path
s=Carriage();x=s.x.copy();E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam)
rng=np.random.default_rng(854);tests=[]
for i in range(3):
 v=rng.normal(size=x.size);v/=np.linalg.norm(v);rows=[]
 for eps in (1e-6,1e-7,1e-8):
  ep,gp=s.evaluate(x+eps*v);cp=s.constraints().copy();Jp=s.jacobian()
  em,gm=s.evaluate(x-eps*v);cm=s.constraints().copy();Jm=s.jacobian()
  fd=(gp+Jp.T@lam-gm-Jm.T@lam)/(2*eps);hv=np.r_[op.matvec(v[:-1]),0.]
  rows.append(dict(eps=eps,energy_gradient_error=abs((ep-em)/(2*eps)-g@v)/max(abs(g@v),1e-15),constraint_error=np.linalg.norm((cp-cm)/(2*eps)-J@v)/max(np.linalg.norm(J@v),1e-15),Hessian_error=np.linalg.norm(fd-hv)/max(np.linalg.norm(hv),1e-15)))
 tests.append(rows)
s.evaluate(x);p0=s.c.points().copy();oldheight=s.c.fixed_axis_height;s.c.fixed_axis_height-=.005;e1,g1=s.evaluate(x);p1=s.c.points().copy();translated_op=RodHessianOperator(s.c,lam)
v=rng.normal(size=x.size-1);v/=np.linalg.norm(v)
translation=dict(energy_relative_error=abs(e1-E)/max(abs(E),1e-30),gradient_relative_error=np.linalg.norm(g1-g)/max(np.linalg.norm(g),1e-30),Hessian_relative_error=np.linalg.norm(translated_op.matvec(v)-op.matvec(v))/max(np.linalg.norm(op.matvec(v)),1e-30),zero_plane_Hessian=translated_op.surface.nnz==0,zero_plane_gradient=np.max(abs(s.c.last['backing_gradient_N']))==0,zero_plane_energy=s.c.last['backing_contact_energy_J']==0)
path=certify_path(p0,p1,s.c.radius,s.c.ref,plane_bounds=None);backed_path=certify_path(p0,p1,s.c.radius,s.c.ref,plane_bounds=(0.,.02));s.c.fixed_axis_height=oldheight;s.evaluate(x)
backed=Carriage(plane_bounds=(0.,.02));Eb,gb=backed.evaluate(x);backed_op=RodHessianOperator(backed.c,lam)
out=dict(scope='Suspended spinning boundary source check. Not equilibrium, process calibration, or wool realism.',derivative_checks=tests,rigid_translation=translation,path_suspended=path,path_backed_control=backed_path,inherited_state_backed_energy=Eb,inherited_state_suspended_energy=E,initial_suspended_force_residual=float(abs(pg).max()),backed_plane_reaction_N=backed.c.last['planes'],backed_plane_Hessian_nonzeros=backed_op.surface.nnz,passed=all(min(r['energy_gradient_error'] for r in t)<1e-6 and min(r['constraint_error'] for r in t)<1e-6 and min(r['Hessian_error'] for r in t)<1e-4 for t in tests) and all(translation[k]<1e-7 for k in ('energy_relative_error','gradient_relative_error','Hessian_relative_error')) and translation['zero_plane_Hessian'] and translation['zero_plane_gradient'] and translation['zero_plane_energy'] and path['passed'] and not backed_path['passed'])
(P/'receipts/stage_boundary_tests.json').write_text(json.dumps(out,indent=2,default=lambda a:a.item()));print(json.dumps(out,indent=2,default=lambda a:a.item()));assert out['passed']
