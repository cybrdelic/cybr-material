import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from incremental_envelope import IncrementalEnvelope
from directional_candidate import synthetic_candidate
P=Path(__file__).resolve().parents[1];old=synthetic_candidate();new=IncrementalEnvelope(old.tension,old.compression,old.G,old.axes,old.minimum_stretch,old.maximum_stretch);rng=np.random.default_rng(741);regular=[]
for k in range(15):
 H=rng.normal(size=(3,3))*.025;wa,pa,_=old.evaluate(np.eye(3)+H);wb,pb,_=new.evaluate_incremental(H);v=rng.normal(size=(3,3));v/=np.linalg.norm(v);h=1e-6;fd=(new.evaluate_incremental(H+h*v)[0]-new.evaluate_incremental(H-h*v)[0])/(2*h);regular.append(dict(relative_energy=float(abs(wb-wa)/max(abs(wa),1e-30)),relative_stress=float(np.linalg.norm(pb-pa)/max(np.linalg.norm(pa),1e-30)),directional_gradient_error_Pa=float(abs(fd-np.sum(pb*v)))))
small=[]
for eps in (1e-3,1e-6,1e-9,1e-12):
 H=np.diag([eps,0.,0.]);w,p,d=new.evaluate_incremental(H);expected=.5*old.tension[0].E*eps**2
 shear=np.zeros((3,3));shear[0,1]=eps;ws,ps,ds=new.evaluate_incremental(shear);angular_expected=.5*old.G*np.log1p(eps*eps)
 small.append(dict(eps=eps,axial_energy_relative_error=float(abs(w-expected)/expected),pure_extension_angular_energy_J_m3=float(d['angular_energy_J_m3']),shear_angular_energy_relative_error=float(abs(ds['angular_energy_J_m3']-angular_expected)/angular_expected),shear_stress_scaled=ps[0,1]/(old.G*eps)))
Q,_=np.linalg.qr(rng.normal(size=(3,3)));Q[:,0]*=np.linalg.det(Q);H=rng.normal(size=(3,3))*.02;w,p,_=new.evaluate_incremental(H);wr,pr,_=new.evaluate(Q@(np.eye(3)+H));obj=max(abs(w-wr)/1e6,np.max(abs(pr-Q@p))/1e6)
covariant=[]
for eps in (1e-6,1e-9,1e-12):
 H=np.diag([eps,-.5*eps,.2*eps]);H[0,1]=.3*eps;w,p,_=new.evaluate_incremental(H)
 rotated=IncrementalEnvelope(old.tension,old.compression,old.G,Q,old.minimum_stretch,old.maximum_stretch)
 wr,pr,_=rotated.evaluate_incremental(Q@H@Q.T)
 covariant.append(dict(eps=eps,relative_energy=float(abs(wr-w)/max(abs(w),1e-30)),relative_stress=float(np.linalg.norm(pr-Q@p@Q.T)/max(np.linalg.norm(p),1e-30))))
assert max(v['relative_energy'] for v in covariant)<1e-12 and max(v['relative_stress'] for v in covariant)<1e-12
r=dict(scope='Numerical evaluation correction only; same uncalibrated material form and parameter fixture',regular_parity=regular,small_strain=small,small_strain_frame_covariance=covariant,objectivity_error_over_modulus=float(obj),passed=bool(max(c['relative_energy'] for c in regular)<1e-10 and max(c['relative_stress'] for c in regular)<1e-10 and max(c['directional_gradient_error_Pa'] for c in regular)<.01 and max(c['axial_energy_relative_error'] for c in small)<1e-12 and max(c['shear_angular_energy_relative_error'] for c in small)<1e-12 and obj<1e-12))
(P/'receipts/incremental_evaluation.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));assert r['passed']
