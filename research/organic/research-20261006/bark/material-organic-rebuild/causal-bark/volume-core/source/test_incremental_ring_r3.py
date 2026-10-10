from pathlib import Path
import json,hashlib,numpy as np
from scipy.spatial.transform import Rotation
from incremental_ring_r3 import IncrementalRing
R=Path(__file__).resolve().parents[1];m=IncrementalRing(4,0.);rng=np.random.default_rng(9261);u=rng.normal(0,4e-7,m.q.shape);E,g,op,rep=m.evaluate_displacement(u);Ea,ga,opa,repa=m.evaluate(m.origin+u);comparison={'energy_relative':abs(E-Ea)/max(abs(E),1e-20),'gradient_relative':float(np.linalg.norm(g-ga)/np.linalg.norm(ga)),'opening_absolute_m':float(np.max(abs(op-opa)))}
assert comparison['energy_relative']<1e-8 and comparison['gradient_relative']<1e-7
# Difference quotient at resolved field scale, complete nonlinear operator.
v=rng.normal(size=u.shape);v/=np.linalg.norm(v);h=2e-10
fd=(m.evaluate_displacement(u+h*v)[0]-m.evaluate_displacement(u-h*v)[0])/(2*h);gd=float(np.sum(v*g));err=abs(fd-gd)/max(abs(fd),abs(gd),1e-20);assert err<1e-6
# Removing spatial birth prestrain must change stored energy: this origin is
# not a stress-free metric. This catches the tempting F=I+grad(du) bug.
Ebirth,gb,_,_=m.evaluate_displacement(np.zeros_like(u));H0=m.H0.copy();m.H0[:]=0.;Ewrong,_,_,_=m.evaluate_displacement(np.zeros_like(u));m.H0[:]=H0;assert abs(Ebirth-Ewrong)>1e-8
# A rigid rotation moves both prestrain and current surface normals. Remove
# gravity before assessing objectivity, because the laboratory field rotates
# neither with the body nor its potential datum.
Q=Rotation.from_rotvec([.4,-.3,.7]).as_matrix();qr=(m.origin+u)@Q.T;ur=qr-m.origin;Er,gr,_,_=m.evaluate_displacement(ur);Ein=E+float(np.sum(m.gravity*(m.origin+u)));Erin=Er+float(np.sum(m.gravity*qr));gin=g+m.gravity;grin=gr+m.gravity;obj=max(abs(Ein-Erin)/abs(Ein),float(np.linalg.norm(grin-gin@Q.T)/np.linalg.norm(gin)));assert obj<1e-7
# Finite-difference Hessian at much smaller increments than absolute q can
# reliably retain. The derivative should converge then round off, not disappear.
from ring_system_r2 import RingSystem
hsmall=2e-13;gp=m.evaluate_displacement(u+hsmall*v)[1];gm=m.evaluate_displacement(u-hsmall*v)[1];hvp=(gp-gm)/(2*hsmall);assert np.linalg.norm(hvp)>1.
out={'all_pass':True,'absolute_coordinate_comparison':comparison,'directional_energy_gradient_relative':err,'complete_internal_objectivity_relative':obj,'retained_birth_energy_J':Ebirth,'incorrect_reset_metric_energy_J':Ewrong,'initial_interface_roundoff_m':m.reference_interface_gap_roundoff_m,'small_increment_m':hsmall,'resolved_Hv_norm_N_per_m':float(np.linalg.norm(hvp)),'sources':{name:hashlib.sha256((R/'source'/name).read_bytes()).hexdigest() for name in ('incremental_ring_r3.py','batch_volume_incremental.cpp','cohesion_incremental.cpp')},'scope':'Equivalent-law incremental coordinates and gradient/objectivity checks. No loading-path or morphology qualification.'}
(R/'receipts/incremental_ring_r3_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
