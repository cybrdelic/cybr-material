from pathlib import Path
import numpy as np,json,time
from ring_system import RingSystem
R=Path(__file__).resolve().parents[1];m=RingSystem(8,.001);checks=[];rng=np.random.default_rng(156);v=rng.normal(size=m.q.shape);v[~m.free_mask]=0.;v*=1e-4
# The reference linear operator uses each material cell's stress-free state.
# Check bulk block Hessians separately; the cohesive tangent is already linear
# for intact isotropic kt=kn bonds, regardless of current surface normal.
x=np.array([c.X for c in m.cells]);eps=1e-5;Gp=m.bulk.evaluate(x+eps*v)[1];Gm=m.bulk.evaluate(x-eps*v)[1];fd=(Gp-Gm)/(2*eps);exact=np.einsum('cij,cj->ci',m.K,v.reshape(m.N,54)).reshape(v.shape);err=np.linalg.norm(fd-exact)/np.linalg.norm(exact);checks.append({'case':'bulk_reference_tangent','relative_error':err,'pass':err<1e-6})
Gp=m.bonds.evaluate(m.q+eps*v)[1];Gm=m.bonds.evaluate(m.q-eps*v)[1];fd=(Gp-Gm)/(2*eps);actual=m.linear_operator(v)-exact;err=np.linalg.norm(fd-actual)/np.linalg.norm(actual);checks.append({'case':'intact_cohesive_tangent','relative_error':err,'pass':err<1e-6})
Minvv=m.acceleration(-np.einsum('cij,cja->cia',m.M,v));err=np.linalg.norm(Minvv-v)/np.linalg.norm(v);checks.append({'case':'free_inertia_inverse','relative_error':err,'pass':err<1e-12});omega=m.maximum_reference_frequency();E,g,op,report=m.evaluate();checks.append({'case':'closed_ring_initial_state','minimum_det_F':report['minimum_det_F'],'cohesive_opening_m':float(op.max()),'reference_frequency_per_s':omega,'pass':report['minimum_det_F']>0 and op.max()<1e-14 and np.isfinite(omega)});start=time.time()
for _ in range(30):m.evaluate()
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'cells':m.N,'free_DOFs':len(m.free),'seconds_per_complete_force':(time.time()-start)/30,'scope':'Complete radial-load ring operator and consistent inertia; evolution and fracture validation pending.'};(R/'receipts/ring_system_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
