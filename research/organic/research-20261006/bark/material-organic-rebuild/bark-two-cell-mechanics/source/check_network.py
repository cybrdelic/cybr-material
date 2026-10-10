import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json, hashlib
import numpy as np
from two_cell_network import TwoCellNetwork

P=Path(__file__).resolve().parents[1]
n=TwoCellNetwork(24,8)
rng=np.random.default_rng(80841)
q=n.q0+rng.normal(size=n.ndof)*1e-4
e,g,H,parts=n.evaluate(q,True)
v=rng.normal(size=n.ndof);v/=np.linalg.norm(v)
checks=[]
for h in (1e-5,1e-6,1e-7):
    ep,gp,_=n.evaluate(q+h*v);em,gm,_=n.evaluate(q-h*v)
    checks.append(dict(step=h,gradient_error=abs((ep-em)/(2*h)-g@v)/max(abs(g@v),1.),
                       Hessian_error=float(np.linalg.norm((gp-gm)/(2*h)-H@v)/np.linalg.norm(H@v))))
a=.83;R=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
qr=q.copy();qr[:2*n.nn]=(q[:2*n.nn].reshape(-1,2)@R.T+[.31,-.27]).ravel();qr[2*n.nn:]+=a
er,gr,Hr,_=n.evaluate(qr,True)
T=np.eye(n.ndof)
for i in range(n.nn): T[2*i:2*i+2,2*i:2*i+2]=R
x=q[:2*n.nn].reshape(-1,2);gx=g[:2*n.nn].reshape(-1,2)
f_identity=np.sum(gx,axis=0)
m_identity=float(np.sum(x[:,0]*gx[:,1]-x[:,1]*gx[:,0])+np.sum(g[2*n.nn:]))
e0,g0,H0,_=n.evaluate(n.q0,True)
closure=[]
for boundary in n.cells:
    for (wi,d),(wj,dj) in zip(boundary,boundary[1:]+boundary[:1]):
        i=n.walls[wi]['nodes'][::d][-1];j=n.walls[wj]['nodes'][::dj][0]
        closure.append(bool(i==j))
# The single shared-wall node array is used in opposite directions by both cells.
shared=[d for b in n.cells for w,d in b if w==1]
receipt=dict(source_snapshot_sha256=hashlib.sha256((P/'source/corrugated_wall_frozen.py').read_bytes()).hexdigest(),
    reference_energy_J=e0*n.unit,reference_gradient_max=float(np.max(abs(g0))),
    derivative_checks=checks,
    energy_objectivity_error=abs(er-e)/max(abs(e),1.),
    gradient_objectivity_error=float(np.linalg.norm(gr-T@g)/np.linalg.norm(g)),
    Hessian_objectivity_error=float(np.linalg.norm(Hr-T@H@T.T)/np.linalg.norm(H)),
    Hessian_symmetry_error=float(np.linalg.norm(H-H.T)/np.linalg.norm(H)),
    translation_identity_error=float(np.linalg.norm(f_identity)),
    rotation_identity_error=abs(m_identity),
    all_cell_boundary_ids_closed=all(closure),shared_wall_used_once_opposite_orientation=shared==[1,-1],
    junction_end_terms=14,independent_junction_rotations=6,
    reference_free_minimum_eigenvalue=float(np.linalg.eigvalsh(H0[np.ix_(n.free,n.free)])[0]))
receipt['passed']=bool(e0*n.unit<1e-30 and np.max(abs(g0))<1e-9
    and min(c['gradient_error'] for c in checks)<1e-7
    and min(c['Hessian_error'] for c in checks)<1e-7
    and receipt['energy_objectivity_error']<1e-9
    and receipt['gradient_objectivity_error']<1e-9
    and receipt['Hessian_objectivity_error']<1e-9
    and receipt['Hessian_symmetry_error']<1e-12
    and receipt['translation_identity_error']<1e-9
    and receipt['rotation_identity_error']<1e-9 and all(closure) and shared==[1,-1]
    and receipt['reference_free_minimum_eigenvalue']>0)
(P/'receipts/primitive_checks.json').write_text(json.dumps(receipt,indent=2))
(P/'data/topology_24.json').write_text(json.dumps(n.topology(),indent=2))
print(json.dumps(receipt,indent=2));assert receipt['passed']
