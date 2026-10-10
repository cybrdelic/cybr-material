from pathlib import Path
import json,numpy as np
from rod_pullback import RodCoordinates,log,mv
from test_rod_coupling import fixture
P=Path(__file__).resolve().parents[1];f=RodCoordinates();q=fixture(f);g=f.unpack(q);s=np.array([0.,1e-8,.000375,.0012-1e-10,.0012,.0012+1e-10,.00239999,.0024]);D=f.surface(g,1,s)[0];rng=np.random.default_rng(553);rows=[]
for _ in range(3):
 v=rng.normal(size=len(q));v/=np.linalg.norm(v);_,omega=f.velocity(g,v,1,s);steps=[]
 for h in [1e-4,1e-5,1e-6]:
  plus=f.surface(f.unpack(q+h*v),1,s)[0];minus=f.surface(f.unpack(q-h*v),1,s)[0];numeric=mv(D,log(np.swapaxes(minus,-1,-2)@plus)/(2*h));steps.append(dict(h=h,relative_error=float(np.linalg.norm(numeric-omega)/np.linalg.norm(omega))))
 rows.append(steps)
r=dict(scope='Independent finite-rotation check of the analytic geodesic surface-spin Jacobian, including centreline nodes and endpoints.',rows=rows,passed=all(min(s['relative_error'] for s in row)<1e-7 for row in rows));(P/'receipts/surface_spin.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
