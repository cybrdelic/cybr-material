from pathlib import Path
import json,time,numpy as np
from prepare_three import ThreeCarriage
from carriage import RodHessianOperator,stability
P=Path(__file__).resolve().parents[1];start=time.monotonic();s=ThreeCarriage();z=np.load(P/'data/three_preparation_terminal.npz');x=z['x'];E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);stable=stability(op,s,J,S);bridge=s.bridge;geom=bridge.geometry(bridge.pack(s.c.origins,__import__('rod_pullback').exp(s.c.w)));ne,ng,_=bridge.native(geom['q']);r=dict(force_residual=float(abs(pg).max()),native_energy_relative=float(abs(ne/geom['normal_energy_J']-1)),native_contact_force_relative=float(np.linalg.norm(ng-geom['normal_gp'])/np.linalg.norm(ng)),stability_native=stable)
rng=np.random.default_rng(257);checks=[]
for _ in range(3):
 v=rng.normal(size=x.shape);v-=J.T@np.linalg.solve(S,J@v);v/=np.linalg.norm(v);expected=np.r_[op.matvec(v[:-1]),0.];samples=[]
 for h in [1e-5,1e-6,1e-7]:
  _,gp=s.evaluate(x+h*v);rp=gp+s.jacobian().T@lam;_,gm=s.evaluate(x-h*v);rm=gm+s.jacobian().T@lam;numeric=(rp-rm)/(2*h);samples.append(dict(h=h,native_tangent_relative=float(np.linalg.norm(numeric-expected)/np.linalg.norm(numeric))))
 checks.append(samples)
s.evaluate(x);r.update(tangent_checks=checks,elapsed_s=time.monotonic()-start);(P/'receipts/preparation_diagnosis.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
