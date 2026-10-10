import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from corrugated_wall import CorrugatedWall
P=Path(__file__).resolve().parents[1];w=CorrugatedWall(24);rng=np.random.default_rng(774);x=w.X+rng.normal(size=w.X.shape)*1e-4;e,g,H,p=w.evaluate(x,hessian=True);v=rng.normal(size=x.shape);v/=np.linalg.norm(v);rows=[]
for eps in (1e-5,1e-6,1e-7):
 ep,gp,_=w.evaluate(x+eps*v);em,gm,_=w.evaluate(x-eps*v);rows.append(dict(eps=eps,energy_gradient_error=abs((ep-em)/(2*eps)-np.sum(g*v))/max(abs(np.sum(g*v)),1.),Hessian_error=np.linalg.norm(((gp-gm)/(2*eps)).ravel()-H@v.ravel())/np.linalg.norm(H@v.ravel())))
a=.77;R=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]]);er,gr,hr,pr=w.evaluate(x@R.T+[.3,-.7],rotation_rad=a,hessian=True)
out=dict(derivatives=rows,objectivity_energy_error=abs(er-e)/max(abs(e),1.),objectivity_force_error=np.linalg.norm(gr-g@R.T)/np.linalg.norm(g),stress_free_energy_J=w.evaluate(w.X)[0]*w.unit,passed=bool(min(r['energy_gradient_error'] for r in rows)<1e-7 and min(r['Hessian_error'] for r in rows)<1e-7 and abs(er-e)<1e-9))
(P/'receipts/primitive_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
