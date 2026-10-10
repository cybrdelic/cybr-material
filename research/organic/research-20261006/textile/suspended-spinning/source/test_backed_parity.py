import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import subprocess,sys,json,numpy as np
from pathlib import Path
P=Path(__file__).resolve().parents[1]
code='''import sys,numpy as np
sys.path.insert(0,sys.argv[1])
from carriage import Carriage,RodHessianOperator
s=Carriage(plane_bounds=(0.,.02)) if sys.argv[3]=='new' else Carriage()
x=s.x.copy();E,g=s.evaluate(x);pg,lam,J,S=s.projection(g);op=RodHessianOperator(s.c,lam);rng=np.random.default_rng(721);V=rng.normal(size=(3,len(x)-1));V/=np.linalg.norm(V,axis=1)[:,None];H=np.array([op.matvec(v) for v in V]);np.savez(sys.argv[2],energy=E,gradient=g,projected_gradient=pg,roots=s.constraints(),points=s.c.points(),Hessian_products=H)
'''
for label,source in [('old',P.parent/'spinning-carriage/source'),('new',P/'source')]:
 subprocess.run([sys.executable,'-c',code,str(source),str(P/'data'/f'backed_parity_{label}.npz'),label],check=True,timeout=30)
a=np.load(P/'data/backed_parity_old.npz');b=np.load(P/'data/backed_parity_new.npz');rows={k:dict(maximum_absolute_error=float(np.max(np.abs(a[k]-b[k]))),relative_error=float(np.linalg.norm(a[k]-b[k])/max(np.linalg.norm(a[k]),1e-30))) for k in a.files};out=dict(tests=rows,passed=all(r['relative_error']<1e-12 for r in rows.values()),scope='Old backed model and explicitly backed suspended-stage code agree at identical material state and inputs.')
(P/'receipts/backed_parity.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
