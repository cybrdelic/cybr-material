import os,resource,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
from pathlib import Path
import json,numpy as np
from corrugated_wall import CorrugatedWall,Parameters
P=Path(__file__).resolve().parents[1];start=time.monotonic();runs=[]
for segments in (24,48,96):
 for direction,targets in (('compression',np.linspace(0,-.28,15)),('tension',np.linspace(0,.133,15))):
  w=CorrugatedWall(segments);x=w.X.copy();states=[];history=[]
  for strain in targets:
   x,r=w.solve(strain,x,160);history.append(r);states.append(x*w.L)
   print('WALL',segments,direction,strain,r['passed'],r['dimensionless_force_residual'],r['max_axial_strain'],time.monotonic()-start,flush=True)
   if not r['passed']:break
  np.savez_compressed(P/f'data/wall_{segments}_{direction}.npz',reference_m=w.X*w.L,positions_m=states,engineering_strains=[h['engineering_strain'] for h in history],parameters=list(vars(w.p).values()))
  summary=dict(segments=segments,direction=direction,all_passed=all(r['passed'] for r in history),history=history,elapsed_s=time.monotonic()-start)
  (P/f'receipts/wall_{segments}_{direction}.json').write_text(json.dumps(summary,indent=2));runs.append(summary)
  if time.monotonic()-start>120:raise TimeoutError('Wall mechanism qualification exceeded declared120s budget')
(P/'receipts/wall_controls.json').write_text(json.dumps(dict(scope='Resolved elastic corrugated wall mechanism; not whole cork constitutive or fracture qualification',runs=runs,wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024),indent=2))
