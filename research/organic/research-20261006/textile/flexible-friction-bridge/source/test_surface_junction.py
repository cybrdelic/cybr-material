"""Point-couple migration across a rod node: force-map continuity control."""
from pathlib import Path
import json,numpy as np
from rod_pullback import RodCoordinates
from rod_pullback_piecewise_control import RodCoordinates as Piecewise
from test_rod_coupling import fixture
P=Path(__file__).resolve().parents[1];results=[]
for cls in [Piecewise,RodCoordinates]:
 f=cls();g=f.unpack(fixture(f));rows=[]
 for eps in [1e-7,1e-8,1e-9,1e-10]:
  loads=[]
  for side in [-1,1]:
   gp=np.zeros_like(g['p']);gm=np.zeros((f.n,f.M,3));f.point_load(g,1,np.array([.0012+side*eps]),np.zeros((1,3)),np.array([[1e-8,2e-8,-.5e-8]]),gp,gm);loads.append(f.pullback(g,gp,gm))
  rows.append(dict(eps_m=eps,force_map_jump_relative=float(np.linalg.norm(loads[1]-loads[0])/max(np.linalg.norm(loads[0]),np.linalg.norm(loads[1])))))
 results.append(dict(implementation=cls.__module__,rows=rows))
r=dict(scope='Pure surface couple at a changing material coordinate; exact rigid stock and centreline are unchanged.',results=results,piecewise_failed=results[0]['rows'][-1]['force_map_jump_relative']>.1,continuous_passed=results[1]['rows'][-1]['force_map_jump_relative']<1e-5);(P/'receipts/surface_junction.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
