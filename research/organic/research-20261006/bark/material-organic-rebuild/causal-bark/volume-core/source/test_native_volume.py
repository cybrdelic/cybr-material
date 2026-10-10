from pathlib import Path
import numpy as np,json,time,hashlib
from scipy.spatial.transform import Rotation
from ring_geometry import build
from batch_volume import BatchPrisms
from native_volume import NativePrisms
R=Path(__file__).resolve().parents[1];h,cells,q,meta,tri=build(16);a=BatchPrisms(cells);b=NativePrisms(cells);checks=[];rng=np.random.default_rng(761)
for name,x in [('grown',q),('rotated',q@Rotation.from_rotvec([1.,.2,-.5]).as_matrix().T+[.4,-.3,.2]),('reference',np.array([c.X for c in cells]))]:
 E,g,v=a.evaluate(x);En,gn,vn=b.evaluate(x);error=abs(En-E)/max(E,1e-22);gerr=np.linalg.norm(gn-g)/max(np.linalg.norm(g),1e-12);checks.append({'case':name,'energy_relative_or_floor_error':error,'gradient_relative_or_floor_error':gerr,'absolute_energy_difference_J':abs(En-E),'pass':(error<1e-9 or abs(En-E)<1e-23) and (gerr<1e-8 or np.linalg.norm(gn-g)<1e-11)})
start=time.time()
for _ in range(30):b.evaluate(q)
native=(time.time()-start)/30;start=time.time()
for _ in range(5):a.evaluate(q)
python=(time.time()-start)/5
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'cells':len(cells),'DOFs':q.size,'native_seconds_per_force':native,'numpy_seconds_per_force':python,'speedup':python/native,'source_sha256':hashlib.sha256((R/'source/batch_volume.cpp').read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256((R/'source/libvolume_native.so').read_bytes()).hexdigest(),'scope':'Same material law and geometry with numerical small-strain energy stabilization. Full ring assembly and evolution still pending.'};(R/'receipts/native_volume_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
