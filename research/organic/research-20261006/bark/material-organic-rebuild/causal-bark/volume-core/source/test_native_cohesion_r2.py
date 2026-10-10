from pathlib import Path
import json,hashlib,numpy as np
from ring_geometry import build
from ring_interfaces import build_interfaces
from native_cohesion import NativeCohesion as Frozen
from native_cohesion_r2 import NativeCohesion
R=Path(__file__).resolve().parents[1]
h,cells,q,meta,tris=build(4);ops=build_interfaces(cells,q,meta,4);N=len(ops['area']);checks=[]
def rejects(label,fn):
 try:fn();ok=False
 except (ValueError,TypeError,OverflowError):ok=True
 checks.append({'case':label,'pass':ok});assert ok,label
old=Frozen(len(cells),dict(ops,strength=np.full(N,np.nan)));old_accepts_nan=bool(np.isnan(old.delta0).all())
checks.append({'case':'Frozen wrapper reproduces unchecked NaN strength','pass':old_accepts_nan});assert old_accepts_nan
for key in ('strength','Gc'):
 for value in (np.full(N,np.nan),np.full(N,np.inf),np.zeros(N),np.full(N,-1.),np.ones(N+1),1.):
  rejects(key+' invalid '+str(np.shape(value)),lambda key=key,value=value:NativeCohesion(len(cells),dict(ops,**{key:value})))
rejects('fractional cell count',lambda:NativeCohesion(len(cells)+.5,ops))
p=ops['pairs'].astype(float);p[0,0]+=.5
rejects('fractional material indices',lambda:NativeCohesion(len(cells),dict(ops,pairs=p)))
m=NativeCohesion(len(cells),ops);base=Frozen(len(cells),ops)
for value in (np.full(N,np.nan),np.full(N,-1.),np.zeros(N+1),0.):rejects('invalid evolving opening',lambda value=value:m.evolve(value))
rejects('inconsistent restored damage',lambda:m.restore(np.zeros(N),np.ones(N)*.2))
assert not m.damage.any() and not m.maximum_opening.any()
# Owned material arrays must survive mutations to the caller's input.
a0=m.ops['area'].copy();ops['area']*=2
checks.append({'case':'owned operator arrays do not alias caller','pass':bool(np.array_equal(a0,m.ops['area']))})
# Actual law: work minus recoverable energy equals Gc at full failure;
# unloading and closure retain the maximum-opening state.
for ratio in (0.,.2,.8,1.,1.5):
 k=m.deltaf*ratio;m.evolve(k);base.evolve(k)
 assert np.array_equal(m.damage,base.damage) and np.allclose(m.dissipation(),base.dissipation())
 e,g,o=m.evaluate(q);eb,gb,ob=base.evaluate(q)
 # Base material area aliases the caller in this reproducible old behavior;
 # compare against the matching original area explicitly.
 base.ops['area']=a0.copy();eb,gb,ob=base.evaluate(q)
 assert np.allclose(e,eb,rtol=1e-13,atol=1e-22) and np.allclose(g,gb,rtol=1e-13,atol=1e-22)
D=m.dissipation();target=float(a0@m.Gc);assert np.isclose(D,target,rtol=1e-14)
before=(m.damage.copy(),m.maximum_opening.copy());m.evolve(np.zeros(N));assert np.array_equal(before[0],m.damage) and np.array_equal(before[1],m.maximum_opening)
checks.extend([{'case':'same valid traction/energy law','pass':True},{'case':'full failure dissipates reference-area times Gc','pass':True},{'case':'unloading irreversibility','pass':True}])
out={'all_pass':all(x['pass'] for x in checks),'cases':checks,'failure_dissipation_J':D,'expected_Gc_area_J':target,'source_sha256':hashlib.sha256(Path(__file__).with_name('native_cohesion_r2.py').read_bytes()).hexdigest(),'scope':'Input/history integrity and law-preserving wrapper. Not physical calibration or native fracture evolution qualification.'}
(R/'receipts/native_cohesion_r2_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
