import sys,math,json
from pathlib import Path
from mathutils import Vector
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'source'),str(ROOT/'inputs/CYBR_structures_portable_review/source')]
from expansion.geometry import Tubes as LegacyTubes
from transport_tubes import TransportTubes
points=[(t-.5,0,2*math.sin(math.pi*t)) for t in [i/64 for i in range(65)]]
tangents=[]
for i,p in enumerate(points):
 tangents.append((Vector(points[min(i+1,len(points)-1)])-Vector(points[max(i-1,0)])).normalized())
tangent_dots=[a.dot(b) for a,b in zip(tangents,tangents[1:])]
rows=[]
for typ in (LegacyTubes,TransportTubes):
 obj=typ();obj.path(points,.03,6)
 frames=[(Vector(obj.v[i*6])-Vector(p)).normalized() for i,p in enumerate(points)]
 dots=[a.dot(b) for a,b in zip(frames,frames[1:])]
 maxraderr=max(abs((Vector(obj.v[i*6+k])-Vector(p)).length-.03) for i,p in enumerate(points) for k in range(6))
 row={'builder':typ.__name__,'minimum_adjacent_frame_dot':min(dots),'abrupt_half_turns':sum(d<-.8 for d in dots),'max_radius_error':maxraderr};rows.append(row)
 print(row,flush=True)
 if typ is TransportTubes:
  assert min(d-td for d,td in zip(dots,tangent_dots))>-1e-5, 'Frame turns faster than the planar centerline tangent'
  assert maxraderr<1e-6
assert rows[0]['abrupt_half_turns']>=1
(ROOT/'receipts/tube_frame_regression.json').write_text(json.dumps({'pass':True,'test':'Vertical loop tangent crosses legacy0.92 threshold twice','results':rows},indent=2));print(json.dumps(rows),flush=True)
