"""A partner-fibre ID alone is insufficient for disjoint contact components."""
from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parents[1]
def components(point,vertices,material_arc,radius):
 intervals=[]
 for a,b,sa,sb in zip(vertices[:-1],vertices[1:],material_arc[:-1],material_arc[1:]):
  v=b-a;d=a-point;A=v@v;B=d@v;C=d@d-radius**2;disc=B*B-A*C
  if disc<=0:continue
  lo=max(0.,(-B-np.sqrt(disc))/A);hi=min(1.,(-B+np.sqrt(disc))/A)
  if hi>lo:intervals.append([float(sa+lo*(sb-sa)),float(sa+hi*(sb-sa))])
 merged=[]
 for a,b in intervals:
  if merged and a<=merged[-1][1]+1e-14:merged[-1][1]=max(merged[-1][1],b)
  else:merged.append([a,b])
 return merged
R=32e-6;straight=np.array([[-1e-3,0,0],[0,0,0],[1e-3,0,0]]);single=components(np.array([0,0,31e-6]),straight,np.array([0,1e-3,2e-3]),R)
u=np.array([[-15e-6,0,0],[-15e-6,500e-6,0],[15e-6,500e-6,0],[15e-6,0,0]]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(u,axis=0),axis=1))];multiple=components(np.array([0,20e-6,0]),u,arc,R)
out=dict(straight_support_across_segment_boundary=single,folded_polyline_disjoint_supports=multiple,passed=len(single)==1 and len(multiple)==2,meaning='Merge adjacent segment supports in material arclength. More than one connected support requires separate tracked component IDs; partner identity alone is insufficient. This is a topology counterexample, not a qualified rod shape.');(P/'receipts/component_guard.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
