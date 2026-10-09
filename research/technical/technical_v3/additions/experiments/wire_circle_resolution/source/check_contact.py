"""Exact centerline-segment clearance gives a conservative tube-surface bound."""
import bpy,numpy as np,math,json,itertools
from pathlib import Path
R=Path(__file__).resolve().parents[1];bpy.ops.wm.open_mainfile(filepath=str(R/"scenes/wire_circle32_native16.blend"))
me=bpy.context.scene.objects["Lattice / actual strand network"].data;a=np.empty(len(me.vertices)*3,"f4");me.vertices.foreach_get("co",a);q=a.reshape(28,157,32,3).astype(float)
c=q.mean(2);radius=float(np.max(np.linalg.norm(q-c[:,:,None],axis=-1)))
def pointseg(p,a,b):
 v=b-a;t=np.clip(np.dot(p-a,v)/np.dot(v,v),0,1);return float(np.dot(p-a-t*v,p-a-t*v))
def segments(a,b,c,d):
 u=b-a;v=d-c;w=a-c;A=np.dot(u,u);B=np.dot(u,v);C=np.dot(v,v);D=np.dot(u,w);E=np.dot(v,w);den=A*C-B*B
 values=[pointseg(a,c,d),pointseg(b,c,d),pointseg(c,a,b),pointseg(d,a,b)]
 if den>1e-32:
  s=(B*E-C*D)/den;t=(A*E-B*D)/den
  if 0<=s<=1 and 0<=t<=1:values.append(float(np.dot(w+s*u-t*v,w+s*u-t*v)))
 return min(values)
best=float(np.dot(c[0,0]-c[14,0],c[0,0]-c[14,0]));pair=None;tested=0
for i,j in itertools.combinations(range(28),2):
 x,y=c[i],c[j];al=np.minimum(x[:-1],x[1:]);ah=np.maximum(x[:-1],x[1:]);bl=np.minimum(y[:-1],y[1:]);bh=np.maximum(y[:-1],y[1:])
 sep=np.maximum(np.maximum(al[:,None]-bh[None],bl[None]-ah[:,None]),0);lower=np.sum(sep*sep,axis=-1)
 for k,l in np.argwhere(lower<=best+1e-20):
  d=segments(x[k],x[k+1],y[l],y[l+1]);tested+=1
  if d<best:best=d;pair=[i,j,int(k),int(l)]
clearance=math.sqrt(best)-2*radius
assert clearance>0,clearance
report=dict(all_strand_pairs=378,exact_segment_pairs_tested=tested,min_centerline_segment_distance_m=math.sqrt(best),maximum_vertex_radius_m=radius,conservative_surface_clearance_m=clearance,closest_pair=pair,passed=True,method="AABB lower bounds then exact3D segment distances. Every tube triangle lies inside its centerline-segment capsule; positive separation certifies no inter-strand intersections.")
(R/"receipts/contact.json").write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

