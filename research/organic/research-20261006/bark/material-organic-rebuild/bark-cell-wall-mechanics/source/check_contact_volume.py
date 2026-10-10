from pathlib import Path
import json,numpy as np
from corrugated_wall import CorrugatedWall
P=Path(__file__).resolve().parents[1]
def distance(a,b,c,d):
 u=b-a;v=d-c;cross=lambda u,v:u[0]*v[1]-u[1]*v[0];den=cross(u,v)
 if abs(den)>1e-30:
  r=c-a;t=cross(r,v)/den;s=cross(r,u)/den
  if 0<=t<=1 and 0<=s<=1:return 0.
 def point(p,a,b):
  e=b-a;t=np.clip((p-a)@e/(e@e),0,1);return np.linalg.norm(p-(a+t*e))
 return min(point(a,c,d),point(b,c,d),point(c,a,b),point(d,a,b))
out=[]
for kind in ('compression','tension'):
 z=np.load(P/f'data/refined_wall_192_{kind}.npz');w=CorrugatedWall(192);ref=w.X*w.L;arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(ref,axis=0),axis=1))];mid=(arc[:-1]+arc[1:])/2;stock=np.diff(arc)
 for index in (0,len(z['positions_m'])-1):
  points=z['positions_m'][index];minimum=np.inf;pair=None
  for i in range(w.n):
   for j in range(i+2,w.n):
    if mid[j]-mid[i]<2*w.p.wall_thickness_m:continue
    dd=distance(points[i],points[i+1],points[j],points[j+1])
    if dd<minimum:minimum=dd;pair=(i,j)
  length=np.linalg.norm(np.diff(points,axis=0),axis=1);density=w.p.wall_density_kg_m3*stock/length;current_segment_volume=w.A*length;mass=float(density@current_segment_volume);d=np.diff(points,axis=0);theta=np.unwrap(np.arctan2(d[:,1],d[:,0]));curvature=np.diff(theta)/((length[1:]+length[:-1])/2)
  out.append(dict(control=kind,engineering_strain=float(z['engineering_strains'][index]),minimum_nonlocal_centerline_distance_m=float(minimum),wall_thickness_m=w.p.wall_thickness_m,nonlocal_pair=pair,material_arc_exclusion_m=2*w.p.wall_thickness_m,max_half_thickness_curvature=float(np.max(abs(curvature))*w.p.wall_thickness_m/2),material_mass_relative_error=abs(mass-w.mass_kg)/w.mass_kg,current_volume_m3=float(current_segment_volume.sum()),reference_volume_m3=w.volume_m3,current_density_range_kg_m3=[float(density.min()),float(density.max())],passed=bool(minimum>w.p.wall_thickness_m and np.max(abs(curvature))*w.p.wall_thickness_m/2<1 and abs(mass-w.mass_kg)/w.mass_kg<1e-12)))
r=dict(scope='Post-solve absence-of-contact gate, not an interacting cell-wall contact solver. Fixed strip section with density advected from reference stock; material compressibility not calibrated.',checks=out,passed=all(r['passed'] for r in out));(P/'receipts/contact_volume_check.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));assert r['passed']
