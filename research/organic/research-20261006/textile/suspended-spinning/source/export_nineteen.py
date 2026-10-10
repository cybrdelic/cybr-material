from pathlib import Path
import numpy as np,json,hashlib
from scipy.spatial import ConvexHull
P=Path(__file__).resolve().parents[1];state=P/'data/yarn_19_qualified.npz';z=np.load(state);h=json.loads((P/'receipts/yarn_19_qualified.json').read_text());assert h['passed'] and h['stability']['passed'];p=z['points'];r=float(z['radius_m']);rows=[];phi=np.arange(64)*2*np.pi/64;circle=np.c_[np.cos(phi),np.sin(phi)]
for x in [.0004,.0008,.0012,.0016,.002]:
 sections=[]
 for i in range(len(p)):
  for j in range(p.shape[1]-1):
   a,b=p[i,j:j+2]
   if a[0]<=x<b[0]:
    u=(x-a[0])/(b[0]-a[0]);center=(a+u*(b-a))[1:];t=b-a;t/=np.linalg.norm(t);Q=np.eye(2)-np.outer(t[1:],t[1:]);v,V=np.linalg.eigh(Q);shape=V@np.diag(1/np.sqrt(v))@V.T;sections.append((i,center,r*circle@shape.T+center,np.pi*r*r/abs(t[0])));
 for label in ['all','clamped_core']:
  keep=[s for s in sections if label=='all' or z['endpoint_fixed'][s[0]]];cloud=np.concatenate([s[2] for s in keep]);area=ConvexHull(cloud).volume;rows.append(dict(x_m=x,subset=label,crossing_fibres=len(keep),summed_fibre_section_area_m2=float(sum(s[3] for s in keep)),sampled_ellipse_envelope_area_m2=float(area),packing_fraction=float(sum(s[3] for s in keep)/area)))
meta=dict(fibres=len(p),stock_length_per_fibre_m=.0024,radius_m=r,state=str(state),state_sha256=hashlib.sha256(state.read_bytes()).hexdigest(),normal_contact=dict(max_penetration_m=h['capsule_penetration_m']),maximum_length_drift=h['relative_length_error'],root_error_m=h['root_constraint_error_m'],force_residual=h['maximum_scaled_KKT_residual'],stability=h['stability'],packing_sections=rows,scope='19 actual fibres in a stable nonlinear normal-contact equilibrium with prescribed intrinsic crimp and coherent hexagonal roots; 18 right roots clamped and one true free tip. No repeated seven-fibre substitution, no opaque core, no simulated spinning or finite forming-friction history, and spatial refinement remains open. This is a short yarn construction/material diagnostic, not a carpet panel.',render_promotion_scope='May show the qualified force/contact state; not a selected carpet material.');(P/'receipts/nineteen_geometry.json').write_text(json.dumps(meta,indent=2));print(json.dumps({k:v for k,v in meta.items() if k!='stability'},indent=2))
