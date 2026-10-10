from pathlib import Path
import json,numpy as np,hashlib
from scipy.spatial import ConvexHull
P=Path(__file__).resolve().parents[1]
old=P.parent/'spinning-carriage/ramp/batch3';new=P/'ramp03';paths=[(old/'data/accepted_004.npz',old/'receipts/accepted_004.json'),(new/'data/accepted_011.npz',new/'receipts/accepted_011.json')]
phi=np.linspace(0,2*np.pi,128,endpoint=False);circle=np.c_[np.cos(phi),np.sin(phi)];rows=[];zs=[]
for state,receipt in paths:
 z=np.load(state);h=json.loads(receipt.read_text());p=z['points'];r=float(z['radius_m']);cross=[];x=.5*h['carriage_m']
 for i in np.flatnonzero(z['endpoint_fixed']):
  for a,b in zip(p[i,:-1],p[i,1:]):
   if a[0]<=x<b[0]:
    u=(x-a[0])/(b[0]-a[0]);center=(a+u*(b-a))[1:];t=(b-a)/np.linalg.norm(b-a);e,V=np.linalg.eigh(np.eye(2)-np.outer(t[1:],t[1:]));ellipse=center+r*circle@(V@np.diag(1/np.sqrt(e))@V.T).T;cross.append((ellipse,np.pi*r*r/abs(t[0])))
 area=ConvexHull(np.concatenate([c[0] for c in cross])).volume;rows.append(dict(state=str(state),state_sha256=hashlib.sha256(state.read_bytes()).hexdigest(),bearing_angle_rad=h['twist_rad'],turns_per_m=h['twist_turns_per_m'],omega_R=h['dimensionless_root_twist'],carriage_m=h['carriage_m'],torque_Nm=h['torque_conjugate_Nm'],stored_energy_J=h['stored_energy_J'],midspan_core_packing_fraction=float(sum(c[1] for c in cross)/area),midspan_envelope_m2=float(area),minimum_fibre_surface_z_m=float(p[:,:,2].min()-r),free_tip_world_m=p[-1,-1].tolist(),force_residual=h['force_residual'],stability=h['stability']['minimum_generalized_eigenvalue']));zs.append(p)
delta=zs[1]-zs[0];out=dict(backed=rows[0],suspended=rows[1],maximum_material_labelled_displacement_m=float(np.linalg.norm(delta,axis=2).max()),rms_material_labelled_displacement_m=float(np.sqrt(np.mean(np.sum(delta**2,axis=2)))),same_bearing_angle=rows[0]['bearing_angle_rad']==rows[1]['bearing_angle_rad'],scope='Same nineteen stocks, radius, intrinsic crimp, root guide axes, force control and bearing angle. The suspended stage has independently re-equilibrated without support planes, so preparation and mechanical boundary differ. This is not an exact render A/B or a finite-friction spinning comparison.')
(P/'receipts/backed_suspended_same_angle.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
