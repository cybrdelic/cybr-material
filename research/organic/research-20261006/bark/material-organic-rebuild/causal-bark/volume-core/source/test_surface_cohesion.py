from pathlib import Path
import numpy as np,json
from scipy.spatial.transform import Rotation
from surface_cohesion import evaluate
R=Path(__file__).resolve().parents[1];rng=np.random.default_rng(989);I=np.eye(3);J=np.zeros((3,18));T1=J.copy();T2=J.copy();weights=np.array([.2,.35,.45])
for i in range(3):J[:,3*i:3*i+3]=-weights[i]*I;J[:,9+3*i:12+3*i]=weights[i]*I
for offset in (0,9):T1[:,offset:offset+3]=-.5*I;T1[:,offset+3:offset+6]=.5*I;T2[:,offset:offset+3]=-.5*I;T2[:,offset+6:offset+9]=.5*I
base=np.array([[0,0,0],[.002,0,0],[0,.0015,0.]]);checks=[];errs=[];obj=[];force=[];torque=[]
for sign in (-1,1):
 x=np.vstack((base,base+[.00004,.00002,sign*.00003]));x+=rng.normal(size=x.shape)*2e-6;q=x.ravel();E,g,d=evaluate(q,J,T1,T2,1e-6,1.6e7,4e7,.7);direction=rng.normal(size=18);direction/=np.linalg.norm(direction);eps=1e-8;fd=(evaluate(q+eps*direction,J,T1,T2,1e-6,1.6e7,4e7,.7)[0]-evaluate(q-eps*direction,J,T1,T2,1e-6,1.6e7,4e7,.7)[0])/(2*eps);errs.append(abs(fd-g@direction)/abs(fd));Q=Rotation.from_rotvec([.9,-.7,1.4]).as_matrix();qr=(x@Q.T+[.01,-.02,.03]).ravel();Er,gr,_=evaluate(qr,J,T1,T2,1e-6,1.6e7,4e7,.7);obj.append(max(abs(Er-E)/E,np.linalg.norm(gr.reshape(-1,3)-g.reshape(-1,3)@Q.T)/np.linalg.norm(g)));force.append(float(np.linalg.norm(g.reshape(-1,3).sum(0))));torque.append(float(np.linalg.norm(np.cross(x,g.reshape(-1,3)).sum(0))))
checks=[{'case':'opening_and_closing_exact_gradients','maximum_relative_error':max(errs),'pass':max(errs)<1e-7},{'case':'superposed_finite_rotation_objectivity','maximum_relative_error':max(obj),'pass':max(obj)<1e-10},{'case':'internal_force_and_moment_balance','force_error_N':max(force),'moment_error_N_m':max(torque),'pass':max(force)<1e-15 and max(torque)<1e-18}]
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'scope':'One finite-surface cohesive quadrature point. Does not provide non-neighbor collision detection or full fractured-body contact.'};(R/'receipts/surface_cohesion_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
