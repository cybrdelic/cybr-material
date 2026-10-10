from pathlib import Path
import numpy as np,json,time
from scipy.spatial.transform import Rotation
from ring_geometry import build
from ring_interfaces import build_interfaces
from native_cohesion import NativeCohesion
from surface_cohesion import evaluate
R=Path(__file__).resolve().parents[1];h,cells,q,meta,triangles=build(8);ops=build_interfaces(cells,q,meta,8);m=NativeCohesion(len(cells),ops);E0,g0,d0=m.evaluate(q);checks=[{'case':'all_born_interfaces_coincident','maximum_initial_opening_m':float(d0.max()),'energy_J':E0,'pass':float(d0.max())<1e-14}];rng=np.random.default_rng(932);x=q+rng.normal(size=q.shape)*1e-6;m.damage=rng.uniform(0,.8,m.N);En,gn,dn=m.evaluate(x);g=np.zeros_like(x);E=0;d=[]
for p,(i,j) in enumerate(ops['pairs']):
 J=np.kron(ops['J'][p][None,:],np.eye(3));T1=np.kron(ops['T1'][p][None,:],np.eye(3));T2=np.kron(ops['T2'][p][None,:],np.eye(3));qq=np.r_[x[i].ravel(),x[j].ravel()];e,gg,info=evaluate(qq,J,T1,T2,ops['area'][p],ops['kt'][p],ops['kn'][p],m.damage[p],ops['normal_sign'][p]);E+=e;g[i]+=gg[:54].reshape(18,3);g[j]+=gg[54:].reshape(18,3);d.append(info['equivalent_opening_m'])
eerr=abs(En-E)/E;gerr=np.linalg.norm(gn-g)/np.linalg.norm(g);checks.append({'case':'native_matches_scalar_surface_law','energy_relative_error':eerr,'gradient_relative_error':gerr,'pass':eerr<1e-10 and gerr<1e-9})
Q=Rotation.from_rotvec([.8,-.4,1.2]).as_matrix();Er,gr,_=m.evaluate(x@Q.T+[.2,-.3,.4]);obj=max(abs(Er-En)/En,np.linalg.norm(gr-gn@Q.T)/np.linalg.norm(gn));checks.append({'case':'complete_interface_objectivity','maximum_relative_error':obj,'pass':obj<1e-9});net=np.linalg.norm(gn.sum((0,1)));tor=np.linalg.norm(np.cross(x,gn).sum((0,1)));checks.append({'case':'assembled_internal_force_moment','force_N':float(net),'moment_N_m':float(tor),'pass':net<1e-9 and tor<1e-10})
start=time.time()
for _ in range(40):m.evaluate(x)
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'cells':len(cells),'cohesive_points':m.N,'seconds_per_cohesive_evaluation':(time.time()-start)/40,'scope':'Closed-ring finite cohesive assembly with tissue-continuity test stiffness. No ring evolution or material appearance qualification.'};(R/'receipts/ring_interface_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
