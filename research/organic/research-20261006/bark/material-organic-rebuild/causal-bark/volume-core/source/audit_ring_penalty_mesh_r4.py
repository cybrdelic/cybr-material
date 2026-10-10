"""Same intact birth history, spatial and numerical-compliance refinement.
Only small receipts are saved. This does not establish fracture-path convergence.
"""
from pathlib import Path
import json,time,resource,numpy as np
from scipy.sparse.linalg import splu
from metric_ring_r4 import MetricRing
from prepare_ring_growth_state_r3 import reference_stiffness
from source_receipt_r3 import hashes
R=Path(__file__).resolve().parents[1]
def solve(sectors,eta):
 m=MetricRing(sectors,continuity_compliance=eta);K=reference_stiffness(m);kf=K[m.free][:,m.free].tocsc();solver=splu(kf);u=np.zeros_like(m.q);trace=[];t0=time.time()
 for it in range(40):
  E,g,op,rep=m.evaluate_displacement(u);norm=float(np.linalg.norm(g[m.free_mask]));support=float(np.linalg.norm(g[~m.free_mask]));ratio=norm/max(support,1e-20);trace.append({'iteration':it,'force_ratio':ratio,'energy_J':E})
  if np.max(op/m.bonds.delta0)>=1:raise RuntimeError('Refinement left the intact regime')
  if ratio<1e-9:break
  d=np.zeros(u.size);d[m.free]=solver.solve(-g.ravel()[m.free]);d=d.reshape(u.shape);s=1.;descent=float(np.sum(d*g))
  if descent>=0:raise RuntimeError('Invalid metric descent')
  for ls in range(30):
   Et,gt,ot,rt=m.evaluate_displacement(u+s*d);roundoff=32*np.finfo(float).eps*max(abs(E),abs(Et),1e-20)
   if Et<E+.0001*s*descent or (Et<E+roundoff and np.linalg.norm(gt[m.free_mask])<norm):u+=s*d;break
   s*=.5
  else:raise RuntimeError('No accepted refinement equilibrium step')
 else:raise RuntimeError('Equilibrium residual not converged')
 q=m.origin+u;rad=np.sqrt(q[:,:,0]**2+q[:,:,2]**2);rad0=np.sqrt(m.origin[:,:,0]**2+m.origin[:,:,2]**2);old=np.arange(m.N)>=2*sectors
 return {'sectors':sectors,'continuity_compliance':eta,'energy_J':E,'force_ratio':ratio,'maximum_opening_over_initiation':float(np.max(op/m.bonds.delta0)),'mean_old_radial_displacement_m':float(np.mean((rad-rad0)[old])),'maximum_gap_m':float(np.max(op)),'material_mass_kg':float(m.bulk.mass.sum()),'elapsed_s':time.time()-t0,'iterations':it,'trace':trace}
if __name__=='__main__':
 dep=hashes();start=time.time();rows=[]
 for n in (8,16,32):
  for eta in (.02,.01,.005):
   row=solve(n,eta);rows.append(row);print('CASE',json.dumps({k:v for k,v in row.items() if k!='trace'}),flush=True)
 out={'completed':True,'cases':rows,'dependencies_unchanged':dep==hashes(),'source_sha256':dep,'elapsed_s':time.time()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'scope':'Intact birth-state mesh/compliance study. No fracture morphology, tissue growth loading, or material calibration qualification.'};(R/'receipts/ring_penalty_mesh_r4.json').write_text(json.dumps(out,indent=2));print('RESULT',json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','cases')}),flush=True)
