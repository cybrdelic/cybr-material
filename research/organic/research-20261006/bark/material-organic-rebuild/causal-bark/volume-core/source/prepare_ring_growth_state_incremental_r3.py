"""Equilibrate the intact, minimally prestressed birth state before a growth ramp.
Uses the true nonlinear energy/gradient and a reference-tangent search metric.
There is no cohesive history evolution here: damage must remain identically zero,
and maximum opening must remain below every damage-initiation threshold.
"""
from pathlib import Path
import json,time,hashlib,resource
import numpy as np
from scipy.sparse import coo_matrix,block_diag,kron,eye
from scipy.sparse.linalg import splu
from incremental_ring_r3 import IncrementalRing as RingSystem
from source_receipt_r3 import hashes
R=Path(__file__).resolve().parents[1]
def reference_stiffness(m):
 rows=[];cols=[];vals=[]
 for p,(c1,c2) in enumerate(m.bonds.ops['pairs']):
  jc=m.bonds.ops['J'][p];idx=np.r_[c1*18+np.arange(18),c2*18+np.arange(18)];active=np.flatnonzero(abs(jc)>1e-15);i=idx[active];a=jc[active];rows.extend(np.repeat(i,len(i)));cols.extend(np.tile(i,len(i)));vals.extend((m.bonds.ops['area'][p]*m.bonds.ops['kt'][p]*np.outer(a,a)).ravel())
 C=coo_matrix((vals,(rows,cols)),shape=(m.N*18,m.N*18)).tocsr();K=block_diag(list(m.K),format='csr')+kron(C,eye(3),format='csr');return K

def run(label='ring8_intact_birth_equilibrium_incremental_verified_r3',sectors=8):
 dest=R/'data'/label
 if dest.exists():raise FileExistsError(dest)
 m=RingSystem(sectors,0.);dep=hashes();start=time.time();K=reference_stiffness(m);Kf=K[m.free][:,m.free].tocsc();factor=splu(Kf);q=np.zeros_like(m.q);trace=[];ok=False
 for it in range(101):
  E,g,opening,rep=m.evaluate_displacement(q);gf=g.ravel()[m.free];norm=float(np.linalg.norm(gf));support=float(np.linalg.norm(g.ravel()[m.fixed]));ratio=norm/max(support,1e-12);threshold=float(np.max(opening/m.bonds.delta0));trace.append({'iteration':it,'stored_energy_J':E,'free_force_N':norm,'support_force_N':support,'free_support_ratio':ratio,'maximum_opening_over_initiation':threshold})
  if threshold>=1:raise RuntimeError('Intact-state preparation reached damage initiation')
  if ratio<1e-9:ok=True;break
  p=np.zeros(q.size);p[m.free]=factor.solve(-gf);p=p.reshape(q.shape);descent=float(np.sum(g*p))
  if descent>=0:raise RuntimeError('Reference metric is not a descent direction')
  step=1.
  for ls in range(30):
   trial=q+step*p
   try:Et,gt,ot,rt=m.evaluate_displacement(trial)
   except ValueError:Et=float('inf')
   # At the final floating-point floor, require a falling true gradient as well
   # as an energy tolerance derived from uncancelled local work.
   roundoff=32*np.finfo(float).eps*max(abs(E),abs(Et),1e-20)
   if Et<=E+.0001*step*descent or (Et<=E+roundoff and np.linalg.norm(gt.ravel()[m.free])<norm):q=trial;break
   step*=.5
  else:raise RuntimeError('No energy-descent preparation step')
  trace[-1]['step_fraction']=step
 if not ok:
  np.savez_compressed(R/'receipts'/(label+'_failure_state.npz'),q=q);fullE,fullg,_,_=m.evaluate_displacement(q+p);fail={'completed':False,'reason':'Preparation did not converge','full_step_energy_difference_J':fullE-E,'full_step_gradient_N':float(np.linalg.norm(fullg.ravel()[m.free])),'predicted_energy_difference_J':descent,'maximum_increment_m':float(np.max(np.abs(p))),'trace':trace,'elapsed_s':time.time()-start};(R/'receipts'/(label+'_failure.json')).write_text(json.dumps(fail,indent=2));print(json.dumps(trace[-3:],indent=2));raise RuntimeError('Preparation did not converge')
 m.bonds.evolve(opening);assert not m.bonds.damage.any();assert dep==hashes();dest.mkdir(parents=True);np.savez_compressed(dest/'state.npz',q=m.origin+q,displacement=q,velocity=np.zeros_like(q),damage=m.bonds.damage,maximum_opening=m.bonds.maximum_opening,q_reference_support=m.q0)
 out={'revision':'BIRTH_EQUILIBRIUM_INCREMENTAL_R3','completed':ok,'sectors':sectors,'stem_increment_m':0.,'final':trace[-1],'trace':trace,'elapsed_s':time.time()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'source_sha256':dep,'dependencies_unchanged':True,'state_sha256':hashlib.sha256((dest/'state.npz').read_bytes()).hexdigest(),'scope':'Stationary intact initial state. No fracture avalanche, damping history, or growth work is inferred from this static preparation.'};(dest/'receipt.json').write_text(json.dumps(out,indent=2));(R/'receipts'/(label+'.json')).write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ('trace','source_sha256')},indent=2));return out
if __name__=='__main__':run()
