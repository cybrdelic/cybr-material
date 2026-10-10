"""Finite-volume/cohesive relaxation of an explicitly specified growth state.
The prescribed bottom face stays fixed. Initial stored energy is accounted,
with real inertia, cohesive dissipation, kinetic energy and explicit damping heat.
This does not claim to supply the biochemical work of growing the initial state.
"""
from pathlib import Path
import numpy as np,json,time,sys,hashlib
from scipy.linalg import eigh
from birth_pair import BirthPair
R=Path(__file__).resolve().parents[1]

def run(factor=.4,duration=.0002,label='pair_dt04_probe'):
 dest=R/'data'/label
 if dest.exists():raise FileExistsError('Frozen state exists')
 m=BirthPair();start=time.time();q=m.q.copy();q0=q.copy();E0,g,openings,volumes=m.evaluate(q);H=np.zeros((len(m.free),len(m.free)));eps=1e-8
 for j,index in enumerate(m.free):
  a=q.copy();b=q.copy();a[index]+=eps;b[index]-=eps;H[:,j]=(m.evaluate(a)[1][m.free]-m.evaluate(b)[1][m.free])/(2*eps)
 symmetry=np.linalg.norm(H-H.T)/np.linalg.norm(H);H=(H+H.T)/2;lam=eigh(H,m.Mf,eigvals_only=True);omega=np.sqrt(max(lam));n=int(np.ceil(duration*omega/factor));dt=duration/n;v=np.zeros(len(m.free));rate=20000.;tau=np.exp(-rate*dt/2);heat=0.;maximum_error=0.;trace=[];print('START',json.dumps({'steps':n,'dt_s':dt,'omega_max':omega,'tangent_symmetry_error':symmetry,'minimum_initial_eigenvalue':float(lam.min()),'initial_energy_J':E0}),flush=True)
 def kinetic(v):return float(.5*v@m.Mf@v)
 for i in range(n):
  K=kinetic(v);v*=tau;heat+=K-kinetic(v);vh=v-.5*dt*(m.Minv@g[m.free]);q[m.free]+=dt*vh;_,_,opening,_=m.evaluate(q);m.evolve_damage(opening);E,g,opening,volumes=m.evaluate(q);v=vh-.5*dt*(m.Minv@g[m.free]);K=kinetic(v);v*=tau;heat+=K-kinetic(v);K=kinetic(v);D=m.dissipation();balance=E+K+D+heat-E0;maximum_error=max(maximum_error,abs(balance))
  if i%max(1,n//100)==0 or i==n-1:trace.append({'time_s':(i+1)*dt,'stored_J':E,'kinetic_J':K,'cohesive_dissipation_J':D,'heat_J':heat,'energy_error_J':balance,'free_force_norm_N':float(np.linalg.norm(g[m.free])),'support_force_norm_N':float(np.linalg.norm(g[m.fixed])),'maximum_damage':float(m.damage.max()),'minimum_det_F':min(c['minimum_det_F'] for c in volumes)})
  if i and i%2000==0:print('PROGRESS',i,n,'seconds',time.time()-start,'relative_energy_error',balance/E0,flush=True)
 out={'completed':True,'initial_material_metadata':m.material_metadata,'dt_s':dt,'duration_s':duration,'damping_rate_per_s':rate,'no_mass_scaling':True,'maximum_relative_energy_error':maximum_error/E0,'tangent_symmetry_error':symmetry,'omega_max_initial_per_s':float(omega),'final':trace[-1],'trace':trace,'elapsed_s':time.time()-start,'scope':'Relaxation of the specified initial grown configuration. Young/old material reference and mass are retained. Growth injection work, multi-cell fracture morphology and measured cork constitutive response remain unqualified.','source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'source').glob('*.py')}};dest.mkdir(parents=True);np.savez_compressed(dest/'state.npz',q=q,q_initial=q0,velocity_free=v,damage=m.damage,maximum_opening=m.maximum_opening);out['state_sha256']=hashlib.sha256((dest/'state.npz').read_bytes()).hexdigest();(dest/'receipt.json').write_text(json.dumps(out,indent=2));(R/'receipts'/(label+'.json')).write_text(json.dumps(out,indent=2));print('RESULT',json.dumps({k:v for k,v in out.items() if k not in ['trace','source_hashes','initial_material_metadata']}),flush=True);return out
if __name__=='__main__':run(float(sys.argv[1]) if len(sys.argv)>1 else .4,float(sys.argv[2]) if len(sys.argv)>2 else .0002,sys.argv[3] if len(sys.argv)>3 else 'pair_dt04_probe')
