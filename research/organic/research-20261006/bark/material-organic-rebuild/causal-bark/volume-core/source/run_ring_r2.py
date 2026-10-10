"""Energy-accounted relaxation of retained curved cohorts in a closed ring.
Growth is an explicit prescribed initial condition, not a solved biochemical
source. All clamps hold actual inner-surface positions; no tangential foundation.
"""
from pathlib import Path
import numpy as np,json,time,hashlib,sys,resource
from ring_system_r2 import RingSystem
R=Path(__file__).resolve().parents[1]
def hashes():
 out={}
 for mod in list(sys.modules.values()):
  fn=getattr(mod,'__file__',None)
  if fn:
   p=Path(fn).resolve()
   if p.is_file() and p.is_relative_to(R.parent) and p.suffix in ('.py','.so'):out[str(p.relative_to(R.parent))]=hashlib.sha256(p.read_bytes()).hexdigest()
 for name in ('libvolume_native.so','libcohesion_native.so','batch_volume.cpp','cohesion_batch.cpp'):
  p=R/'source'/name;out['volume-core/source/'+name]=hashlib.sha256(p.read_bytes()).hexdigest()
 return out

def run(label='ring8_growth6mm_probe',sectors=8,growth=.006,factor=.35,duration=.0005):
 dest=R/'data'/label
 if dest.exists():raise FileExistsError('Frozen run exists')
 m=RingSystem(sectors,growth);source=hashes();omega=m.maximum_reference_frequency();n=int(np.ceil(duration*omega/factor));dt=duration/n;rate=20000.;tau=np.exp(-rate*dt/2);q=m.q.copy();v=np.zeros_like(q);E,g,op,report=m.evaluate(q);E0=E;heat=0.;maxerr=0.;start=time.time();trace=[];failure=None;print('START',json.dumps({'cells':m.N,'DOFs':q.size,'dt_s':dt,'steps':n,'initial_stored_energy_J':E0,'reference_omega_max_per_s':omega}),flush=True)
 try:
  for i in range(n):
   K=m.kinetic(v);v*=tau;heat+=K-m.kinetic(v);vh=v+.5*dt*m.acceleration(g);q+=dt*vh;_,_,op,_=m.evaluate(q);m.bonds.evolve(op);E,g,op,report=m.evaluate(q);v=vh+.5*dt*m.acceleration(g);K=m.kinetic(v);v*=tau;heat+=K-m.kinetic(v);K=m.kinetic(v);D=m.bonds.dissipation();error=E+K+D+heat-E0;maxerr=max(maxerr,abs(error))
   if i%max(1,n//150)==0 or i==n-1:
    damage=m.bonds.damage;kind=m.ops['kind'];layer=np.array([m.metadata[int(a)]['layer'] for a,b in m.ops['pairs']]);row={'step':i+1,'time_s':(i+1)*dt,'stored_J':E,'kinetic_J':K,'cohesive_dissipation_J':D,'damping_heat_J':heat,'balance_error_J':error,'minimum_det_F':report['minimum_det_F'],'young_fracture_points':int(np.sum((damage>.01)&(kind=='fracture')&(layer==0))),'old_fracture_points':int(np.sum((damage>.01)&(kind=='fracture')&(layer==1))),'birth_interface_damage_points':int(np.sum((damage>.01)&(kind=='birth_interface'))),'fully_failed_points':int(np.sum(damage>1-1e-12)),'free_force_N':float(np.linalg.norm(g.ravel()[m.free])),'support_force_N':float(np.linalg.norm(g.ravel()[m.fixed]))};trace.append(row)
   if i and i%3000==0:print('PROGRESS',i,n,'elapsed',time.time()-start,'energy_error',error/E0,flush=True)
 except Exception as exc:failure=type(exc).__name__+': '+str(exc)
 unchanged=source==hashes();out={'revision':'CLOSED_RADIAL_BIRTH_VOLUME_R2_EXPERIMENT','completed':failure is None and i==n-1,'failure':failure,'cells':m.N,'sectors':sectors,'stem_increment_m':growth,'dt_s':dt,'duration_s':duration,'mass_scaling':False,'damping_rate_per_s':rate,'maximum_relative_energy_error':maxerr/max(abs(E0),1e-20),'initial_energy_J':E0,'material_mass_kg':float(m.bulk.mass.sum()),'final':trace[-1] if trace else None,'trace':trace,'elapsed_s':time.time()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'dependencies_unchanged':unchanged,'source_sha256':source,'scope':'Prescribed initial growth of material-born, curved, closed cohorts. Uncalibrated isotropic constitutive and damping inputs; no natural-surface realism or multi-season formation qualification.'};dest.mkdir(parents=True);np.savez_compressed(dest/'state.npz',q=q,q_initial=m.q0,velocity=v,damage=m.bonds.damage,maximum_opening=m.bonds.maximum_opening,pairs=m.ops['pairs'],kind=m.ops['kind'],area=m.ops['area']);out['state_sha256']=hashlib.sha256((dest/'state.npz').read_bytes()).hexdigest();(dest/'receipt.json').write_text(json.dumps(out,indent=2));(R/'receipts'/(label+'.json')).write_text(json.dumps(out,indent=2));print('RESULT',json.dumps({k:v for k,v in out.items() if k not in ('trace','source_sha256')}),flush=True)
 if not unchanged:raise RuntimeError('Source changed during execution')
 return out
if __name__=='__main__':run()
