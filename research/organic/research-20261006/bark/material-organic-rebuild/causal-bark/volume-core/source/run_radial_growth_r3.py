"""Gradual radial support growth with full moving-Dirichlet energy accounting.

Physical time is deliberately accelerated as a numerical fixture. The supplied
mass-proportional damping is uncalibrated. No biological growth-rate or mature
cork morphology is inferred. Birth metrics and mass remain fixed during this
loading interval; new material deposition has its own unimplemented gate.
"""
from pathlib import Path
import json,time,hashlib,resource,numpy as np
from incremental_ring_r3 import IncrementalRing
from source_receipt_r3 import hashes
R=Path(__file__).resolve().parents[1]

def boundary(t,ramp,growth,direction):
 if t<=0:s=ds=dds=0.
 elif t>=ramp:s=1.;ds=dds=0.
 else:
  x=t/ramp;s=x*x*x*(10+x*(-15+6*x));ds=30*x*x*(1-x)**2/ramp;dds=60*x*(1-x)*(1-2*x)/ramp**2
 return growth*s*direction,growth*ds*direction,growth*dds*direction

def run(label='ring8_radial_growth_1mm_r3_dt035',growth=.001,ramp=.0005,hold=.0005,factor=.35,alpha=20000.,initial='ring8_intact_birth_equilibrium_incremental_verified_r3'):
 if not np.isfinite([growth,ramp,hold,factor,alpha]).all() or growth<0 or ramp<=0 or hold<0 or factor<=0 or alpha<0:raise ValueError('Invalid growth loading parameters')
 dest=R/'data'/label
 if dest.exists():raise FileExistsError('Frozen run exists: '+str(dest))
 init=R/'data'/initial;initial_receipt=json.loads((init/'receipt.json').read_text());m=IncrementalRing(initial_receipt['sectors'],0.);s=np.load(init/'state.npz');ih=hashlib.sha256((init/'state.npz').read_bytes()).hexdigest()
 if ih!=initial_receipt['state_sha256']:raise ValueError('Initial state source digest mismatch')
 u=s['displacement'].copy();v=s['velocity'].copy();m.bonds.restore(s['maximum_opening'],s['damage']);assert not np.any(u[~m.free_mask]);direction=np.zeros_like(u);rad=np.sqrt(m.origin[:,:,0]**2+m.origin[:,:,2]**2);direction[:,:,0]=m.origin[:,:,0]/rad;direction[:,:,2]=m.origin[:,:,2]/rad;direction[m.free_mask]=0.
 dep=hashes();omega=m.maximum_reference_frequency();duration=ramp+hold;n=int(np.ceil(duration*omega/factor));dt=duration/n;decay=np.exp(-alpha*dt/2);gain=dt/2 if alpha==0 else -np.expm1(-alpha*dt/2)/alpha
 E,g,op,rep=m.evaluate_displacement(u);uc,vc,ac=boundary(0,ramp,growth,direction);v[~m.free_mask]=vc[~m.free_mask];F=m.moving_force(g,vc,ac,alpha);a,res=m.moving_acceleration_reaction(g,v,ac,alpha);power=float(np.sum(v[~m.free_mask]*res[~m.free_mask]));heatpower=alpha*float(np.sum(v*m.mass_product(v)));E0=E+m.kinetic(v)+m.bonds.dissipation();work=heat=0.;maxerr=0.;trace=[];failure=None;start=time.time();dest.mkdir(parents=True)
 config={'revision':'GRADUAL_RADIAL_BIRTH_RING_R3','growth_m':growth,'ramp_s':ramp,'hold_s':hold,'factor':factor,'alpha_per_s':alpha,'dt_s':dt,'steps':n,'mass_scaling':False,'initial_state_sha256':ih,'initial_state':str(init),'source_sha256':dep};(dest/'configuration.json').write_text(json.dumps(config,indent=2));print('START',json.dumps({k:v for k,v in config.items() if k!='source_sha256'}),flush=True)
 def row(i,E,g,rep,error):
  kind=m.ops['kind'];layer=np.array([m.metadata[int(k)]['layer'] for k,j in m.ops['pairs']]);d=m.bonds.damage
  return {'step':i,'time_s':i*dt,'stored_J':E,'kinetic_J':m.kinetic(v),'fracture_dissipation_J':m.bonds.dissipation(),'damping_heat_J':heat,'support_work_J':work,'balance_error_J':error,'minimum_det_F':rep['minimum_det_F'],'free_gradient_N':float(np.linalg.norm(g[m.free_mask])),'free_dynamic_residual_N':float(np.linalg.norm(res[m.free_mask])),'true_support_reaction_N':float(np.linalg.norm(res[~m.free_mask])),'young_damage_points':int(np.sum((d>.01)&(kind=='fracture')&(layer==0))),'old_damage_points':int(np.sum((d>.01)&(kind=='fracture')&(layer==1))),'birth_interface_damage_points':int(np.sum((d>.01)&(kind=='birth_interface'))),'fully_failed_points':int(np.sum(d>=1-1e-12))}
 trace.append(row(0,E,g,rep,0.));i=-1
 try:
  for i in range(n):
   # Audited symmetric exponential kick/drift/kick with exact prescribed data.
   vh=decay*np.where(m.free_mask,v,0.)+gain*F;u+=dt*vh;uc,vc,ac=boundary((i+1)*dt,ramp,growth,direction);u[~m.free_mask]=uc[~m.free_mask]
   E,g,op,rep=m.evaluate_displacement(u);old_damage=m.bonds.damage.copy();m.bonds.evolve(op)
   if not np.array_equal(old_damage,m.bonds.damage):E,g,op,rep=m.evaluate_displacement(u)
   F=m.moving_force(g,vc,ac,alpha);v=decay*vh+gain*F;v[~m.free_mask]=vc[~m.free_mask];a,res=m.moving_acceleration_reaction(g,v,ac,alpha);newpower=float(np.sum(v[~m.free_mask]*res[~m.free_mask]));newheatpower=alpha*float(np.sum(v*m.mass_product(v)));work+=.5*dt*(power+newpower);heat+=.5*dt*(heatpower+newheatpower);power,heatpower=newpower,newheatpower;K=m.kinetic(v);D=m.bonds.dissipation();error=E+K+D+heat-E0-work;maxerr=max(maxerr,abs(error))
   if i%max(1,n//150)==0 or i==n-1:trace.append(row(i+1,E,g,rep,error))
   if i and i%6000==0:print('PROGRESS',i,n,'elapsed_s',round(time.time()-start,2),'balance_error_J',error,flush=True)
 except Exception as exc:failure=type(exc).__name__+': '+str(exc)
 unchanged=dep==hashes();np.savez_compressed(dest/'state.npz',displacement=u,q=m.origin+u,velocity=v,damage=m.bonds.damage,maximum_opening=m.bonds.maximum_opening,pairs=m.ops['pairs'],kind=m.ops['kind'],area=m.ops['area']);out=dict(config,completed=failure is None and i==n-1,failure=failure,initial_total_energy_J=E0,maximum_balance_error_J=maxerr,maximum_balance_error_relative_to_input=maxerr/max(abs(E0),abs(work),1e-20),final=trace[-1],trace=trace,elapsed_s=time.time()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,dependencies_unchanged=unchanged,state_sha256=hashlib.sha256((dest/'state.npz').read_bytes()).hexdigest(),scope='Moving radial support of retained birth-metric tissue, with full consistent-mass inertia and support work. Accelerated uncalibrated material/damping fixture. No new tissue deposition or morphology qualification.');(dest/'receipt.json').write_text(json.dumps(out,indent=2));(R/'receipts'/(label+'.json')).write_text(json.dumps(out,indent=2));print('RESULT',json.dumps({k:v for k,v in out.items() if k not in ('trace','source_sha256')}),flush=True)
 if not unchanged:raise RuntimeError('Source/dependencies changed during run')
 if failure:raise RuntimeError(failure)
 return out
if __name__=='__main__':run()
