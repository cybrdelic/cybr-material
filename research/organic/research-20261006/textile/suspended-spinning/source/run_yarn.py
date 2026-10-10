import resource
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
from pathlib import Path
import sys,json,numpy as np
from scipy.optimize import brentq
from newton_globalized import MatrixFreeBundle,log,root_blocks,roots_tv,mv
from line_contact import audit
P=Path(__file__).resolve().parents[1];N=int(sys.argv[1]) if len(sys.argv)>1 else 19;amp=40e-6;label=f'yarn_{N}_128';c=MatrixFreeBundle(N,128,crimp_amplitude_m=amp,crimp_wavelength_m=.0012);d=.00238;c.clamp_twist_rad=0;c.controls(d,0);t=np.linspace(0,1,2049)
def points(A):
 p=np.c_[d*t,np.zeros_like(t),A*np.sin(np.pi*t)**2];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))];ss=np.linspace(0,arc[-1],c.M+1);return np.column_stack([np.interp(ss,arc,p[:,i]) for i in range(3)])
A=brentq(lambda a:np.linalg.norm(np.diff(points(a),axis=0),axis=1).sum()-.0024,0,.0004,xtol=1e-14);p=points(A);T=np.diff(p,axis=0);T/=np.linalg.norm(T,axis=1)[:,None];D1=np.c_[T[:,2],np.zeros(c.M),-T[:,0]];R=np.stack([D1,np.cross(T,D1),T],axis=-1);c.w[:]=log(R)[None,:,:]
for _ in range(6):
 err=c.endpoint_constraint();J,iv=root_blocks(c);c.w[:]=(c.w.ravel()-roots_tv(J,mv(iv,err.reshape(c.N,3)))).reshape(c.w.shape)
init=dict(root_error_m=float(abs(c.endpoint_constraint()).max()*.0024),contact=audit(c.points(),c.radius,c.ref,self_contact=True));print('INITIAL',json.dumps(init),flush=True);assert init['contact']['max_penetration_m']<1e-6;initial=c.points().copy()
try:
 c.checkpoint_path=P/'data'/(label+'_accepted_checkpoint.npz');row=c.solve(d,0,wall_s=120,maxiter=240);row['initialization']=init;row['scope']='Bounded actual-fibre short yarn with coherent hex root packing and prescribed intrinsic crimp; N-1 right clamps and one true free right endpoint. Nonlinear normal equilibrium; finite forming friction history remains absent. No appearance or wool material qualification.';row['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024;np.savez_compressed(P/'data'/(label+'.npz'),w=c.w,points=c.points(),initial_points=initial,material_arc_m=c.rod.s,endpoint_fixed=c.endpoint_fixed,intrinsic_curvature_per_m=np.array([r.k0 for r in c.rods]),roots=np.stack([c.origins,c.targets],axis=1))
except Exception as e:row=dict(passed=False,error=type(e).__name__+': '+str(e))
(P/'receipts'/(label+'.json')).write_text(json.dumps(row,indent=2));print('RESULT',json.dumps({k:v for k,v in row.items() if k not in ['solver_iteration_trace','root_positions_m','root_forces_world_N','clamp_torques_world_Nm']}),flush=True)
