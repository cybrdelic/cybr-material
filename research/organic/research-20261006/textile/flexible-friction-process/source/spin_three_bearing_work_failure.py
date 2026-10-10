"""Bounded nonsymmetric frictional residual solve; immutable accepted history."""
import time,json,hashlib,resource,numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator,gmres,splu
from mixed_friction import P,BRIDGE,mixed_metric
from prepare_three import ThreeCarriage
from carriage import RodHessianOperator,matrix,certify_path
from rod_pullback import exp,right_jacobian,mv
class SpinThree(ThreeCarriage):
 def __init__(self):
  self.history=None;super().__init__(128)
 def evaluate_rods(self,a):
  if self.history is None:return super().evaluate_rods(a)
  c=self.c;c.w[:]=np.asarray(a).reshape(c.w.shape);R=exp(c.w);q=self.bridge.pack(c.origins,R);t=self.bridge.trial(q,self.history);self.last_trial=t;b=t['state']['g'];g=t['residual'].reshape(3,-1)[:,3:].reshape(c.w.shape).copy();E=t['normal_energy_J']+t['tangential_energy_J'];Eb=0.;torques=np.zeros((c.N,2,3))
  for i in range(c.N):
   e,gi=c.rods[i].energy_gradient(R[i],c.frames_for(i),c.boundary_k0[i]);E+=e;Eb+=e;g[i]+=mv(np.swapaxes(right_jacobian(c.w[i]),-1,-2),gi);_,_,eg=c.rods[i].boundary_terms(R[i],c.frames_for(i),c.boundary_k0[i]);torques[i]=[c.ends[j]@eg[j] for j in range(2)]
  gp=b['normal_gp']+t['friction_gp'];c.last=dict(energy_J=float(E),bending_energy_J=float(Eb),pair_contact_energy_J=t['normal_energy_J'],tangential_energy_J=t['tangential_energy_J'],backing_contact_energy_J=0.,position_gradient_N=gp,backing_gradient_N=np.zeros_like(gp),contact=dict(active_samples=int(t['active_set'].sum()),law='symmetric material-pair normal plus two-sided friction'),planes=dict(top_force_N=0.,bottom_force_N=0.,top_penetration_m=0.,bottom_penetration_m=0.,boundary_mode='suspended'),clamp_torques_world_Nm=torques,points=b['p']);return E/self.unit,g.ravel()/self.unit

def history_hash(h):return hashlib.sha256(h['q'].tobytes()+h['elastic'].tobytes()).hexdigest()
def solve_first(theta=1e-4,wall_s=90):
 start=time.monotonic();deadline=start+wall_s;s=SpinThree();z=np.load(BRIDGE/'data/three_refined128_equilibrium.npz');s.x=z['x'].copy();oldrow=s.record(s.x,True,[]);oldpoints=s.c.points().copy();s.history=s.bridge.state(s.bridge.pack(s.c.origins,exp(s.c.w)));initial_hash=history_hash(s.history);x=s.predictor_boundary(theta);s.theta=theta;_,g=s.evaluate(x);pg,lam,J,S=s.projection(g);nf=int(s.fixed.sum());y=np.r_[x,lam[s.fixed]];nx=len(x);trace=[];switches=[];evaluations=0;status='running';(P/'data').mkdir(exist_ok=True)
 def residual(y):
  nonlocal evaluations
  if time.monotonic()>deadline:raise TimeoutError('Bounded frictional corrector')
  evaluations+=1;_,g=s.evaluate(y[:nx]);J=s.jacobian();la=np.zeros(9);la[s.fixed]=y[nx:];return np.r_[g+J.T@la,s.constraints()[s.fixed]]
 def merit(r):return np.linalg.norm(np.r_[r[:nx],48*r[nx:]])
 try:
  for it in range(16):
   R=residual(y);base_t=s.last_trial;pg,projlam,J,S=s.projection(s.evaluate(y[:nx])[1]);con=s.constraints();row=dict(iteration=it,force_residual=float(abs(pg).max()),root_error_m=float(abs(con).max()*s.L),merit=float(merit(R)),elapsed_s=time.monotonic()-start);trace.append(row);print('SPIN',row,flush=True)
   if row['force_residual']<2e-5 and row['root_error_m']<1e-9:status='equilibrium';break
   lamfull=np.zeros(9);lamfull[s.fixed]=y[nx:];op=RodHessianOperator(s.c,lamfull);A,_=matrix(op,s,positive=True);H,B,stats=mixed_metric(s.bridge,base_t);HC=H.tocoo();A=A+coo_matrix((HC.data/s.unit,(HC.row,HC.col)),shape=A.shape).tocsc();factor=splu(A,permc_spec='COLAMD');n=op.size
   def inv(v):
    rhs=np.r_[v[:n],np.zeros(n),v[n],np.zeros(n),v[nx:]];sol=factor.solve(rhs);sol+=factor.solve(rhs-A@sol);return np.r_[sol[:n],sol[2*n],sol[3*n+1:]]
   base=y.copy();baseR=R.copy();active=base_t['active_set'];slip=base_t['sliding_set'];proj=base_t['partner_segments'];baseJ=s.jacobian().copy()
   def jvp(v):
    h=1e-7/max(np.max(abs(v[:nx])),1e-30);probe=base.copy();probe[:nx]+=h*v[:nx];rp=residual(probe);t=s.last_trial;switches.append(dict(active=int(np.count_nonzero(t['active_set']!=active)),slip=int(np.count_nonzero(t['sliding_set']!=slip)),projection=int(np.count_nonzero(t['partner_segments']!=proj))));d=(rp-baseR)/h;d[:nx]+=baseJ[s.fixed,:].T@v[nx:];return d
   counter=[];delta,info=gmres(LinearOperator((len(y),len(y)),matvec=jvp,dtype=float),-R,M=LinearOperator((len(y),len(y)),matvec=inv,dtype=float),rtol=1e-5,atol=0,restart=30,maxiter=2,callback=lambda e:counter.append(float(e)),callback_type='pr_norm');row.update(gmres_info=int(info),linear_iterations=len(counter),last_linear_residual=counter[-1] if counter else None,mixed_stats=stats);residual(base);p0=s.c.points().copy();angle=np.linalg.norm(delta[:n].reshape(s.c.w.shape),axis=2).max();scale=min(1.,.05/max(angle,1e-30));accepted=False
   for ls in range(18):
    try:
     candidate=base+scale*delta;rt=residual(candidate)
     if merit(rt)<(1-1e-4*scale)*merit(R):
      cert=certify_path(p0,s.c.points(),s.c.radius,s.c.ref,plane_bounds=None)
      if cert['passed']:accepted=True;break
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):pass
    scale*=.5
   row.update(step_scale=scale,line_search_trials=ls+1,accepted=accepted)
   if not accepted:status='no_residual_decreasing_contact_safe_step';break
   row['path_contact']=cert;y=candidate;np.savez_compressed(P/'data/first_spin_pending.npz',y=y,theta=theta,history_q=s.history['q'],history_elastic=s.history['elastic'],history_identity=s.history['identity']);(P/'receipts/first_spin_progress.json').write_text(json.dumps(dict(status='running',trace=trace,evaluations=evaluations),indent=2))
 except TimeoutError:status='bounded_pause'
 except Exception as e:status='error:'+type(e).__name__+':'+str(e)
 # Terminal inspection is allowed outside the corrector timer; no further solve.
 terminal=s.record(y[:nx],status=='equilibrium',trace);t=s.last_trial;W=.5*(oldrow['torque_conjugate_Nm']+terminal['torque_conjugate_Nm'])*theta+.5*(oldrow['tension_N']+terminal['tension_N'])*(terminal['carriage_m']-oldrow['carriage_m']);dE=terminal['stored_energy_J']-oldrow['stored_energy_J'];gap=W-dE-t['friction_heat_J']-t['numerical_loss_J'];den=max(abs(W),abs(dE)+t['friction_heat_J']+t['numerical_loss_J'],1e-18);path=certify_path(oldpoints,s.c.points(),s.c.radius,s.c.ref,plane_bounds=None);passed=status=='equilibrium' and abs(gap)/den<.02 and terminal['relative_stock_error']<1e-12 and terminal['capsule_penetration_m']<1e-6 and path['passed'] and history_hash(s.history)==initial_hash;terminal.update(status=status,passed=bool(passed),force_equilibrium=bool(status=='equilibrium'),work_J=W,stored_increment_J=dE,friction_heat_J=t['friction_heat_J'],return_map_numerical_loss_J=t['numerical_loss_J'],work_gap_J=gap,work_relative=abs(gap)/den,history_committed=bool(passed),old_history_unchanged=history_hash(s.history)==initial_hash,derivative_switch_counts=dict(active=sum(r['active'] for r in switches),slip=sum(r['slip'] for r in switches),projection=sum(r['projection'] for r in switches)),evaluations=evaluations,wall_s=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,endpoint_path_contact=path,scope='Internal first flexible frictional spin increment; no visible material promotion. Numerical directional tangent with switches explicitly recorded.');np.savez_compressed(P/'data/first_spin_terminal.npz',y=y,theta=theta,q=t['state']['q'],elastic=t['state']['elastic'],identity=t['state']['identity'],qualified=passed);(P/'receipts/first_spin_terminal.json').write_text(json.dumps(terminal,indent=2));print(json.dumps({k:v for k,v in terminal.items() if k not in ('trace','root_forces_N','root_torques_Nm')},indent=2));return terminal
if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));solve_first()
