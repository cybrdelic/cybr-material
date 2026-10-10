"""Small interacting whole-fibre loop with physical root clamps and backing.
No postwarping, surface coating, opaque yarn core, or feed-reserve substitution.
The imposed formation controls move and rotate two rigid clamp cross-sections.
"""
from pathlib import Path
import sys,time
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from kirchhoff3d import *
from segment_contact_sine import evaluate as pair_energy
from line_contact import planes,audit

class LoopBundle:
 def __init__(self,fibres=7,segments=128,contact_order=16,seed_sign=1,crimp_amplitude_m=40e-6,crimp_wavelength_m=.0012,natural_twist_per_m=0.,plane_bounds=(0.,.02)):
  if fibres not in (3,7,19,37):raise ValueError('Only the bounded 3/7/19/37 fibre diagnostic is implemented.')
  self.plane_bounds=plane_bounds
  if plane_bounds is not None and (len(plane_bounds)!=2 or not np.isfinite(plane_bounds).all() or plane_bounds[1]<=plane_bounds[0]):raise ValueError("Invalid boundary planes")
  self.N=fibres;self.M=segments;self.endpoint_fixed=np.ones(fibres,dtype=bool);self.endpoint_fixed[-1]=False;self.rod=KirchhoffRod(np.linspace(0,.0024,segments+1));self.radius=self.rod.mat.radius_m;self.unit=self.rod.mat.B[0]/.0024;self.contact_order=contact_order;self.kc=1e8;self.ref=np.tile(self.rod.ds,(fibres,1));self.history=[];self.states=[]
  if fibres==3:
   a=2.15*self.radius/np.sqrt(3);theta=np.arange(3)*2*np.pi/3;self.offset=np.c_[a*np.cos(theta),a*np.sin(theta),np.zeros(3)]
  elif fibres==7:
   theta=np.arange(6)*np.pi/3;self.offset=np.vstack([np.zeros(3),np.c_[2.15*self.radius*np.cos(theta),2.15*self.radius*np.sin(theta),np.zeros(6)]])
  else:
   rings={19:2,37:3}[fibres];cells=[(a,b) for a in range(-rings,rings+1) for b in range(-rings,rings+1) if max(abs(a),abs(b),abs(a+b))<=rings];cells.sort(key=lambda ab:(max(abs(ab[0]),abs(ab[1]),abs(sum(ab))),np.arctan2(np.sqrt(3)*ab[1]/2,ab[0]+ab[1]/2)%(2*np.pi)));cells=np.array(cells);self.offset=2.15*self.radius*np.c_[cells[:,0]+cells[:,1]/2,cells[:,1]*np.sqrt(3)/2,np.zeros(fibres)]
   assert self.offset.shape==(fibres,3)

  wave=2*np.pi/crimp_wavelength_m;amplitude=crimp_amplitude_m*wave*wave;phases=np.arctan2(self.offset[:,1],self.offset[:,0]);self.rods=[];self.boundary_k0=[]
  for phase in phases:
   def intrinsic(s):return np.c_[amplitude*np.sin(wave*s+phase),np.zeros_like(s),np.full_like(s,natural_twist_per_m)]
   self.rods.append(KirchhoffRod(self.rod.s,self.rod.mat,intrinsic(self.rod.s[1:-1])));self.boundary_k0.append(intrinsic(np.array([0.,.0024])))
  self.boundary_k0=np.array(self.boundary_k0);self.clamp_twist_rad=0.;self.manufacturing_parameters=dict(intrinsic_crimp_amplitude_m=crimp_amplitude_m,intrinsic_crimp_wavelength_m=crimp_wavelength_m,intrinsic_material_twist_per_m=natural_twist_per_m,crimp_phase_rule='clamp cross-section polar material label',root_packing_pitch_over_fibre_radius=2.15,rest_state_classification='prescribed/calibrated intrinsic curvature, not simulated fibre manufacture')
  self.w=np.zeros((fibres,segments,3));self.w[:,:,1]=np.pi/2
  # Declared tiny non-equilibrium 3D imperfection chooses a bifurcation branch.
  # This initial imperfection is separate from the prescribed intrinsic crimp above.
  self.w[:,:,0]=seed_sign*1e-4*np.cos(np.arange(fibres))[:,None]*np.sin(2*np.pi*(np.arange(segments)+.5)/segments)[None,:]
  self.seed_sign=seed_sign
 def frames_for(self,i):return (self.ends[0],self.ends[1] if self.endpoint_fixed[i] else None)
 def endpoint_constraint(self):return (((self.points()[:,-1]-self.targets)/.0024)*self.endpoint_fixed[:,None]).ravel()
 def controls(self,d,alpha):
  self.ends=(exp(np.array([0.,np.pi/2-alpha,0.]))@exp(np.array([0.,0.,-self.clamp_twist_rad/2])),exp(np.array([0.,np.pi/2+alpha,0.]))@exp(np.array([0.,0.,self.clamp_twist_rad/2])));start_offset=self.offset@self.ends[0].T;end_offset=self.offset@self.ends[1].T;height=self.radius-min(start_offset[:,2].min(),end_offset[:,2].min(),0.);self.origins=start_offset+[0,0,height];self.targets=end_offset+[d,0,height]
  return np.stack([self.origins,self.targets],axis=1)
 def points(self):return np.array([self.rod.points(exp(self.w[i]),self.origins[i]) for i in range(self.N)])
 def endpoint_jacobian(self):
  R=exp(self.w);J=-self.rod.ds[None,:,None,None]*(R@skew(np.array([0.,0.,1.]))@right_jacobian(self.w));full=np.zeros((self.N*3,self.N*self.M*3))
  for i in range(self.N):full[3*i:3*i+3,i*self.M*3:(i+1)*self.M*3]=np.transpose(J[i],(1,0,2)).reshape(3,-1)/.0024
  full[~np.repeat(self.endpoint_fixed,3)]=0.
  return full
 def boundary_contact(self,p):
  if self.plane_bounds is None:
   return 0.,np.zeros_like(p),dict(top_force_N=0.,bottom_force_N=0.,top_penetration_m=0.,bottom_penetration_m=0.,boundary_mode="suspended"),np.zeros(p.shape[:2])
  bottom,top=self.plane_bounds
  return planes(p,self.radius,self.kc,top,bottom,reference_lengths=self.ref,current_measure=False)
 def evaluate(self,a):
  self.w[:]=np.asarray(a).reshape(self.w.shape);R=exp(self.w);p=self.points();Ec,gc,cs=pair_energy(p,self.radius,self.ref,self.kc,order=self.contact_order,self_contact=True);Ep,gp,ps,_=self.boundary_contact(p);E=Ec+Ep;g=np.zeros_like(self.w);endtorque=np.zeros((self.N,2,3));Eb=0.
  for i in range(self.N):
   e,gi=self.rods[i].energy_gradient(R[i],self.frames_for(i),self.boundary_k0[i]);E+=e;Eb+=e;gi+=self.rod.pullback_position_gradient(R[i],(gc+gp)[i]);g[i]=mv(np.swapaxes(right_jacobian(self.w[i]),-1,-2),gi);_,_,eg=self.rods[i].boundary_terms(R[i],self.frames_for(i),self.boundary_k0[i]);endtorque[i]=[self.ends[j]@eg[j] for j in range(2)]
  self.last=dict(energy_J=float(E),bending_energy_J=float(Eb),pair_contact_energy_J=Ec,backing_contact_energy_J=Ep,position_gradient_N=gc+gp,backing_gradient_N=gp,contact=cs,planes=ps,clamp_torques_world_Nm=endtorque,points=p)
  return E/self.unit,g.ravel()/self.unit
 def solve(self,d,alpha,wall_s=120,maxiter=1800):
  roots=self.controls(d,alpha);start=time.monotonic();iterations=[0];last_candidate=[self.w.ravel().copy()]
  def callback(a):
   iterations[0]+=1;last_candidate[0]=a.copy()
   if time.monotonic()-start>wall_s:raise TimeoutError('Bounded interacting-loop solve budget exceeded')
  def con(a):self.w[:]=a.reshape(self.w.shape);return ((self.points()[:,-1]-self.targets)/.0024).ravel()
  def jac(a):self.w[:]=a.reshape(self.w.shape);return self.endpoint_jacobian()
  original=self.w.ravel().copy();failure=None
  try:r=minimize(self.evaluate,original,jac=True,method='SLSQP',constraints=[dict(type='eq',fun=con,jac=jac)],callback=callback,options=dict(maxiter=maxiter,ftol=2e-11,disp=False));candidate=r.x;message=str(r.message);success=bool(r.success);nit=int(r.nit)
  except (TimeoutError,FloatingPointError,ValueError) as e:candidate=last_candidate[0];message=str(e);success=False;nit=iterations[0];failure=type(e).__name__
  return self.record(candidate,d,alpha,roots,start,success,message,nit,failure)
 def record(self,candidate,d,alpha,roots,start,success,message,nit,failure=None):
  E,g=self.evaluate(candidate);constraint=self.endpoint_constraint();J=self.endpoint_jacobian();lag=np.linalg.lstsq(J.T,-g,rcond=None)[0].reshape(self.N,3);kkt=g+J.T@lag.ravel();p=self.last['points'];ca=audit(p,self.radius,self.ref,self_contact=True);root_force=np.stack([self.last['position_gradient_N'].sum(axis=1)+lag*self.unit/.0024,-lag*self.unit/.0024],axis=1);mom=(np.cross(roots,root_force)+self.last['clamp_torques_world_Nm']).sum(axis=(0,1))-np.cross(p,self.last['backing_gradient_N']).sum(axis=(0,1));balance=root_force.sum(axis=(0,1))-self.last['backing_gradient_N'].sum(axis=(0,1));lengths=np.linalg.norm(np.diff(p,axis=1),axis=2);lengtherr=float(np.max(abs(lengths-self.ref)/self.ref));endcurv=max(self.rod.diagnostics(exp(self.w[i]))['max_radius_curvature'] for i in range(self.N));plane=max(self.last['planes']['bottom_penetration_m'],self.last['planes']['top_penetration_m']);residual=float(abs(kkt).max());passed=bool(success and max(abs(constraint))*.0024<1e-9 and residual<2e-5 and ca['max_penetration_m']<1e-6 and plane<1e-6 and lengtherr<1e-12 and endcurv<.1)
  row=dict(endpoint_fixed=self.endpoint_fixed.tolist(),manufacturing_parameters=self.manufacturing_parameters,applied_clamp_twist_rad=self.clamp_twist_rad,fibres=self.N,segments=self.M,root_distance_m=d,clamp_half_turn_rad=alpha,passed=passed,optimizer_success=success,optimizer_message=message,failure=failure,iterations=nit,wall_s=time.monotonic()-start,energy_J=self.last['energy_J'],bending_energy_J=self.last['bending_energy_J'],pair_contact_energy_J=self.last['pair_contact_energy_J'],backing_contact_energy_J=self.last['backing_contact_energy_J'],maximum_scaled_KKT_residual=residual,root_constraint_error_m=float(abs(constraint).max()*.0024),capsule_penetration_m=ca['max_penetration_m'],backing_penetration_m=plane,relative_length_error=lengtherr,maximum_radius_curvature=endcurv,total_mass_kg=float(self.N*np.pi*self.radius**2*.0024*self.rod.mat.density_kg_m3),total_reference_material_volume_m3=float(self.N*np.pi*self.radius**2*.0024),root_positions_m=roots.tolist(),root_forces_world_N=root_force.tolist(),clamp_torques_world_Nm=self.last['clamp_torques_world_Nm'].tolist(),force_balance_N=float(np.linalg.norm(balance)),moment_balance_Nm=float(np.linalg.norm(mom)),contact=self.last['contact'],max_height_m=float(p[:,:,2].max()),contact_order=self.contact_order,initial_imperfection_sign=self.seed_sign,render_promotion_allowed=False)
  self.history.append(row);self.states.append(self.w.copy());return row
