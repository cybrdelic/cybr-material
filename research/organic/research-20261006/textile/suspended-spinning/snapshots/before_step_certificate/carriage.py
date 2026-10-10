"""One explicit tension-controlled axial carriage, fixed material stock.
The roots rotate about fixed bearing axes. d is a solved coordinate; geometry
is integrated from unit directors and is never rescaled to make it fit.
"""
from pathlib import Path
import sys,time,numpy as np
from scipy.sparse import bmat,block_diag,coo_matrix,diags,eye,kron
from scipy.sparse.linalg import splu,LinearOperator,eigsh
P=Path(__file__).resolve().parents[1];BASE=P.parent/'short-yarn-scaling';sys.path.insert(0,str(P/'source'))
from loop_bundle import LoopBundle,exp,log,skew,right_jacobian,mv
from operators import RodHessianOperator
from rod_bands import sparse_bands
from contact_blocks import UPPER
from contact_path import certify_path
from line_contact import audit

class MovingRoots(LoopBundle):
 def controls(self,d,alpha):
  roots=super().controls(d,alpha)
  if hasattr(self,'fixed_axis_height'):
   shift=self.fixed_axis_height-self.origins[0,2];self.origins[:,2]+=shift;self.targets[:,2]+=shift
  return np.stack([self.origins,self.targets],axis=1)

class Carriage:
 def __init__(self,plane_bounds=None):
  import json
  self.c=MovingRoots(19,128,plane_bounds=plane_bounds);self.L=.0024;self.d=.00238;self.theta=0.;self.c.controls(self.d,0);self.c.fixed_axis_height=float(self.c.origins[0,2]);z=np.load(BASE/'data/yarn_19_qualified.npz');self.c.w[:]=z['w'];self.n=self.c.w.size;h=json.loads((BASE/'receipts/yarn_19_qualified.json').read_text());self.tension=float(np.array(h['root_forces_world_N'])[:,1,0].sum());self.x=np.r_[self.c.w.ravel(),self.d/self.L];self.fixed=np.repeat(self.c.endpoint_fixed,3);self.hd=np.zeros(3*self.c.N);self.hd[::3]=-self.c.endpoint_fixed.astype(float);self.unit=self.c.unit
 def evaluate(self,x):
  x=np.asarray(x);self.d=float(x[-1]*self.L)
  if not .001<self.d<self.L:raise ValueError('Carriage outside declared short-stock domain')
  self.c.clamp_twist_rad=self.theta;self.c.controls(self.d,0);E,g=self.c.evaluate(x[:-1]);return E-self.tension*self.d/self.unit,np.r_[g,-self.tension*self.L/self.unit]
 def constraints(self):return self.c.endpoint_constraint()
 def jacobian(self):return np.column_stack([self.c.endpoint_jacobian(),self.hd])
 def projection(self,g):
  J=self.jacobian();S=J@J.T;S[~self.fixed,~self.fixed]=1.;lam=np.linalg.solve(S,-J@g);return g+J.T@lam,lam,J,S
 def predictor(self,delta_theta):
  c=self.c;R=exp(c.w);p=c.points();sm=(c.rod.s[:-1]+c.rod.s[1:])/2;phi=(sm/self.L-.5)*delta_theta;Rx=exp(np.c_[phi,np.zeros_like(phi),np.zeros_like(phi)]);relative=(p[:,:-1]+p[:,1:])/2;relative[:,:,1]-=0.;relative[:,:,2]-=c.fixed_axis_height;radial=mv(Rx[None],relative);t=mv(Rx[None],R[:,:,:,2])+delta_theta/self.L*np.cross(np.array([1.,0,0]),radial);t/=np.linalg.norm(t,axis=2)[:,:,None];a=mv(Rx[None],R[:,:,:,0]);a-=t*np.sum(t*a,axis=2)[:,:,None];a/=np.linalg.norm(a,axis=2)[:,:,None];RR=np.stack([a,np.cross(t,a),t],axis=-1);return np.r_[log(RR).ravel(),self.x[-1]]
 def predictor_boundary(self,delta_theta):
  # A local extension of changed root boundary data is only an initial guess.
  # Its unit tangents retain stock; the physical equilibrium decides the state.
  c=self.c;R=exp(c.w);p=c.points();origins=c.origins.copy();targets=c.targets.copy();old=c.clamp_twist_rad;c.clamp_twist_rad=self.theta+delta_theta;c.controls(self.d,0);left=c.origins-origins;right=(c.targets-targets)*c.endpoint_fixed[:,None];c.clamp_twist_rad=old;c.controls(self.d,0)
  arc=c.rod.s;span=.2*self.L
  def fade(u):
   u=np.clip(u,0,1);return 1-3*u*u+2*u*u*u
  wl=fade(arc/span);wr=fade((self.L-arc)/span);guess=p+left[:,None,:]*wl[None,:,None]+right[:,None,:]*wr[None,:,None];t=np.diff(guess,axis=1);t/=np.linalg.norm(t,axis=2)[:,:,None];sm=(arc[:-1]+arc[1:])/2;phi=delta_theta*(-.5*fade(sm/span)[None,:]+.5*c.endpoint_fixed[:,None]*fade((self.L-sm)/span)[None,:]);Rx=exp(np.stack([phi,np.zeros_like(phi),np.zeros_like(phi)],axis=-1));a=mv(Rx,R[:,:,:,0]);a-=t*np.sum(t*a,axis=2)[:,:,None];a/=np.linalg.norm(a,axis=2)[:,:,None];RR=np.stack([a,np.cross(t,a),t],axis=-1);return np.r_[log(RR).ravel(),self.x[-1]]
 def record(self,x,passed,trace,stability=None):
  potential,g=self.evaluate(x);pg,lam,J,S=self.projection(g);c=self.c;p=c.points();length=np.linalg.norm(np.diff(p,axis=1),axis=2);ca=audit(p,c.radius,c.ref,self_contact=True);root_force=np.stack([c.last['position_gradient_N'].sum(axis=1)+lam.reshape(c.N,3)*self.unit/self.L,-lam.reshape(c.N,3)*self.unit/self.L],axis=1);ends=np.array(c.ends);roots=np.stack([c.origins,c.targets],axis=1);torques=c.last['clamp_torques_world_Nm'];M=0.
  for end,sign in [(0,-.5),(1,.5)]:
   omega=np.array([sign,0,0]);cent=np.array([0 if end==0 else self.d,0,c.fixed_axis_height]);velocity=np.cross(omega,roots[:,end]-cent);M+=float(np.sum(root_force[:,end]*velocity)+np.sum(torques[:,end]*omega))
  self.x=x.copy();return dict(passed=bool(passed),force_residual=float(abs(pg).max()),root_error_m=float(abs(self.constraints()).max()*self.L),stored_energy_J=c.last['energy_J'],total_potential_J=potential*self.unit,carriage_m=self.d,tension_N=self.tension,actual_axial_reaction_N=float(root_force[:,1,0].sum()),twist_rad=self.theta,twist_turns_per_m=self.theta/(2*np.pi*self.d),dimensionless_root_twist=self.theta*np.max(np.linalg.norm(c.offset[:,:2],axis=1))/self.d,nominal_helix_angle_degrees=float(np.degrees(np.arctan(self.theta*np.max(np.linalg.norm(c.offset[:,:2],axis=1))/self.d))),torque_conjugate_Nm=M,total_mass_kg=float(c.N*np.pi*c.radius**2*self.L*c.rod.mat.density_kg_m3),effective_tex=float(c.N*np.pi*c.radius**2*self.L*c.rod.mat.density_kg_m3/self.d*1e6),relative_stock_error=float(abs(length-c.ref).max()/c.ref.min()),capsule_penetration_m=ca['max_penetration_m'],stability=stability,trace=trace,boundary_planes_m=c.plane_bounds,plane_contact=c.last['planes'],scope='Normal-only quasistatic tension/twist process with fixed stock. No finite sliding-friction history or wool calibration.',root_forces_N=root_force.tolist(),root_torques_Nm=torques.tolist())

def matrix(op,system,shift=0.,positive=False):
 c=op.c;N,M,K=c.N,c.M,c.M+1;n=op.size;C=op.contact;ai=C.nodes[:,UPPER[:,0]];bi=C.nodes[:,UPPER[:,1]];BB=C.gn if positive else C.full+shift*C.gn;rows=np.broadcast_to(ai[:,:,None,None]*3+np.arange(3)[None,None,:,None],BB.shape).ravel();cols=np.broadcast_to(bi[:,:,None,None]*3+np.arange(3)[None,None,None,:],BB.shape).ravel();off=UPPER[:,0]!=UPPER[:,1];Bo=BB[:,off];r2=np.broadcast_to(bi[:,off,None,None]*3+np.arange(3)[None,None,None,:],Bo.shape).ravel();c2=np.broadcast_to(ai[:,off,None,None]*3+np.arange(3)[None,None,:,None],Bo.shape).ravel();H=coo_matrix((np.r_[BB.ravel(),Bo.ravel()],(np.r_[rows,r2],np.r_[cols,c2])),shape=(N*K*3,N*K*3)).tocsr()+(1+shift)*op.surface;keep=np.ones(N*K*3,dtype=bool);keep[np.arange(N)[:,None]*K*3+np.arange(3)[None,:]]=False;Hp=H[keep][:,keep]*(op.L**2/op.unit);Hq=block_diag(op.positive_rods,format='csc') if positive else block_diag([sparse_bands(op.diag[i]+op.geom[i]+shift*op.gndiag[i],op.off[i]+shift*op.gnoff[i])/op.unit for i in range(N)],format='csc');D=block_diag(list(op.D.reshape(-1,3,3)/op.L),format='csc');diff=kron(eye(N),kron(diags([np.ones(M),-np.ones(M-1)],[0,-1],shape=(M,M)),eye(3)),format='csc');ids=np.array([(i*M+M-1)*3+a for i in range(N) if c.endpoint_fixed[i] for a in range(3)]);Jp=coo_matrix((np.ones(len(ids)),(np.arange(len(ids)),ids)),shape=(len(ids),n)).tocsc();hd=coo_matrix(system.hd[system.fixed,None]);zero=coo_matrix((1,1));A=bmat([[Hq,None,None,-D.T,None],[None,Hp,None,diff.T,Jp.T],[None,None,zero,None,hd.T],[-D,diff,None,None,None],[None,Jp,hd,None,None]],format='csc');return A,splu(A,permc_spec='COLAMD')

def direction(op,system,g,con,shift):
 A,factor=matrix(op,system,shift);n=op.size;rhs=np.r_[-g[:-1],np.zeros(n),-g[-1],np.zeros(n),-con[system.fixed]];sol=factor.solve(rhs)
 for _ in range(2):sol+=factor.solve(rhs-A@sol)
 dx=np.r_[sol[:n],sol[2*n]];lam=np.zeros(3*op.N);lam[system.fixed]=sol[3*n+1:];J=system.jacobian();physical=np.r_[op.matvec(dx[:-1])+shift*op.matvec(dx[:-1],True),0.];res=np.r_[physical+J.T@lam+g,J@dx+con];return dx,dict(relative_linear_residual=float(np.linalg.norm(res)/max(np.linalg.norm(np.r_[g,con]),1e-30)),factor_nonzeros=factor.L.nnz+factor.U.nnz,max_abs_multiplier=float(abs(lam).max()))

def stability(op,system,J,S):
 n=op.size;size=n+1
 def project(v):return v-J.T@np.linalg.solve(S,J@v)
 def apply(v,positive=False):
  p=project(v);return project(np.r_[op.matvec(p[:-1],positive),0.])+v-p
 Apos,factor=matrix(op,system,positive=True)
 def inverse(v):
  p=project(v);rhs=np.r_[p[:-1],np.zeros(n),p[-1],np.zeros(n),np.zeros(system.fixed.sum())];sol=factor.solve(rhs)
  for _ in range(2):sol+=factor.solve(rhs-Apos@sol)
  return np.r_[sol[:n],sol[2*n]]+v-p
 A=LinearOperator((size,size),matvec=apply,dtype=float);B=LinearOperator((size,size),matvec=lambda v:apply(v,True),dtype=float);BI=LinearOperator((size,size),matvec=inverse,dtype=float);rng=np.random.default_rng(76);test=rng.normal(size=size);ie=float(np.linalg.norm(B@inverse(test)-test)/np.linalg.norm(test));ev,v=eigsh(A,M=B,Minv=BI,k=2,which='SA',tol=1e-8,maxiter=600,ncv=20,v0=rng.normal(size=size));rr=[float(np.linalg.norm(A@v[:,i]-ev[i]*(B@v[:,i]))/np.linalg.norm(B@v[:,i])) for i in range(2)];return dict(minimum_generalized_eigenvalue=float(ev[0]),eigenvalues=ev.tolist(),eigen_residuals=rr,inverse_error=ie,passed=bool(ev[0]>1e-6 and max(rr)<1e-6 and ie<1e-7))

def solve(system,target_theta,wall_s=90,checkpoint=None):
 start=time.monotonic();initial=system.x.copy();system.evaluate(initial);oldtheta=system.theta;oldpoints=system.c.points().copy();x=system.predictor_boundary(target_theta-oldtheta);system.theta=target_theta;E,g=system.evaluate(x);predictor_cert=certify_path(oldpoints,system.c.points(),system.c.radius,system.c.ref,plane_bounds=system.c.plane_bounds)
 if not predictor_cert['passed']:system.theta=oldtheta;system.evaluate(initial);return system.record(initial,False,[dict(predictor_contact_failure=predictor_cert)])
 shift=.05;trace=[];stable=None;success=False
 for it in range(150):
  if time.monotonic()-start>wall_s:break
  E,g=system.evaluate(x);pg,lam,J,S=system.projection(g);con=system.constraints();op=RodHessianOperator(system.c,lam);res=float(abs(pg).max());row=dict(iteration=it,force_residual=res,carriage_m=system.d,theta=system.theta,stored_J=system.c.last['energy_J']);trace.append(row);print('CARRIAGE',it,'res',res,'d',system.d,'seconds',time.monotonic()-start,flush=True)
  if res<2e-5 and abs(con).max()*system.L<1e-9:
   stable=stability(op,system,J,S);success=stable['passed'];break
  p0=system.c.points().copy()
  for attempt in range(12):
   dx,linear=direction(op,system,g,con,shift);rho=max(1.,1.5*max(float(abs(lam).max()),linear['max_abs_multiplier'])+1);pred=float(pg@(dx-J.T@np.linalg.solve(S,J@dx))+lam@con-rho*np.linalg.norm(con,1))
   if pred<0 and linear['relative_linear_residual']<1e-7:break
   shift=max(.001,shift*2)
  else:row['failed_direction']=linear;break
  angle=np.linalg.norm(dx[:-1].reshape(system.c.w.shape),axis=2).max();scale=min(1.,.1/max(angle,1e-30));accepted=False
  for ls in range(24):
   trial=x+scale*dx
   try:
    for _ in range(5):
     system.evaluate(trial);ct=system.constraints();remainder=ct-(1-scale)*con
     if abs(remainder).max()<1e-12:break
     _,_,Jt,St=system.projection(np.zeros_like(g));trial-=Jt.T@np.linalg.solve(St,remainder)
    et,gt=system.evaluate(trial)
    if et+rho*np.linalg.norm(system.constraints(),1)<=E+rho*np.linalg.norm(con,1)+1e-4*scale*pred+1e-13:
     cert=certify_path(p0,system.c.points(),system.c.radius,system.c.ref,plane_bounds=system.c.plane_bounds)
     if cert['passed']:accepted=True;break
   except (ValueError,FloatingPointError,np.linalg.LinAlgError):pass
   scale*=.5
  row.update(metric_shift=shift,linear=linear,step_scale=scale,line_search_steps=ls+1)
  if not accepted:row['failure']='No energy-decreasing feasible step';break
  x=trial;x[:-1]=log(exp(x[:-1].reshape(system.c.w.shape))).ravel();shift=shift*2 if ls>2 else (shift*.5 if shift>1e-5 else 0.)
  if checkpoint is not None:checkpoint(x,system.theta,trace)
 r=system.record(x,success,trace,stable);r.update(wall_s=time.monotonic()-start,predictor_contact= predictor_cert,endpoint_path_contact=certify_path(oldpoints,system.c.points(),system.c.radius,system.c.ref,plane_bounds=system.c.plane_bounds));return r
