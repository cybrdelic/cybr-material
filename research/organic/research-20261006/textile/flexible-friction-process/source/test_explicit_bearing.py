import json,numpy as np
from spin_three import SpinThree,P,BRIDGE
from rod_pullback import exp
s=SpinThree();prep=np.load(BRIDGE/'data/three_refined128_equilibrium.npz');s.x=prep['x'];s.evaluate(s.x);s.history=s.bridge.state(s.bridge.pack(s.c.origins,exp(s.c.w)));z=np.load(P/'data/first_spin_terminal_centroid_revision.npz');x=z['y'][:len(s.x)];s.theta=float(z['theta']);_,g=s.evaluate(x);_,lam,J,S=s.projection(g);base=s.record(x,True,[]);pbase=s.c.points().copy();fric=s.last_trial['friction_gp'].copy();theta=s.theta;fd=[]
for h in [1e-3,1e-4,1e-5]:
 values=[]
 for sign in [-1,1]:
  s.theta=theta+sign*h;s.evaluate(x);E=s.c.last['energy_J']-s.last_trial['tangential_energy_J'];values.append(E+s.unit*(lam@s.constraints())+np.sum(fric*(s.c.points()-pbase)))
 numeric=(values[1]-values[0])/(2*h);fd.append(dict(h=h,derivative_Nm=float(numeric),reaction_Nm=base['torque_conjugate_Nm'],relative_error=abs(numeric/base['torque_conjugate_Nm']-1)))
s.theta=theta;s.evaluate(x);offset=s.c.offset.copy();flags=s.c.endpoint_fixed.copy();root=np.stack([s.c.origins,s.c.targets],axis=1).copy();perm=np.array([2,0,1]);s.c.offset=offset[perm];s.c.endpoint_fixed=flags[perm];s.c.controls(s.d,0);perr=float(abs(np.stack([s.c.origins,s.c.targets],axis=1)-root[perm]).max());s.c.offset=offset;s.c.endpoint_fixed=flags;s.c.controls(s.d,0);axis=s.c.bearing_axis_point.copy();translation=np.array([.0002,-.0001,.0003]);s.c.bearing_axis_point=axis+translation;s.c.fixed_axis_height=s.c.bearing_axis_point[2];s.c.controls(s.d,0);terr=float(abs(np.stack([s.c.origins,s.c.targets],axis=1)-root-translation).max());F=np.array(base['root_forces_N']);T=np.array(base['root_torques_Nm'])
def torque(points,force,couple,left):
 total=0.
 for end,sign in [(0,-.5),(1,.5)]:
  omega=np.array([sign,0.,0.]);centre=left+np.array([end*s.d,0.,0.]);total+=np.sum(force[:,end]*np.cross(omega,points[:,end]-centre))+np.sum(couple[:,end]*omega)
 return float(total)
t0=torque(root,F,T,axis);tp=torque(root[perm],F[perm],T[perm],axis);tt=torque(root+translation,F,T,axis+translation);r=dict(scope='Frozen uncommitted candidate: exact bearing motion and actuator virtual work. Friction forces are frozen in the finite-difference work functional; no derivative of friction stored energy is substituted for traction.',axis_left_m=axis.tolist(),permutation_position_error_m=perr,translation_position_error_m=terr,permutation_torque_error_Nm=tp-t0,translation_torque_error_Nm=tt-t0,finite_difference_work=fd,passed=perr<1e-18 and terr<1e-18 and abs(tp-t0)<1e-18 and abs(tt-t0)<1e-18 and min(a['relative_error'] for a in fd)<1e-5);(P/'receipts/explicit_bearing_audit.json').write_text(json.dumps(r,indent=2,default=lambda v:v.item()));print(json.dumps(r,indent=2,default=lambda v:v.item()))
