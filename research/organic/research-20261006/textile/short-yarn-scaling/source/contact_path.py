"""Conservative continuous contact gate for a linear polygonal homotopy.
This is a numerical continuation safeguard, not an inertial forming law. The
homotopy retains material node correspondence and exact fixed roots; intermediate
polygons need not satisfy inextensibility. Its existence prevents axis crossing
between accepted equilibria. A Lipschitz bound certifies every time interval,
including between audit samples; uncertain intervals are rejected, never passed.
"""
import numpy as np
from line_contact import audit

def certify_path(p0,p1,radius,ref,tolerance_m=1e-6,max_audits=4097,bottom=0.,top=.02):
 p0=np.asarray(p0,dtype=float);p1=np.asarray(p1,dtype=float)
 if p0.shape!=p1.shape or p0.ndim!=3 or p0.shape[2]!=3 or not np.isfinite([*p0.ravel(),*p1.ravel(),radius,tolerance_m]).all() or radius<=0 or not 0<tolerance_m<radius:raise ValueError('Invalid continuation geometry or penetration tolerance.')
 # One common translation cancels out of every pair distance, tightening the bound.
 displacement=p1-p0;relative=displacement-displacement.mean(axis=(0,1));speed=float(np.linalg.norm(relative,axis=-1).max());surface=max(0.,radius-float(min(p0[:,:,2].min(),p1[:,:,2].min()))+bottom,float(max(p0[:,:,2].max(),p1[:,:,2].max()))+radius-top);calls=0;peak=0.;certified=0.;deepest=0;stack=[(0.,1.,0)]
 if surface>tolerance_m:return dict(passed=False,reason='Plane penetration along linear path',audits=0,peak_sample_penetration_m=surface)
 while stack:
  lo,hi,depth=stack.pop();deepest=max(deepest,depth)
  if calls>=max_audits:return dict(passed=False,reason='Continuous-contact certificate budget exhausted',audits=calls,peak_sample_penetration_m=peak)
  t=(lo+hi)*.5;p=(1-t)*p0+t*p1
  try:pen=audit(p,radius,ref,self_contact=True)['max_penetration_m']
  except ValueError:return dict(passed=False,reason='Degenerate intermediate material segment',audits=calls,peak_sample_penetration_m=peak)
  calls+=1;peak=max(peak,pen)
  if pen>tolerance_m:return dict(passed=False,reason='Intermediate capsule penetration exceeds threshold',audits=calls,at_parameter=t,peak_sample_penetration_m=peak)
  # Each segment is a convex combination of its endpoint velocities. Pair
  # distance is 2*speed-Lipschitz; any t is at most (hi-lo)/2 from this sample.
  bound=pen+speed*(hi-lo)
  if bound<=tolerance_m:certified=max(certified,bound);continue
  if depth>=24:return dict(passed=False,reason='Continuous-contact recursion exhausted',audits=calls,peak_sample_penetration_m=peak)
  stack.extend([(lo,t,depth+1),(t,hi,depth+1)])
 return dict(passed=True,reason='All continuous-time intervals bounded',audits=calls,maximum_proved_penetration_bound_m=certified,peak_sample_penetration_m=peak,plane_penetration_m=surface,maximum_relative_node_motion_m=speed,deepest_subdivision=deepest)
