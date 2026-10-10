import json,numpy as np
from spin_three import SpinThree,P
from rod_pullback import log
s=SpinThree();rows=[]
for theta in [0.,.027,.2]:
 s.c.clamp_twist_rad=theta;s.c.controls(s.d,0);p=np.stack([s.c.origins,s.c.targets],axis=1).copy();ends=np.array(s.c.ends);h=1e-6;parts=[];rots=[]
 for sign in [-1,1]:
  s.c.clamp_twist_rad=theta+sign*h;s.c.controls(s.d,0);parts.append(np.stack([s.c.origins,s.c.targets],axis=1).copy());rots.append(np.array(s.c.ends))
 numeric=(parts[1]-parts[0])/(2*h);expected=np.zeros_like(numeric);spins=[]
 for end,sign in [(0,-.5),(1,.5)]:
  centre=np.array([0. if end==0 else s.d,0.,s.c.fixed_axis_height]);expected[:,end]=np.cross([sign,0.,0.],p[:,end]-centre);spins.append(float(np.linalg.norm(ends[end]@log(rots[0][end].T@rots[1][end])/(2*h)-np.array([sign,0,0]))))
 rows.append(dict(theta=theta,position_velocity_error_m=float(abs(numeric-expected).max()),angular_velocity_error=max(spins)))
k=100.;du=1e-6;U=.5*k*du*du;F=k*du;Wt=.5*F*du;We=F*du;A=.5*k*du*du;r=dict(bearing_checks=rows,stick_work=dict(stored_J=U,trapezoid_work_J=Wt,endpoint_work_J=We,endpoint_quadrature_A_J=A,correct_trapezoid_gap_J=Wt-U,incorrect_double_count_gap_J=Wt-U-A,correct_endpoint_gap_J=We-U-A),passed=max(z['position_velocity_error_m'] for z in rows)<1e-12 and max(z['angular_velocity_error'] for z in rows)<1e-8);(P/'receipts/bearing_work_audit.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
