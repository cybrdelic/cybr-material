"""Static curved partner: changed closest material label must not create spin."""
from pathlib import Path
import json,numpy as np
from flexible_contact import FlexibleContact
from rod_pullback import S,log,mv,skew
from test_rod_coupling import fixture
P=Path(__file__).resolve().parents[1];f=FlexibleContact(order=9);q0=fixture(f);old=f.state(q0);q=q0.copy();delta=np.array([0,0,20e-6]);q.reshape(f.n,-1)[0,:3]+=delta/S;t=f.trial(q,old);errors=[];crossings=0;spurious=[]
for a,b in zip(old['g']['sides'],t['state']['g']['sides']):
 ids=np.flatnonzero(b['N']>0)
 if not len(ids):continue
 crossings+=int(np.count_nonzero((a['segment'][ids]!=b['segment'][ids])&(a['N'][ids]>0)));spurious.append(float(np.linalg.norm(log(b['Rb'][ids]@np.swapaxes(a['Rb'][ids],-1,-2)),axis=1).max(initial=0)))
for s in t['surface']:
 side=f.pairs.index((s['i'],s['j']));a=old['g']['sides'][side];b=t['state']['g']['sides'][side];ids=np.flatnonzero(b['N']>0);motion=delta*(int(s['i']==0)-int(s['j']==0));n0=a['n'][ids];n1=b['n'][ids];cross=np.cross(n0,n1);c=np.sum(n0*n1,axis=1);K=skew(cross);Q=np.eye(3)+K+K@K/(1+c)[:,None,None];v1=motion-n1*(n1@motion)[:,None];v0=motion-n0*(n0@motion)[:,None];expected=.5*(v1+mv(Q,v0));errors.append(float(np.linalg.norm(expected-s['du'],axis=1).max(initial=0)))
r=dict(scope='Finite axial translation of a flexible curved rod against static curved partners. Analytic material translation with normal transport; this tests migration kinematics, not quasistatic equilibrium or piecewise-junction force convergence.',active_segment_crossings=crossings,max_surface_increment_error_m=max(errors),incorrect_closest_frame_rotation_rad=max(spurious),correct_static_material_rotation_rad=0.,passed=crossings>0 and max(errors)<1e-14);(P/'receipts/material_spin_migration.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
