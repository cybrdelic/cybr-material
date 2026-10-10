import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from moving_contact import quadrature,geometry,RADIUS,KN
P=Path(__file__).resolve().parents[1];penetrations=np.r_[0.,np.geomspace(1e-10,.5e-6,80)];u,ww=np.polynomial.legendre.leggauss(200);u=(u+1)/2;ww=ww/2;refs=[]
for pen in penetrations:
 R=2*RADIUS-pen;a=np.sqrt(max(0.,(2*RADIUS)**2-R*R));dist=np.sqrt(R*R+(a*u)**2);delta=np.maximum(2*RADIUS-dist,0.);refs.append([2*KN*a*np.sum(ww*delta*R/dist),KN*a*np.sum(ww*delta**2)])
refs=np.array(refs);peak=refs.max(axis=0);orders=[]
for q in [5,9,17,33,65]:
 s,w,_=quadrature(q);values=[]
 for pen in penetrations:
  t=np.arccos(np.clip(-pen/(.5e-6),-1,1))/(2*np.pi);g=geometry(s,w,t,'birth_loss');R=2*RADIUS-pen;values.append([float(np.sum(g['N']*R/np.linalg.norm(g['pa']-g['pb'],axis=1))),g['Un']])
 values=np.array(values);absolute=np.max(abs(values-refs),axis=0)/peak;relative=[]
 for j in range(2):
  mask=refs[:,j]>=.01*peak[j];relative.append(float(np.max(abs(values[mask,j]-refs[mask,j])/refs[mask,j])))
 orders.append(dict(order=q,peak_normalized_maximum_force_error=float(absolute[0]),peak_normalized_maximum_energy_error=float(absolute[1]),maximum_relative_force_error_above_one_percent_peak=relative[0],maximum_relative_energy_error_above_one_percent_peak=relative[1],passed=bool(max(*absolute,*relative)<.02)))
out=dict(declared_tolerance=.02,near_zero_reporting='Also require2% pointwise error wherever each observable exceeds1% of its own maximum; below that, retain2% peak-normalized absolute gate.',reference='200-point Gauss integration over exact analytic support of the declared crossed straight-fibre normal energy.',penetration_range_m=[float(penetrations[0]),float(penetrations[-1])],load_samples=len(penetrations),orders=orders);(P/'receipts/normal_sweep.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
