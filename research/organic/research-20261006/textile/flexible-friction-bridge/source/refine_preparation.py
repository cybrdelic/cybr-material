from pathlib import Path
import numpy as np
from scipy.spatial.transform import Slerp,Rotation
from prepare_three import ThreeCarriage,P
from prepare_three_spectral import prepare
from rod_pullback import exp,log
s=ThreeCarriage(128);z=np.load(P/'data/three_spectral_equilibrium.npz');oldR=exp(z['w']);oldsm=(np.arange(64)+.5)*s.L/64;sm=(np.arange(128)+.5)*s.L/128;frames=[]
for i in range(3):
 right=s.c.ends[1] if s.c.endpoint_fixed[i] else oldR[i,-1];rr=np.concatenate([s.c.ends[0][None],oldR[i],right[None]]);frames.append(Slerp(np.r_[0.,oldsm,s.L],Rotation.from_matrix(rr))(sm).as_matrix())
s.x=np.r_[log(np.array(frames)).ravel(),z['x'][-1]];s.theta=0.;prepare(s,s.x,90,prefix='three_refined128')
