"""Support-clipped quadrature of the same additive material-pair contact energy.
Weights are fixed material arclength; compact-support boundary terms vanish
because the potential is zero there. Returns integral forces, not the derivative
of finite-order quadrature-node placement; verify convergence independently.
"""
from pathlib import Path
import ctypes,numpy as np
P=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib=ctypes.CDLL(str(Path(__file__).with_name('_segment_contact_sine.so')));fun=lib.segment_pair_contact_sine;fun.argtypes=[ctypes.c_int,ctypes.c_int,P,P,P,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double,P,P];fun.restype=ctypes.c_double

def evaluate(positions,radii,reference_lengths,kc,order=16,self_contact=False,exclusion_radii=4.):
 p=np.ascontiguousarray(positions,dtype=np.float64);ref=np.ascontiguousarray(reference_lengths,dtype=np.float64)
 if p.ndim!=3 or p.shape[2]!=3 or not np.isfinite(p).all():raise ValueError('Expected finite (fibres,nodes,3) points.')
 N,K,_=p.shape;r=np.full(N,radii,dtype=np.float64) if np.isscalar(radii) else np.ascontiguousarray(radii,dtype=np.float64)
 if not 1<=N<=700 or not 2<=K<=1025 or N*(K-1)>50000 or r.shape!=(N,) or ref.shape!=(N,K-1) or not np.isfinite(r).all() or not np.isfinite(ref).all() or np.min(r)<1e-10 or np.min(ref)<=0 or np.any(np.linalg.norm(np.diff(p,axis=1),axis=2)<=1e-14):raise ValueError('Reference material/radius array contract failed.')
 if not np.isfinite([kc,exclusion_radii]).all() or kc<0 or not isinstance(order,int) or not 4<=order<=48 or exclusion_radii<0:raise ValueError('Invalid bounded coefficient/quadrature.')
 g=np.zeros_like(p);stats=np.zeros(7);E=fun(N,K,p,r,ref,float(kc),order,int(self_contact),float(exclusion_radii),g,stats)
 if not np.isfinite(E) or not np.isfinite(g).all():raise FloatingPointError('Singular state or segment-pair resource limit reached.')
 return float(E),g,dict(max_sample_overlap_m=float(stats[0]),active_segment_pairs=int(stats[1]),quadrature_point_pairs=int(stats[2]),segments=int(stats[3]),sum_pair_force_magnitudes_N=float(stats[4]),measure='fixed reference material arclength',quadrature_order=order)
