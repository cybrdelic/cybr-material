"""Validated-shape Python boundary for the SI contact integral kernel.
U_ij = 1/2*(integral_i psi(distance(x_i,curve_j)) ds_i + integral_j ...),
psi(d)=kc/2*max(ri+rj-d,0)^2. No inflated section radii enter the law.
"""
from pathlib import Path
import ctypes
import numpy as np
P=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS');lib=ctypes.CDLL(str(Path(__file__).with_name('_line_contact_native.so')))
f=lib.line_contact;f.argtypes=[ctypes.c_int,ctypes.c_int,P,P,P,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_int,P,P];f.restype=ctypes.c_double
pl=lib.plane_contact;pl.argtypes=[ctypes.c_int,ctypes.c_int,P,P,P,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,P,P,P];pl.restype=ctypes.c_double

def checked(positions,radii,reference_lengths,current_measure,self_contact):
 p=np.ascontiguousarray(positions,dtype=np.float64)
 if p.ndim!=3 or p.shape[2]!=3:raise ValueError('positions must have exact shape (curves,nodes,3).')
 N,K,_=p.shape
 if not (1<=N<=700 and 2<=K<=1025) or not np.isfinite(p).all():raise ValueError('positions outside bounded finite shape contract.')
 r=np.full(N,radii,dtype=np.float64) if np.isscalar(radii) else np.ascontiguousarray(radii,dtype=np.float64)
 if r.shape!=(N,) or not np.isfinite(r).all() or np.any(r<=0):raise ValueError('One positive finite radius per curve is required.')
 lengths=np.linalg.norm(np.diff(p,axis=1),axis=2)
 if np.any(lengths<=1e-14):raise ValueError('Zero/degenerate segment refused.')
 if reference_lengths is None:
  if self_contact or not current_measure:raise ValueError('Fixed material reference lengths are required for this mode.')
  ref=np.ascontiguousarray(lengths)
 else:ref=np.ascontiguousarray(reference_lengths,dtype=np.float64)
 if ref.shape!=(N,K-1) or not np.isfinite(ref).all() or np.any(ref<=0):raise ValueError('reference_lengths must match (curves,nodes-1), positive and finite.')
 assert p.shape==(N,K,3) and r.shape==(N,) and ref.shape==(N,K-1)
 return p,r,ref

def contact(positions,radii,kc,reference_lengths=None,current_measure=True,self_contact=False,exclusion_radii=4.,order=3,grid=2):
 p,r,ref=checked(positions,radii,reference_lengths,current_measure,self_contact);N,K,_=p.shape
 if not np.isfinite(kc) or kc<0 or order not in [1,2,3,4] or grid not in [0,1,2]:raise ValueError('Invalid contact coefficient/quadrature/broad-phase mode.')
 if not np.isfinite(exclusion_radii) or exclusion_radii<0:raise ValueError('Invalid nonlocal exclusion.')
 g=np.zeros_like(p);s=np.zeros(7);e=f(N,K,p,r,ref,float(kc),order,int(current_measure),int(self_contact),float(exclusion_radii),grid,g,s)
 if not np.isfinite(e) or not np.isfinite(g).all():raise FloatingPointError('Nonfinite native contact result.')
 return float(e),g,dict(max_quadrature_penetration_m=float(s[1]),active_samples=int(s[2]),undefined_normal_samples=int(s[3]),grid_used=bool(s[4]),grid_entries=int(s[5]),integrated_normal_force_N=float(s[6]))

def planes(positions,radii,kc,top,bottom=0.,reference_lengths=None,current_measure=True):
 p,r,ref=checked(positions,radii,reference_lengths,current_measure,False);N,K,_=p.shape
 if not np.isfinite([kc,top,bottom]).all() or kc<0 or top<=bottom:raise ValueError('Invalid plane coefficient or bounds.')
 g=np.zeros_like(p);s=np.zeros(5);node=np.zeros((N,K));e=pl(N,K,p,r,ref,float(kc),float(top),float(bottom),int(current_measure),g,s,node)
 if not np.isfinite(e) or not np.isfinite(g).all():raise FloatingPointError('Nonfinite plane result.')
 return float(e),g,dict(top_force_N=float(s[0]),bottom_force_N=float(s[1]),top_penetration_m=float(s[2]),bottom_penetration_m=float(s[3])),node

au=lib.capsule_audit;au.argtypes=[ctypes.c_int,ctypes.c_int,P,P,P,ctypes.c_int,ctypes.c_double,P]
def audit(positions,radii,reference_lengths=None,self_contact=False,exclusion_radii=4.):
 p,r,ref=checked(positions,radii,reference_lengths,True,self_contact);N,K,_=p.shape;s=np.zeros(3)
 au(N,K,p,r,ref,int(self_contact),float(exclusion_radii),s)
 return dict(max_penetration_m=float(s[0]),minimum_tested_distance_m=(float(s[1]) if s[1]<1e90 else None),penetrating_segment_pairs=int(s[2]))

fm=lib.line_contact_modal;fm.argtypes=f.argtypes+[ctypes.c_int,P,P];fm.restype=ctypes.c_double
pm=lib.plane_contact_modal;pm.argtypes=pl.argtypes+[ctypes.c_int,P,P];pm.restype=ctypes.c_double

def modal_system(positions,radii,kc,top,phi,reference_lengths=None,bottom=0.,order=3):
 """Exact energy/position gradient plus a positive Gauss–Newton modal approximation.
 The approximate Hessian is a numerical search metric, not an altered energy law.
 Modes move x and z only, with nodal basis phi. Shape contracts precede native calls.
 """
 p,r,ref=checked(positions,radii,reference_lengths,True,False);N,K,_=p.shape;phi=np.ascontiguousarray(phi,dtype=np.float64)
 if phi.ndim!=2 or phi.shape[0]!=K or not 1<=phi.shape[1]<=4 or not np.isfinite(phi).all():raise ValueError('Invalid modal basis.')
 M=phi.shape[1];d=N*M*2
 if d>2048 or not np.isfinite([kc,top,bottom]).all() or kc<0 or top<=bottom or order not in (1,2,3,4):raise ValueError('Modal solver shape/parameter admission refused.')
 g=np.zeros_like(p);gp=np.zeros_like(p);H=np.zeros((d,d));Hp=np.zeros_like(H);s=np.zeros(7);sp=np.zeros(5);node=np.zeros((N,K))
 Ec=fm(N,K,p,r,ref,float(kc),order,1,0,4.,2,g,s,M,phi,H)
 Ep=pm(N,K,p,r,ref,float(kc),float(top),float(bottom),1,gp,sp,node,M,phi,Hp)
 if not np.isfinite(Ec+Ep) or not np.isfinite(H).all() or not np.isfinite(Hp).all() or not np.isfinite(g+gp).all():raise FloatingPointError('Invalid modal native result.')
 cs=dict(max_quadrature_penetration_m=float(s[1]),active_samples=int(s[2]),undefined_normal_samples=int(s[3]),integrated_normal_force_N=float(s[6]))
 ps=dict(top_force_N=float(sp[0]),bottom_force_N=float(sp[1]),top_penetration_m=float(sp[2]),bottom_penetration_m=float(sp[3]))
 return float(Ec+Ep),g+gp,H+Hp,cs,ps,node
