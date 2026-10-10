"""Independent analytic derivative of the exact rod geometric force map.
Frozen physical nodal forces and edge couples isolate geometric derivatives
from the numerical contact law tangent. Nonzero surface-spin couples included.
"""
from pathlib import Path
import json,numpy as np
from rod_pullback import RodCoordinates,skew,right_jacobian,mv
from test_rod_coupling import fixture
P=Path(__file__).resolve().parents[1]
def dJ(w,v):
 t2=np.sum(w*w,axis=-1);t=np.sqrt(np.maximum(t2,1e-30));dot=np.sum(w*v,axis=-1);small=t2<1e-4;A=np.where(small,.5-t2/24+t2*t2/720,(1-np.cos(t))/t2);B=np.where(small,1/6-t2/120+t2*t2/5040,(t-np.sin(t))/(t*t2));dA=np.where(small,(-1/12+t2/180-t2*t2/6720)*dot,(t*np.sin(t)-2*(1-np.cos(t)))*dot/t**4);dB=np.where(small,(-1/60+t2/1260-t2*t2/60480)*dot,(t*(1-np.cos(t))-3*(t-np.sin(t)))*dot/t**5);K=skew(w);dK=skew(v);return -dA[...,None,None]*K-A[...,None,None]*dK+dB[...,None,None]*(K@K)+B[...,None,None]*(dK@K+K@dK)
f=RodCoordinates();q=fixture(f);g=f.unpack(q);rng=np.random.default_rng(381);F=rng.normal(size=g['p'].shape)*1e-4;M=rng.normal(size=(f.n,f.M,3))*1e-8;rows=[]
for _ in range(3):
 v=rng.normal(size=q.shape);v/=np.linalg.norm(v);a=v.reshape(f.n,-1)[:,3:].reshape(f.n,f.M,3);J=right_jacobian(g['w']);body_spin=mv(J,a);out=np.zeros((f.n,3+3*f.M))
 for i in range(f.n):
  cum=np.cumsum(F[i,::-1],axis=0)[::-1][1:];bodyF=mv(np.swapaxes(g['R'][i],-1,-2),cum);bodyM=mv(np.swapaxes(g['R'][i],-1,-2),M[i]);b=f.rod.ds[:,None]*np.cross([0,0,1],bodyF)+bodyM;db=f.rod.ds[:,None]*np.cross([0,0,1],-np.cross(body_spin[i],bodyF))-np.cross(body_spin[i],bodyM);out[i,3:]=(mv(np.swapaxes(dJ(g['w'][i],a[i]),-1,-2),b)+mv(np.swapaxes(J[i],-1,-2),db)).ravel()
 exact=out.ravel();samples=[]
 for h in [1e-4,1e-5,1e-6]:
  numeric=(f.pullback(f.unpack(q+h*v),F,M)-f.pullback(f.unpack(q-h*v),F,M))/(2*h);samples.append(dict(h=h,error=float(np.linalg.norm(numeric-exact)/np.linalg.norm(exact))))
 rows.append(samples)
r=dict(scope='Analytic derivative of node/director geometric pullback under frozen physical forces and surface-spin couples. Independent of finite-difference contact Jacobian construction.',rows=rows,passed=all(min(a['error'] for a in row)<1e-7 for row in rows));(P/'receipts/pullback_analytic_derivative.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
