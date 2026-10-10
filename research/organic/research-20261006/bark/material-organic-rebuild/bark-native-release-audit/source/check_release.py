"""Exact-tie production native bulk release accounting at the frozen small load.
No bulk-only J claim is made for finite interface penalties.
"""
import os,resource
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
from pathlib import Path
import sys,json,time,hashlib,numpy as np
P=Path(__file__).resolve().parents[1];BROOT=P.parent/'bark-notch-reference/production-bridge'
sys.path.insert(0,str(BROOT/'source'))
from bridge import Bridge

def fields(model,state):
 u=model.lift(state['u']);H=model.H0+np.einsum('cia,cqib->cqab',u-u.mean(1,keepdims=True),model.bulk.grad);F=np.eye(3)+H;invT=np.linalg.inv(F).swapaxes(-2,-1)
 A=H+H.swapaxes(-2,-1)+H.swapaxes(-2,-1)@H
 term=A@A;diff=.5*np.trace(term,axis1=-2,axis2=-1)
 for k in range(3,25):
  term=term@A;diff+=((-1)**k/k)*np.trace(term,axis1=-2,axis2=-1)
 logJ=.5*(np.trace(A,axis1=-2,axis2=-1)-diff)
 W=.5*model.bulk.mu[:,None]*diff+.5*model.bulk.lam[:,None]*logJ**2
 left=H+H.swapaxes(-2,-1)+H@H.swapaxes(-2,-1)
 stress=model.bulk.mu[:,None,None,None]*(left@invT)+model.bulk.lam[:,None,None,None]*logJ[...,None,None]*invT
 g=np.einsum('cqib,cqab,cq->cia',model.bulk.grad,stress,model.bulk.weights);g-=g.mean(1,keepdims=True)
 native=state['full_nodal_bulk_gradient'];energy=float(np.sum(W*model.bulk.weights))
 checks=dict(field_energy_vs_native_relative=abs(energy-state['bulk_energy_half_J'])/abs(energy),field_gradient_vs_native_relative=float(np.linalg.norm(g-native)/np.linalg.norm(native)),interface_energy_fraction=abs(state['interface_energy_half_J'])/state['energy_half_J'])
 return H,W,stress,checks

def domain(model,state,radii):
 H,W,stress,checks=fields(model,state);rows=[]
 for ro in radii:
  ri=.5*ro;r=np.linalg.norm(model.X[:,:,:2]-[model.a,0.],axis=2);z=np.clip((r-ri)/(ro-ri),0,1);q=.5*(1+np.cos(np.pi*z))
  gradq=np.einsum('ci,cqib->cqb',q,model.bulk.grad)
  # Material-reference J: [P_ij * u_i,1 - W delta_1j] q_,j.
  flux=np.einsum('cqij,cqi->cqj',stress,H[:,:,:,0]);flux[:,:,0]-=W
  J=2*float(np.sum(np.einsum('cqi,cqi->cq',flux,gradq)*model.bulk.weights))/model.t
  rows.append(dict(outer_radius_m=ro,J_J_m2=J))
 return rows,checks

def run(n):
 cfg=json.loads((BROOT/'PROTOCOL.json').read_text());start=time.monotonic();m=Bridge(n,cfg);delta=cfg['half_grip_displacement_m'];base_a=m.a;states={};rows=[]
 for offset in (0,-1,1,-2,2):
  m.a=base_a+offset*m.h;m.prepare(None);st=m.solve(delta);states[offset]=st
  result=dict(offset_elements=offset,crack_m=m.a,energy_full_J=st['energy_full_J'],P_N=st['P_N'],compliance_m_N=st['compliance_m_N'],residual=st['free_residual_relative'])
  if offset==0:domains,parity=domain(m,st,(.006,.009,.012));result['domain_J']=domains;result['field_parity']=parity
  rows.append(result);print('NATIVE_RELEASE',n,offset,result['P_N'],time.monotonic()-start,flush=True)
  if time.monotonic()-start>115:raise TimeoutError('Declared115s bound exceeded')
 G=[]
 for step in (1,2):
  da=2*step*m.h;lo,hi=states[-step],states[step]
  gd=-(hi['energy_full_J']-lo['energy_full_J'])/(da*m.t);gf=states[0]['P_N']**2*(hi['compliance_m_N']-lo['compliance_m_N'])/(2*da*m.t)
  G.append(dict(step_elements=step,G_fixed_grip_J_m2=gd,G_compliance_proxy_J_m2=gf,scope='Compliance expression is the small-load linear proxy; native energy is used for the actual matched-grip release.'))
 j=np.array([r['J_J_m2'] for r in domains]);g=G[0]['G_fixed_grip_J_m2'];out=dict(mesh=n,prisms=m.N,rows=rows,release_derivatives=G,J_vs_G_max_relative=float(np.max(abs(j-g))/g),J_domain_spread=float(np.ptp(j)/j.mean()),release_step_sensitivity=abs(G[1]['G_fixed_grip_J_m2']-g)/g,wall_s=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,passed=bool(np.all(j>0) and g>0 and np.max(abs(j-g))/g<.03 and np.ptp(j)/j.mean()<.02 and abs(G[1]['G_fixed_grip_J_m2']-g)/g<.03 and parity['field_energy_vs_native_relative']<1e-10 and parity['field_gradient_vs_native_relative']<1e-10 and parity['interface_energy_fraction']<1e-15),scope='Small-load exact-tie constrained production response only; no finite-penalty interface J, damage, crack propagation or bark calibration.')
 (P/f'receipts/native_release_{n}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))

if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 16)
