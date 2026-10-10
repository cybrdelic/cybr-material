import json,time,resource,numpy as np
from mixed_friction import P,BRIDGE,mixed_metric
from prepare_three import ThreeCarriage
from rod_pullback import exp,L
start=time.monotonic();s=ThreeCarriage(128);z=np.load(BRIDGE/'data/three_refined128_equilibrium.npz');s.evaluate(z['x']);f=s.bridge;q=f.pack(s.c.origins,exp(s.c.w));old=f.state(q);t=f.trial(q,old);H,B,stats=mixed_metric(f,t);g=t['state']['g'];rng=np.random.default_rng(26);errors=[];energies=[]
for _ in range(3):
 v=rng.normal(size=(f.n,f.M,3));v/=np.linalg.norm(v);qv=np.c_[np.zeros((f.n,3)),v.reshape(f.n,-1)].ravel();node_velocity=np.array([f.velocity(g,qv,i,f.rod.s[1:])[0] for i in range(f.n)]);mixed=np.r_[v.ravel(),node_velocity.ravel()/L];expected=[]
 for a in t['surface']:
  vi,oi=f.velocity(g,qv,a['i'],a['s']);vj,oj=f.velocity(g,qv,a['j'],a['sb']);expected.append((vi+np.cross(oi,a['armA'])-vj-np.cross(oj,a['armB'])).ravel())
 expected=np.concatenate(expected);errors.append(float(np.linalg.norm(B@mixed-expected)/max(np.linalg.norm(expected),1e-30)));energies.append(float(mixed@(H@mixed)))
asym=H-H.T;ratio=float(np.linalg.norm(asym.data)/max(np.linalg.norm(H.data),1e-30));checks=dict(velocity=max(errors)<1e-8,positive=min(energies)>=-1e-12*max(energies),symmetric=ratio<1e-12);r=dict(scope='Mixed-coordinate positive friction preconditioner only; full Coulomb residual and nonsymmetric tangent still required for the actual solve.',checks={k:bool(v) for k,v in checks.items()},passed=bool(all(checks.values())),velocity_relative=errors,quadratic_values=energies,symmetry_relative=ratio,stats=stats,wall_s=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024);(P/'receipts/mixed_metric.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
