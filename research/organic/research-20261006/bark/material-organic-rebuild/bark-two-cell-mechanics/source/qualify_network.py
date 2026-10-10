import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib
import numpy as np
from two_cell_network import TwoCellNetwork
from contact_gate import check_contact,segment_distances,outside_disk

P=Path(__file__).resolve().parents[1]
checks=[]
for a,b,c,d,target in [([0,0],[1,0],[.5,-1],[.5,1],0),
    ([0,0],[1,0],[0,2],[1,2],2),([0,0],[1,0],[.5,0],[2,0],0),
    ([0,0],[1,0],[2,0],[3,0],1)]:
    actual=float(segment_distances(*[np.array(v,dtype=float) for v in (a,b,c,d)]))
    checks.append(abs(actual-target)<1e-12)
pieces=outside_disk(np.array([[0.,0.],[2.,0.]]),np.zeros(2),1.)
checks.append(len(pieces)==1 and np.allclose(pieces[0][0],[1.,0.]))
rows=[]
for ns,steps in [(24,10),(48,10),(96,10),(48,20),(48,40)]:
    name=f'network_{ns}_{steps}_5pct'
    data=np.load(P/f'data/{name}.npz')
    old=json.loads((P/f'receipts/{name}.json').read_text())
    n=TwoCellNetwork(ns,ns//3)
    states=[]
    for q,load in zip(data['states'],old['states']):
        e,g,H,_=n.evaluate(q,True);f=n.free
        pos=q[:2*n.nn].reshape(-1,2);support=np.zeros_like(g);support[n.fixed]=g[n.fixed]
        forces=support[:2*n.nn].reshape(-1,2)*n.unit/n.L
        force=float(np.linalg.norm(np.sum(forces,axis=0)))
        moment=float(abs(np.sum((pos[:,0]*forces[:,1]-pos[:,1]*forces[:,0])*n.L)))
        free_res=float(np.max(abs(g[f])));phi_res=float(np.max(abs(g[2*n.nn:]))*n.unit)
        eig=float(np.linalg.eigvalsh(H[np.ix_(f,f)])[0])
        contact=check_contact(n,q)
        control=max(float(abs(pos[0,0])),float(np.max(abs(pos[:3,1]))),
                    float(np.max(abs(pos[3:6,1]-(1+load['engineering_strain'])))))
        passed=bool(free_res<2e-8 and eig>0 and force<1e-13 and moment<1e-18
                    and phi_res<1e-18 and control<1e-12 and contact['passed'])
        states.append(dict(engineering_strain=load['engineering_strain'],energy_J=e*n.unit,
            residual=free_res,min_free_eigenvalue=eig,support_force_balance_N=force,
            support_moment_balance_Nm=moment,max_junction_moment_residual_Nm=phi_res,
            contact=contact,control_error=control,passed=passed))
    reaction=np.array([s['radial_reaction_N'] for s in old['states']])
    displacements=data['engineering_strains']*n.L
    work=float(np.trapezoid(reaction,displacements))
    energy=states[-1]['energy_J']
    rows.append(dict(name=name,radial_segments=ns,steps=steps,states=states,passed=all(s['passed'] for s in states),
        final_energy_J=energy,final_reaction_N=float(reaction[-1]),work_J=work,
        work_energy_relative_gap=abs(work-energy)/energy))
    # Replace stale descriptive topology/contact fields after independently
    # re-evaluating each saved equilibrium with the current code.
    old['topology']=n.topology()
    for state,check in zip(old['states'],states):state['contact_gate']=check['contact']
    (P/f'receipts/{name}.json').write_text(json.dumps(old,indent=2))
mesh=[]
for coarse,fine in zip(rows[:2],rows[1:3]):
    mesh.append(dict(coarse=coarse['radial_segments'],fine=fine['radial_segments'],
        final_energy_relative_difference=abs(coarse['final_energy_J']/fine['final_energy_J']-1),
        final_reaction_relative_difference=abs(coarse['final_reaction_N']/fine['final_reaction_N']-1)))
workrows=[rows[1],rows[3],rows[4]]
work_gaps=[r['work_energy_relative_gap'] for r in workrows]
end_energy_spread=(max(r['final_energy_J'] for r in workrows)-min(r['final_energy_J'] for r in workrows))/workrows[-1]['final_energy_J']
n=TwoCellNetwork(24,8);bad=n.q0.copy();bad[:2*n.nn].reshape(-1,2)[:,1]*=.01
rejects_collapsed=not check_contact(n,bad)['passed']
out=dict(scope='Qualified only for this two-cell topology and the sampled 0 to -5% radial shortening path.',
    wall_modulus_Pa=1e9,modulus_calibrated=False,fracture_implemented=False,
    contact_solver_implemented=False,contact_distance_unit_tests_passed=bool(all(checks)),
    rejects_deliberately_collapsed_network=rejects_collapsed,
    mesh_refinement=mesh,work_refinement=[dict(steps=r['steps'],relative_gap=r['work_energy_relative_gap']) for r in workrows],
    increment_refinement_final_energy_relative_spread=end_energy_spread,
    paths=rows,
    acceptance=dict(force_residual_dimensionless=2e-8,positive_free_Hessian=True,
        support_force_balance_N=1e-13,support_moment_balance_Nm=1e-18,
        last_mesh_energy_difference=0.001,last_mesh_reaction_difference=0.001,
        last_work_relative_gap=0.001,increment_final_energy_spread=1e-8),
    passed=bool(all(r['passed'] for r in rows) and all(checks) and rejects_collapsed
        and mesh[-1]['final_energy_relative_difference']<.001
        and mesh[-1]['final_reaction_relative_difference']<.001
        and work_gaps[0]>work_gaps[1]>work_gaps[2] and work_gaps[-1]<.001
        and end_energy_spread<1e-8))
(P/'receipts/network_qualification.json').write_text(json.dumps(out,indent=2))
manifest={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((P/'source').glob('*.py'))}
(P/'receipts/source_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='paths'},indent=2));assert out['passed']
