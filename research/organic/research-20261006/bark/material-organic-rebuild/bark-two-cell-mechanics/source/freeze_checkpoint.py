"""Recheck completed controls, create per-state atomic files, and hash deliverables."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,time
import numpy as np
from scipy.linalg import eigh
from two_cell_network import TwoCellNetwork
from zero_shear_network import ZeroAverageShearNetwork
from contact_gate import check_contact
from atomic_io import save_npz,save_json

P=Path(__file__).resolve().parents[1];start=time.monotonic();branches=[]
for prefix,cls in [('',TwoCellNetwork),('zero_shear_',ZeroAverageShearNetwork)]:
    name=prefix+'contact_limit_96'
    record=json.loads((P/f'receipts/{name}.json').read_text())
    assert record['terminal'] and record['reached_requested_target']
    data=np.load(P/f'data/{name}.npz');n=cls(96,32);checks=[]
    folder=P/'data/checkpoints'/('zero_average_shear' if prefix else 'free_shear')
    folder.mkdir(parents=True,exist_ok=True)
    for i,(q,strain) in enumerate(zip(data['states'],data['engineering_strains'])):
        e,g,H,_=n.evaluate(q,True);support=n.support_gradient(g)
        pos=q[:2*n.nn].reshape(-1,2);force=support[:2*n.nn].reshape(-1,2)*n.unit/n.L
        imbalance=float(np.linalg.norm(np.sum(force,axis=0)))
        moment=abs(float(np.sum((pos[:,0]*force[:,1]-pos[:,1]*force[:,0])*n.L)))
        residual=float(np.max(abs(n.reduced_gradient(g))))
        eig=float(eigh(n.reduced_hessian(H),eigvals_only=True,subset_by_index=[0,0])[0])
        contact=check_contact(n,q)
        shear=float((pos[3:6,0].mean()-pos[:3,0].mean())*n.L)
        passed=bool(residual<2e-8 and eig>0 and imbalance<1e-13 and moment<1e-18
            and contact['passed'] and (not prefix or abs(shear)<1e-15))
        checks.append(dict(engineering_strain=float(strain),energy_J=e*n.unit,projected_residual=residual,
            min_constrained_eigenvalue=eig,support_force_balance_N=imbalance,
            support_moment_balance_Nm=moment,mean_transverse_shear_m=shear,contact=contact,passed=passed))
        save_npz(folder/f'state_{i:03d}.npz',state=q,engineering_strain=strain,
            positions_m=pos*n.L,junction_rotations_rad=q[2*n.nn:])
        save_json(folder/f'state_{i:03d}.json',checks[-1])
    endpoint=[json.loads((P/f'receipts/{prefix}endpoint_48_28pct.json').read_text())['state'],
              record['accepted_states'][-1],
              json.loads((P/f'receipts/{prefix}endpoint_192_28pct.json').read_text())['state']]
    differences=[]
    for coarse,fine,ns in zip(endpoint[:-1],endpoint[1:],[48,96]):
        row=dict(coarse=ns,fine=2*ns,energy_relative_difference=abs(coarse['energy_J']/fine['energy_J']-1),
            radial_reaction_relative_difference=abs(coarse['radial_reaction_N']/fine['radial_reaction_N']-1))
        if prefix:row['shear_reaction_relative_difference']=abs(coarse['generalized_shear_reaction_N']/fine['generalized_shear_reaction_N']-1)
        differences.append(row)
    work=float(np.trapezoid([s['radial_reaction_N'] for s in record['accepted_states']],data['engineering_strains']*n.L))
    gap=abs(work-checks[-1]['energy_J'])/checks[-1]['energy_J']
    last=differences[-1]
    branch_pass=bool(all(c['passed'] for c in checks) and all(e['passed'] for e in endpoint)
        and last['energy_relative_difference']<.001 and last['radial_reaction_relative_difference']<.001
        and last.get('shear_reaction_relative_difference',0)<.005 and gap<.001)
    branches.append(dict(control='zero_average_shear' if prefix else 'free_shear',states=checks,
        endpoint_refinement=differences,work_energy_relative_gap=gap,
        passed=branch_pass,atomic_checkpoint_directory=str(folder)))
qualification=dict(scope='Two-cell planar radial elastic mechanism; both declared controls 0 to -28%, no contact encountered within sampled path',
    wall_modulus_Pa=1e9,modulus_calibrated=False,circumferential_fracture_calibrated=False,
    unresolved_junction_radius_m=2.5e-6,contact_forces_implemented=False,
    increment=.005,branches=branches,passed=all(b['passed'] for b in branches),
    elapsed_seconds=time.monotonic()-start)
save_json(P/'receipts/final_qualification.json',qualification)
assert qualification['passed']
files={}
for folder in ['source','data','receipts','scenes']:
    for f in sorted((P/folder).rglob('*')):
        if not f.is_file() or f.suffix in ['.pyc','.log','.tmp'] or f.name in ['final_manifest.json','source_manifest.json']:continue
        files[str(f.relative_to(P))]=dict(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
files['README.md']=dict(bytes=(P/'README.md').stat().st_size,sha256=hashlib.sha256((P/'README.md').read_bytes()).hexdigest())
manifest=dict(title='Closed two-cell radial mechanism checkpoint',passed=True,
    qualification='receipts/final_qualification.json',render_handoff='receipts/render_handoff.json',
    source_wall_snapshot_sha256='edf2978d97ac4ae34240528e1fbf1dadfe0b738cbf1dcf6cddad99e349839c40',
    controls=['free_shear','zero_average_shear'],terminal_radial_strain=-.28,
    no_render_performed=True,files=files)
save_json(P/'receipts/final_manifest.json',manifest)
save_json(P/'receipts/source_manifest.json',{k:v['sha256'] for k,v in files.items() if k.startswith('source/')})
print(json.dumps(dict(passed=True,manifest=str(P/'receipts/final_manifest.json'),
    branch_summaries=[{k:v for k,v in b.items() if k!='states'} for b in branches],
    file_count=len(files),elapsed_seconds=time.monotonic()-start),indent=2))
