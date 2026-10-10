"""Advance the qualified network until the first failed equilibrium/contact gate.

Every accepted state is checkpointed. Never seed from a rejected state.
The first rejected trial is evidence of the model limit, not a usable specimen.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import argparse,json,time
import numpy as np
from two_cell_network import TwoCellNetwork
from contact_gate import check_contact
from atomic_io import save_npz,save_json

P=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--target',type=float,default=-.28)
ap.add_argument('--increment',type=float,default=.005)
ap.add_argument('--seed',default='network_96_10_5pct')
ap.add_argument('--radial-segments',type=int,default=96)
ap.add_argument('--resume',action='store_true')
ap.add_argument('--zero-shear',action='store_true')
args=ap.parse_args()
if args.zero_shear:
    from zero_shear_network import ZeroAverageShearNetwork
    network_class=ZeroAverageShearNetwork
    if args.seed=='network_96_10_5pct':args.seed='zero_shear_'+args.seed
else:network_class=TwoCellNetwork
n=network_class(args.radial_segments,args.radial_segments//3)
prefix='zero_shear_' if args.zero_shear else ''
name=prefix+f'contact_limit_{args.radial_segments}'
receipt_path=P/f'receipts/{name}.json';data_path=P/f'data/{name}.npz'
if args.resume and data_path.exists():
    prior=json.loads(receipt_path.read_text())
    if prior['terminal']:raise SystemExit('Already terminal; do not continue a rejected branch')
    data=np.load(data_path);states=list(data['states']);strains=list(data['engineering_strains'])
    records=prior['accepted_states']
else:
    data=np.load(P/f'data/{args.seed}.npz');states=list(data['states']);strains=list(data['engineering_strains'])
    records=json.loads((P/f'receipts/{args.seed}.json').read_text())['states']
q=states[-1];current=float(strains[-1]);rejected=None;terminal=False;start=time.monotonic()
while current>args.target+1e-12:
    candidate_strain=max(args.target,round(current-args.increment,12))
    candidate,record=n.solve(candidate_strain,q)
    contact=check_contact(n,candidate);record['contact_gate']=contact
    equilibrium=record['passed']
    record['passed']=bool(equilibrium and contact['passed'])
    if record['passed']:
        q=candidate;current=candidate_strain;states.append(q.copy());strains.append(current);records.append(record)
        save_npz(P/f'data/{prefix}largest_admissible_{args.radial_segments}.npz',
            state=q,positions_m=q[:2*n.nn].reshape(-1,2)*n.L,junction_rotations_rad=q[2*n.nn:],
            engineering_strain=current,reference_positions_m=n.X*n.L)
    else:
        rejected=record;terminal=True
        save_npz(P/f'data/{prefix}first_rejected_{args.radial_segments}.npz',state=candidate,
            positions_m=candidate[:2*n.nn].reshape(-1,2)*n.L,engineering_strain=candidate_strain)
    reached_target=current<=args.target+1e-12
    terminal=terminal or reached_target
    out=dict(topology=n.topology(),accepted_states=records,largest_accepted_strain=current,
        first_rejected_trial=rejected,terminal=terminal,reached_requested_target=reached_target,
        increment=args.increment,requested_target=args.target,
        stop_reason=('target_reached' if reached_target else
            ('contact_gate_failed' if equilibrium else 'equilibrium_gate_failed') if rejected else 'running'),
        qualified_accepted_states=all(r['passed'] for r in records),elapsed_seconds_this_process=time.monotonic()-start)
    save_npz(data_path,states=np.array(states),engineering_strains=np.array(strains),
        positions_m=np.array(states)[:,:2*n.nn].reshape(-1,n.nn,2)*n.L,
        junction_rotations_rad=np.array(states)[:,2*n.nn:],reference_positions_m=n.X*n.L)
    save_json(receipt_path,out)
    print(json.dumps(dict(strain=candidate_strain,accepted=record['passed'],iterations=record['iterations'],
        free_residual=record['dimensionless_residual'],energy_J=record['energy_J'],
        radial_reaction_N=record['radial_reaction_N'],contact=contact)),flush=True)
    if terminal:break
    if time.monotonic()-start>95:
        print('Safe checkpoint before process budget; use --resume',flush=True);break
print(json.dumps(dict(largest_accepted_strain=current,terminal=terminal,receipt=str(receipt_path))),flush=True)
