import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import argparse,json,time
import numpy as np
from two_cell_network import TwoCellNetwork
from contact_gate import check_contact
from atomic_io import save_npz,save_json

P=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--radial-segments',type=int,default=24)
ap.add_argument('--steps',type=int,default=10);ap.add_argument('--strain',type=float,default=-.05)
ap.add_argument('--zero-shear',action='store_true')
args=ap.parse_args()
if args.zero_shear:
    from zero_shear_network import ZeroAverageShearNetwork
    network_class=ZeroAverageShearNetwork
else:network_class=TwoCellNetwork
n=network_class(args.radial_segments,args.radial_segments//3)
name=('zero_shear_' if args.zero_shear else '')+f'network_{args.radial_segments}_{args.steps}_{abs(args.strain)*100:g}pct'
states=[];receipts=[];q=n.q0.copy();start=time.monotonic()
for strain in np.linspace(0,args.strain,args.steps+1):
    q,row=n.solve(float(strain),q)
    row['contact_gate']=check_contact(n,q)
    row['passed']=bool(row['passed'] and row['contact_gate']['passed'])
    states.append(q.copy());receipts.append(row)
    data=dict(states=np.array(states),positions_m=np.array(states)[:,:2*n.nn].reshape(-1,n.nn,2)*n.L,
              junction_rotations_rad=np.array(states)[:,2*n.nn:],reference_positions_m=n.X*n.L,
              engineering_strains=np.array([r['engineering_strain'] for r in receipts]))
    save_npz(P/f'data/{name}.npz',**data)
    out=dict(scope='Two closed planar cells, elastic strip network; post-solve contact admissibility; no full bark law',
             topology=n.topology(),states=receipts,passed=all(r['passed'] for r in receipts),
             elapsed_seconds=time.monotonic()-start)
    save_json(P/f'receipts/{name}.json',out)
    print(json.dumps({k:v for k,v in row.items() if k not in ('trace','contact_gate')}),flush=True)
    if not row['passed']:raise SystemExit('A load-state gate failed; checkpoint retained')
print(json.dumps(dict(passed=out['passed'],elapsed_seconds=out['elapsed_seconds'],path=str(P/f'receipts/{name}.json'))),flush=True)
