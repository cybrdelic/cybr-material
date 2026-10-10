"""Compare a fixed-strain endpoint on another material-coordinate mesh."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import argparse,json,time
import numpy as np
from two_cell_network import TwoCellNetwork
from contact_gate import check_contact

P=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--radial-segments',type=int,required=True)
ap.add_argument('--zero-shear',action='store_true')
args=ap.parse_args();start=time.monotonic()
if args.zero_shear:
    from zero_shear_network import ZeroAverageShearNetwork
    network_class=ZeroAverageShearNetwork
else:network_class=TwoCellNetwork
prefix='zero_shear_' if args.zero_shear else ''
source=network_class(96,32);target=network_class(args.radial_segments,args.radial_segments//3)
data=np.load(P/f'data/{prefix}largest_admissible_96.npz');old=data['state'];strain=float(data['engineering_strain'])
q=target.q0.copy();xx=q[:2*target.nn].reshape(-1,2);ox=old[:2*source.nn].reshape(-1,2)
xx[:6]=ox[:6];q[2*target.nn:]=old[2*source.nn:]
for sw,tw in zip(source.walls,target.walls):
    us=np.linspace(0,1,len(sw['nodes']));ut=np.linspace(0,1,len(tw['nodes']))
    for dim in range(2):xx[tw['nodes'],dim]=np.interp(ut,us,ox[sw['nodes'],dim])
solved,row=target.solve(strain,q,maxiter=100)
row['contact_gate']=check_contact(target,solved)
row['passed']=bool(row['passed'] and row['contact_gate']['passed'])
name=prefix+f'endpoint_{args.radial_segments}_{abs(strain)*100:g}pct'
np.savez_compressed(P/f'data/{name}.npz',state=solved,
    positions_m=solved[:2*target.nn].reshape(-1,2)*target.L,
    junction_rotations_rad=solved[2*target.nn:],engineering_strain=strain,
    reference_positions_m=target.X*target.L)
out=dict(topology=target.topology(),state=row,elapsed_seconds=time.monotonic()-start,
         initialization='Interpolation in each wall material coordinate from the accepted 96-segment endpoint, then fully coupled equilibrium',
         passed=row['passed'])
(P/f'receipts/{name}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in row.items() if k!='trace'},indent=2));assert out['passed']
