"""Exchange accepted histories at actual reversal/unload controls, with replay."""
from pathlib import Path
import json,numpy as np
from symmetric_spans import SymmetricSpans,swap,reorder
from coupled_spans import S
from driven_cycle import digest
P=Path(__file__).resolve().parents[1];n=128;f=SymmetricSpans();r=json.loads((P/'receipts/driven_cycle_128.json').read_text());checks=[]
for leg in [1,2,3]:
 a=np.load(P/'receipts'/f'driven_cycle_128_leg_{leg}.npz');old=f.state(a['q'],a['elastic_a'],a['elastic_b']);row=r['rows'][leg*n];q=old['q'].copy();q[0]=row['x_m']/S;q[2]=row['normal_z_m']/S;t=f.trial(q,old);s=f.trial(swap(q),f.exchange_state(old));w=reorder(s['wrenches']);force=float(np.linalg.norm(t['wrenches']-w)/max(np.linalg.norm(t['wrenches']),1e-20));heat=abs(t['friction_heat_J']-s['friction_heat_J']);work=abs(t['tangential_endpoint_work_J']-s['tangential_endpoint_work_J']);replay=np.linalg.norm(t['wrenches']-np.array(row['normal_wrenches'])-np.array(row['friction_wrenches']));checks.append(dict(control_leg=leg,history_digest_matches=digest(old)==row['accepted_history_digest'],exchange_force_relative=force,exchange_heat_error_J=heat,exchange_work_error_J=work,replay_wrench_error=float(replay)))
result=dict(scope='Actual committed reversal and unload histories, exchanged with bodies; replay of next accepted control.',checks=checks,passed=all(z['history_digest_matches'] and z['exchange_force_relative']<1e-10 and z['exchange_heat_error_J']<1e-20 and z['exchange_work_error_J']<1e-20 and z['replay_wrench_error']<1e-12 for z in checks));(P/'receipts/driven_exchange.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
