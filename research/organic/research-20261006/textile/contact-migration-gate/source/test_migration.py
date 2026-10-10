import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time,numpy as np,argparse
from moving_contact import run,quadrature,geometry
P=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,default=9);args=ap.parse_args();start=time.monotonic();results={}
for case in ['rolling','sliding','birth_loss']:
 for n in ([128] if case!='birth_loss' else [64,128,256]):
  r=run(case,n,order=args.order);results[f'{case}_{n}']=r;print(case,n,'time',time.monotonic()-start,'balance',r['relative_work_balance_residual'],flush=True)
# Objectivity uses the same prescribed relative motion under a common finite,
# time-dependent rigid rotation and translation, including material velocities.
objective=run('birth_loss',128,order=args.order,superpose=True);plain=results['birth_loss_128'];objective_errors={k:abs(objective[k]-plain[k])/max(abs(plain[k]),1e-30) for k in ['friction_heat_J','numerical_loss_J','normal_work_J','tangential_endpoint_work_J']}
# Normal-force quadrature control for the explicit unique-projection law.
normal=[]
for q in [5,9,17,33,65]:
 s,w,_=quadrature(q);g=geometry(s,w,.5,'birth_loss');normal.append(dict(order=q,force_N=float(g['N'].sum()),energy_J=g['Un']))
ref=normal[-1];qerrors={str(a['order']):dict(force=abs(a['force_N']-ref['force_N'])/ref['force_N'],energy=abs(a['energy_J']-ref['energy_J'])/ref['energy_J']) for a in normal[:-1]}
ss,ww,_=quadrature(args.order);peak_normal_energy=geometry(ss,ww,.5,'birth_loss')['Un'];normal_objective_abs=abs(objective['normal_work_J']-plain['normal_work_J']);normal_objective_scaled=normal_objective_abs/peak_normal_energy
rolling=results['rolling_128'];slide=results['sliding_128'];last=results['birth_loss_256'];checks=dict(rolling_no_material_slip=rolling['mean_integrated_material_slip_m']<1e-12,rolling_no_heat=rolling['friction_heat_J']<1e-22,rolling_memory_preserved=abs(rolling['final_tangential_energy_J']-rolling['initial_tangential_energy_J'])/rolling['initial_tangential_energy_J']<1e-8,rolling_negative_control_detected=rolling['mean_wrong_geometric_point_slip_m']>1e-6,sliding_material_slip=abs(slide['mean_integrated_material_slip_m']-20e-6)<1e-10,sliding_geometric_point_wrongly_stationary=slide['mean_wrong_geometric_point_slip_m']<1e-12,boundary_crossed_without_reset=all(results[k]['active_segment_boundary_crossings']>0 and results[k]['history_resets_at_segment_boundary']==0 for k in results),birth_loss_complete=last['total_births']>0 and last['maximum_births_per_label']==1 and last['final_active_contacts']==0,work_gate=last['relative_work_balance_residual']<.02 and last['numerical_loss_fraction']<.02,normal_work_refines=last['normal_loading_work_relative_error']<results['birth_loss_64']['normal_loading_work_relative_error'],normal_work_objectivity=normal_objective_scaled<1e-7,objectivity_heat_work=max(objective_errors[k] for k in ['friction_heat_J','numerical_loss_J','tangential_endpoint_work_J'])<1e-7,normal_quadrature_full_loading_range_pass=next(x for x in json.loads((P/'receipts/normal_sweep.json').read_text())['orders'] if x['order']==args.order)['passed'])
checks={k:bool(v) for k,v in checks.items()}
out=dict(passed=all(checks.values()),checks=checks,results=results,superposed_motion=objective,objectivity_relative_errors=objective_errors,normal_work_objectivity_absolute_J=normal_objective_abs,normal_work_objectivity_over_peak_energy=normal_objective_scaled,normal_quadrature=normal,normal_quadrature_relative_errors=qerrors,wall_s=time.monotonic()-start,scope='Prescribed two-straight-fibre unique-projection moving-contact gate. Does not replace the double material-pair force law or qualify multi-branch loop contact, full rods, or spinning history.');(P/'receipts'/f'migration_tests_q{args.order}.json').write_text(json.dumps(out,indent=2));print('SUMMARY',json.dumps(dict(passed=out['passed'],checks=checks,objectivity=objective_errors,quadrature=qerrors,wall_s=out['wall_s'])),flush=True)
