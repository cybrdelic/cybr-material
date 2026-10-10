import os,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from corrugated_wall import CorrugatedWall
P=Path(__file__).resolve().parents[1];start=time.monotonic();runs={}
for segments,steps in ((96,28),(192,28)):
 for name,target in (('compression',-.28),('tension',.133)):
  w=CorrugatedWall(segments);x=w.X.copy();history=[];states=[]
  for strain in np.linspace(0,target,steps+1):
   x,r=w.solve(strain,x,160);history.append(r);states.append(x*w.L)
   if not r['passed']:break
  strain=np.array([r['engineering_strain'] for r in history]);force=np.array([r['end_reaction_N'] for r in history]);work=float(np.trapezoid(force,w.L*strain));stored=history[-1]['energy_J'];gap=abs(work-stored)/max(abs(stored),1e-30)
  r=dict(segments=segments,loading_steps=steps,history=history,work_J=work,stored_energy_J=stored,work_gap_fraction=gap,passed_equilibrium=all(r['passed'] for r in history),scope='Elastic wall bending/stretches, no tissue plasticity, fracture, or gas pressure')
  key=f'{segments}_{name}';runs[key]=r;(P/f'receipts/refined_wall_{key}.json').write_text(json.dumps(r,indent=2));np.savez_compressed(P/f'data/refined_wall_{key}.npz',reference_m=w.X*w.L,positions_m=states,engineering_strains=strain)
checks=[]
for name in ('compression','tension'):
 a=runs['96_'+name];b=runs['192_'+name]
 assert a['passed_equilibrium'] and b['passed_equilibrium'], 'Cannot label unfinished equilibrium endpoints as refinement'
 assert abs(a['history'][-1]['engineering_strain']-b['history'][-1]['engineering_strain'])<1e-13, 'Refinement endpoints must share the same target strain'
 energy=abs(a['stored_energy_J']-b['stored_energy_J'])/abs(b['stored_energy_J']);reaction=abs(a['history'][-1]['end_reaction_N']-b['history'][-1]['end_reaction_N'])/abs(b['history'][-1]['end_reaction_N']);coarse=json.loads((P/f'receipts/wall_96_{name}.json').read_text());h=coarse['history'];work14=float(np.trapezoid([r['end_reaction_N'] for r in h],43e-6*np.array([r['engineering_strain'] for r in h])));gap14=abs(work14-h[-1]['energy_J'])/abs(h[-1]['energy_J']);checks.append(dict(control=name,energy_refinement_fraction=energy,reaction_refinement_fraction=reaction,work_gap14=gap14,work_gap28=a['work_gap_fraction'],work_error_reduction=gap14/a['work_gap_fraction']))
out=dict(status='Wall-scale mechanism qualification only',checks=checks,wall_s=time.monotonic()-start,passed=all(r['passed_equilibrium'] for r in runs.values()) and all(c['energy_refinement_fraction']<.003 and c['reaction_refinement_fraction']<.003 and c['work_gap28']<.001 for c in checks))
(P/'receipts/refinement_qualification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
