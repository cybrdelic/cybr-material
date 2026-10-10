import os,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from corrugated_wall import CorrugatedWall
P=Path(__file__).resolve().parents[1];start=time.monotonic();runs=[]
for steps in (28,56,112):
 w=CorrugatedWall(192);x=w.X.copy();history=[];states=[]
 for strain in np.linspace(0,-.28,steps+1):
  x,r=w.solve(strain,x,200);history.append(r);states.append(x*w.L)
  if not r['passed']:break
 reached=len(history)==steps+1 and all(r['passed'] for r in history)
 eps=np.array([r['engineering_strain'] for r in history]);force=np.array([r['end_reaction_N'] for r in history]);work=float(np.trapezoid(force,w.L*eps));stored=history[-1]['energy_J'];drops=[t.get('energy_drop',0.)*w.unit for r in history for t in r['trace'] if t.get('negative_curvature_step')]
 report=dict(steps=steps,reached_target=reached,stored_J=stored,work_J=work,relative_work_gap=abs(work-stored)/max(abs(stored),1e-30),negative_curvature_steps=len(drops),numerical_descent_energy_drops_J=drops,warning='Numerical relaxation descent is not physical heat. This is a quasistatic elastic branch; work convergence is a check, not a simulated transient energy balance.',history=history)
 (P/f'receipts/loading_work_{steps}.json').write_text(json.dumps(report,indent=2));np.savez_compressed(P/f'data/loading_work_{steps}.npz',reference_m=w.X*w.L,positions_m=states,engineering_strains=eps);runs.append({k:v for k,v in report.items() if k!='history'});print('WORK',steps,reached,report['relative_work_gap'],time.monotonic()-start,flush=True)
 if time.monotonic()-start>90:raise TimeoutError('Loading-work refinement exceeded90s budget')
out=dict(scope='Quasistatic branch work check, not dynamic snap prediction',runs=runs,passed=all(r['reached_target'] for r in runs) and runs[-1]['relative_work_gap']<1e-4,wall_s=time.monotonic()-start);(P/'receipts/loading_work_qualification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
