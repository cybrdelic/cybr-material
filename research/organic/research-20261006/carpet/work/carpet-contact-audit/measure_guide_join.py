"""Read-only evaluated guide diagnostic; no geometry or process change."""
from pathlib import Path
import sys,hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent
UNIT=ROOT.parent/'carpet-rounded-yarn-unit'
sys.path.insert(0,str(UNIT/'revisions/r2_section_feasibility/source'))
from advance_probe import Field
DATA=ROOT.parent/'carpet-r5-packed-core/receipts/baseline_loops.npz'
f=Field(968,np.load(DATA));arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(f.guide,axis=0),axis=1))];join=float(arc[512]);rows=[]
for delta in [-10,-5,-2,-1,-.5,0,.5,1,2,5,10]:
 s=join+delta*1e-6;v=f.spline(s,1);a=f.spline(s,2);signed=float(np.cross(v,a)@f.normal/np.linalg.norm(v)**3)
 rows.append({'distance_from_join_um':delta,'signed_kappa_times_R':signed*f.R,'parameter_speed':float(np.linalg.norm(v))})
report={'loop':968,'join_parameter_m':join,'interpretation':'The evaluated dense cubic guide is C2, but concentrates a sign-reversing curvature change and small overshoot into a few micrometres at the authored C1 leg/circle join. This is a geometry-source diagnostic, not proof that every contact failure has this cause.','scope':'Original guide and saved geometry unchanged.','inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [UNIT/'arrays/rounded_unit_968.npz',DATA]},'rows':rows}
(ROOT/'guide_join_diagnostic.json').write_text(json.dumps(report,indent=2));print(json.dumps({'receipt':str(ROOT/'guide_join_diagnostic.json'),'join_m':join}))
