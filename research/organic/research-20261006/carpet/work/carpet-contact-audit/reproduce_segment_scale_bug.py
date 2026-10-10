"""Analytic regression for the frozen r2 distance routine; no construction changes."""
from pathlib import Path
import hashlib, importlib.util, json, sys
import numpy as np
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'carpet-rounded-yarn-unit/revisions/r2_section_feasibility/source/advance_probe.py'
sys.path.insert(0,str(SOURCE.parent))
spec=importlib.util.spec_from_file_location('frozen_advance',SOURCE)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
a=np.array([[-1.,0.,0.]]);b=-a;c=np.array([[0.,-1.,0.]]);d=-c
rows=[]
for scale in [1.,1e-3,1e-6,1e-7,1e-8,1e-9]:
 distance=float(module.exact_segments(a*scale,b*scale,c*scale,d*scale)[0])
 rows.append({'scale_m':scale,'computed_distance_m':distance,'expected_distance_m':0.,'error_over_scale':distance/scale})
report={'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'case':'Two perpendicular finite segments crossing exactly at their midpoints. Uniformly scale every coordinate. The minimum distance is always zero.','results':rows,'scale_invariance_pass':all(r['computed_distance_m']==0 for r in rows),'cause':'The interior closest-point branch compares a determinant with units m^4 against a fixed 1e-30 threshold. It drops the actual crossing at small scales and returns endpoint distance.','scope':'Demonstrates a contact-gate defect. It does not change geometry, qualify the advancement model, or establish a material appearance result.'}
(ROOT/'segment_scale_bug.json').write_text(json.dumps(report,indent=2))
assert not report['scale_invariance_pass'], 'Frozen failure was not reproduced; investigate source identity.'
print(json.dumps({'failure_reproduced':True,'source_sha256':report['source_sha256'],'receipt':str(ROOT/'segment_scale_bug.json')}))
