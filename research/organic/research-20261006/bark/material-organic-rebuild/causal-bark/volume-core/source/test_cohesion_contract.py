from pathlib import Path
import numpy as np,json
from surface_cohesion import evaluate
R=Path(__file__).resolve().parents[1];q=np.array([1.,1.,1.]);J=np.eye(3);a=np.diag([1.,0.,0.]);b=np.diag([0.,1.,0.]);base={'q':q,'Jjump':J,'Jt1':a,'Jt2':b,'area':1.,'kt':1.,'kn':1.,'damage':.5,'normal_sign':1.};checks=[]
for key,values in [('normal_sign',[0.,2.,float('nan')]),('area',[float('nan'),-1]),('kt',[float('nan'),0]),('kn',[float('inf'),0]),('damage',[float('nan'),2]),('q',[np.array([1.,np.nan,1.])]),('Jjump',[np.zeros((2,3))])]:
 for v in values:
  case=base.copy();case[key]=v
  try:evaluate(**case);rejected=False
  except ValueError:rejected=True
  checks.append({'parameter':key,'value':str(v),'pass':rejected})
out={'all_pass':all(c['pass'] for c in checks),'checks':checks};(R/'receipts/cohesion_contract_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['all_pass']
