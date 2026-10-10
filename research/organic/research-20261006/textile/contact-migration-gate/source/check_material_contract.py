from pathlib import Path
import hashlib,json,numpy as np
from moving_contact import quadrature,geometry,controls,L
P=Path(__file__).resolve().parents[1];s,w,labels=quadrature(65);identity=hashlib.sha256(labels.tobytes()+s.tobytes()+w.tobytes()).hexdigest();drift=0.;projection=0.
for case in ['rolling','sliding','birth_loss']:
 for t in np.linspace(0,1,17):
  for rigid in [False,True]:
   g=geometry(s,w,t,case,rigid);length=np.linalg.norm(np.diff(g['partner_nodes'],axis=0),axis=1);drift=max(drift,float(np.max(abs(length-L/2)/(L/2))));projection=max(projection,float(np.max(abs(g['sb']-controls(t,case)[4]))))
s2,w2,l2=quadrature(33);other=hashlib.sha256(l2.tobytes()+s2.tobytes()+w2.tobytes()).hexdigest();out=dict(reference_weight_sum_m=float(w.sum()),reference_span_m=L,maximum_partner_segment_relative_length_drift=drift,maximum_closest_material_coordinate_error_m=projection,material_identity_sha256=identity,changed_quadrature_has_different_identity=identity!=other,passed=bool(abs(w.sum()-L)<1e-16 and drift<1e-12 and projection<1e-15 and identity!=other));(P/'receipts/material_contract.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert out['passed']
