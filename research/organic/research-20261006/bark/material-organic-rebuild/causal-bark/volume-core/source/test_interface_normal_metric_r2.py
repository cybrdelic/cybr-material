from pathlib import Path
import json,numpy as np
from ring_geometry import build
from ring_interfaces_r2 import build_interfaces
R=Path(__file__).resolve().parents[1];rows=[];by_sectors={}
for sectors in (4,8,16):
 h,c,q,m,t=build(sectors,stem_increment_m=0.);ops=build_interfaces(c,q,m,sectors);eta=.01;check=np.max(abs(ops['kt']*eta*ops['dual_normal_length_m']/ops['effective_penalty_modulus_Pa']-1));assert check<2e-15
 pairs=ops['pairs'];kind=ops['kind'];circ=np.array([kind[p]=='fracture' and m[int(i)]['sector']!=m[int(j)]['sector'] for p,(i,j) in enumerate(pairs)]);birth=kind=='birth_interface';vals=ops['kt'][circ];span=2*np.pi*.10515/sectors;expected_length=2*span/3
 old=circ & np.array([m[int(i)]['cohort_id']=='older_outer' for i,j in pairs]);err=float(np.max(abs(ops['dual_normal_length_m'][old]-expected_length)));assert err<1e-14
 by_sectors[sectors]=float(np.mean(ops['kt'][old]));assert np.allclose(ops['dual_normal_length_m'][birth],.5*(.00025+.0003),rtol=1e-13)
 stiff=build_interfaces(c,q,m,sectors,compliance=.005);assert np.allclose(stiff['kt'],2*ops['kt'])
 rows.append({'sectors':sectors,'old_circumferential_face_dual_length_m':float(np.mean(ops['dual_normal_length_m'][old])),'old_circumferential_penalty_Pa_per_m':by_sectors[sectors],'metric_identity_error':float(check),'analytic_dual_length_error_m':err,'minimum_penalty_Pa_per_m':float(ops['kt'].min()),'maximum_penalty_Pa_per_m':float(ops['kt'].max())})
assert np.isclose(by_sectors[8],2*by_sectors[4],rtol=1e-13) and np.isclose(by_sectors[16],2*by_sectors[8],rtol=1e-13)
out={'all_pass':True,'checks':rows,'scope':'Material-normal scale and dimensional refinement of numerical continuity. Not full ring mesh/fracture convergence or measured anatomical cohesion.'};(R/'receipts/interface_normal_metric_r2_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
