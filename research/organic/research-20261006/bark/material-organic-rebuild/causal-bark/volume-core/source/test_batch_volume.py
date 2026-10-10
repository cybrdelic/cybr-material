from pathlib import Path
import numpy as np,json,time,resource
from batch_volume import BatchPrisms
from ring_geometry import build
R=Path(__file__).resolve().parents[1];h,cells,q,meta,triangles=build(16);b=BatchPrisms(cells);start=time.time();E,g,r=b.evaluate(q);batchtime=time.time()-start;esc=[];gsc=[]
for c,x in zip(cells,q):
 e,d,_=c.evaluate(x);esc.append(e);gsc.append(d)
Eg=abs(E-sum(esc))/E;G=np.linalg.norm(g-np.array(gsc))/np.linalg.norm(g);mass=sum(c.mass for c in cells);expected=sum(float(c.material_mass_kg.sum()) for c in h.cohorts);checks=[{'case':'exact_scalar_batch_law','energy_relative_error':Eg,'gradient_relative_error':G,'pass':Eg<1e-13 and G<1e-13},{'case':'closed_ring_retains_all_cohort_mass','material_mass_kg':mass,'expected_mass_kg':expected,'relative_error':abs(mass-expected)/expected,'pass':abs(mass-expected)/expected<1e-13}]
# The full angular chart closes geometrically; there is no artificial open cut.
for layer in range(2):
 first=q[layer*32+1];last=q[layer*32+30];a=first[[0,2,6,8,12,14]];z=last[[1,2,7,8,13,14]];gap=float(np.max(np.linalg.norm(a-z,axis=1)));checks.append({'case':f'periodic_geometric_seam_layer_{layer}','gap_m':gap,'pass':gap<1e-14})
# Young material was born in its present geometry; old material retains prior size.
energy=np.array(esc);checks.append({'case':'birth_age_loads_outer_cohort','young_energy_J':float(energy[:32].sum()),'old_energy_J':float(energy[32:].sum()),'pass':energy[:32].sum()<1e-22 and energy[32:].sum()>1e-5})
out={'all_pass':all(c['pass'] for c in checks),'checks':checks,'cells':len(cells),'position_DOFs':q.size,'one_batch_force_evaluation_s':batchtime,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'scope':'Closed-ring geometry and exact batched constitutive evaluation. No ring equilibrium/fracture claim.'};(R/'receipts/batch_volume_tests.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item()));print(json.dumps(out,indent=2,default=lambda v:v.item()));assert out['all_pass']
