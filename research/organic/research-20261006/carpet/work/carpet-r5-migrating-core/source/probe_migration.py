from pathlib import Path
import json,time,resource,numpy as np
from migrating_construction import *
start=time.time();data=np.load(R1/'receipts/baseline_loops.npz');rows=[]
for index in [968,955,356]:
    t0=time.time();color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));source=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi].copy()
    p,rr,construction=construct(index,data,source)
    np.save(ROOT/f'arrays/probe_loop_{index}.npy',p)
    a=evaluate_catmull(source,3);b=evaluate_catmull(p,3);ra=np.full(a.shape[:2],rr[0,0],dtype='f4')
    before=old.clearance(a,ra,k=20);after=old.clearance(b,ra,k=20)
    wb=retained_proximity(a,float(rr[0,0]),data,index);wa=retained_proximity(b,float(rr[0,0]),data,index)
    src=before['regions']['exposed_pile'];dst=after['regions']['exposed_pile']
    accepted=dst['negative_segment_fraction']<=src['negative_segment_fraction']+1e-6 and dst['max_curvature_times_radius']<=max(1,src['max_curvature_times_radius'])*1.05 and wa['exposed_potential_penetration_fraction']<=wb['exposed_potential_penetration_fraction']+1e-6
    row={'construction':construction,'r1_body_clearance':before,'candidate_body_clearance':after,'r1_retained_wrap_proximity':wb,'candidate_retained_wrap_proximity':wa,'numerical_gate_pass':bool(accepted),'seconds':time.time()-t0}
    rows.append(row);print(json.dumps({'loop':index,'pass':bool(accepted),'seconds':row['seconds'],'r1_overlap_fraction':src['negative_segment_fraction'],'candidate_overlap_fraction':dst['negative_segment_fraction'],'r1_max_curvature_radius':src['max_curvature_times_radius'],'candidate_max_curvature_radius':dst['max_curvature_times_radius'],'r1_wrap_potential_fraction':wb['exposed_potential_penetration_fraction'],'candidate_wrap_potential_fraction':wa['exposed_potential_penetration_fraction'],'migration_median_um':construction['radial_excursion_median_um']}),flush=True)
report={'scope':'One bounded exact-r5-phase, migrating packing construction on three actual central loops. No scene built or rendered.','all_tests_hold':all(x['numerical_gate_pass'] for x in rows),'rows':rows,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'sample_count_per_fibre':430,'final_check':'Catmull-Rom evaluated threefold after final interpolation/fairing; finite3D segment distances, not section-only gaps.','limits':'Body clearance is a bounded nearest-segment check, not a global collision certificate. Retained-wrap capsule test is conservative; exact polygon-surface contact would be required before acceptance.'}
(ROOT/'receipts/migration_probe.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
