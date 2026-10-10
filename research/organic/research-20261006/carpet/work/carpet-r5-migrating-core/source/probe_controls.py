from pathlib import Path
import json,time,resource,numpy as np
from migrating_construction import *
start=time.time();data=np.load(R1/'receipts/baseline_loops.npz');rows=[]
cases=[('exact_no_fairing_old_ends','exact',0.,'r1_endpoints'),('exact_no_fairing_continuous_ends','exact',0.,'continuous_new_roots'),('arclength_no_fairing_continuous_ends','arclength_budget',0.,'continuous_new_roots')]
for index in [968,955,356]:
    color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));source=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi].copy();a=evaluate_catmull(source,3);r=float(data['radius'][index].mean()*.045);ra=np.full(a.shape[:2],r,dtype='f4');before=old.clearance(a,ra,k=20);wb=retained_proximity(a,r,data,index)
    for name,phase,fairing,end in cases:
        t0=time.time()
        try:
            p,rr,construction=construct(index,data,source,phase,fairing,end)
            b=evaluate_catmull(p,3);after=old.clearance(b,ra,k=20);wa=retained_proximity(b,r,data,index)
            src=before['regions']['exposed_pile'];dst=after['regions']['exposed_pile'];accepted=dst['negative_segment_fraction']<=src['negative_segment_fraction']+1e-6 and dst['max_curvature_times_radius']<=max(1,src['max_curvature_times_radius'])*1.05 and wa['exposed_potential_penetration_fraction']<=wb['exposed_potential_penetration_fraction']+1e-6
            row={'case':name,'loop':index,'construction':construction,'r1_body_clearance':before,'candidate_body_clearance':after,'r1_retained_wrap_proximity':wb,'candidate_retained_wrap_proximity':wa,'numerical_gate_pass':bool(accepted),'seconds':time.time()-t0}
            short={'case':name,'loop':index,'pass':bool(accepted),'r1_exposed_overlap':src['negative_segment_fraction'],'candidate_exposed_overlap':dst['negative_segment_fraction'],'candidate_buried_overlap':after['regions']['backing_or_root_join']['negative_segment_fraction'],'candidate_curvature_radius':dst['max_curvature_times_radius'],'wrap_potential_overlap':wa['exposed_potential_penetration_fraction'],'seconds':row['seconds']}
        except Exception as error:
            row={'case':name,'loop':index,'numerical_gate_pass':False,'error':str(error)};short=row
        rows.append(row);print(json.dumps(short),flush=True)
        (ROOT/'receipts/controlled_migration_probe.json').write_text(json.dumps({'scope':'Confounder controls and one arclength-budget follow-up; no source scene or render.','rows':rows,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024},indent=2))
