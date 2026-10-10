"""Reconstructed UNVERIFIED proposal. Do not run before dependency confirmation.

The pre-reset patch/run had no completion receipt. This version adds the strict
kappa*radius<1 and final2micrometre envelope gates requested after that draft.
"""
import json,time,resource,numpy as np
from migrating_construction import *
from correct_3d import correct,validate_pair,curvature_radius_upper_bound,final_envelope_excess

def run():
    start=time.time();pair=validate_pair();data=np.load(R1/'receipts/baseline_loops.npz');rows=[]
    print(json.dumps({'pair_validation':pair}),flush=True)
    for index in [968,955,356]:
        t0=time.time();color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));source=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi].copy()
        p,rr,construction=construct(index,data,source,'arclength_budget',0.,'continuous_new_roots');r=float(rr[0,0]);center=data['center'][index];R=float(data['radius'][index].mean());phase=data['phase'][index][None]+np.arange(6)[:,None]*np.pi/3;wraps=center[None]+.99*R*(np.cos(phase)[:,:,None]*data['frame'][index,None,:,0]+np.sin(phase)[:,:,None]*data['frame'][index,None,:,1])
        q,correction=correct(p,r,center,R,wraps);np.save(ROOT/f'arrays/corrected_3d_loop_{index}.npy',q)
        aa=evaluate_catmull(source,3);bb=evaluate_catmull(p,3);cc=evaluate_catmull(q,3);rad=np.full(aa.shape[:2],r,dtype='f4')
        baseline=old.clearance(aa,rad,k=20);before=old.clearance(bb,rad,k=20);after=old.clearance(cc,rad,k=20);wr0=retained_proximity(aa,r,data,index);wr1=retained_proximity(cc,r,data,index)
        old_region=baseline['regions']['exposed_pile'];a=after['regions']['exposed_pile']
        curvature=curvature_radius_upper_bound(q,r,subdivisions=16)
        fine=evaluate_catmull(q,6);envelope=final_envelope_excess(fine,center,r,R)
        ends=q[:,[0,-1]];buried=bool(np.all(abs(ends[:,:,:2])<.0298) and np.all((ends[:,:,2]<-.0004)&(ends[:,:,2]>-.0009)))
        gates={'body_overlap_not_increased':bool(a['negative_segment_fraction']<=old_region['negative_segment_fraction']+1e-6),'strict_kappa_radius_below_one':curvature['strict_nonfold_pass'],'wrap_proximity_not_worse':bool(wr1['exposed_potential_penetration_fraction']<=wr0['exposed_potential_penetration_fraction']+1e-6),'evaluated_envelope_expansion_at_most_2um':bool(envelope<=2e-6),'loop_centroid_drift_at_most_3um':bool(correction['loop_centroid_drift_um']<=3),'all_ends_buried':buried}
        gates['buried_body_overlap_not_increased']=bool(after['regions']['backing_or_root_join']['negative_segment_fraction']<=baseline['regions']['backing_or_root_join']['negative_segment_fraction']+1e-6)
        # A possible pass gets the denser finite-segment check before acceptance.
        finer=None
        if all(gates.values()):
            rf=np.full(fine.shape[:2],r,dtype='f4');ref_fine=evaluate_catmull(source,6)
            fine_before=old.clearance(ref_fine,rf,k=32);fine_after=old.clearance(fine,rf,k=32)
            gates['finer_body_overlap_not_increased']=bool(fine_after['regions']['exposed_pile']['negative_segment_fraction']<=fine_before['regions']['exposed_pile']['negative_segment_fraction']+1e-6)
            finer={'r1':fine_before,'candidate':fine_after}
        row={'loop':index,'construction':construction,'correction':correction,'r1_clearance':baseline,'before_3d_clearance':before,'after_3d_clearance':after,'r1_wrap_proximity':wr0,'after_wrap_proximity':wr1,'curvature_certificate':curvature,'evaluated_envelope_excess_um':envelope*1e6,'gates':gates,'finer_check':finer,'pass':all(gates.values()),'seconds':time.time()-t0}
        rows.append(row);print(json.dumps({'loop':index,'pass':row['pass'],'seconds':row['seconds'],'gates':gates,'after_overlap':a['negative_segment_fraction'],'curvature_upper_bound':curvature['maximum_kappa_radius_upper_bound'],'envelope_excess_um':envelope*1e6,'maximum_control_displacement_um':correction['maximum_control_displacement_um']}),flush=True)
        (ROOT/'receipts/fixed_budget_3d_probe.json').write_text(json.dumps({'scope':'Three-loop fixed-budget geometric correction only. No scene, field run, render, mechanical solve or post-contact smoothing.','reconstruction_status':'New post-recovery execution of reconstructed proposal, after parent verified final r1 array hashes.','terminal':len(rows)==3,'all_loops_pass':len(rows)==3 and all(row['pass'] for row in rows),'envelope_check_limit':'Final Catmull positions at six subdivisions per native span; a sampled surface-envelope test, not a continuous envelope certificate.','pair_validation':pair,'rows':rows,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024},indent=2))

if __name__=='__main__':run()
