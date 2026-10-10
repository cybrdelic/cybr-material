"""One three-loop qualification probe; no scene or render."""
import json,time,resource,hashlib
import numpy as np
from scipy.spatial import cKDTree
from construct_unit import *

def old_wraps(data,index):
    c=data['center'][index];R=float(data['radius'][index].mean());phase=data['phase'][index][None]+np.arange(6)[:,None]*np.pi/3
    return c[None]+.99*R*(np.cos(phase)[:,:,None]*data['frame'][index,None,:,0]+np.sin(phase)[:,:,None]*data['frame'][index,None,:,1])

def pair_proximity(paths,r,obstacles,ro):
    pa=paths[:,:-1].reshape(-1,3).astype('f8');pb=paths[:,1:].reshape(-1,3).astype('f8');mid=(pa+pb)*.5
    a=obstacles[:,:-1].reshape(-1,3);b=obstacles[:,1:].reshape(-1,3);tree=cKDTree((a+b)*.5);gaps=[]
    for start in range(0,len(mid),4096):
        ids=tree.query(mid[start:start+4096],k=min(24,len(a)),workers=1)[1]
        dist=migration.old.segment_distance(pa[start:start+4096,None],pb[start:start+4096,None],a[ids],b[ids]).min(1)
        gaps.append(dist-r-ro)
    gaps=np.concatenate(gaps);exposed=mid[:,2]>.00013
    return {'exposed_min_gap_um':float(gaps[exposed].min()*1e6),'exposed_negative_segment_fraction':float(np.mean(gaps[exposed]<0)),'buried_min_gap_um':float(gaps[~exposed].min()*1e6),'buried_negative_segment_fraction':float(np.mean(gaps[~exposed]<0)),'method':'Finite-segment distances against 24 nearest obstacle segment midpoints; circular tube radii. This is a finite-sampling proximity check.'}

def neighbor_check(body,wraps,r,data,index):
    centers=data['center'];rootmid=centers[:,[0,-1]].mean(1);indices=np.flatnonzero(np.linalg.norm(rootmid-rootmid[index],axis=1)<.0035);indices=indices[indices!=index]
    a=centers[indices,:-1].reshape(-1,3);b=centers[indices,1:].reshape(-1,3);radius=((data['radius'][indices,:-1]+data['radius'][indices,1:])*.5).ravel()
    # Add the six retained strands around each neighbouring original core.
    ww=np.concatenate([old_wraps(data,int(i)) for i in indices]);a=np.concatenate((a,ww[:,:-1].reshape(-1,3)));b=np.concatenate((b,ww[:,1:].reshape(-1,3)));radius=np.r_[radius,np.full(len(ww)*48,16e-6)]
    tree=cKDTree((a+b)*.5);allgap=[]
    for paths,rad in [(body,r),(wraps,16e-6)]:
        pa=paths[:,:-1].reshape(-1,3).astype('f8');pb=paths[:,1:].reshape(-1,3).astype('f8');mid=(pa+pb)*.5;visible=mid[:,2]>.00013
        pa=pa[visible];pb=pb[visible];mid=mid[visible]
        for start in range(0,len(mid),2048):
            ids=tree.query(mid[start:start+2048],k=48,workers=1)[1]
            distance=migration.old.segment_distance(pa[start:start+2048,None],pb[start:start+2048,None],a[ids],b[ids]);gap=(distance-rad-radius[ids]).min(1);allgap.append(gap)
    gap=np.concatenate(allgap)
    return {'neighbor_indices':indices.tolist(),'exposed_min_conservative_gap_um':float(gap.min()*1e6),'exposed_potential_overlap_fraction':float(np.mean(gap<0)),'method':'Finite unit segments versus circumscribed original-core capsules and six retained-wrap capsules, 48 midpoint neighbours. Negative gaps are possible original polygon interactions, not exact polygon intersection depths.'}

def run():
    start=time.time();pair=contact.validate_pair();data=np.load(R1/'receipts/baseline_loops.npz');rows=[]
    print(json.dumps({'pair_validation':pair}),flush=True)
    for index in (968,955,356):
        t0=time.time();color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));source=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi].copy()
        p,wraps,center,dense,r,R,construction=construct(index,data,source)
        q,correction=correct(p,wraps,center,r,R)
        np.savez_compressed(ROOT/f'arrays/rounded_unit_{index}.npz',body=q,wraps=wraps,guide=dense,body_radius=r,wrap_radius=16e-6)
        reference=migration.evaluate_catmull(source,3);before=migration.evaluate_catmull(p,3);after=migration.evaluate_catmull(q,3);ww=migration.evaluate_catmull(wraps,3);ow=old_wraps(data,index)
        rad=np.full(after.shape[:2],r,dtype='f4')
        baseline=migration.old.clearance(reference,rad,k=20);pre=migration.old.clearance(before,rad,k=20);post=migration.old.clearance(after,rad,k=20)
        original_wrap=pair_proximity(reference,r,ow,16e-6);new_wrap=pair_proximity(after,r,ww,16e-6)
        wrap_self=migration.old.clearance(ww,np.full(ww.shape[:2],16e-6,dtype='f4'),k=20)
        body_curvature=contact.curvature_radius_upper_bound(q,r,16);wrap_curvature=contact.curvature_radius_upper_bound(wraps,16e-6,16)
        fine=migration.evaluate_catmull(q,6);envelope=contact.final_envelope_excess(fine,center,r,R)
        neighboring_before=neighbor_check(reference,ow,r,data,index);neighboring_after=neighbor_check(after,ww,r,data,index)
        ends=np.concatenate((q[:,[0,-1]],wraps[:,[0,-1]]));buried=bool(np.all(ends[:,:,2]<-.0004)&np.all(ends[:,:,2]>-.0009)&np.all(abs(ends[:,:,:2])<.0298))
        old_bounds_min=np.minimum((reference-r).min((0,1)),(ow-16e-6).min((0,1)));old_bounds_max=np.maximum((reference+r).max((0,1)),(ow+16e-6).max((0,1)))
        new_bounds_min=np.minimum((after-r).min((0,1)),(ww-16e-6).min((0,1)));new_bounds_max=np.maximum((after+r).max((0,1)),(ww+16e-6).max((0,1)))
        gates={'body_overlap_fraction_not_increased':all(post['regions'][key]['negative_segment_fraction']<=baseline['regions'][key]['negative_segment_fraction']+1e-6 for key in baseline['regions']),'body_maximum_penetration_not_increased':post['min_gap_um']>=baseline['min_gap_um']-.001,'body_strict_kappa_radius_below_one':body_curvature['strict_nonfold_pass'],'wrap_strict_kappa_radius_below_one':wrap_curvature['strict_nonfold_pass'],'wrap_body_overlap_not_increased':new_wrap['exposed_negative_segment_fraction']<=original_wrap['exposed_negative_segment_fraction']+1e-6,'wrap_wrap_no_overlap':wrap_self['negative_segment_fraction']==0,'new_guide_body_envelope_expansion_at_most_2um':envelope<=2e-6,'neighbor_potential_contacts_not_increased':neighboring_after['exposed_potential_overlap_fraction']<=neighboring_before['exposed_potential_overlap_fraction']+1e-6,'all_ends_buried':buried,'guide_roots_apex_exact':construction['roots_exact'] and construction['apex_exact'],'correction_centroid_drift_at_most_3um':correction['loop_centroid_drift_um']<=3}
        row={'loop':index,'construction':construction,'correction':correction,'r1_clearance':baseline,'before_3d_clearance':pre,'after_3d_clearance':post,'original_wrap_proximity':original_wrap,'new_wrap_proximity':new_wrap,'wrap_self_clearance':wrap_self,'body_curvature_certificate':body_curvature,'wrap_curvature_certificate':wrap_curvature,'sampled_new_guide_envelope_excess_um':envelope*1e6,'neighboring_original':neighboring_before,'neighboring_corrected':neighboring_after,'bounds_change_m':{'minimum_delta':(new_bounds_min-old_bounds_min).tolist(),'maximum_delta':(new_bounds_max-old_bounds_max).tolist(),'extent_delta':((new_bounds_max-new_bounds_min)-(old_bounds_max-old_bounds_min)).tolist()},'gates':gates,'pass':all(gates.values()),'seconds':time.time()-t0}
        rows.append(row);print(json.dumps({'loop':index,'pass':row['pass'],'seconds':row['seconds'],'exposed_overlap':post['regions']['exposed_pile']['negative_segment_fraction'],'minimum_gap_um':post['min_gap_um'],'body_curvature_bound':body_curvature['maximum_kappa_radius_upper_bound'],'wrap_curvature_bound':wrap_curvature['maximum_kappa_radius_upper_bound'],'envelope_excess_um':envelope*1e6,'neighbor_potential_contact_fraction':neighboring_after['exposed_potential_overlap_fraction'],'gates':gates}),flush=True)
        (ROOT/'receipts/rounded_unit_probe.json').write_text(json.dumps({'scope':'Three-loop coherent rounded guide/body/wrap geometric construction. No scene, render or full field.','pair_validation':pair,'terminal':len(rows)==3,'all_loops_pass':len(rows)==3 and all(x['pass'] for x in rows),'rows':rows,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024},indent=2))

if __name__=='__main__':run()
