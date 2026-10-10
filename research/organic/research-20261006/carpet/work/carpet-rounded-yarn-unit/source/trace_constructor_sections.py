"""One authorized deterministic constructor replay; no 3D correction."""
import sys,json,time,resource,hashlib
import numpy as np
from scipy.spatial import cKDTree
import construct_unit as construction
from construct_unit import ROOT,R1,migration

def metric_table(xy,curvature,rate,r):
    gaps=[]
    for k,x in enumerate(xy):
        d=x[:,None]-x[None];mid=(x[:,None]+x[None])*.5;v=rate*np.stack((-mid[:,:,1],mid[:,:,0]),axis=-1);axial=1-mid@curvature[k]
        metric=d-v*(np.sum(v*d,axis=-1)/(axial*axial+np.sum(v*v,axis=-1)))[:,:,None]
        gap=np.sqrt(np.maximum(np.sum(d*metric,axis=-1),0))-2*r
        np.fill_diagonal(gap,np.inf);gaps.append(gap)
    return np.asarray(gaps)

def collisions(paths,r,metric,section,kind):
    a=paths[:,:-1].reshape(-1,3).astype('f8');b=paths[:,1:].reshape(-1,3).astype('f8');mid=(a+b)*.5
    segments=paths.shape[1]-1;length=np.linalg.norm(b-a,axis=1);tree=cKDTree(mid)
    affected=np.zeros(len(a),bool);worst=np.inf;counts={'same_original_span':0,'adjacent_original_span':0,'two_spans_apart':0,'three_or_more_spans_apart':0};tested_unique=set();physical_clear=set();margin_clear=set();cubic_only=set();incomplete=0;worst_pair=None
    for start in range(0,len(a),1024):
        distances,near=tree.query(mid[start:start+1024],k=64,workers=1)
        safe_horizon=2*r+(length[start:start+1024]+length.max())*.5
        incomplete+=int((distances[:,-1]<=safe_horizon).sum())
        first=np.repeat(np.arange(start,min(start+1024,len(a))),64);second=near.ravel();valid=(first//segments!=second//segments)
        first=first[valid];second=second[valid]
        gap=migration.old.segment_distance(a[first],b[first],a[second],b[second])-2*r
        if len(gap):
            location=int(np.argmin(gap))
            if gap[location]<worst:worst=float(gap[location]);worst_pair=[int(first[location]//segments),int(second[location]//segments),int(first[location]%segments),int(second[location]%segments)]
        negative=gap<0;np.logical_or.at(affected,first[negative],True)
        for first_id,second_id in zip(first[negative],second[negative]):
            fi=int(first_id//segments);fj=int(second_id//segments)
            ki=min(int((first_id%segments+.5)*96/segments),95);kj=min(int((second_id%segments+.5)*96/segments),95)
            if fi>fj:fi,fj,ki,kj=fj,fi,kj,ki
            key=(fi,fj,ki,kj)
            if key in tested_unique:continue
            tested_unique.add(key);delta=abs(ki-kj)
            counts[['same_original_span','adjacent_original_span','two_spans_apart','three_or_more_spans_apart'][min(delta,3)]]+=1
            values=metric[[ki,ki+1,kj,kj+1],fi,fj]
            if values.min()>=0:physical_clear.add(key)
            if values.min()>=2e-6:margin_clear.add(key)
            if kind=='cubic':
                straight_gap=float(migration.old.segment_distance(section[ki,fi],section[ki+1,fi],section[kj,fj],section[kj+1,fj])-2*r)
                if straight_gap>=0:cubic_only.add(key)
    visible=mid[:,2]>.00013
    return {'paths':len(paths),'segments_per_path':segments,'minimum_gap_um':worst*1e6,'affected_exposed_segment_fraction':float(affected[visible].mean()),'affected_buried_segment_fraction':float(affected[~visible].mean()),'unique_original_span_pair_contacts':len(tested_unique),'contact_locality':counts,'contacts_with_all_four_section_metric_gaps_nonnegative':len(physical_clear),'contacts_with_all_four_section_metric_gaps_at_least_2um':len(margin_clear),'cubic_contacts_whose_same_original_straight_span_pair_is_clear':len(cubic_only),'worst_pair_fibres_then_sample_spans':worst_pair,'query_k':64,'segments_without_complete_midpoint_horizon_certificate':incomplete,'total_segments':len(a),'method':'Finite segment distances. Pair classifications use the two original section endpoints bracketing each sampled segment. Horizon certificate checks whether all possible contact midpoints fit inside the 64-neighbour query.'}

def run():
    start=time.time();index=968;data=np.load(R1/'receipts/baseline_loops.npz');color=int(data['color'][index]);wi=int(np.count_nonzero(data['color'][:index]==color));source=np.load(R1/f'arrays/colour_{color}_positions.npy',mmap_mode='r')[wi].copy();captured={}
    original_hash=hashlib.sha256((ROOT/'source/construct_unit.py').read_bytes()).hexdigest()
    def capture(frame,event,arg):
        if event=='return' and frame.f_code is construction.construct.__code__:
            for key in ('xy','section','s','center','curvature','phase','body','sample','rate','p'):
                value=frame.f_locals[key];captured[key]=value.copy() if isinstance(value,np.ndarray) else value
    sys.setprofile(capture)
    try:p,wraps,center,dense,r,R,details=construction.construct(index,data,source)
    finally:sys.setprofile(None)
    assert original_hash==hashlib.sha256((ROOT/'source/construct_unit.py').read_bytes()).hexdigest()
    np.savez_compressed(ROOT/'arrays/constructor_trace_968.npz',**captured,wraps=wraps)
    xy=captured['xy'];metric=metric_table(xy,captured['curvature'],captured['rate'],r);i,j=np.triu_indices(187,1);section_gaps=metric[:,i,j]
    actual=np.linalg.norm(xy[:,:,None]-xy[:,None,:],axis=-1)-2*r;actual=actual[:,i,j]
    rows=[]
    for k in range(len(xy)):
        rows.append({'section':k,'guide_z_m':float(captured['center'][k,2]),'declared_metric_min_physical_gap_um':float(section_gaps[k].min()*1e6),'declared_metric_negative_pairs':int((section_gaps[k]<0).sum()),'declared_metric_pairs_missing_2um_margin':int((section_gaps[k]<2e-6).sum()),'euclidean_section_min_physical_gap_um':float(actual[k].min()*1e6),'euclidean_section_negative_pairs':int((actual[k]<0).sum())})
    straight=collisions(captured['section'].transpose(1,0,2),r,metric,captured['section'],'straight');print(json.dumps({'straight':straight}),flush=True)
    cubic=collisions(captured['body'],r,metric,captured['section'],'cubic');print(json.dumps({'cubic':cubic}),flush=True)
    replay=migration.evaluate_catmull(p,3);old_report=migration.old.clearance(replay,np.full(replay.shape[:2],r,dtype='f4'),k=20);expected=json.loads((ROOT/'receipts/rounded_unit_probe.json').read_text())['rows'][0]['before_3d_clearance']
    step=np.linalg.norm(np.diff(xy,axis=0),axis=-1);visible=captured['center'][:,2]>.00013
    result={'loop':index,'scope':'Single authorized deterministic replay of unchanged section constructor. No 16-pass 3D correction, parameter change or scene.','constructor_source_sha256':original_hash,'precorrection_position_sha256':hashlib.sha256(p.tobytes()).hexdigest(),'saved_precorrection_hash_available':False,'replayed_precorrection_clearance_matches_prior_receipt_exactly':old_report==expected,'section_spacing_um':float(np.diff(captured['s'])[0]*1e6),'fibre_diameter_um':2*r*1e6,'maximum_transverse_step_um':float(step.max()*1e6),'median_transverse_step_um':float(np.median(step)*1e6),'transverse_steps_over_one_diameter':int((step>2*r).sum()),'total_transverse_steps':step.size,'sections':{'count':len(xy),'visible_count':int(visible.sum()),'visible_sections_failing_declared_metric_physical_clearance':int((section_gaps[visible].min(1)<0).sum()),'visible_sections_missing_declared_2um_margin':int((section_gaps[visible].min(1)<2e-6).sum()),'minimum_declared_metric_physical_gap_um':float(section_gaps.min()*1e6),'visible_sections_with_euclidean_overlap':int((actual[visible].min(1)<0).sum()),'minimum_euclidean_physical_gap_um':float(actual.min()*1e6),'rows':rows},'straight_connections':straight,'cubic_interpolation':cubic,'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    (ROOT/'receipts/constructor_section_trace_968.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ('sections','straight_connections','cubic_interpolation')}|{'sections_summary':{k:v for k,v in result['sections'].items() if k!='rows'}}))

if __name__=='__main__':run()
