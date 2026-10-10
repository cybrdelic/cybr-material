from pathlib import Path
import numpy as np,json,time,resource
from adaptive_keys import *
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
source=np.load(ROOT/'receipts/baseline_loops.npz');colors=source['color'];records=[];start=time.time()
for idx in [230,0,266,672,1006,400,900,1400]:
    t0=time.time();color=int(colors[idx]);wi=int((colors[:idx]==color).sum())
    data=np.load(ROOT/f'arrays/colour_{color}_positions.npy',mmap_mode='r');p=data[wi].copy();del data
    rr=float(np.load(ROOT/f'arrays/colour_{color}_radius_per_loop.npy')[wi]);r=np.full(len(p),rr);n=p.shape[1]
    vv=np.diff(p,axis=1);ll=np.linalg.norm(vv,axis=2);tan=vv/ll[:,:,None];kr=np.linalg.norm(np.diff(tan,axis=1),axis=2)/((ll[:,:-1]+ll[:,1:])*.5)*rr
    protected=set(range(6))|set(range(n-6,n))
    for k in np.flatnonzero((kr>.12).any(0))+1:protected.update(range(max(0,k-2),min(n,k+3)))
    plan=sorted(protected|set(range(0,n,4))|{n-1});history=[]
    for attempt in range(7):
        tested_plan=plan.copy()
        q,pos,surf,ang,deg,detail=candidate_and_bound(p,r,tuple(plan),details=True)
        history.append({'keys':len(plan),'position_bound_um':float(pos.max()*1e6),'surface_bound_um':float(surf.max()*1e6),'failed_segments':len(detail['failed_candidate_spans'])})
        if pos.max()<=POSITION_LIMIT and surf.max()<=SURFACE_LIMIT and not deg.any():break
        extra=set()
        for j in detail['failed_candidate_spans']:
            if plan[j+1]-plan[j]>1:extra.add((plan[j]+plan[j+1])//2)
            else:
                for jj in [j-1,j+1]:
                    if 0<=jj<len(plan)-1 and plan[jj+1]-plan[jj]>1:extra.add((plan[jj]+plan[jj+1])//2)
        if not extra:plan=list(range(n))
        else:plan=sorted(set(plan)|extra)
    plan=tested_plan
    assert len(plan)==q.shape[1]
    failed=(pos>POSITION_LIMIT)|(surf>SURFACE_LIMIT)|deg
    counts=np.where(failed,n,len(plan)).astype('i4')
    pos[failed]=0;surf[failed]=0;ang[failed]=0
    assert all(np.array_equal(q[j,[0,-1]],p[j,[0,-1]]) for j in range(len(p)))
    record={'loop':idx,'fibres':187,'mean_keys_per_fibre':float(counts.mean()),'reduced_common_key_count':len(plan),'exact_fallback_fibres':int(failed.sum()),'key_count_histogram':{str(k):int((counts==k).sum()) for k in np.unique(counts)},'key_plan':plan,'position_bound_um':float(pos.max()*1e6),'swept_circle_bound_um':float(surf.max()*1e6),'maximum_tangent_bound_degrees':float(np.degrees(ang.max())),'history':history,'seconds':time.time()-t0,'endpoints_bitwise_equal':True}
    records.append(record);print(json.dumps(record),flush=True)
average=np.mean([x['mean_keys_per_fibre'] for x in records]);wall=time.time()-start
report={'scope':'Bounded eight-loop certificate sample only; shared adaptive key plan per187-fibre loop; no full-field run, scene build or render','position_limit_um':POSITION_LIMIT*1e6,'swept_circle_limit_um':SURFACE_LIMIT*1e6,'macro_pixel_pitch_um':.008/960*1e6,'records':records,'average_keys':float(average),'estimated_total_keys':int(round(average*329868)),'reduction_fraction':float(1-average/144),'sample_wall_seconds':wall,'projected_python_full_field_seconds':wall/len(records)*1764,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
(HERE/'shared_plan_probe.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='records'}),flush=True)
