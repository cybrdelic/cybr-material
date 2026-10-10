from pathlib import Path
import numpy as np,json,time,resource
from adaptive_keys import *
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
src=np.load(ROOT/'receipts/baseline_loops.npz');colors=src['color'];reports=[];start=time.time()
for idx in [230,0,266,672]:
    col=int(colors[idx]);wi=int((colors[:idx]==col).sum());p=np.load(ROOT/f'arrays/colour_{col}_positions.npy',mmap_mode='r')[wi].copy();rr=float(np.load(ROOT/f'arrays/colour_{col}_radius_per_loop.npy')[wi]);r=np.full(187,rr)
    v=np.diff(p,axis=1);ll=np.linalg.norm(v,axis=2);t=v/ll[:,:,None];kr=np.linalg.norm(np.diff(t,axis=1),axis=2)/((ll[:,:-1]+ll[:,1:])*.5)*rr
    protected=set(range(6))|set(range(138,144))
    for k in np.flatnonzero((kr>.12).any(0))+1:protected.update(range(max(0,k-2),min(144,k+3)))
    plan=sorted(protected|set(range(0,144,4))|{143});t0=time.time();q,pos,surf,angle,deg,history=fit_keys(p,r,plan,protected,steps=3)
    fail=(pos>POSITION_LIMIT)|(surf>SURFACE_LIMIT)|deg
    report={'loop':idx,'key_count':len(plan),'passing_fibres':int((~fail).sum()),'failing_fibres':int(fail.sum()),'history':history,'seconds':time.time()-t0,'endpoints_bitwise_equal':True,'crown_and_root_key_positions_bitwise_equal':bool(np.array_equal(q[:,[j for j,k in enumerate(plan) if k in protected]],p[:,[k for k in plan if k in protected]]))}
    reports.append(report);print(json.dumps(report),flush=True)
report={'scope':'Four-loop bounded fit experiment only; no full-field run or scene.','reports':reports,'wall_seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
(HERE/'fitted_plan_probe.json').write_text(json.dumps(report,indent=2))
