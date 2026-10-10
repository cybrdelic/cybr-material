from pathlib import Path
import numpy as np,json,time,resource
from adaptive_keys import *
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
source=np.load(ROOT/'receipts/baseline_loops.npz');colors=source['color']
indices=[230,0,266,672,1006,400,900,1400]
report=[];start=time.time();allcounts=[];allpos=[];allsurf=[];allangles=[]
for idx in indices:
    color=int(colors[idx]);within=int(np.count_nonzero(colors[:idx]==color))
    data=np.load(ROOT/f'arrays/colour_{color}_positions.npy',mmap_mode='r');p=data[within].copy();del data
    rr=float(np.load(ROOT/f'arrays/colour_{color}_radius_per_loop.npy')[within]);r=np.full(len(p),rr)
    out,counts,pos,surf,angle,degen=reduce_batch(p,r)
    assert all(np.array_equal(q[[0,-1]],p[i,[0,-1]]) for i,q in enumerate(out))
    item={'original_loop':idx,'fibres':len(p),'original_keys':int(p.shape[0]*p.shape[1]),'reduced_keys':int(counts.sum()),'key_counts_histogram':{str(x):int(np.count_nonzero(counts==x)) for x in np.unique(counts)},'maximum_continuous_position_bound_um':float(pos.max()*1e6),'maximum_swept_circle_bound_um':float(surf.max()*1e6),'maximum_tangent_angle_bound_degrees':float(np.degrees(angle.max())),'fibres_with_degenerate_rejected_attempt':int(degen.sum()),'endpoints_bitwise_equal':True,'radii_unchanged':True}
    report.append(item);allcounts.extend(counts.tolist());allpos.extend(pos.tolist());allsurf.extend(surf.tolist());allangles.extend(angle.tolist());print(json.dumps(item),flush=True)
out={'scope':'Eight complete frozen r1 loops. No construction changes, no crown guard fix, no Blender scene or render.','certification':'Every overlap of original and candidate cubic spans is compared by restricted cubic Bezier control differences. Tangent cones use quadratic-derivative cross/dot Bernstein bounds. Surface bound pairs radius circles by minimal rotation.','position_limit_m':POSITION_LIMIT,'swept_circle_limit_m':SURFACE_LIMIT,'macro_width_m':.008,'macro_pixels':960,'surface_bound_pixels':SURFACE_LIMIT/(.008/960),'loops':report,'average_keys_per_fibre':float(np.mean(allcounts)),'original_keys_per_fibre':144,'estimated_full_panel_keys':int(round(np.mean(allcounts)*329868)),'key_reduction_fraction':float(1-np.mean(allcounts)/144),'seconds':time.time()-start,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
(HERE/'adaptive_key_probe.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='loops'}),flush=True)
