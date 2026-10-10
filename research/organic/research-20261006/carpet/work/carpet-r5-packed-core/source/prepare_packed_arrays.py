"""Stream the full selected panel's constructed fibre points to bounded arrays."""
from pathlib import Path
import json,time,resource,hashlib
import numpy as np
import packed_construction as model
ROOT=Path(__file__).resolve().parents[1]
d=np.load(ROOT/'receipts/baseline_loops.npz');outdir=ROOT/'arrays';outdir.mkdir(exist_ok=True)
start=time.time();report={'construction':'Geometric constructed yarn, not solved manufacturing mechanics. No volume or stock-length conservation claim.','groups':[],'endpoints_z_m':[],'bounds_min_m':[1.,1.,1.],'bounds_max_m':[-1.,-1.,-1.]}
for color in range(4):
    ids=np.flatnonzero(d['color']==color);path=outdir/f'colour_{color}_positions.npy'
    data=np.lib.format.open_memmap(path,mode='w+',dtype='f4',shape=(len(ids),model.FIBRES,model.POINTS,3))
    radii=[];length=[]
    for j,index in enumerate(ids):
        p,rr=model.construct(d['center'][index],d['radius'][index],int(index));data[j]=p;radii.append(rr[0,0]);length.extend(np.linalg.norm(np.diff(p,axis=1),axis=2).sum(1).tolist())
        report['endpoints_z_m'].extend([float(p[:,0,2].min()),float(p[:,0,2].max()),float(p[:,-1,2].min()),float(p[:,-1,2].max())])
        report['bounds_min_m']=np.minimum(report['bounds_min_m'],p.min(axis=(0,1))).tolist();report['bounds_max_m']=np.maximum(report['bounds_max_m'],p.max(axis=(0,1))).tolist()
        if j%100==0:print(f'color {color}, loop {j}/{len(ids)}, seconds {time.time()-start:.1f}',flush=True)
    data.flush();del data
    np.save(outdir/f'colour_{color}_radius_per_loop.npy',np.array(radii,dtype='f4'))
    report['groups'].append({'color':color,'loops':len(ids),'fibre_count':len(ids)*model.FIBRES,'points':len(ids)*model.FIBRES*model.POINTS,'loop_indices':ids.tolist(),'radius_min_m':float(min(radii)),'radius_max_m':float(max(radii)),'fibre_length_min_m':float(min(length)),'fibre_length_max_m':float(max(length)),'stock_length_m':float(sum(length))})
    print(f'color {color} complete',flush=True)
endz=report.pop('endpoints_z_m');report['endpoint_z_range_m']=[min(endz),max(endz)]
report['seconds']=time.time()-start;report['rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024;report['fibres_per_loop']=model.FIBRES;report['samples_per_fibre']=model.POINTS
(ROOT/'receipts/full_construction.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='groups'}),flush=True)
