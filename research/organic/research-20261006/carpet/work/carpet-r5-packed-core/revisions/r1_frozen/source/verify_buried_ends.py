"""Apply the final below-backing edge tuck to staged arrays and verify all ends."""
from pathlib import Path
import numpy as np,json
from packed_construction import smooth
ROOT=Path(__file__).resolve().parents[1]
count=0;tucked=0;lo=np.full(3,np.inf);hi=np.full(3,-np.inf)
for color in range(4):
    data=np.load(ROOT/f'arrays/colour_{color}_positions.npy',mmap_mode='r+')
    for i in range(len(data)):
        p=data[i];ends=p[:,[0,-1],:];tucked+=int(np.any(np.abs(ends[:,:,:2])>.02972,axis=2).sum())
        weight=smooth((-p[:,:,2]-.00010)/.00032);target=np.clip(p[:,:,:2],-.02972,.02972);p[:,:,:2]+=weight[:,:,None]*(target-p[:,:,:2])
        ends=p[:,[0,-1],:];assert np.all(abs(ends[:,:,:2])<.029721);assert np.all((ends[:,:,2]<-.00047)&(ends[:,:,2]>-.0009))
        lo=np.minimum(lo,ends.min(axis=(0,1)));hi=np.maximum(hi,ends.max(axis=(0,1)));count+=ends.shape[0]*2
    data.flush();del data
report={'all_body_endpoints_verified':count,'edge_endpoints_tucked':tucked,'endpoints_min_m':lo.tolist(),'endpoints_max_m':hi.tolist(),'backing_nominal_bounds_m':[[-.03,-.03,-.001],[.03,.03,.0001]],'method':'Below-surface smooth inward tuck to ±29.72mm. All visible loop anchors remain unchanged; no new crown ends.'}
(ROOT/'receipts/all_endpoints_buried.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
