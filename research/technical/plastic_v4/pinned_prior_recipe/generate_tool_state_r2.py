"""Native4096 physical-scale rounded tool impressions; no straight cell tessellation.
Empirical blast/etch tool-relief and replication approximation; no particle/melt/heat solver.
"""
import numpy as np,json,hashlib
from pathlib import Path
from scipy.ndimage import gaussian_filter
from PIL import Image
R=Path('/workspace/shared/material-native4k/plastic-detail');D=R/'maps'/'r2';D.mkdir(exist_ok=True)
N=4096;WIDTH=.08;PIX=WIDTH/N;rng=np.random.default_rng(21052027)
tool=np.zeros((N,N),np.float32)
count=95000
for k in range(count):
 cx,cy=rng.uniform(-.0003,WIDTH+.0003,2);radius=rng.uniform(.00012,.000245);depth=rng.uniform(.000020,.000044)
 ix0=max(0,int((cx-radius)/PIX));ix1=min(N,int((cx+radius)/PIX)+2);iy0=max(0,int((cy-radius)/PIX));iy1=min(N,int((cy+radius)/PIX)+2)
 if ix0>=ix1 or iy0>=iy1:continue
 x=((np.arange(ix0,ix1,dtype=np.float32)+.5)*PIX-cx)/radius;y=((np.arange(iy0,iy1,dtype=np.float32)+.5)*PIX-cy)/radius
 rr=x[None,:]**2+y[:,None]**2
 # Rounded finite-footprint impression. Overlapping tool removal uses the depth envelope,
 # followed by light etch/polish rounding, rather than arbitrary image noise or cell borders.
 cap=depth*np.maximum(0,1-rr)**1.75
 patch=tool[iy0:iy1,ix0:ix1];np.maximum(patch,cap,out=patch)
tool=gaussian_filter(tool,1.05,mode='reflect')
replication=.90;height=(tool-float(tool.mean()))*replication
state=tool/float(tool.max())
rough=np.clip(.49+.030*(float(tool.mean())-tool)/float(tool.std()),.445,.535)
np.save(D/'tool_state_native4096.npy',state);np.save(D/'replicated_relief_native4096_m.npy',height)
Image.fromarray(np.rint(np.flipud(state)*65535).astype(np.uint16)).save(D/'replicated_relief_native4096.png')
Image.fromarray(np.rint(np.flipud(rough)*65535).astype(np.uint16)).save(D/'tool_roughness_native4096.png')
gy,gx=np.gradient(height,PIX);sl=np.sqrt(gx*gx+gy*gy)
r={'revision':'r2','native_size':[N,N],'extent_m':[WIDTH,WIDTH],'texel_pitch_m':PIX,'seed':21052027,'method':'95000 overlapping rounded finite-footprint tool removal impressions; max-depth envelope followed by20.5um Gaussian etch/polish rounding. Empirical geometric process model, not particle/flow/heat simulation.','impression_count':count,'impression_radius_range_m':[.00012,.000245],'tool_impression_depth_range_m':[.000020,.000044],'replication_ratio':replication,'replicated_relief_peak_to_valley_m':float(height.max()-height.min()),'relief_min_m':float(height.min()),'relief_max_m':float(height.max()),'roughness_range':[float(rough.min()),float(rough.max())],'roughness_mean':float(rough.mean()),'surface_slope_percentiles':[float(v) for v in np.percentile(sl,[0,50,95,99,100])],'no_upsampling':True,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.iterdir() if p.is_file()}}
(R/'receipts'/'native_tool_state_r2.json').write_text(json.dumps(r,indent=2)); print(json.dumps(r),flush=True)
