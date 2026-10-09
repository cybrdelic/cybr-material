"""Original physical-scale seated mineral fragments; no image input.

Deterministic jittered cells own angular finite fragments. Height, optical
exposure and mineral identity are outputs of the same geometric construction.
This is an authored reduced surface model, not a particle-packing simulation.
"""
import numpy as np

def hash01(x,y,seed):
    # uint64 arithmetic avoids signed overflow; repeat cells across tile edges.
    v=(x.astype('u8')*np.uint64(374761393)+y.astype('u8')*np.uint64(668265263)+np.uint64(seed))
    v=(v^(v>>np.uint64(13)))*np.uint64(1274126177)
    v=v^(v>>np.uint64(16))
    return (v&np.uint64(0xffffff)).astype('f4')/16777216.

def bonded_fragment_field(n,tile,seed):
    if not isinstance(n,int) or n<4 or not np.isfinite(tile) or tile<=0 or not isinstance(seed,int) or seed<0:
        raise ValueError('Expected integer resolution>=4, finite positive tile and nonnegative seed')
    cells=round(tile/.00034)
    if cells<3:raise ValueError('Tile must span at least three fragment cells')
    pitch=tile/cells
    x=(np.arange(n,dtype='f4')+.5)/n*tile
    ix=np.floor(x/pitch).astype('i4');iy=ix.copy();fx=x[None,:];fy=x[:,None]
    height=np.zeros((n,n),'f4');coverage=height.copy();ids=np.zeros((n,n),'u1');radius_sum=0.
    palette=np.array([[.245,.253,.25],[.37,.36,.325],[.165,.18,.178],[.47,.458,.403],[.285,.30,.29],[.52,.515,.47]],'f4')
    for oy in [-1,0,1]:
      cy=iy[:,None]+oy
      for ox in [-1,0,1]:
        cx=ix[None,:]+ox;a=hash01(cx%cells,cy%cells,seed);b=hash01(cx%cells,cy%cells,seed+41);c=hash01(cx%cells,cy%cells,seed+79)
        centerx=(cx+.5+(a-.5)*.56)*pitch;centery=(cy+.5+(b-.5)*.56)*pitch
        dx=fx-centerx;dy=fy-centery;angle=c*np.float32(6.283185307179586);co=np.cos(angle);si=np.sin(angle);u=dx*co+dy*si;v=(-dx*si+dy*co)/(.72+.25*b)
        radius=.000105+.000068*a
        rho=np.maximum.reduce([np.abs(u),np.abs(.5*u+.8660254*v),np.abs(.5*u-.8660254*v)])/radius
        edge=(1-rho)*radius;cover=np.clip(edge/(tile/n)+.5,0,1)
        peak=.000060+.000055*c
        # Truncated angular roof: finite bevels and tilted planar faces.
        plane=peak+u*(a-.5)*.16+v*(b-.5)*.16
        h=np.maximum(0,np.minimum(plane,edge*.85))*cover
        winner=h>height;height=np.maximum(height,h);coverage=np.maximum(coverage,cover);ids=np.where(winner,np.minimum((c*6).astype('u1'),5),ids)
    rgb=palette[ids];mean=float(height.mean())
    return {'height_m':height,'coverage':coverage,'mineral_rgb':rgb,'mean_height_m':mean,'report':{'cell_pitch_m':pitch,'cells_per_axis':cells,'nominal_apothem_m':[.000105,.000173],'nominal_peak_m':[.000060,.000115],'coverage_fraction':float((coverage>.5).mean()),'mean_height_m':mean,'maximum_height_m':float(height.max()),'seed':seed,'state':'Shared geometric fragment height, exposure coverage and mineral identity','packing_qualified':False}}
