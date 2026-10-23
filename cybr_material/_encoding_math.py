"""Unchanged historical metric normals and seam diagnostics."""
import numpy as np

def normal_from_height(h,s,n):
    if s.get('tiling')=='finite-cut':
        dy,dx=np.gradient(h)
        dx*=n*s['height_scale_m']/s['tile_m']
        dy*=n*s['height_scale_m']/s.get('tile_y_m',s['tile_m'])
    else:
        dx=(np.roll(h,-1,1)-np.roll(h,1,1))*(n/2)*s['height_scale_m']/s['tile_m']
        dy=(np.roll(h,-1,0)-np.roll(h,1,0))*(n/2)*s['height_scale_m']/s['tile_m']
    v=np.stack((-dx,dy,np.ones_like(dx)),axis=-1)
    v/=np.linalg.norm(v,axis=-1,keepdims=True)
    return .5+.5*v

def seam(a):
    edge=np.concatenate((abs(a[0]-a[-1]).ravel(),abs(a[:,0]-a[:,-1]).ravel()))
    inside=np.concatenate((abs(a[1]-a[0]).ravel(),abs(a[:,1]-a[:,0]).ravel()))
    return dict(wrap_mean=float(edge.mean()),neighbor_mean=float(inside.mean()),wrap_max=float(edge.max()),neighbor_max=float(inside.max()))
