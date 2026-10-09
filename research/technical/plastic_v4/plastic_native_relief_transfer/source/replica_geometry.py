"""Reversible r2 geometry adapter; no Blender import and no random geometry."""
import numpy as np
WIDTH=.08

def sample_native(field,xy):
    ny,nx=field.shape
    q=np.clip((np.asarray(xy)/WIDTH+.5)*[nx,ny]-.5,[0.,0.],[nx-1,ny-1])
    ix=np.floor(q[:,0]).astype(int);iy=np.floor(q[:,1]).astype(int)
    fx=q[:,0]-ix;fy=q[:,1]-iy
    return (field[iy,ix]*(1-fx)*(1-fy)+field[iy,np.minimum(ix+1,nx-1)]*fx*(1-fy)+field[np.minimum(iy+1,ny-1),ix]*(1-fx)*fy+field[np.minimum(iy+1,ny-1),np.minimum(ix+1,nx-1)]*fx*fy)

def construction_normal(base):
    q=np.asarray(base,dtype=np.float64).copy()
    q[:,:2]/=(1-np.tan(np.deg2rad(.8))*np.abs(q[:,2]-.004)/.04)[:,None]
    q[:,2]-=.004
    core=np.clip(q,[-.037,-.037,-.001],[.037,.037,.001])
    dv=q-core;length=np.linalg.norm(dv,axis=1)
    return np.divide(dv,length[:,None],out=np.zeros_like(dv),where=length[:,None]>1e-12)

def remove_r2_relief(original,old_height,iterations=8):
    """Invert precisely the documented r2 normal*height*nz**6 displacement.

    No points below the upper rounded rim move. Side parting geometry is left
    intact. Return base points, construction normals, and reconstruction error.
    """
    original=np.asarray(original,dtype=np.float64)
    base=original.copy();affected=original[:,2]>.005-5e-5
    pts=base[affected].copy();source=original[affected]
    for _ in range(iterations):
        n=construction_normal(pts);amp=np.maximum(n[:,2],0.)**6
        h=sample_native(old_height,pts[:,:2])
        pts=source-n*(h*amp)[:,None]
    n=construction_normal(pts);amp=np.maximum(n[:,2],0.)**6
    reconstructed=pts+n*(sample_native(old_height,pts[:,:2])*amp)[:,None]
    error=float(np.max(np.linalg.norm(reconstructed-source,axis=1))) if len(pts) else 0.
    base[affected]=pts
    normals=np.zeros_like(base);normals[affected]=n
    return base,normals,error

def replace_relief(original,old_height,new_height):
    base,n,error=remove_r2_relief(original,old_height)
    amp=np.maximum(n[:,2],0.)**6
    new=base+n*(sample_native(new_height,base[:,:2])*amp)[:,None]
    return new,error
