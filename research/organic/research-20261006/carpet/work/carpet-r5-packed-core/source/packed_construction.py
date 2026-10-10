"""Compact constructed yarn, not solved manufacturing mechanics.

Packed disk-cut hexagonal cross-sections follow the frozen r5 loop paths.
Coherent twist uses physical arclength and pauses at the tight crown. Folded
inner crown lobes are removed and their joins are rounded by a fixed physical
arclength fairing. Correlated intrinsic crimp has independent small components.
All fibres remain continuous between buried endpoints. No surrogate core remains.
Residual local overlap is measured and explicitly retained as a construction limit.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d

ROOT=Path(__file__).resolve().parents[1]
PACK=np.array([(q+.5*r,np.sqrt(3)*.5*r) for q in range(-8,9) for r in range(-8,9) if (q+.5*r)**2+(.75*r*r)<=49.000001],dtype='f8')
FIBRES=len(PACK)
POINTS=144
REL_RADIUS=.045
REL_PITCH=.130
TWIST_TURNS=.35

def trim_inner_crown(p,axis,shear=0.,fillet_radius=35e-6):
    """Discard a folded inner parallel-curve lobe, then round its actual join.

    This is a geometric yarn-construction correction inside the selected outline,
    not a relaxed mechanical state or a claim about manufacturing strain.
    """
    from scipy.optimize import brentq
    x=p@axis-shear*p[:,2];z=p[:,2];dx=np.diff(x)
    bad=np.flatnonzero((dx<0)&(np.arange(len(dx))>30)&(np.arange(len(dx))<len(dx)-30))
    if not len(bad):return p
    lo=int(bad[0]);hi=int(bad[-1])+1
    lx=x[:lo+1];lz=z[:lo+1];rx=x[hi:];rz=z[hi:]
    if np.any(np.diff(lx)<=0) or np.any(np.diff(rx)<=0):return p
    low=max(lx[0],rx[0]);high=min(lx[-1],rx[-1])
    if low>=high:return p
    fun=lambda xx:np.interp(xx,lx,lz)-np.interp(xx,rx,rz)
    if fun(low)*fun(high)>0:return p
    xx=brentq(fun,low,high,xtol=1e-14)
    il=int(np.searchsorted(lx,xx));ir=int(np.searchsorted(rx,xx))+hi
    tl=(xx-x[il-1])/(x[il]-x[il-1]);tr=(xx-x[ir-1])/(x[ir]-x[ir-1])
    cross=(p[il-1]*(1-tl)+p[il]*tl+p[ir-1]*(1-tr)+p[ir]*tr)*.5
    # Keep the exact non-folded boundary. A common arclength fairing below rounds
    # this join; independent finite-radius fillets can cross a neighbouring level.
    return np.concatenate((p[:il],cross[None],p[ir:]))

def unit(x):return x/np.maximum(np.linalg.norm(x,axis=-1,keepdims=True),1e-20)
def smooth(t):
    t=np.clip(t,0,1);return t*t*(3-2*t)

def construct(center,radius,index,points=POINTS):
    # Initial samples are dense, followed by each fibre's own arclength sampling.
    t=np.linspace(0,1,513);base_t=np.linspace(0,1,49)
    cs=CubicSpline(base_t,center,axis=0);c=cs(t);tangent=unit(cs(t,1))
    axis=unit(center[-1]-center[0]);axis[2]=0;axis=unit(axis)
    binormal=np.array([-axis[1],axis[0],0.]);up=np.array([0.,0.,1.])
    height=center[24,2]-center[0,2];f=(c[:,2]-center[0,2])/height
    fan=(2*t-1)[:,None]*axis[None]+f[:,None]*up[None]
    perpendicular=fan-tangent*np.sum(fan*tangent,axis=1)[:,None]
    # B is parallel-transported from the loop's plane with a small sway correction.
    b=binormal[None]-tangent*np.sum(binormal[None]*tangent,axis=1)[:,None];b=unit(b)
    normal=unit(np.cross(tangent,b))
    r=float(np.mean(radius));rng=np.random.default_rng(241006+index*104729)
    phase=rng.uniform(0,2*np.pi)
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(c,axis=0),axis=1))]
    # Uniform physical arclength rate away from the crown, easing to zero where
    # inner offsets need the local rounding correction.
    weight=np.maximum(smooth((.30-t)/.08),smooth((t-.70)/.08))
    wds=(weight[1:]+weight[:-1])*.5*np.diff(s);twist_arc=np.r_[0,np.cumsum(wds)]
    twist=phase+2*np.pi*TWIST_TURNS*twist_arc/twist_arc[-1]
    xy=PACK*REL_PITCH*r
    aa=xy[:,0,None]*np.cos(twist)[None]-xy[:,1,None]*np.sin(twist)[None]
    bb=xy[:,0,None]*np.sin(twist)[None]+xy[:,1,None]*np.cos(twist)[None]
    p=c[None]+aa[:,:,None]*normal[None]+bb[:,:,None]*b[None]
    phase_f=rng.uniform(0,2*np.pi,FIBRES)[:,None]
    window=np.sin(np.pi*t)**.6
    common=5e-6*np.sin(2*np.pi*s/.00063+phase)
    individual=1.2e-6*np.sin(2*np.pi*s[None]/.00047+phase_f)
    cr=(common[None]+individual)*window[None]
    # Extend both ends into existing backing. The material is cut below its surface.
    ends=np.stack((p[:,0],p[:,-1]),axis=1).copy()
    endz=rng.uniform(-.00062,-.00048,FIBRES)
    left=ends[:,0].copy();left[:,2]=endz
    right=ends[:,1].copy();right[:,2]=endz-rng.uniform(0,.00004,FIBRES)
    out=np.empty((FIBRES,points,3),dtype='f4')
    for j in range(FIBRES):
        shear=np.dot(center[24]-(center[0]+center[-1])*.5,axis)/height
        path=trim_inner_crown(p[j],axis,shear)
        path=np.concatenate((left[j,None],path,right[j,None]),axis=0)
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        ss=np.linspace(0,arc[-1],points)
        out[j]=np.stack([np.interp(ss,arc,path[:,k]) for k in range(3)],axis=1)
        win=np.maximum(np.sin(np.pi*ss/ss[-1]),0)**.6
        cr=(5e-6*np.sin(2*np.pi*ss/.00063+phase)+1.2e-6*np.sin(2*np.pi*ss/.00047+phase_f[j,0]))*win
        out[j]+=cr[:,None]*binormal[None]
    rr=np.full((FIBRES,points),REL_RADIUS*r,dtype='f4')
    # One bounded arclength fairing removes the discretized fillet join's tiny
    # hooks. It is geometry fairing, with no force, strain, friction or time solve.
    faired=np.empty_like(out)
    for j in range(FIBRES):
        step=float(np.linalg.norm(np.diff(out[j],axis=0),axis=1).mean())
        faired[j]=gaussian_filter1d(out[j],45e-6/max(step,1e-9),axis=0,mode='nearest')
    faired[:,0]=out[:,0];faired[:,-1]=out[:,-1]
    out=faired
    # A few original edge roots sit beyond the nominal board. Below its surface,
    # tuck their body ends inward; the visible loop and inherited halo stay put.
    weight=smooth((-out[:,:,2]-.00010)/.00032)
    target=np.clip(out[:,:,:2],-.02972,.02972)
    out[:,:,:2]+=weight[:,:,None]*(target-out[:,:,:2])
    return out,rr

def segment_distance(a,b,c,d):
    """Robust closest distance of corresponding 3D segment pairs."""
    u=b-a;v=d-c;w=a-c
    aa=np.sum(u*u,axis=-1);bb=np.sum(u*v,axis=-1);cc=np.sum(v*v,axis=-1)
    dd=np.sum(u*w,axis=-1);ee=np.sum(v*w,axis=-1);det=aa*cc-bb*bb
    s=np.divide(bb*ee-cc*dd,det,out=np.zeros_like(det),where=det>1e-32);s=np.clip(s,0,1)
    t=np.clip((bb*s+ee)/np.maximum(cc,1e-32),0,1)
    s=np.clip((bb*t-dd)/np.maximum(aa,1e-32),0,1)
    return np.linalg.norm(w+s[...,None]*u-t[...,None]*v,axis=-1)

def clearance(points,radii,k=20):
    # Segment midpoint neighbours, then actual finite-segment distance. Every
    # tested segment gets a spatial search, not same-ring or 2D spacing alone.
    a=points[:,:-1].reshape(-1,3).astype('f8');b=points[:,1:].reshape(-1,3).astype('f8');mid=(a+b)*.5
    owner=np.repeat(np.arange(len(points)),points.shape[1]-1)
    rad=((radii[:,:-1]+radii[:,1:])*.5).reshape(-1)
    tree=cKDTree(mid);ds,ids=tree.query(mid,k=k,workers=1)
    valid=owner[:,None]!=owner[ids]
    dist=segment_distance(a[:,None],b[:,None],a[ids],b[ids]);gap=dist-rad[:,None]-rad[ids];gap[~valid]=np.inf
    loc=np.unravel_index(np.argmin(gap),gap.shape)
    vv=np.diff(points,axis=1);ll=np.linalg.norm(vv,axis=2);tan=unit(vv)
    bend=np.linalg.norm(np.diff(tan,axis=1),axis=2)/np.maximum((ll[:,:-1]+ll[:,1:])*.5,1e-20)
    closest=gap.min(1)
    split={}
    for name,mask in [('exposed_pile',mid[:,2]>.00013),('backing_or_root_join',mid[:,2]<=.00013)]:
        kmask=points[:,1:-1,2]>.00013 if name=='exposed_pile' else points[:,1:-1,2]<=.00013
        split[name]={'min_gap_um':float(closest[mask].min()*1e6),'negative_segment_fraction':float(np.mean(closest[mask]<0)),'max_curvature_times_radius':float((bend*radii[:,1:-1])[kmask].max())}
    return {'min_gap_um':float(gap[loc]*1e6),'gap_quantiles_um':list(np.quantile(closest,[0,.001,.01,.1,.5,1])*1e6),'negative_segment_fraction':float(np.mean(closest<0)),'max_curvature_times_radius':float(np.max(bend*radii[:,1:-1])),'max_segment_um':float(ll.max()*1e6),'worst_pair':[int(owner[loc[0]]),int(owner[ids[loc]])],'neighbor_k':k,'regions':split}

if __name__=='__main__':
    import json,time
    data=np.load(ROOT/'receipts/baseline_loops.npz');c=data['center'];r=data['radius'];h=c[:,24,2]-c[:,0,2]
    span=np.linalg.norm(c[:,-1]-c[:,0],axis=1)
    risk=(span-2*r[:,24])**2/(.82*np.pi*np.pi*(h-r[:,24]))
    indices=list(dict.fromkeys([int(np.argmin(risk)),int(np.argmax(risk)),int(np.argmin(h)),int(np.argmax(h)),0,400,900,1400]))
    result=[];start=time.time()
    for i in indices:
        p,rr=construct(c[i],r[i],i);rep={'index':i,**clearance(p,rr)};result.append(rep);print(rep,flush=True)
    report={'fibres_per_loop':FIBRES,'points_per_fibre':POINTS,'relative_radius':REL_RADIUS,'relative_pitch':REL_PITCH,'twist_turns':TWIST_TURNS,'tested':result,'seconds':time.time()-start}
    (ROOT/'receipts/packing_preflight.json').write_text(json.dumps(report,indent=2))
