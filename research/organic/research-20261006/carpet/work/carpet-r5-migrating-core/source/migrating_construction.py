"""Fixed-budget migrating packing test; geometric constraints, not mechanics."""
from pathlib import Path
import sys,math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
R1=ROOT.parent/'carpet-r5-packed-core'
sys.path.insert(0,str(R1/'revisions/r1_frozen/source'))
import packed_construction as old
SECTIONS=81
PASSES=24
SEED_PASSES=48
KEYS=144

def unit(a):return a/np.maximum(np.linalg.norm(a,axis=-1,keepdims=True),1e-30)
def project(x,radius,core_radius,phase_rate,passes,wrapper_relative_phase=0.):
    """Bounded sectionwise separation in a local helical metric and envelope."""
    cap=.91*core_radius;minimum=2*radius+2e-6
    wa=np.arange(6)*np.pi/3+wrapper_relative_phase
    wraps=.99*core_radius*np.stack((np.cos(wa),np.sin(wa)),axis=1)
    horizon=minimum*np.sqrt(1+(phase_rate*cap)**2)+2e-6
    for _ in range(passes):
        pairs=cKDTree(x).query_pairs(horizon,output_type='ndarray')
        delta=np.zeros_like(x);weight=np.zeros(len(x))
        if len(pairs):
            i=pairs[:,0];j=pairs[:,1];d=x[i]-x[j];middle=(x[i]+x[j])*.5
            v=phase_rate*np.stack((-middle[:,1],middle[:,0]),axis=1)
            vd=np.sum(v*d,axis=1);vv=1+np.sum(v*v,axis=1)
            metric=d-v*(vd/vv)[:,None]
            distance=np.sqrt(np.maximum(np.sum(d*metric,axis=1),1e-24))
            active=distance<minimum
            i=i[active];j=j[active];metric=metric[active];distance=distance[active]
            if len(i):
                g=metric/distance[:,None];gg=np.maximum(np.sum(g*g,axis=1),.02)
                correction=.52*(minimum-distance)[:,None]*g/gg[:,None]
                length=np.linalg.norm(correction,axis=1);correction*=np.minimum(1,7e-6/np.maximum(length,1e-30))[:,None]
                np.add.at(delta,i,correction);np.add.at(delta,j,-correction);np.add.at(weight,i,1);np.add.at(weight,j,1)
        # These are conservative sectionwise exclusions of the retained wraps.
        d=x[:,None]-wraps[None];length=np.linalg.norm(d,axis=-1);over=radius+16e-6+2e-6-length
        correction=np.maximum(over,0)[:,:,None]*d/np.maximum(length[:,:,None],1e-20)
        delta+=correction.sum(1);weight+=(over>0).sum(1)
        x+=delta/np.maximum(weight[:,None],1)
        norm=np.linalg.norm(x,axis=1);x*=np.minimum(1,cap/np.maximum(norm,1e-20))[:,None]
    return x

def construct(index,data,r1_points,phase_mode='exact',fairing_m=0.,endpoint_mode='continuous_new_roots'):
    center=data['center'][index].astype('f8');R=float(data['radius'][index].mean());r=.045*R
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(center,axis=0),axis=1))]
    s=np.linspace(0,arc[-1],SECTIONS);tt=np.interp(s,arc,np.linspace(0,1,49))
    cs=CubicSpline(np.linspace(0,1,49),center,axis=0);c=cs(tt);tangent=unit(cs(tt,1))
    frame=data['frame'][index]
    a=np.stack([np.interp(s,arc,frame[:,0,k]) for k in range(3)],axis=1)
    a=unit(a-tangent*np.sum(a*tangent,axis=1)[:,None]);b=unit(np.cross(tangent,a))
    inherited_phase=np.interp(s,arc,data['phase'][index])
    if phase_mode=='exact':phase=inherited_phase.copy()
    else:phase=inherited_phase[0]+(inherited_phase[-1]-inherited_phase[0])*s/s[-1]
    rate=np.gradient(phase,s)
    rho=np.linalg.norm(old.PACK*.13*R,axis=1)
    projected_fraction=.045**2*np.sqrt(1+(rate[:,None]*rho[None])**2).sum(1)
    rng=np.random.default_rng(241016+index*104729)
    seed=old.PACK*.13*R+rng.normal(0,19e-6,(187,2))
    seed=project(seed,r,R,float(rate[0]),SEED_PASSES,float(inherited_phase[0]-phase[0]))
    amplitude=rng.uniform(30e-6,60e-6,(187,2));wavelength=rng.uniform(.0006,.0012,(187,2));offset=rng.uniform(0,2*np.pi,(187,2))
    window=np.sin(np.pi*s/s[-1])**2
    targets=seed[None]+window[:,None,None]*amplitude[None]*(np.sin(2*np.pi*s[:,None,None]/wavelength[None]+offset[None])-np.sin(offset)[None])*.5
    xy=np.empty_like(targets);xy[0]=seed
    for k in range(1,SECTIONS):
        trial=xy[k-1]+targets[k]-targets[k-1]
        xy[k]=project(trial,r,R,float(rate[k]),PASSES,float(inherited_phase[k]-phase[k]))
    rotated=np.stack((xy[:,:,0]*np.cos(phase)[:,None]-xy[:,:,1]*np.sin(phase)[:,None],xy[:,:,0]*np.sin(phase)[:,None]+xy[:,:,1]*np.cos(phase)[:,None]),axis=-1)
    section_points=c[:,None]+rotated[:,:,0,None]*a[:,None]+rotated[:,:,1,None]*b[:,None]
    dense_s=np.linspace(0,s[-1],513);dense=CubicSpline(s,section_points,axis=0)(dense_s).transpose(1,0,2)
    out=np.empty((187,KEYS,3),dtype='f4')
    for i in range(187):
        if endpoint_mode=='r1_endpoints':left=r1_points[i,0].copy();right=r1_points[i,-1].copy()
        else:
            first=unit(dense[i,1]-dense[i,0]);last=unit(dense[i,-1]-dense[i,-2])
            if first[2]<=.05 or last[2]>=-.05:raise RuntimeError('New root tangent cannot extend continuously into backing')
            left=dense[i,0]+first*((float(r1_points[i,0,2])-dense[i,0,2])/first[2])
            right=dense[i,-1]+last*((float(r1_points[i,-1,2])-dense[i,-1,2])/last[2])
        path=np.concatenate((left[None],dense[i],right[None]),axis=0)
        length=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        ss=np.linspace(0,length[-1],KEYS);p=np.stack([np.interp(ss,length,path[:,j]) for j in range(3)],axis=1)
        if fairing_m:p=gaussian_filter1d(p,fairing_m/(length[-1]/(KEYS-1)),axis=0,mode='nearest')
        p[0]=left;p[-1]=right
        out[i]=p
    ends=out[:,[0,-1]]
    assert np.all(abs(ends[:,:,:2])<.0298) and np.all((ends[:,:,2]<-.0004)&(ends[:,:,2]>-.0009))
    radial=np.linalg.norm(xy,axis=-1);rank=np.argsort(np.argsort(radial,axis=1),axis=1)
    mobile=(radial.max(0)-radial.min(0))>26e-6
    report={'loop_index':int(index),'fibres':187,'radius_m':r,'sections':SECTIONS,'projection_passes_per_section':PASSES,'initial_projection_passes':SEED_PASSES,'phase_turns':float((phase[-1]-phase[0])/(2*np.pi)),'phase_mode':phase_mode,'phase_source':'Measured r5 total turns and endpoint orientation, in exact original transported frame. Arclength mode redistributes phase uniformly over physical length; retained wrappers keep their exact inherited phase. No0.35-turn fallback.','max_phase_rate_rad_per_m':float(abs(rate).max()),'maximum_projected_area_fraction_estimate':float(projected_fraction.max()),'amplitude_controls_m':[30e-6,60e-6],'wavelength_controls_m':[.0006,.0012],'rank_range_median':float(np.median(rank.max(0)-rank.min(0))),'radial_excursion_median_um':float(np.median(radial.max(0)-radial.min(0))*1e6),'fibres_migrating_more_than_one_diameter':int(mobile.sum()),'fibre_ends_bitwise_equal':True,'construction_limit':'Fixed-budget sectionwise projection. Neither manufacture dynamics nor a contact/force solve. Final three-dimensional checks determine acceptance.'}
    report.update({'fairing_m':fairing_m,'endpoint_mode':endpoint_mode,'fibre_ends_bitwise_equal':bool(np.array_equal(ends,r1_points[:,[0,-1]])),'all_new_ends_buried':True,'endpoint_bounds_min_m':ends.min(axis=(0,1)).tolist(),'endpoint_bounds_max_m':ends.max(axis=(0,1)).tolist()})
    return out,np.full((187,KEYS),r,dtype='f4'),report

def evaluate_catmull(p,factor=3):
    prev=np.concatenate((p[:,:1],p[:,:-1]),axis=1);nex=np.concatenate((p[:,1:],p[:,-1:]),axis=1);tangent=(nex-prev)*.5
    t=np.arange(factor)/factor;t=t[None,None,:,None]
    q=(2*t**3-3*t*t+1)*p[:,:-1,None]+(t**3-2*t*t+t)*tangent[:,:-1,None]+(-2*t**3+3*t*t)*p[:,1:,None]+(t**3-t*t)*tangent[:,1:,None]
    return np.concatenate((q.reshape(len(p),-1,3),p[:,-1:]),axis=1).astype('f4')

def retained_proximity(p,r,data,index):
    center=data['center'][index];R=data['radius'][index].mean();a=data['frame'][index,:,0];b=data['frame'][index,:,1]
    theta=data['phase'][index][None]+np.arange(6)[:,None]*np.pi/3
    wraps=center[None]+.99*R*(np.cos(theta)[:,:,None]*a[None]+np.sin(theta)[:,:,None]*b[None])
    c=wraps[:,:-1].reshape(-1,3);d=wraps[:,1:].reshape(-1,3)
    pa=p[:,:-1].reshape(-1,3).astype('f8');pb=p[:,1:].reshape(-1,3).astype('f8');mid=(pa+pb)*.5
    ids=cKDTree((c+d)*.5).query(mid,k=12,workers=1)[1]
    dist=old.segment_distance(pa[:,None],pb[:,None],c[ids],d[ids]).min(1)
    exposed=mid[:,2]>.00013;gap=dist-r-16e-6
    return {'exposed_min_conservative_gap_um':float(gap[exposed].min()*1e6),'exposed_potential_penetration_fraction':float(np.mean(gap[exposed]<-1e-6)),'exposed_contact_band_fraction':float(np.mean(abs(gap[exposed])<=1e-6)),'method':'3D finite body segments against all six exact-phase recovered retained strand centerlines, with16µm circumscribed tube radius. Negative gaps are conservative penetration candidates, not exact polygon penetration depths.','retained_geometry_changed':False}
