"""Bounded rounded yarn unit: geometric construction, not manufacturing dynamics."""
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
R1=ROOT.parent/'carpet-r5-packed-core'
sys.path.insert(0,str(ROOT.parent/'carpet-r5-migrating-core/source'))
sys.path.insert(0,str(ROOT.parent/'carpet-crown-feasibility/source'))
import migrating_construction as migration
import correct_3d as contact
from measure_guides_v2 import guide
unit=migration.unit
SECTIONS=97
KEYS=144
PASSES=24

def project(x,r,R,rate,curvature,passes):
    cap=.91*R;minimum=2*r+2e-6
    angles=np.arange(6)*np.pi/3
    wraps=.99*R*np.stack((np.cos(angles),np.sin(angles)),axis=1)
    # Curved-guide compression changes the axial tangent component.
    horizon=minimum*np.sqrt(1+(rate*cap/max(1-np.linalg.norm(curvature)*cap,.12))**2)+2e-6
    for _ in range(passes):
        pairs=cKDTree(x).query_pairs(horizon,output_type='ndarray');delta=np.zeros_like(x);weight=np.zeros(len(x))
        if len(pairs):
            i,j=pairs.T;d=x[i]-x[j];middle=(x[i]+x[j])*.5
            v=rate*np.stack((-middle[:,1],middle[:,0]),axis=1)
            axial=1-middle@curvature
            metric=d-v*(np.sum(v*d,axis=1)/(axial*axial+np.sum(v*v,axis=1)))[:,None]
            distance=np.sqrt(np.maximum(np.sum(d*metric,axis=1),1e-24));active=distance<minimum
            i=i[active];j=j[active];metric=metric[active];distance=distance[active]
            if len(i):
                gradient=metric/distance[:,None];square=np.maximum(np.sum(gradient*gradient,axis=1),.02)
                correction=.52*(minimum-distance)[:,None]*gradient/square[:,None]
                length=np.linalg.norm(correction,axis=1);correction*=np.minimum(1,7e-6/np.maximum(length,1e-30))[:,None]
                np.add.at(delta,i,correction);np.add.at(delta,j,-correction);np.add.at(weight,i,1);np.add.at(weight,j,1)
        d=x[:,None]-wraps[None];length=np.linalg.norm(d,axis=-1);depth=r+16e-6+2e-6-length
        correction=np.maximum(depth,0)[:,:,None]*d/np.maximum(length[:,:,None],1e-20)
        delta+=correction.sum(1);weight+=(depth>0).sum(1)
        x+=delta/np.maximum(weight[:,None],1)
        length=np.linalg.norm(x,axis=1);x*=np.minimum(1,cap/np.maximum(length,1e-20))[:,None]
    return x

def resample_with_roots(paths,left_z,right_z):
    out=[]
    for i,path in enumerate(paths):
        first=unit(path[1]-path[0]);last=unit(path[-1]-path[-2])
        if first[2]<=.05 or last[2]>=-.05:raise RuntimeError('Root tangent cannot extend into backing')
        left=path[0]+first*((left_z[i]-path[0,2])/first[2]);right=path[-1]+last*((right_z[i]-path[-1,2])/last[2])
        p=np.concatenate((left[None],path,right[None]));s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
        sample=np.linspace(0,s[-1],KEYS);q=np.stack([np.interp(sample,s,p[:,j]) for j in range(3)],axis=1);q[0]=left;q[-1]=right;out.append(q)
    return np.asarray(out,dtype='f4')

def construct(index,data,original):
    c=data['center'][index];R=float(data['radius'][index].mean());r=.045*R
    dense,description=guide(c,R);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(dense,axis=0),axis=1))]
    spline=CubicSpline(arc,dense,axis=0);s=np.linspace(0,arc[-1],SECTIONS);center=spline(s);tangent=unit(spline(s,1))
    axis=unit(c[-1]-c[0]);normal=unit(np.cross(axis,c[24]-(c[0]+c[-1])*.5));in_plane=unit(np.cross(normal[None],tangent))
    old_a=data['frame'][index,0,0];alpha=np.arctan2(old_a@normal,old_a@in_plane[0])
    a=in_plane*np.cos(alpha)+normal*np.sin(alpha);b=unit(np.cross(tangent,a))
    phase0=float(data['phase'][index,0]);turns=float((data['phase'][index,-1]-phase0)/(2*np.pi));rate=2*np.pi*turns/arc[-1];phase=phase0+rate*s
    e1=np.cos(phase)[:,None]*a+np.sin(phase)[:,None]*b;e2=-np.sin(phase)[:,None]*a+np.cos(phase)[:,None]*b
    velocity=spline(s,1);acceleration=spline(s,2);speed=np.linalg.norm(velocity,axis=1)
    bend=(acceleration-tangent*np.sum(acceleration*tangent,axis=1)[:,None])/speed[:,None]**2
    curvature=np.stack((np.sum(bend*e1,axis=1),np.sum(bend*e2,axis=1)),axis=1)
    rng=np.random.default_rng(241016+index*104729);seed=migration.old.PACK*.13*R+rng.normal(0,19e-6,(187,2))
    seed=project(seed,r,R,rate,curvature[0],48)
    amplitude=rng.uniform(30e-6,60e-6,(187,2));wavelength=rng.uniform(.0006,.0012,(187,2));offset=rng.uniform(0,2*np.pi,(187,2))
    window=np.sin(np.pi*s/s[-1])**2
    target=seed[None]+window[:,None,None]*amplitude[None]*(np.sin(2*np.pi*s[:,None,None]/wavelength[None]+offset[None])-np.sin(offset)[None])*.5
    xy=np.empty_like(target);xy[0]=seed
    for k in range(1,SECTIONS):xy[k]=project(xy[k-1]+target[k]-target[k-1],r,R,rate,curvature[k],PASSES)
    section=center[:,None]+xy[:,:,0,None]*e1[:,None]+xy[:,:,1,None]*e2[:,None]
    sample=np.linspace(0,s[-1],1025);body=CubicSpline(s,section,axis=0)(sample).transpose(1,0,2)
    p=resample_with_roots(body,original[:,0,2],original[:,-1,2])
    # Six wraps share the guide, physical-arclength phase and radius budget.
    cc=spline(sample);tt=unit(spline(sample,1));aa=unit(np.cross(normal[None],tt))*np.cos(alpha)+normal*np.sin(alpha);bb=unit(np.cross(tt,aa))
    theta=phase0+rate*sample[None]+np.arange(6)[:,None]*np.pi/3
    wraps=cc[None]+.99*R*(np.cos(theta)[:,:,None]*aa[None]+np.sin(theta)[:,:,None]*bb[None])
    wraps=resample_with_roots(wraps,np.full(6,float(np.median(original[:,0,2]))),np.full(6,float(np.median(original[:,-1,2]))))
    compact=spline(np.linspace(0,arc[-1],257));radial=np.linalg.norm(xy,axis=-1)
    report={'loop':index,**description,'roots_exact':bool(np.array_equal(dense[[0,-1]],c[[0,-1]])),'apex_exact':bool(np.array_equal(dense[1024],c[24])),'body_fibres':187,'body_radius_m':r,'wraps':6,'wrap_radius_m':16e-6,'phase_turns':turns,'phase_rate_rad_per_m':rate,'guide_length_m':float(arc[-1]),'sections':SECTIONS,'section_projection_passes':PASSES,'seed_passes':48,'radial_excursion_median_um':float(np.median(radial.max(0)-radial.min(0))*1e6),'native_keys':KEYS,'new_guide_envelope_body_limit_m':R+2e-6,'construction':'Orthogonal circular crown, transported in-plane frame, shared physical-arclength twist, two-axis independent migration, curvature-aware section spacing. Geometric construction, not manufacturing simulation.'}
    return p,wraps,compact,dense,r,R,report

def confine(points,reference,center,r,R):
    out=points.copy();flat=out.reshape(-1,3);ids=np.flatnonzero(flat[:,2]>.00013)
    a=center[:-1];v=np.diff(center,axis=0);square=np.sum(v*v,axis=1)
    for start in range(0,len(ids),512):
        selected=ids[start:start+512];q=flat[selected];t=np.clip(np.sum((q[:,None]-a[None])*v[None],axis=2)/square[None],0,1)
        nearest=a[None]+t[:,:,None]*v[None];dist=np.linalg.norm(q[:,None]-nearest,axis=2);which=np.argmin(dist,axis=1);nearest=nearest[np.arange(len(q)),which];d=q-nearest;length=np.linalg.norm(d,axis=1)
        flat[selected]=nearest+d*np.minimum(1,(R-r)/np.maximum(length,1e-20))[:,None]
    displacement=out-reference;norm=np.linalg.norm(displacement,axis=2);out=reference+displacement*np.minimum(1,contact.TOTAL/np.maximum(norm,1e-20))[:,:,None]
    buried=reference[:,:,2]<.00007;out[:,:,2][buried]=np.minimum(out[:,:,2][buried],.00008)
    out[:,:2]=reference[:,:2];out[:,-2:]=reference[:,-2:]
    return out

def correct(p,wraps,center,r,R):
    # Same correction budget as the rejected preceding experiment; only the
    # guide/wrap unit and curvature-aware section field changed here.
    contact.confine=confine
    return contact.correct(p,r,center,R,migration.evaluate_catmull(wraps,3))
