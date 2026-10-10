"""Reconstructed, unexecuted 3D constraint proposal after executor reset.

The pre-reset write/run was unconfirmed. This is geometric construction, never
force/friction/time integration or a manufacturing calibration.
"""
import numpy as np
from scipy.spatial import cKDTree
from migrating_construction import evaluate_catmull

PASSES=16
PER_PASS=2.5e-6
TOTAL=30e-6
CONTACT_MARGIN=2e-6
EVAL_FACTOR=2

def closest(a,b,c,d):
    u=b-a;v=d-c;w=a-c
    aa=np.sum(u*u,axis=-1);bb=np.sum(u*v,axis=-1);cc=np.sum(v*v,axis=-1);dd=np.sum(u*w,axis=-1);ee=np.sum(v*w,axis=-1);den=aa*cc-bb*bb
    s=np.divide(bb*ee-cc*dd,den,out=np.zeros_like(den),where=den>1e-30);s=np.clip(s,0,1)
    t=np.clip((bb*s+ee)/np.maximum(cc,1e-30),0,1);s=np.clip((bb*t-dd)/np.maximum(aa,1e-30),0,1);t=np.clip((bb*s+ee)/np.maximum(cc,1e-30),0,1)
    p=a+s[...,None]*u;q=c+t[...,None]*v
    return p,q,s,t

def pair_correction(a,b,c,d,minimum,fixed=False):
    p,q,s,t=closest(a,b,c,d);v=p-q;distance=np.linalg.norm(v,axis=-1);normal=v/np.maximum(distance[...,None],1e-20);depth=np.maximum(minimum-distance,0)
    scale=1. if fixed else .5
    return normal*(depth*scale)[...,None],np.zeros_like(normal) if fixed else -normal*(depth*.5)[...,None]

def body_contacts(p,r):
    a=p[:,:-1].reshape(-1,3).astype('f8');b=p[:,1:].reshape(-1,3).astype('f8');mid=(a+b)*.5;segments=p.shape[1]-1
    near=cKDTree(mid).query(mid,k=20,workers=1)[1]
    first=np.repeat(np.arange(len(mid)),20);second=near.ravel()
    valid=(first<second)&(first//segments!=second//segments)
    pairs=np.unique(np.stack((first[valid],second[valid]),axis=1),axis=0)
    i=pairs[:,0];j=pairs[:,1];x,y,s,t=closest(a[i],b[i],a[j],b[j]);difference=x-y;distance=np.linalg.norm(difference,axis=1)
    active=distance<2*r+CONTACT_MARGIN
    return i[active],j[active],s[active],t[active],difference[active]/np.maximum(distance[active,None],1e-20),(2*r+CONTACT_MARGIN-distance[active])

def obstacle_contacts(p,r,c,d):
    a=p[:,:-1].reshape(-1,3).astype('f8');b=p[:,1:].reshape(-1,3).astype('f8');mid=(a+b)*.5
    near=cKDTree((c+d)*.5).query(mid,k=12,workers=1)[1];i=np.repeat(np.arange(len(mid)),12);j=near.ravel()
    x,y,s,t=closest(a[i],b[i],c[j],d[j]);difference=x-y;distance=np.linalg.norm(difference,axis=1)
    active=distance<r+16e-6+CONTACT_MARGIN
    return i[active],s[active],difference[active]/np.maximum(distance[active,None],1e-20),(r+16e-6+CONTACT_MARGIN-distance[active])

def weights(segment,parameter,keys):
    eval_segments=(keys-1)*EVAL_FACTOR;owner=segment//eval_segments
    value=(segment%eval_segments+parameter)/EVAL_FACTOR;index=np.minimum(value.astype('i4'),keys-2);t=value-index
    ids=owner[:,None]*keys+np.clip(index[:,None]+np.arange(-1,3)[None],0,keys-1)
    basis=np.stack((-.5*t+t*t-.5*t**3,1-2.5*t*t+1.5*t**3,.5*t+2*t*t-1.5*t**3,-.5*t*t+.5*t**3),axis=1)
    return ids,basis

def add(delta,weight,ids,basis,move):
    for k in range(4):
        np.add.at(delta,ids[:,k],basis[:,k,None]*move)
        np.add.at(weight,ids[:,k],abs(basis[:,k]))

def confine(points,reference,center,r,R):
    out=points.copy();flat=out.reshape(-1,3);exposed=flat[:,2]>.00013
    q=flat[exposed];a=center[:-1];v=np.diff(center,axis=0);vv=np.sum(v*v,axis=1)
    t=np.clip(np.sum((q[:,None]-a[None])*v[None],axis=2)/vv[None],0,1);nearest=a[None]+t[:,:,None]*v[None]
    difference=q[:,None]-nearest;dist=np.linalg.norm(difference,axis=2);ids=np.argmin(dist,axis=1);n=nearest[np.arange(len(q)),ids];d=q-n;length=np.linalg.norm(d,axis=1)
    flat[exposed]=n+d*np.minimum(1,(R-r)/np.maximum(length,1e-20))[:,None]
    displacement=out-reference;norm=np.linalg.norm(displacement,axis=2);out=reference+displacement*np.minimum(1,TOTAL/np.maximum(norm,1e-20))[:,:,None]
    buried=reference[:,:,2]<.00007;out[:,:,2][buried]=np.minimum(out[:,:,2][buried],.00008)
    out[:,:2]=reference[:,:2];out[:,-2:]=reference[:,-2:]
    return out

def correct(points,r,center,R,wraps):
    reference=points.astype('f8');out=reference.copy();keys=out.shape[1];history=[]
    c=wraps[:,:-1].reshape(-1,3).copy();d=wraps[:,1:].reshape(-1,3).copy();c0=c.copy();d0=d.copy()
    for step in range(PASSES):
        evaluated=evaluate_catmull(out,EVAL_FACTOR);delta=np.zeros((out.shape[0]*keys,3));weight=np.zeros(len(delta))
        i,j,s,t,normal,depth=body_contacts(evaluated,r)
        ia,ba=weights(i,s,keys);ib,bb=weights(j,t,keys);den=np.maximum(np.sum(ba*ba+bb*bb,axis=1),.2);move=.85*normal*(depth/den)[:,None]
        add(delta,weight,ia,ba,move);add(delta,weight,ib,bb,-move)
        body_count=len(i);body_depth=float(depth.max()) if len(depth) else 0
        del i,j,s,t,normal,depth,ia,ba,ib,bb,den,move
        i,s,normal,depth=obstacle_contacts(evaluated,r,c,d);ids,basis=weights(i,s,keys);den=np.maximum(np.sum(basis*basis,axis=1),.2);move=.85*normal*(depth/den)[:,None]
        add(delta,weight,ids,basis,move);obstacle_count=len(i);obstacle_depth=float(depth.max()) if len(depth) else 0
        delta=(delta/np.maximum(weight[:,None],1)).reshape(out.shape)
        # Regularization is inside each contact pass; no post-contact smoothing.
        smoothed=delta.copy();smoothed[:,1:-1]=.15*delta[:,:-2]+.70*delta[:,1:-1]+.15*delta[:,2:]
        length=np.linalg.norm(smoothed,axis=2);smoothed*=np.minimum(1,PER_PASS/np.maximum(length,1e-20))[:,:,None]
        proposed=confine(out+smoothed,reference,center,r,R)
        # The envelope projection must not silently exceed the per-pass cap.
        actual_step=proposed-out;step_length=np.linalg.norm(actual_step,axis=2)
        out+=actual_step*np.minimum(1,PER_PASS/np.maximum(step_length,1e-20))[:,:,None]
        out[:,:2]=reference[:,:2];out[:,-2:]=reference[:,-2:]
        history.append({'pass':step+1,'body_contact_constraints':body_count,'wrap_contact_constraints':obstacle_count,'max_body_deficit_um':body_depth*1e6,'max_wrap_deficit_um':obstacle_depth*1e6,'maximum_control_displacement_um':float(np.linalg.norm(out-reference,axis=2).max()*1e6)})
    assert np.array_equal(c,c0) and np.array_equal(d,d0)
    assert np.array_equal(out[:,[0,-1]],reference[:,[0,-1]])
    assert np.linalg.norm(out-reference,axis=2).max()<=TOTAL+1e-12
    return out.astype('f4'),{'passes':PASSES,'per_pass_control_cap_um':PER_PASS*1e6,'total_control_cap_um':TOTAL*1e6,'contact_margin_um':CONTACT_MARGIN*1e6,'projection_evaluation_factor':EVAL_FACTOR,'maximum_control_displacement_um':float(np.linalg.norm(out-reference,axis=2).max()*1e6),'loop_centroid_drift_um':float(np.linalg.norm(out.mean((0,1))-reference.mean((0,1)))*1e6),'fixed_wrap_nodes_unchanged':True,'buried_endpoints_unchanged':True,'post_contact_smoothing':False,'history':history}

def validate_pair():
    a=np.array([[-1.,0.,0.]]);b=np.array([[1.,0.,0.]]);c=np.array([[-1.,1.,0.]]);d=np.array([[1.,1.,0.]])
    da,db=pair_correction(a,b,c,d,2.);pa,pb,_,_=closest(a+da,b+da,c+db,d+db);assert np.allclose(np.linalg.norm(pa-pb,axis=1),2.)
    da,db=pair_correction(a,b,c,d,2.,fixed=True);assert np.array_equal(db,np.zeros_like(db));pa,pb,_,_=closest(a+da,b+da,c,d);assert np.allclose(np.linalg.norm(pa-pb,axis=1),2.)
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);theta=.71;cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    rot=np.eye(3)*np.cos(theta)+(1-np.cos(theta))*np.outer(axis,axis)+np.sin(theta)*cross;shift=np.array([3.,-2.,5.])
    ea,eb=pair_correction(a@rot.T+shift,b@rot.T+shift,c@rot.T+shift,d@rot.T+shift,2.,fixed=True)
    error=float(np.linalg.norm(ea-da@rot.T));assert error<1e-12 and np.array_equal(eb,np.zeros_like(eb))
    return {'isolated_pair_final_distance':2.,'required_distance':2.,'fixed_obstacle_displacement':0.,'rigid_transform_delta_error':error}

def curvature_radius_upper_bound(points,radius,subdivisions=16):
    """Conservative cubic Bernstein bound; strict <1 is a non-fold gate.

    Bounds the whole native curve, including between sampled positions. Failure
    to establish a positive speed lower bound is an unproven/failing result.
    """
    p=points.astype('f8');previous=np.concatenate((p[:,:1],p[:,:-1]),axis=1);following=np.concatenate((p[:,1:],p[:,-1:]),axis=1);tangent=(following-previous)*.5
    bez=np.stack((p[:,:-1],p[:,:-1]+tangent[:,:-1]/3,p[:,1:]-tangent[:,1:]/3,p[:,1:]),axis=2)
    def basis(t):return np.array([(1-t)**3,3*(1-t)**2*t,3*(1-t)*t*t,t**3])
    def derivative(t):return np.array([-3*(1-t)**2,3*(1-t)**2-6*(1-t)*t,6*(1-t)*t-3*t*t,3*t*t])
    maximum=0.;unproven=0
    for part in range(subdivisions):
        a=part/subdivisions;b=(part+1)/subdivisions;h=(b-a)/3
        transform=np.stack((basis(a),basis(a)+h*derivative(a),basis(b)-h*derivative(b),basis(b)))
        restricted=np.einsum('ab,cfbd->cfad',transform,bez,optimize=True)
        d1=3*np.diff(restricted,axis=2);d2=2*np.diff(d1,axis=2)
        mid=(d1[:,:,0]+2*d1[:,:,1]+d1[:,:,2])*.25;ref=mid/np.maximum(np.linalg.norm(mid,axis=-1)[:,:,None],1e-30)
        lower=np.sum(d1*ref[:,:,None],axis=-1).min(axis=-1);numerator=np.zeros_like(lower)
        comb2=(1,2,1);comb3=(1,3,3,1)
        for k in range(4):
            cross=np.zeros_like(mid)
            for i in range(3):
                j=k-i
                if 0<=j<2:cross+=(comb2[i]/comb3[k])*np.cross(d1[:,:,i],d2[:,:,j])
            numerator=np.maximum(numerator,np.linalg.norm(cross,axis=-1))
        good=lower>1e-14;unproven+=int((~good).sum())
        if good.any():maximum=max(maximum,float((radius*numerator[good]/lower[good]**3).max()))
    return {'maximum_kappa_radius_upper_bound':maximum,'unproven_subspans':unproven,'subdivisions_per_native_span':subdivisions,'strict_nonfold_pass':bool(unproven==0 and maximum<1.),'interpretation':'Geometric non-fold sufficient bound, not wool bending strain or mechanical calibration.'}

def final_envelope_excess(points,center,radius,core_radius):
    """Check evaluated fibre surfaces against the original circular envelope."""
    p=points.reshape(-1,3).astype('f8');p=p[p[:,2]>.00013]
    a=center[:-1];v=np.diff(center,axis=0);vv=np.sum(v*v,axis=1);maximum=-np.inf
    for start in range(0,len(p),1024):
        q=p[start:start+1024];t=np.clip(np.sum((q[:,None]-a[None])*v[None],axis=2)/vv[None],0,1)
        distance=np.linalg.norm(q[:,None]-a[None]-t[:,:,None]*v[None],axis=2).min(1)
        maximum=max(maximum,float((distance+radius-core_radius).max()))
    return maximum
