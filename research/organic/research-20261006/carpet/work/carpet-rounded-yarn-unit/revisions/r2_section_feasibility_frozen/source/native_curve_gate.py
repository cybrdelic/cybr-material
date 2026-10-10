"""Fail-closed continuous Catmull contact bound for a future qualified field.

This gate is not executed when swept advancement has already failed. It uses
Bezier convex-hull chord errors, so denser sampling alone cannot create a pass.
"""
import time
import numpy as np
from scipy.spatial import cKDTree
from advance_probe import exact_segments,MINIMUM,old

def bounded_chords(curves,tolerance=0.01e-6,max_depth=10):
    p=np.asarray(curves,dtype='f8');previous=np.concatenate((p[:,:1],p[:,:-1]),axis=1);following=np.concatenate((p[:,1:],p[:,-1:]),axis=1);tangent=(following-previous)*.5
    bez=np.stack((p[:,:-1],p[:,:-1]+tangent[:,:-1]/3,p[:,1:]-tangent[:,1:]/3,p[:,1:]),axis=2).reshape(-1,4,3)
    owners=np.repeat(np.arange(len(p)),p.shape[1]-1);out_a=[];out_b=[];out_error=[];out_owner=[]
    for depth in range(max_depth+1):
        chord=bez[:,3]-bez[:,0];error=np.maximum(np.linalg.norm(bez[:,1]-bez[:,0]-chord/3,axis=1),np.linalg.norm(bez[:,2]-bez[:,0]-2*chord/3,axis=1))
        accept=error<=tolerance;out_a.append(bez[accept,0]);out_b.append(bez[accept,3]);out_error.append(error[accept]);out_owner.append(owners[accept]);bez=bez[~accept];owners=owners[~accept]
        if not len(bez):return np.concatenate(out_a),np.concatenate(out_b),np.concatenate(out_error),np.concatenate(out_owner)
        if depth==max_depth:raise RuntimeError('Native cubic chord-error budget exceeded; no qualification')
        first=(bez[:,:-1]+bez[:,1:])*.5;second=(first[:,:-1]+first[:,1:])*.5;middle=(second[:,0]+second[:,1])*.5
        left=np.stack((bez[:,0],first[:,0],second[:,0],middle),axis=1);right=np.stack((middle,second[:,1],first[:,2],bez[:,3]),axis=1)
        bez=np.concatenate((left,right));owners=np.tile(owners,2)

def qualify_native(body,wraps,body_radius,wall_seconds=60):
    start=time.monotonic();curves=np.concatenate((body,wraps));radii=np.r_[np.full(len(body),body_radius),np.full(len(wraps),16e-6)]
    a,b,error,owner=bounded_chords(curves);radius=radii[owner];mid=(a+b)*.5;half=np.linalg.norm(b-a,axis=1)*.5;tree=cKDTree(mid);minimum=np.inf;pairs=0
    for start_index in range(0,len(a),256):
        if time.monotonic()-start>wall_seconds:return {'qualified':False,'reason':'native_contact_wall_budget_exceeded','complete':False}
        end=min(start_index+256,len(a));horizon=radius[start_index:end]+error[start_index:end]+half[start_index:end]+float((radius+error+half).max())+MINIMUM+1e-12
        lists=tree.query_ball_point(mid[start_index:end],horizon,workers=1);i=np.repeat(np.arange(start_index,end),[len(x) for x in lists]);j=np.concatenate([np.asarray(x,dtype='i8') for x in lists]) if len(i) else np.empty(0,dtype='i8')
        keep=(i<j)&(owner[i]!=owner[j]);i=i[keep];j=j[keep]
        if len(i):
            lower=exact_segments(a[i],b[i],a[j],b[j])-radius[i]-radius[j]-error[i]-error[j];minimum=min(minimum,float(lower.min()));pairs+=len(i)
    body_curvature=old.contact.curvature_radius_upper_bound(body,body_radius,16);wrap_curvature=old.contact.curvature_radius_upper_bound(wraps,16e-6,16)
    return {'qualified':bool(minimum>=MINIMUM and body_curvature['strict_nonfold_pass'] and wrap_curvature['strict_nonfold_pass']),'complete':True,'interfibre_clearance_lower_bound_m':minimum,'minimum_accepted_clearance_m':MINIMUM,'maximum_cubic_chord_error_m':float(error.max()),'chord_segments':len(a),'candidate_pairs':pairs,'body_curvature_certificate':body_curvature,'wrap_curvature_certificate':wrap_curvature,'seconds':time.monotonic()-start,'method':'Complete variable-radius broad phase; exact chord distances minus both continuous cubic chord-error bounds and both physical radii. Different fibre owners only; local non-fold separately requires strict curvature times radius below one.'}
